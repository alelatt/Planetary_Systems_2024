import numpy as np
from matplotlib import pyplot as plt
from decimal import Decimal
from scipy.integrate import trapezoid
from scipy.interpolate import CubicSpline
import h5py
from scipy.optimize import curve_fit
from scipy.stats import binned_statistic

plt.rcParams.update({'font.size': 15})

###### COSTANTI ######
planck = 6.62607015e-34 #J*s
c = 299792458 #m/s
eV = 1.602176634e-19 #J
N_av = 6.02214076e23 #mol-1
GM_E = 3.986004e14 #m3*s-2
Re_E = 6.3781e6 #m
Rp_E = 6.3568e6	#m
R_E = np.cbrt((Re_E**2) * Rp_E)	#m

rho0_H2O = 7.5 #g/m3
m_H2O = 18.015 #g/mol
n0_H2O = rho0_H2O*N_av/m_H2O #m-3
h0_H2O = 2000	#m
E0_H = planck*c/(242.5e-9) #J
m_H = 1.00782503223e-3	#kg/mol
M_H = 2*m_H/N_av	#kg


########## IMPORTING ##########
wave_diss = None
cross_diss = None
with h5py.File('H2O_partial_cross_sections.h5', 'r') as f:
	wave_diss = f['wavelength'][:-1]*1e-9	#m
	cross_diss = f['O'][:-1]*1e-4	#m2

sun2019_data = np.loadtxt("solar_spectrum_2019_10_15.txt")
wave_sun = sun2019_data[:,0]*1e-9 #m
flux_2019 = sun2019_data[:,1]*1e9 #W/m3

sun1989_data = np.loadtxt("solar_spectrum_1989_09_15.txt")
flux_1989 = sun1989_data[:,1]*1e9 #W/m3

unity_data = np.loadtxt('unity_opt_depth.txt')
wave_unity = unity_data[:,0]*1e-9	#m
height_unity = np.round(unity_data[:,1], 1)*1e3 #m


########## BINNING ##########
wave_edges = np.arange(0, np.round(np.max(wave_diss),9) + 1e-9, 1e-9)
wavelength = wave_sun[np.where(wave_sun<np.max(wave_edges))]

wave_unity = wave_unity[np.where(wave_unity<=np.max(wavelength))]
height_unity = height_unity[np.where(wave_unity<=np.max(wavelength))]

if np.any(wavelength-wave_unity != 0):
	print("ERRORE")
	exit()

plt.figure(dpi = 100, layout = 'tight')
plt.plot(wavelength*1e9, height_unity*1e-3, 'k')
plt.xlabel("Wavelength [nm]")
plt.ylabel(r"Height at which $\tau$=1 [km]")
plt.grid()
#plt.show()

binned_diss, out2, out3 = binned_statistic(wave_diss, cross_diss, statistic='mean', bins = wave_edges)
binned_diss[np.where(np.isnan(binned_diss))] = 0
binned_2019, out2, out3 = binned_statistic(wave_sun, flux_2019, statistic='mean', bins = wave_edges)
binned_1989, out2, out3 = binned_statistic(wave_sun, flux_1989, statistic='mean', bins = wave_edges)

print("(F_{1989}-F_{2019})/F_{1989}="+f"{(np.sum(binned_1989)-np.sum(binned_2019))/np.sum(binned_1989):.2f}")
print("(F_{1989}-F_{2019})/F_{1989}(>50nm)="+f"{(np.sum(binned_1989[57:])-np.sum(binned_2019[57:]))/np.sum(binned_1989[57:]):.2f}")
print("(F_{1989}-F_{2019})/F_{1989}(<50nm)="+f"{(np.sum(binned_1989[:57])-np.sum(binned_2019[:57]))/np.sum(binned_1989[:57]):.2f}")

plt.figure(dpi = 100, layout = 'tight')
plt.plot(wavelength*1e9, binned_2019*1e-9, 'k', label = '2019')
plt.plot(wavelength*1e9, binned_1989*1e-9, 'r', label = '1989')
plt.yscale('log')
plt.xlabel("Wavelength [nm]")
plt.ylabel(r"F$_{\lambda, \infty}$ [W m$^{-2}$ nm$^{-1}$]")
plt.legend(loc='lower right')
plt.grid()
#plt.show()


########## ATM DENSITY ##########
atm_data = np.loadtxt('atmos_struct.txt')
altitude_atm = atm_data[:,0]*1e3	#m
dens_atm = atm_data[:,5]	#m-3

spline_atm = CubicSpline(altitude_atm, dens_atm)
altitudes = np.arange(0, 1001, 1)*1e3	#m
dens_atm_spl = spline_atm(altitudes)	#m-3

plt.figure(dpi = 100, layout = 'tight')
plt.plot(altitude_atm*1e-3, dens_atm, 'k.-')
plt.plot(altitudes*1e-3, dens_atm_spl, 'r-')
plt.xlabel("Altitude [km]")
plt.ylabel(r"Atmospheric Density [m$^{-3}$]")
plt.yscale('log')
plt.grid()
#plt.show()

col_dens_atm = np.zeros(len(altitudes))
for i in range(len(altitudes)):
	col_dens_atm[i] = trapezoid(y = dens_atm_spl[i:], x = altitudes[i:])

plt.figure(dpi = 100, layout = 'tight')
plt.plot(altitudes*1e-3, col_dens_atm, 'k-')
plt.xlabel("Altitude [km]")
plt.ylabel(r"Atmospheric Column Density [m$^{-2}$]")
plt.yscale('log')
plt.grid()
#plt.show()

dens_h2o = np.loadtxt('h2o_struct.txt', usecols = 1)	#m-3

########## NUMBER PER ENERGY ##########
number_2019 = np.zeros(len(wavelength))
number_1989 = np.zeros(len(wavelength))
factor = np.zeros(len(wavelength))
integral = np.zeros(len(wavelength))
esc_vel = np.zeros(len(wavelength))
vel = np.zeros(len(wavelength))

for i, wave in enumerate(wavelength):
	index_tau = np.argmin(np.abs(altitudes-height_unity[i]))
	opt_depth = col_dens_atm/col_dens_atm[index_tau]

	if wave == 80.5e-9 or wave == 156.5e-9:
		plt.figure(dpi = 100, layout = 'tight')
		plt.plot(altitudes, opt_depth, 'k-')
		plt.yscale('log')
		plt.xlabel('Altitude [km]')
		plt.ylabel(r'$\tau_{%.0f nm}$' %(wave*1e9))
		plt.grid()
		#plt.show()

	K = (planck*c/wave) - E0_H
	vel[i] = np.sqrt(2*K/M_H)
	esc_vel[i] = np.sqrt(2*GM_E/(R_E+height_unity[i]))

	integral[i] = trapezoid(y = dens_h2o[index_tau:]*np.exp(-opt_depth[index_tau:]), x = altitudes[index_tau:])
	factor[i] = (np.pi/2)*(R_E**2)*binned_diss[i]*wave/(planck*c)

	number_2019[i] = factor[i]*binned_2019[i]*integral[i]	#dN/(dtdl)
	number_1989[i] = factor[i]*binned_1989[i]*integral[i]	#dN/(dtdl)



plt.figure(dpi = 100, layout = 'tight')
plt.plot(wavelength*1e9, factor, 'k.-')
plt.xlabel("Wavelength [nm]")
plt.ylabel(r"$\mathcal{F}_{\lambda}$ [m$^4$ J$^{-1}$]")
plt.yscale('log')
plt.grid()

plt.figure(dpi = 100, layout = 'tight')
plt.plot(wavelength*1e9, integral, 'k.-')
plt.xlabel("Wavelength [nm]")
plt.ylabel(r"$\mathcal{I}_{\lambda}$ [m$^{-2}$]")
plt.yscale('log')
plt.grid()

plt.figure(dpi = 100, layout = 'tight')
plt.plot(wavelength*1e9, binned_diss*1e4, 'k.-')
plt.xlabel("Wavelength [nm]")
plt.ylabel(r"$\sigma_{\lambda}^{O+2H}$ [cm$^2$]")
plt.yscale('log')
plt.grid()
#plt.show()

plt.figure(dpi = 100, layout = 'tight')
plt.plot(wavelength*1e9, (vel-esc_vel)*1e-3, 'k.-')
plt.xlabel("Wavelength [nm]")
plt.ylabel(r"$\Delta$v [km/s]")
plt.yscale('log')
plt.grid()

fig2 = plt.figure(dpi = 100, layout = 'tight')
plt.plot(wavelength*1e9, number_1989*1e-9, 'r-', label = "1989")
plt.plot(wavelength*1e9, number_2019*1e-9, 'k-', label = "2019")
plt.xlabel("Wavelength [nm]")
plt.ylabel(r"dN$_H$/dtd$\lambda$ [s$^{-1}$ nm$^{-1}$]")
plt.yscale('log')
fig2.text(0.22, 0.9, f"Total H produced:\n  1989: {np.sum(number_1989)*M_H*1e-6:.1f} g/s\n  2019: {np.sum(number_2019)*M_H*1e-6:.1f} g/s", ha="left", va="top", bbox=dict(boxstyle='round', fc='blanchedalmond', ec='orange', alpha=1))
plt.legend(loc = "lower right")
plt.grid()
plt.show()

print(f"{np.sum(number_2019)*M_H*1e-9:.2E}, {np.sum(number_1989)*M_H*1e-9:.2E}")