import numpy as np
from matplotlib import pyplot as plt
from decimal import Decimal
from scipy.integrate import trapezoid
from scipy.interpolate import CubicSpline
import h5py
from scipy.optimize import curve_fit
from scipy.stats import binned_statistic

'''
xdim = 664#427
ydim = 469#325
wavedim = 320#300
absdim = 250#250

def round_to_half(array):
	return (np.ceil(array) + np.floor(array))/2


wave, height = np.loadtxt("unity.txt", unpack=True)
wave = round_to_half(wave)
height = ydim - height
height = round_to_half(height)

wave = wave*wavedim/xdim
height = height*absdim/ydim

wavelengths = np.arange(0.5, 200.5, 1)
alt = np.interp(wavelengths, wave, height)

results = np.array([wavelengths, alt]).T
print(results)

plt.figure()
plt.plot(wave, height, 'k')
plt.plot(wavelengths, alt, 'k+')

np.savetxt("unity_opt_depth.txt", results, fmt='%.3E')
'''

N_av = 6.02214076e23 #mol-1
rho0_H2O = 7.5 #g/m3
m_H2O = 18.015 #g/mol
n0_H2O = rho0_H2O*N_av/m_H2O #m-3
h0_H2O = 2	#km


xdim = 1585
ydim = 1122
dens_log_low = -1
dens_log_high = 5
height_log = 2

dens, height = np.loadtxt("h2o.txt", unpack=True)
height = ydim - height

dens = dens*(dens_log_high-dens_log_low)/xdim + dens_log_low
height = height*height_log/ydim

lin = np.arange(2, np.max(10**height), 1)
P = np.interp(np.log10(lin), height, dens)
lin = np.append(np.array([1]), lin)
P = np.append(np.array(dens[0]), P)
dspline = 10**(P)


def fit_func(x, a, b):
	return (a*x) + b

popt, pcov = curve_fit(fit_func, height[-6:], dens[-6:])
errs = np.sqrt(pcov.diagonal())
#print(popt, errs, popt/errs)

lin2 = np.arange(np.max(lin), 1001, 1)
dfit = 10**(fit_func(np.log10(lin2), *popt))

dconst = np.ones(len(lin2))*dspline[-1]

dens = 10**(dens)
height = 10**(height)


plt.plot(lin2, dfit*1e-6, '*b')
plt.plot(lin, dspline*1e-6, '+r')
plt.plot(height, dens*1e-6, '.k')
#plt.xscale('log')
plt.yscale('log')
plt.show()


total_x = np.concatenate((lin, lin2))
total_y = np.concatenate((dspline, dconst))*1e-6
total_y_fit = np.concatenate((dspline, dfit))*1e-6
#total_y[np.where(total_y<2e-6)] = 2e-6


atm_data = np.loadtxt('atmos_struct.txt')
altitude_atm = atm_data[:,0]	#km
dens_atm = atm_data[:,5]	#m-3

spline_atm = CubicSpline(altitude_atm, dens_atm)
dens_atm_spl = spline_atm(total_x)	#m-3

dens_h2o = n0_H2O*np.exp(-total_x/h0_H2O)
ratio = dens_h2o/dens_atm_spl
ratio_corr = dens_h2o/dens_atm_spl
ratio_corr[ratio_corr<2e-6] = 2e-6

plt.plot(total_x[:251], total_y_fit[:251], '+g')
plt.plot(total_x[:251], total_y[:251], '+r')
plt.plot(total_x[:251], ratio[:251], '*b')
plt.plot(total_x[:251], ratio_corr[:251], '.b')
plt.plot(height, dens*1e-6, '.k')
plt.xlabel("Height [km]")
plt.ylabel(r"n$_{H2O}$/n$_{atm}$")
#plt.xscale('log')
plt.yscale('log')
plt.show()

plt.rcParams.update({'font.size': 15})
fig, ax1 = plt.subplots(dpi = 100, layout = 'tight')
ax2 = ax1.twinx()
ax1.plot(total_x, dens_atm_spl, 'k-')
ax2.plot(lin, dspline*1e-6*dens_atm_spl[:len(dspline)], 'k-', label = "Standard Atmosphere")
ax2.plot(lin, dspline*1e-6*dens_atm_spl[:len(dspline)], 'r-', label = "Interpolation")
ax2.plot(lin2, dfit*1e-6*dens_atm_spl[len(dspline):], 'r-.', label = "Fit Extension")
ax1.set_xlabel("Height [km]")
ax1.set_ylabel(r"n$_{atm}$")
ax2.set_ylabel(r"n$_{H_2O}$")
ax1.set_xscale('log')
ax1.set_yscale('log')
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.legend(loc='lower left')
ax1.grid()
ax2.tick_params(axis='y', colors='red')
ax2.spines['right'].set_color('red')
ax2.yaxis.label.set_color('red')
plt.show()
exit()


results = np.array([total_x, total_y_fit*dens_atm_spl]).T
np.savetxt("h2o_struct.txt", results, fmt=('%d','%.3E'))