import numpy as np
from matplotlib import pyplot as plt
from decimal import Decimal
from scipy.integrate import trapezoid
from scipy.interpolate import CubicSpline
import h5py
from scipy.optimize import curve_fit
from scipy.stats import binned_statistic


###### COSTANTI ######
h = 6.62607015e-34 #J*s
c = 299792458 #m/s
eV = 1.602176634e-19 #J
N_av = 6.02214076e23 #mol-1
rho0_H2O = 7.5 #g/m3
m_H2O = 18.015 #g/mol
n0_H2O = rho0_H2O*N_av/m_H2O #m-3
h0_H2O = 2000 #m
'''
E_min = 5.113*eV
nu_min = E_min/h
print(f'Minimum Frequency: {Decimal(nu_min):.2E} Hz')
lam_max = c/nu_min
print(f'Maximum Wavelength: {lam_max*1e9: .2f} nm')

sigma = 1.511e-29

data = np.loadtxt('atmos_struct.txt')
alt = data[:,0]*1e3
dens = data[:,5]

plt.figure()
plt.plot(alt,dens)
plt.yscale('log')
plt.xlabel('Altitude [m]')
plt.ylabel(r'Density [m$^{-3}$]')
plt.show()

tau = np.zeros(len(alt))
for i in range(len(alt)):
	tau[i] = sigma*trapezoid(y=dens[i:], x=alt[i:])

P = CubicSpline(np.flip(tau[:11]), np.flip(alt[:11]))
x = np.linspace(tau[10], tau[0], 1000)
h = P(1)

print(f"h(tau(nu_min)=1) = {h:.0f} m")

plt.figure()
plt.plot(alt, tau,'b.-')
plt.plot(P(x), x, 'r-')
plt.plot(h, 1, 'kx')
plt.yscale('log')
plt.xlabel('Altitude [m]')
plt.ylabel(r'$\tau$')
plt.grid()
plt.show()

n_H2O = n0_H2O*np.exp(-h/h0_H2O)
print(f"n_H2O(h) = {Decimal(n_H2O):.2E} m^-3")
'''

'''
a = np.array([0, 0.1, 0.2, 0.3, 0.5, 0.7, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 2, 2.4, 2.8, 3.5, 4, 4.5, 5])
b = np.random.uniform(0,1,len(a))
print(b)
max_bin = np.max(a[1:]-a[:-1])
print(max_bin)
edges = np.arange(0, np.max(a)+max_bin, max_bin)
bins = (edges[1:]+edges[:-1])/2
print(bins)
binned, out2, out3 = binned_statistic(a, b, statistic='mean', bins = edges)

plt.figure()
plt.plot(a, b)
plt.plot(bins, binned)
plt.show()
'''
'''
def FitFunc(x, a, b, c, d, e):
	return a*(x**4) + b*(x**2) + c*(x**-2) + d*(x**-4) + e

wv, sc = np.loadtxt('scatt.txt', unpack=True)
popt, pcov = curve_fit(FitFunc, wv, sc)
errs = np.sqrt(pcov.diagonal())
print(popt/errs)

plt.figure()
plt.plot(wv, sc)
lin = np.linspace(np.min(wv), np.max(wv), 100)
plt.plot(lin, FitFunc(lin, *popt))
plt.show()
'''
wav = None
diss = None

with h5py.File('H2O_partial_cross_sections.h5', 'r') as f:
	wav = f['wavelength'][:]*1e-9 #m
	diss = f['photodissociation'][:]*1e-4 #m^2
	absorp = f['photoabsorption'][:]*1e-4 #m^2
	ioniz = f['photoionisation'][:]*1e-4 #m^2
	oxy = f['O'][:]*1e-4 #m^2
	hydroxy = f['OH'][:]*1e-4 #m^2

	plt.figure()
	plt.plot(wav*1e9, absorp, 'r', label = "Photoabsorption")
	plt.plot(wav*1e9, ioniz, 'b', label = "Photoionisation")
	plt.plot(wav*1e9, diss, 'k', label = "Photodissociation")
	plt.xlabel("Wavelength [nm]")#
	plt.ylabel(r"$\sigma$ [m$^2$]")
	plt.legend()

	plt.figure()
	plt.plot(wav*1e9, oxy, 'r', lw = 2, label = "O + H2")
	plt.plot(wav*1e9, hydroxy, 'k', label = "H + OH")
	plt.xlabel("Wavelength [nm]")#
	plt.ylabel(r"$\sigma$ [m$^2$]")
	plt.legend()
	plt.show()
exit()
plt.figure()
plt.plot(wav)

plt.figure()
plt.plot(wav*1e9, diss, 'k')
plt.xlabel("Wavelength [nm]")
plt.ylabel(r"$\sigma_{diss}$ [m$^2$]")
plt.show()

sun2019_data = np.loadtxt("solar_spectrum_2019_10_15.txt")
sun2019_wav = sun2019_data[:,0]*1e-9 #m
sun2019_flux = sun2019_data[:,1]*1e-9 #W/(m^3)

sun1989_data = np.loadtxt("solar_spectrum_1989_09_15.txt")
sun1989_wav = sun1989_data[:,0]*1e-9 #m
sun1989_flux = sun1989_data[:,1]*1e-9 #W/(m^3)


plt.figure()
plt.plot(sun2019_wav*1e9, sun2019_flux*1e9, 'k.-', label = '2019')
plt.plot(sun1989_wav*1e9, sun1989_flux*1e9, 'r.-', label = '1989')
plt.xlabel("Wavelength [nm]")
plt.ylabel(r"F$_{\lambda}$ [W/m$^2$nm]")
plt.legend()
plt.show()