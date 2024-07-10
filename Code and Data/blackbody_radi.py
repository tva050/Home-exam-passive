import numpy as np
import matplotlib.pyplot as plt

plt.style.use('ggplot')

# Constants
c = 3e8
h = 6.626e-34
k = 1.38e-23

# temperature's of the different blackbody's

T_0 = np.array([2000, 1800, 1600, 1400, 1200, 1000])
T_1 = np.array([900, 800, 700, 600, 500])

wavelength = np.linspace(0, 15e-6, 1000) # 0 to 15 micrometers

def planck_for(Wavelength, Temperature):
    S = (2*np.pi*h*c**2)/(Wavelength**5)*(1/(np.exp(c*h/(Wavelength*k*Temperature))-1))
    return S

# Plotting
plt.plot(wavelength, planck_for(wavelength, T_0[0]), label='2000 K')
plt.plot(wavelength, planck_for(wavelength, T_0[1]), label='1800 K')
plt.plot(wavelength, planck_for(wavelength, T_0[2]), label='1600 K')
plt.plot(wavelength, planck_for(wavelength, T_0[3]), label='1400 K')
plt.plot(wavelength, planck_for(wavelength, T_0[4]), label='1200 K')
plt.plot(wavelength, planck_for(wavelength, T_0[5]), label='1000 K')
plt.xlabel('Wavelength (m)')
plt.ylabel(r'Spectral radiance $(W m^{-3})$')
plt.legend()
plt.show()

plt.plot(wavelength, planck_for(wavelength, T_1[0]), label='900 K')
plt.plot(wavelength, planck_for(wavelength, T_1[1]), label='800 K')
plt.plot(wavelength, planck_for(wavelength, T_1[2]), label='700 K')
plt.plot(wavelength, planck_for(wavelength, T_1[3]), label='600 K')
plt.plot(wavelength, planck_for(wavelength, T_1[4]), label='500 K')
plt.xlabel('Wavelength (m)')
plt.ylabel(r'Spectral radiance $(W m^{-3})$')
plt.legend()
plt.show()
