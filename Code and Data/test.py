import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

# Read the netCDF file
file_path = r'C:\Users\trym7\OneDrive - UiT Office 365\skole\FYS-3001\Home exam passive\amsr2File_201811120311.nc'
netcdf_file = nc.Dataset(file_path)

# Extract required variables
tb_18_7v = netcdf_file["tb19v"][:]  # Brightness temperature at 18.7GHz
lat = netcdf_file["lat"][:]
lon = netcdf_file["lon"][:]

# Define brightness temperature for open water and 100% sea ice
TO = 273.15  # Example value for open water brightness temperature in Kelvin
Ti = 250.0   # Example value for brightness temperature for 100% sea ice in Kelvin

# Calculate observed brightness temperature Tb
#Tb = tb_18_7v  # Assuming tb_18_7v represents observed brightness temperature

# Calculate observed brightness temperature Tb
Tb = tb_18_7v[0]  # Assuming tb_18_7v represents observed brightness temperature for the first time step

# Estimate sea ice concentration CI using Equation (3)
CI = (Tb - TO) / (Ti - TO)

# Set sea ice concentration range to 0-100%
CI[CI < 0] = 0
CI[CI > 1] = 1

# Plot sea ice concentration on a map
plt.figure(figsize=(10, 8))
plt.contourf(lon, lat, CI, levels=np.linspace(0, 1, 101), cmap='jet')
plt.colorbar(label='Sea Ice Concentration')
plt.title('Sea Ice Concentration')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.grid(True)
plt.show()

plt.show()
