import netCDF4 as nc
import matplotlib.pyplot as plt 
import numpy as np
import cartopy 
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.colors as colors
import matplotlib.colorbar as mcolorbar


path = r'C:\Users\trym7\OneDrive - UiT Office 365\skole\FYS-3001\Home exam passive\amsr2File_201811120311.nc'

netcdf_file = nc.Dataset(path)
print(netcdf_file.variables.keys())

lat = netcdf_file["lat"][:]
lon = netcdf_file["lon"][:]

tb_18_7v = netcdf_file["tb19v"][0]
tb_18_7h = netcdf_file["tb19h"][0]

tb_36_5v = netcdf_file["tb37v"][0]
tb_36_5h = netcdf_file["tb37h"][0]

tb_85_0v = netcdf_file["tb85v"][0]
tb_85_0h = netcdf_file["tb85h"][0]

lat_min, lat_max = 75, 85
lon_min, lon_max = -20, 40

#print(lat.shape)
#print(lon.shape)
#print(tb_18_7v.shape)

def set_to_nan(data, min_val, max_val): # sets the values outside the region to nan
    data = np.where(data < min_val, np.nan, data)
    data = np.where(data > max_val, np.nan, data)
    return data

lat_nan = set_to_nan(lat, lat_min, lat_max)
lon_nan = set_to_nan(lon, lon_min, lon_max)

lat_indices = ~np.isnan(lat_nan)
lon_indices = ~np.isnan(lon_nan)

lat_non_nan = lat[lat_indices]
lon_non_nan = lon[lon_indices]
#print(lat_non_nan.shape)
# makes a mask for the svalbard region, which do not include the nan values,
mask = (lon_nan >= lon_min) & (lon_nan <= lon_max) & (lat_nan >= lat_min) & (lat_nan <= lat_max)  


# extracts the data for the svalbard region
tb_18_7v_svalbard = np.ma.masked_where(~mask, tb_18_7v[:, :])
tb_36_5v_svalbard = np.ma.masked_where(~mask, tb_36_5v[:, :])
tb_85_0v_svalbard = np.ma.masked_where(~mask, tb_85_0v[:, :])

tb_18_7h_svalbard = np.ma.masked_where(~mask, tb_18_7h[:, :])
tb_36_5h_svalbard = np.ma.masked_where(~mask, tb_36_5h[:, :])

# ------------------------------- Task 2 -------------------------------- #
# plotting the data for the three different frequencies
def task_2():
    plt.rcParams["axes.spines.top"]     = False
    plt.rcParams["axes.spines.right"]   = False
    plt.rcParams["axes.spines.left"]    = False
    plt.rcParams["axes.spines.bottom"]  = False
    plt.rcParams["xtick.bottom"]        = False
    plt.rcParams["xtick.labelbottom"]   = False
    plt.rcParams["ytick.left"]          = False
    plt.rcParams["ytick.labelleft"]     = False


    plt.subplot(1, 3, 1)
    plt.imshow(tb_18_7v_svalbard)
    plt.xlim(0, 1609)
    plt.title('18.7GHz')

    plt.subplot(1, 3, 2)
    plt.imshow(tb_36_5v_svalbard)
    plt.xlim(0, 1609)
    plt.title('36.5GHz')

    plt.subplot(1, 3, 3)
    plt.imshow(tb_85_0v_svalbard)
    plt.xlim(0, 1609)
    plt.title('85.0GHz')

    plt.subplots_adjust(bottom=0., right=0.9, top=0.9, wspace=0.1, hspace=0)
    # add a colorbar horizontally at the bottom of all the subplots 
    plt.colorbar(cax=plt.axes([0.14, 0.25, 0.75, 0.03]), orientation='horizontal', label='Brightness Temperature (K)')
    plt.show()

lon_min_geo = np.min(lon)
lon_max_geo = np.max(lon)
lat_min_geo = np.min(lat)
lat_max_geo = np.max(lat)  
#print('The image ranges in longitude from: ', lon_min_geo , 'degrees to ', lon_max_geo, 'degrees.')
#print('The image ranges in latitude from: ', lat_min_geo , 'degrees to ', lat_max_geo, 'degrees.')

# ------------------------------- Task 3 -------------------------------- #

# Constants

Ti = 272    # brightness temperature for 100% sea ice
To = 50     # brightness temperature for open water
 
def sea_ice_concentration(Tb, Ti, To):
    CI = (Tb - To) / (Ti - To)
    CI[CI < 0] = 0
    CI[CI > 1] = 1
    return CI

CI_18_7 = sea_ice_concentration(tb_18_7v_svalbard, Ti, To)
CI_36_5 = sea_ice_concentration(tb_36_5h_svalbard, Ti, To)
CI_85_0 = sea_ice_concentration(tb_85_0v_svalbard, Ti, To)

# calculate the central longitude and latitude of the region of interest
lon_0 = (lon_min_geo + lon_max_geo) / 2
lat_0 = (lat_min_geo + lat_max_geo) / 2
print('Central longitude: ', lon_0)
print('Central latitude: ', lat_0)

orig_proj = ccrs.PlateCarree()
target_proj = ccrs.Stereographic(central_longitude=8.674842834472656, central_latitude=75.75837707519531)

def task_3():
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(1, 1, 1, projection=target_proj)
    ax.set_extent([lon_min_geo-1, lon_max_geo+1, lat_min_geo, lat_max_geo], crs=orig_proj)
    ax.coastlines(resolution='10m')
    ax.add_feature(cfeature.LAND, zorder = 100, edgecolor='black')
    ax.set_xlim(-8.379e+05, 9.029e+05)
    ax.set_ylim(-1.102e+05, 1.142e+06)
    ax.gridlines(draw_labels=True, y_inline=False)

    pcm = ax.pcolormesh(lon, lat, CI_18_7, transform=orig_proj, alpha=1, cmap="viridis") # vmin=0, vmax=1
    #pcm = ax.pcolormesh(lon, lat, CI_85_0, transform=orig_proj, alpha=1, cmap="cool") # vmin=0, vmax=1
    #pcm = ax.pcolormesh(lon, lat, CI_36_5, transform=orig_proj, cmap='viridis', alpha=1) # vmin=0, vmax=1
    cbar = plt.colorbar(pcm, label='Sea Ice Concentration (%)') # ax=ax, shrink=0.6, aspect=30
    cbar.set_ticks([0.6, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90])  # Set the tick locations
    cbar.ax.set_yticklabels(["60%", "65%", "70%", "75%", "80%", "85%", "90%"])  # Set the tick labels
    plt.title('Sea Ice Concentration V-pol 18.7GHz')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.show()




# ------------------------------- Task 4 -------------------------------- #

def task_4():
    
    CI_threshold_water = 0.15
    CI_threshold_new = 0.05
    
    # spectral gradient ratio
    GR = (tb_36_5v_svalbard - tb_18_7v_svalbard) / (tb_36_5v_svalbard + tb_18_7v_svalbard)
    GR_H = (tb_36_5h_svalbard - tb_18_7h_svalbard) / (tb_36_5h_svalbard + tb_18_7h_svalbard)
    print(GR.max(), GR.min(), "GR")
    
    # Polarization gradient ratio
    PR_18 = (tb_18_7v_svalbard - tb_18_7h_svalbard) / (tb_18_7v_svalbard + tb_18_7h_svalbard)
    PR_36 = (tb_36_5v_svalbard - tb_36_5h_svalbard) / (tb_36_5v_svalbard + tb_36_5h_svalbard)
    
    print(PR_18.max(), PR_18.min(), "18")
    print(PR_36.max(), PR_36.min(), "36")
    # Detect the ice edge based on the threshold
    ice_edge_18 = (PR_18 >= 0) & (PR_18 <= 0.20)#PR_threshold
    ice_edge_36 = (PR_36 >= 0) & (PR_36 <= 0.12)#PR_threshold
    ice_edge = ice_edge_18 & ice_edge_36
    
    # Plot the ice edge on a map
    """ plt.contour(lon, lat, ice_edge, levels=[0.5], colors='white', linewidths=1)
    plt.contour(lon, lat, ice_edge_18, levels=[0.5], colors='lime', linewidths=1)
    #plt.contour(lon, lat, ice_edge_36, levels=[0.5], colors='black', linewidths=1)
    plt.pcolormesh(lon, lat, ice_edge, cmap='cool')
    #plt.pcolormesh(lon, lat, PR_18, cmap='cool')
    #plt.pcolormesh(lon, lat, PR_36, cmap='cool')
    plt.colorbar(label='Ice Edge Detection')
    plt.title('Ice Edge Detection')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.show() """
    
    open_water = (GR > -0.051) #& (GR < CI_threshold_water)
    new_ice = (GR > -0.02) & (GR<0.02)   #& (GR <= CI_threshold_new)
    older_ice = (GR>0.04) & (GR < 0.12) #(open_water | new_ice)    
    
    
    open_water_mask = np.ma.masked_where(open_water, GR)
    new_ice_mask = np.ma.masked_where(new_ice, GR)
    older_ice_mask = np.ma.masked_where(~older_ice, GR)
    
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(1, 1, 1, projection=target_proj)
    ax.set_extent([lon_min_geo-1, lon_max_geo+1, lat_min_geo, lat_max_geo], crs=orig_proj)
    ax.coastlines(resolution='10m')
    ax.add_feature(cfeature.LAND, zorder=100, edgecolor='black', facecolor='grey')
    ax.set_xlim(-8.379e+05, 9.029e+05)
    ax.set_ylim(-1.102e+05, 1.142e+06)
    ax.gridlines(draw_labels=True, y_inline=False)

    # Plot the data
    #c1 = plt.pcolormesh(lon, lat, open_water_mask, transform=ccrs.PlateCarree(), cmap='Blues', alpha=1)
    #c2 = plt.pcolormesh(lon, lat, new_ice_mask, transform=ccrs.PlateCarree(), cmap='Greens', alpha=1)
    c3 = plt.pcolormesh(lon, lat, older_ice_mask, transform=ccrs.PlateCarree(), cmap='Oranges', alpha=1)

    plt.contour(lon, lat, ice_edge, levels=[0.5], colors='black', linewidths=1, transform = ccrs.PlateCarree())

    # Create a colorbar for each pcolormesh
    ax1 = fig.add_axes([0.235, 0.05, 0.15, 0.02])
    #cb1 = mcolorbar.ColorbarBase(ax1, cmap=c1.cmap, norm=c1.norm, orientation='horizontal')
    #cb1.set_label("Open water")
    ax1.tick_params(labelbottom=False, labeltop=False, bottom=False, top=False)
    ax2 = fig.add_axes([0.435, 0.05, 0.15, 0.02])
    #cb2 = mcolorbar.ColorbarBase(ax2, cmap=c2.cmap, norm=c2.norm, orientation='horizontal')
    #cb2.set_label("New ice")
    ax2.tick_params(labelbottom=False, labeltop=False, bottom=False, top=False)
    ax3 = fig.add_axes([0.635, 0.05, 0.15, 0.02])
    cb3 = mcolorbar.ColorbarBase(ax3, cmap=c3.cmap, norm=c3.norm, orientation='horizontal')
    cb3.set_label("Multi-year ice")
    ax3.tick_params(labelbottom=False, labeltop=False, bottom=False, top=False)

    fig.suptitle('Sea Ice Classification', fontsize=16, y=1) 
    plt.show()
    

# ------------------------------- Task 5 -------------------------------- #

def task_5():
    # Calculating how large and area is dominated by sea ice below 5cm tickness
    Tb_18_5cm_V = 255 # Example value for brightness temperature for 5cm sea ice in Kelvin 
    Tb_18_5cm_H = 213 # Example value for brightness temperature for 5cm sea ice in Kelvin
    Tb_36_5cm_V = 255 # Example value for brightness temperature for 5cm sea ice in Kelvin
    Tb_36_5cm_H = 205 # Example value for brightness temperature for 5cm sea ice in Kelvin
    
    mask_18_5cm_V = (tb_18_7v_svalbard < Tb_18_5cm_V) & (tb_18_7v_svalbard > 175)
    mask_18_5cm_H = (tb_18_7h_svalbard < Tb_18_5cm_H) & (tb_18_7h_svalbard > 80)
    mask_36_5cm_V = (tb_36_5v_svalbard < Tb_36_5cm_V) & (tb_36_5v_svalbard > 195) 
    mask_36_5cm_H = (tb_36_5h_svalbard < Tb_36_5cm_H) & (tb_36_5h_svalbard > 120) 
    
    mask_5cm_36 = mask_36_5cm_V & mask_36_5cm_H
    mask_5cm_18 = mask_18_5cm_V & mask_18_5cm_H
    mask_5cm = mask_18_5cm_V & mask_18_5cm_H & mask_36_5cm_V & mask_36_5cm_H
    
    
    print('The area dominated by sea ice below 5cm thickness is: ', np.sum(mask_5cm) / mask_5cm.size * 100, '%')
    print("The area dominated by sea ice below 5cm thickness at 36GHz is: ", np.sum(mask_5cm_36)/mask_5cm_36.size, "%")
    print(tb_36_5h_svalbard.max(), tb_36_5h_svalbard.min(), "36.5h")
    print(tb_36_5v_svalbard.max(), tb_36_5v_svalbard.min(), "36.5v")
    
    #tb_18_7_s = np.ma.masked_where((mask_18_5cm_V&mask_18_5cm_H), tb_18_7v_svalbard)
    #tb_36_5_s = np.ma.masked_where(mask_5cm_36, tb_36_5h_svalbard)
    tb18 = np.ma.masked_where(mask_18_5cm_V, tb_18_7v_svalbard)    
    tb36 = np.ma.masked_where(mask_36_5cm_V, tb_36_5v_svalbard)
    
    
    tb_18V = np.ma.masked_where(mask_18_5cm_V, tb_18_7v_svalbard)
    tb_18H = np.ma.masked_where(mask_18_5cm_H, tb_18_7h_svalbard)
    tb_36V = np.ma.masked_where(mask_36_5cm_V, tb_36_5v_svalbard)
    tb_36H = np.ma.masked_where(mask_36_5cm_H, tb_36_5h_svalbard)
    
    # make a mask for sea ice over 5cm thickness
    mask_over5cm_18 = np.ma.masked_where(~mask_5cm, tb_18_7v_svalbard)
    mask_over5cm_36 = np.ma.masked_where(~mask_5cm_36, tb_36_5v_svalbard)
    
    

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(1, 1, 1, projection=target_proj)
    ax.set_extent([lon_min_geo-1, lon_max_geo+1, lat_min_geo, lat_max_geo], crs=orig_proj)
    ax.coastlines(resolution='10m')
    ax.add_feature(cfeature.LAND, zorder=100, edgecolor='black', facecolor='grey')
    ax.set_xlim(-8.379e+05, 9.029e+05)
    ax.set_ylim(-1.102e+05, 1.142e+06)
    ax.gridlines(draw_labels=True, y_inline=False)
    
    c1 = plt.pcolormesh(lon, lat, tb_18V, transform=ccrs.PlateCarree(), cmap='Reds', vmin=0, vmax=213)
    c2 = plt.pcolormesh(lon, lat, tb_18H, transform=ccrs.PlateCarree(), cmap='Reds', vmin=0, vmax=213)
    c3 = plt.pcolormesh(lon, lat, tb_36V, transform=ccrs.PlateCarree(), cmap='Reds', vmin=0, vmax=213)
    c4 = plt.pcolormesh(lon, lat, tb_36H, transform=ccrs.PlateCarree(), cmap='Reds', vmin=0, vmax=213)
    c5 = plt.pcolormesh(lon, lat, mask_over5cm_18, transform=ccrs.PlateCarree(), cmap='Blues_r')
    #c6 = plt.pcolormesh(lon, lat, mask_over5cm_36, transform=ccrs.PlateCarree(), cmap='Blues_r')
    
    
    ax1 = fig.add_axes([0.2, 0.05, 0.15, 0.02])
    cb1 = mcolorbar.ColorbarBase(ax1, cmap=c1.cmap, norm=c1.norm, orientation='horizontal')
    cb1.set_label("<5cm Sea Ice")
    ax1.tick_params(labelbottom=False, labeltop=False, bottom=False, top=False)
    
    
    plt.colorbar(label='Brightness Temperature (K)')
    fig.suptitle('Sea Ice Thickness', fontsize=16, y=1)
    plt.show()


if __name__ == '__main__':
    #task_2()
    task_3()
    #task_4()
    #task_5()
    
    