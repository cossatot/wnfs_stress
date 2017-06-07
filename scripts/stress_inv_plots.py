import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib import colors, gridspec
from tect_stress_functions import *
from stress_plots import rose, scatter_w_marginals
import halfspace.scripts as hss
import halfspace.projections as hsp
#import gdal
import pyproj as pj

husl = colors.ListedColormap(sns.color_palette('husl', 360))

#def get_proj_dict(dem_dataset, epsg):
#    dem_transform = {'upper_left_x': dem_dataset.GetGeoTransform()[0],
#                     'x_res_m' : dem_dataset.GetGeoTransform()[1],
#                     'x_rotation': dem_dataset.GetGeoTransform()[2],
#                     'upper_left_y': dem_dataset.GetGeoTransform()[3],
#                     'y_rotation': dem_dataset.GetGeoTransform()[4],
#                     'y_res_m': dem_dataset.GetGeoTransform()[5],
#                     'n_cols': dem_dataset.RasterXSize,
#                     'n_rows': dem_dataset.RasterYSize}
#    
#    dem_transform['y_res_m'] *= -1
#    dem_transform['east_min'] = dem_transform['upper_left_x']
#
#    dem_transform['east_max'] = (dem_transform['x_res_m'] 
#                                 * dem_transform['n_cols']
#                                 + dem_transform['east_min'])
#    
#    dem_transform['north_max'] = dem_transform['upper_left_y']
#    
#    dem_transform['north_min'] = (-1 * dem_transform['y_res_m'] 
#                                  * dem_transform['n_rows']
#                                  + dem_transform['north_max'])
#
#    wgs84 = pj.Proj(init='epsg:4326')
#    utm = pj.Proj(init='epsg:{}'.format(epsg))
#    
#    dem_transform['lon_min'], dem_transform['lat_min'] = pj.transform(utm, wgs84,
#                                                    dem_transform['east_min'], 
#                                                    dem_transform['north_min'])
#    
#    dem_transform['lon_max'], dem_transform['lat_max'] = pj.transform(utm, wgs84,
#                                                    dem_transform['east_max'], 
#                                                    dem_transform['north_max'])
#    return dem_transform
#
#
#df = gdal.Open('../data/dem/wnfs_dem_utm45_250m_shade.tif')
#data = np.flipud(df.ReadAsArray())
#
#
#proj_d = get_proj_dict(df, 32645)
#
#east_index = np.linspace(proj_d['east_min'], 
#                         proj_d['east_max'], proj_d['n_cols'])
#
#north_index = np.linspace(proj_d['north_min'], 
#                          proj_d['north_max'], proj_d['n_rows'])


T_keeps = pd.read_csv('../results/T_keeps.csv')
T_mle = pd.read_csv('../results/T_best.csv', index_col=0, squeeze=True,
                    header=None)

txx_max = T_mle.txx
tyy_max = T_mle.tyy
txy_max = T_mle.txy
theta_max = T_mle.theta
s1_max = T_mle.s1
s3_max = T_mle.s3


fdf = pd.read_csv('../data/fault_data/fault_df_inv.csv')
fdf = fdf[fdf.depth <-1500]

fdf['tau_rake_misfit'] = hsp.angle_difference(fdf.rake, fdf.tau_rake,
                                              return_abs=False)

fdf['tau_trend'] = trend_from_sd_rake(strike=fdf.strike, dip=fdf.dip,
                                      rake=fdf.tau_rake, aki_richards=True)

fdf['obs_trend'] = trend_from_sd_rake(strike=fdf.strike, dip=fdf.dip,
                                      rake=fdf.rake, aki_richards=True)

fdf['tau_rake_misfit_abs'] = fdf.tau_rake_misfit.abs()

fdf.to_csv('../data/fault_data/fault_df_inv.csv', index=False)

tib = (fdf.fault_name == 'tib_dog')

xx_o, yy_o = trend_to_cart(fdf.obs_trend)

xx_t, yy_t = trend_to_cart(fdf.tau_trend)


gork = fdf[fdf.fault_name == 'galetzka_mainshock']


'''
start plotting
'''

misfit_map = plt.figure(figsize=(16,10))
#plt.pcolormesh(east_index, north_index, data, 
#               cmap='gray')
plt.scatter(fdf.east, fdf.north, c=fdf.tau_rake_misfit.abs(), s=fdf.slip_m*50, 
            lw=0, cmap='viridis')
plt.colorbar()

plt.axis('equal')
plt.gca().set_xlim((-200000.0, 1200000.0))
plt.gca().set_ylim((2900000.0, 3500000.0))
#
#print("xlim: ", plt.gca().get_xlim())
#print("ylim: ", plt.gca().get_ylim())
#
#
#
### Tectonic stresses
#T_fig = plt.figure(figsize=(7,3))
#
#plt.subplot(141, polar=True)
#rose(T_keeps.theta, bidirectional=True, bins=120, lw=0)
#plt.gca().set_yticks([])
#plt.gca().set_xticks(np.radians(np.arange(12) * 30))
#plt.gca().axvline(np.radians(theta_max), color='r')
#plt.gca().axvline(np.radians(theta_max+180), color='r')
#plt.gca().axvline(np.radians(0.074571178026528337), color='g')
#plt.gca().axvline(np.radians(0.074571178026528337+180), color='g')
#plt.title('azimuth of T_max')#, fontsize=24)
#
#plt.subplot(142)
##T_keeps.s1.plot(kind='kde', lw=2)
#plt.hist(T_keeps.s1.values, bins=20, lw=0, normed=True)
#plt.gca().axvline(s1_max, color='r')
#plt.gca().axvline(0.20017105370315322, color='g')
#plt.title('T_max posterior PDF')#, fontsize=24)
#plt.xlabel('T_max, * rho g z')#, fontsize=18)
##plt.gca().tick_params(labelsize=18)
#
#plt.subplot(143)
##T_keeps.s3.plot(kind='kde', lw=2)
#plt.hist(T_keeps.s3.values, bins=20, lw=0, normed=True)
#plt.gca().axvline(s3_max, color='r')
#plt.gca().axvline(-0.00017105370315318745, color='g')
#plt.title('T_min posterior PDF')#, fontsize=24)
#plt.xlabel('T_min, * rho g z')#, fontsize=18)
##plt.gca().tick_params(labelsize=18)
#
#plt.subplot(144)
##T_keeps.s3.plot(kind='kde', lw=2)
#plt.hist((T_keeps.s3 / T_keeps.s1).values, bins=20, lw=0, normed=True)
#plt.gca().axvline(s3_max/s1_max, color='r')
#plt.title('T_min/ T_max posterior PDF')#, fontsize=24)
#plt.xlabel('T_min / T_max')#, fontsize=18)
##plt.gca().tick_params(labelsize=18)





#T_col_fig = plt.figure(figsize=(3.5, 4))
#plt.figure()
#plt.subplot(111, polar=True)
#rose(T_keeps.theta, bidirectional=True, bins=120, lw=0)
#plt.gca().set_yticks([])
#plt.gca().set_xticks(np.radians(np.arange(12) * 30))
#plt.gca().axvline(np.radians(theta_max), color='r')
#plt.gca().axvline(np.radians(theta_max+180), color='r')
#plt.gca().axvline(np.radians(0.074571178026528337), color='g')
#plt.gca().axvline(np.radians(0.074571178026528337+180), color='g')
#plt.title('azimuth of T_max')#, fontsize=24)

#plt.subplot(212)
#plt.figure()

#sns.set(style='white')


gs_ratio = 2
g = (sns.jointplot(x='s1', y='s3', data=T_keeps, 
                   xlim=(0,1), ylim=(-1,1),
                   color='b',
                   alpha=0.01, lw=0, s=10,
                   ratio=gs_ratio,
                   stat_func=None,
                   size=4,
                   space=0.01
                   )
     .set_axis_labels("T_max (* ρ g z)", "T_min (* ρ g z)", fontsize=8))
#g.fig.suptitle('Max and min stress')

gs = gridspec.GridSpec(2, 2, width_ratios=[gs_ratio,1], 
                             height_ratios=[1,gs_ratio])

g.fig.add_subplot(gs[0,1], polar=True)
rose(T_keeps.theta, bidirectional=True, bins=120, lw=0, facecolor='b')
plt.gca().set_yticks([])
plt.gca().set_xticklabels([])
#plt.gca().set_xticks(np.radians(np.arange(12) * 30))
plt.gca().axvline(np.radians(theta_max), color='r', lw=0.5)
plt.gca().axvline(np.radians(theta_max+180), color='r', lw=0.5)
#plt.gca().axvline(np.radians(0.074571178026528337), color='g')
#plt.gca().axvline(np.radians(0.074571178026528337+180), color='g')

g.ax_joint.tick_params(axis='both', labelsize=8)

g.ax_marg_x.axvline(s1_max, lw=0.5)
g.ax_marg_y.axhline(s3_max, lw=0.5)


#plt.savefig('../manuscript/figures/stress_hists.pdf', bbox_inches='tight')

plt.show()

