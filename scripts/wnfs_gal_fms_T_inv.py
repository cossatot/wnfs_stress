import numpy as np
import pandas as pd
import time
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import colors
import seaborn as sns
from scipy.stats import gaussian_kde
import halfspace.scripts as hss
import halfspace.projections as hsp
import sys
sys.path.append('/Users/itchy/research/stress/aux_scripts/')
sys.path.append('/Users/itchy/research/stress/nepal_2015/scripts/')
from tect_stress_functions import *
from stress_plots import rose, scatter_w_marginals
import colormaps as cmaps

#plt.register_cmap(name='viridis', cmap=cmaps.viridis)
#plt.set_cmap(cmaps.viridis)

wnfs = pd.read_csv('../data/fault_data/wnfs_tris.csv')
wnfs = wnfs[wnfs.rake > -180]

wnfs['slip_m'] = 0.1# * wnfs['area_sq_km']

eqs = pd.read_csv('/Users/itchy/research/stress/nepal_2015/data/fault_data/eq_pts_stresses_and.csv')

eqs['fault_name'] = 'fm'

#eqs['slip_m'] *= 4

#rup = pd.read_csv('/Users/itchy/research/stress/nepal_2015/data/fault_data/'
#                 + 'fault_pts_stresses.csv', index_col=0)
#rup = rup[rup.fault_name == 'MHT_rupture']

#rup = pd.read_csv('/Users/itchy/research/stress/nepal_2015/data/fault_data/'
#                  + 'elliott_ng_stresses.csv', index_col=0)

rup = pd.read_csv('../data/fault_data/mht_rup_stresses.csv')

#plt.plot(np.sort(rup.slip_m), range(len(rup.index)) )
#plt.show()

#rup = rup[rup.slip_m > 0.001]

#rup.slip_m

fault_df_inspect = pd.concat((wnfs, 
                              eqs, 
                              rup), ignore_index=True, axis=0) 

#wnfs['slip_m'] = 1.5# * wnfs['area_sq_km']
wnfs_shallow = (wnfs.depth < -1500)
#wnfs.loc[(wnfs.depth > -1500), 'slip_m'] = 0.01

rup['slip_m'] *= 0.01 #scale to cm

fault_df = pd.concat((wnfs, 
                      #eqs, 
                      rup), ignore_index=True, axis=0)

#fault_df = wnfs#[wnfs.fault_name=='tib_dog']

fault_df = wnfs

#fault_df.loc[(fault_df.fault_name != 'galetzka_mainshock


# start tect stress inversion
t_start = time.time()

np.random.seed(69)

n_mcs = int(5e5)

trials_per_loop = 1000

n_loops = int(n_mcs / trials_per_loop)

T_likes = {}

print('doing T inv')
for i in range(n_loops):
    
    first_iter = i * trials_per_loop
    
    T_likes[i] = do_stress_calcs(fault_df, n_trials=trials_per_loop, 
                                 s1_range=(0., 1.), s3_range=(-2.,0.),
                                 theta_range=(-np.pi/2, np.pi/2),
                                 first_iter=first_iter, l_norm=1)
    try:
        if i % (n_loops // 10) == 0:
            print(i * trials_per_loop)
    except ZeroDivisionError:
        pass
    
t_end = time.time()

T_like_df = pd.concat(T_likes.values())

print(T_like_df.likelihood.describe())

T_like_df['rel_likelihood'] = (T_like_df['likelihood'] / 
                               T_like_df.likelihood.max() )

T_most_like = T_like_df[T_like_df.rel_likelihood == 1.]
np.random.seed(1010)

rand_filter = np.random.uniform(size=n_mcs)

T_keeps = T_like_df[T_like_df.rel_likelihood > rand_filter]

T_keep_eigs = T_keeps.apply(pandas_car_to_eigs, axis=1)
T_keeps['s1'] = T_keep_eigs.s1.values
T_keeps['s3'] = T_keep_eigs.s3.values
T_keeps['theta'] = T_keep_eigs.theta.values + 180
T_keeps['theta'][T_keeps['theta'] > 360] -= 360

print('done in {0:.1f} m'.format((t_end - t_start)/60))
del T_likes

T_keeps.to_csv('../results/temp/T_keeps.csv')

### 
'''
Get the multivariate posterior probability density through kernel density
on the posterior samples. The samples are transformed back into the original
sampling strategy (i.e. T1 uniform, T3 = uniform fraction of T1), the density
is calculated, and then re-transformed to find the MLE absolute values.
'''

dens_data = T_keeps[['s1', 's3', 'theta']].values.T
dens_data_s1_s3_ratio = dens_data.copy()
dens_data_s1_s3_ratio[1,:] /= 2.
dens_data_s1_s3_ratio[2,:] /= 180.

kde = gaussian_kde(dens_data_s1_s3_ratio, bw_method=0.1)
density = kde(dens_data_s1_s3_ratio)

#s1_max, s3_max, theta_max = dens_data[:,np.argmax(density)]
T_mle = T_keeps.iloc[np.argmax(density)]

txx_max = T_mle.txx
tyy_max = T_mle.tyy
txy_max = T_mle.txy
theta_max = T_mle.theta
s1_max = T_mle.s1
s3_max = T_mle.s3


#tyy_max = T_most_like['tyy'].values[0]

print('xx max =', txx_max)
print('yy max =', tyy_max)
print('xy max =', txy_max)

#s1_max, s3_max, theta_max = cart_stresses_to_eigs(txx_max, tyy_max, txy_max)

print('t1' , s1_max)
print('t3' , s3_max)
print('theta' , theta_max)

T_mle.to_csv('../results/temp/T_best.csv')


fdf = fault_df_inspect.copy(deep=True)

fdf['txx'] = -fdf.depth * 2700 * 9.81 * txx_max / 1e6
fdf['tyy'] = -fdf.depth * 2700 * 9.81 * tyy_max / 1e6
fdf['txy'] = -fdf.depth * 2700 * 9.81 * txy_max / 1e6

fdf['xx_stress'] += fdf.txx
fdf['yy_stress'] += fdf.tyy
fdf['xy_stress'] += fdf.txy

fdf = hss.resolve_stresses(fdf)
fdf['tau_rake'] = hsp.get_rake_from_shear_components(fdf.tau_ss, fdf.tau_dd)

fdf['tau_rake_misfit'] = hsp.angle_difference(fdf.rake, fdf.tau_rake,
                                              return_abs=False)

fdf['tau_trend'] = trend_from_sd_rake(strike=fdf.strike, dip=fdf.dip,
                                      rake=fdf.tau_rake, aki_richards=True)

fdf['obs_trend'] = trend_from_sd_rake(strike=fdf.strike, dip=fdf.dip,
                                      rake=fdf.rake, aki_richards=True)
tib = (fdf.fault_name == 'tib_dog')

xx_o, yy_o = trend_to_cart(fdf.obs_trend)

xx_t, yy_t = trend_to_cart(fdf.tau_trend)

fdf['tau_rake_misfit_abs'] = np.abs(fdf['tau_rake_misfit'])

fdf.to_csv('../data/fault_data/fault_df_inv.csv', index=False)

#### plotting

misfit_map = plt.figure(figsize=(16,10))
plt.scatter(fdf.east, fdf.north, c=fdf.tau_rake_misfit.abs(), s=fdf.slip_m*50, 
            lw=0, cmap=cmaps.viridis)
plt.colorbar()
plt.axis('equal')

husl = colors.ListedColormap(sns.color_palette('husl', 360))



# tibrikot results map
tib_map = plt.figure(figsize=(16,10))

n_skip = 3

plt.subplot(221)
plt.title('rake')
plt.scatter(fdf[tib].east, fdf[tib].north, c=fdf[tib].rake, lw=0,
            cmap = husl, vmin=-180, vmax=180)
plt.colorbar()

plt.quiver(fdf[tib].east[::n_skip], fdf[tib].north[::n_skip], 
           xx_o[tib][::n_skip], yy_o[tib][::n_skip])

#plt.quiver(fdf[tib].east[::n_skip], fdf[tib].north[::n_skip], 
#           trend_to_cart(trend_from_sd_rake(fdf.strike, fdf.dip, fdf.rake)))

plt.subplot(222)
plt.title('pred. rake')
plt.scatter(fdf[tib].east, fdf[tib].north, c=fdf[tib].tau_rake, lw=0,
            vmin=-180, vmax=180, cmap=husl)
plt.colorbar()

plt.subplot(223)
plt.title('rake misfit')
plt.scatter(fdf[tib].east, fdf[tib].north, c=fdf[tib].tau_rake_misfit.abs(), 
            lw=0, cmap=cmaps.viridis, vmin=0, vmax=180)
plt.colorbar()

plt.subplot(224)
plt.title('pred. trend')
plt.scatter(fdf[tib].east, fdf[tib].north, c=fdf[tib].tau_trend, 
            lw=0, cmap=husl, vmin=0, vmax=359)
plt.colorbar()


#tib = 

plt.quiver(fdf[tib].east[::n_skip], fdf[tib].north[::n_skip], 
           xx_t[tib][::n_skip], yy_t[tib][::n_skip])

## Gurla
gur = (fdf.fault_name == 'gurla')

G_misfit = plt.figure(figsize=(12,12))

plt.scatter(fdf[gur].east, fdf[gur].north, c=fdf[gur].tau_rake_misfit.abs(), 
            lw=0, cmap=cmaps.viridis)

plt.quiver(fdf[gur].east[::n_skip], fdf[gur].north[::n_skip], 
           xx_o[gur][::n_skip], yy_o[gur][::n_skip], color='k')

plt.quiver(fdf[gur].east[::n_skip], fdf[gur].north[::n_skip], 
           xx_t[gur][::n_skip], yy_t[gur][::n_skip], color='r')

plt.axis('equal')

## Gorkha
gal = (fdf.fault_name == 'galetzka_mainshock')

G_misfit = plt.figure(figsize=(12,12))

plt.scatter(fdf[gal].east, fdf[gal].north, c=fdf[gal].tau_rake_misfit.abs(), 
            lw=0, cmap=cmaps.viridis)

plt.quiver(fdf[gal].east[::n_skip], fdf[gal].north[::n_skip], 
           xx_o[gal][::n_skip], yy_o[gal][::n_skip], color='k')

plt.quiver(fdf[gal].east[::n_skip], fdf[gal].north[::n_skip], 
           xx_t[gal][::n_skip], yy_t[gal][::n_skip], color='r')

plt.axis('equal')



## Tectonic stresses
T_fig = plt.figure(figsize=(12,12))

plt.subplot(221, polar=True)
rose(T_keeps.theta, bidirectional=True, bins=120, lw=0)
plt.gca().set_yticks([])
plt.gca().set_xticks(np.radians(np.arange(12) * 30))
plt.gca().axvline(np.radians(theta_max), color='r')
plt.gca().axvline(np.radians(theta_max+180), color='r')
plt.gca().axvline(np.radians(0.074571178026528337), color='g')
plt.gca().axvline(np.radians(0.074571178026528337+180), color='g')
plt.title('azimuth of T_max')#, fontsize=24)

plt.subplot(222)
#T_keeps.s1.plot(kind='kde', lw=2)
plt.hist(T_keeps.s1.values, bins=20, lw=0, normed=True)
plt.gca().axvline(s1_max, color='r')
plt.gca().axvline(0.20017105370315322, color='g')
plt.title('T_max posterior PDF')#, fontsize=24)
plt.xlabel('T_max, * rho g z')#, fontsize=18)
#plt.gca().tick_params(labelsize=18)

plt.subplot(223)
#T_keeps.s3.plot(kind='kde', lw=2)
plt.hist(T_keeps.s3.values, bins=20, lw=0, normed=True)
plt.gca().axvline(s3_max, color='r')
plt.gca().axvline(-0.00017105370315318745, color='g')
plt.title('T_min posterior PDF')#, fontsize=24)
plt.xlabel('T_min, * rho g z')#, fontsize=18)
#plt.gca().tick_params(labelsize=18)

plt.subplot(224)
#T_keeps.s3.plot(kind='kde', lw=2)
plt.hist((T_keeps.s3 / T_keeps.s1).values, bins=20, lw=0, normed=True)
plt.gca().axvline(s3_max/s1_max, color='r')
plt.title('T_min/ T_max posterior PDF')#, fontsize=24)
plt.xlabel('T_min / T_max')#, fontsize=18)
#plt.gca().tick_params(labelsize=18)



#plt.figure()
scatter_w_marginals(T_keeps.s1, T_keeps.s3, #c=T_keeps.theta,
                    bins=20, xlabel='s1', ylabel='s3', alpha=0.1, lw=0,
                    #cmap='viridis', 
                    plot_max=True)


scatter_w_marginals(T_keeps.s1, T_keeps.theta, #c=T_keeps.s3,
                    bins=20, xlabel='s1', ylabel='theta', alpha=0.1, lw=0,
                    #cmap='viridis', 
                    plot_max=True)

scatter_w_marginals(T_keeps.s3, T_keeps.theta, #c=T_keeps.s3,
                    bins=20, xlabel='s3', ylabel='theta', alpha=0.1, lw=0,
                    #cmap='viridis', 
                    plot_max=True)



#fig3d = plt.figure()
#ax3 = fig3d.add_subplot(111, projection='3d')
#ax3.scatter(T_keeps.s1, T_keeps.s3/T_keeps.s1, T_keeps.theta)
#ax3.set_xlabel('s1')
#ax3.set_ylabel('s3/s1')
#ax3.set_zlabel('theta')

#plt.show(block=False)


# print('doing fault params inv')
# 
# 
# mht_df = fault_df[fault_df.fault_name == 'galetzka_mainshock']
# 
# t_start = time.time()
# 
# np.random.seed(69)
# 
# n_mcs = 10000
# 
# trials_per_loop = 100
# 
# n_loops = int(n_mcs / trials_per_loop)
# 
# P_likes = {}
# 
# for i in range(n_loops):
#     
#     first_iter = i * trials_per_loop
#     
#     P_likes[i] = do_fault_prop_calcs(mht_df, T_keeps, n_trials=trials_per_loop,
#                               mu_range=[0., 1.], phi_range=[0., 1.],
#                                    first_iter=first_iter)
#     if i % 10 == 0:
#         print(i * trials_per_loop)
#     
# t_end = time.time()
# 
# P_like_df = pd.concat(P_likes.values())
# 
# n_pts = len(mht_df)
# 
# p_filter = np.random.uniform(size=n_mcs)
# 
# P_keeps = P_like_df[(P_like_df.num_hits / n_pts) > p_filter]
# 
# print('done in', int((t_end - t_start) / 60 ), 'm')
# 
# plt.figure()
# pd.tools.plotting.scatter_matrix(P_keeps, alpha=0.05)

plt.show()

