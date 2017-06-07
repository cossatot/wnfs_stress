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
from tect_stress_functions import *
from stress_plots import rose, scatter_w_marginals
#import colormaps as cmaps


wnfs = pd.read_csv('../data/fault_data/wnfs_tris.csv')
wnfs = wnfs[wnfs.rake > -180]

wnfs['slip_m'] = 0.1# * wnfs['area_sq_km']

eqs = pd.read_csv('/Users/itchy/research/stress/nepal_2015/data/fault_data/'+
                  'eq_pts_stresses_and.csv')

eqs['fault_name'] = 'fm'

rup = pd.read_csv('../data/fault_data/mht_rup_stresses.csv')

fault_df_inspect = pd.concat((wnfs, 
                              eqs, 
                              rup), ignore_index=True, axis=0) 

rup['slip_m'] *= 0.01 #scale to cm

fault_df = wnfs[wnfs.depth < -1500]



# start tect stress inversion
t_start = time.time()

np.random.seed(69)

n_mcs = int(1e6)

trials_per_loop = 1000

n_loops = int(n_mcs / trials_per_loop)

T_likes = {}

print('doing T inv')
for i in range(n_loops):
    
    first_iter = i * trials_per_loop
    
    T_likes[i] = do_stress_calcs(fault_df, n_trials=trials_per_loop, 
                                 s1_range=(0., 1.), s3_range=(-1.,1.),
                                 theta_range=(-np.pi/2, np.pi/2),
                                 topo_stress=True,
                                 first_iter=first_iter, l_norm=1)
    try:
        if i % (n_loops // 10) == 0:
            print(i * trials_per_loop)
    except ZeroDivisionError:
        pass
    
t_end = time.time()

T_like_df = pd.concat(T_likes.values())

T_like_df['rel_likelihood'] = (T_like_df['likelihood'] / 
                               T_like_df.likelihood.max() )

T_most_like = T_like_df[T_like_df.rel_likelihood == 1.]
np.random.seed(1010)

rand_filter = np.random.uniform(size=n_mcs)

T_keeps = T_like_df[T_like_df.rel_likelihood > rand_filter]

T_keeps_eigs = T_keeps.apply(pandas_car_to_eigs, axis=1)
T_keeps['s1'] = T_keeps_eigs.s1.values
T_keeps['s3'] = T_keeps_eigs.s3.values
T_keeps['theta'] = T_keeps_eigs.theta.values + 180
T_keeps['theta'][T_keeps['theta'] > 360] -= 360

print('done in {0:.1f} m'.format((t_end - t_start)/60))
#del T_likes

T_keeps.to_csv('../results/T_keeps.csv')
print(np.degrees(T_keeps.misfit).describe())

### 
'''
Get the multivariate posterior probability density through kernel density
on the posterior samples. The samples are transformed back into the original
sampling strategy (i.e. T1 uniform, T3 = uniform fraction of T1), the density
is calculated, and then re-transformed to find the MLE absolute values.
'''

dens_data = T_keeps[['s1', 's3', 'theta']].values.T
dens_data_s1_s3_ratio = dens_data.copy()
dens_data_s1_s3_ratio[1,:] /= dens_data_s1_s3_ratio[0,:]

kde = gaussian_kde(dens_data_s1_s3_ratio, bw_method=0.1)
density = kde(dens_data_s1_s3_ratio)

T_mle = T_keeps.iloc[np.argmax(density)]

txx_max = T_mle.txx
tyy_max = T_mle.tyy
txy_max = T_mle.txy
theta_max = T_mle.theta
s1_max = T_mle.s1
s3_max = T_mle.s3


print('xx max =', txx_max)
print('yy max =', tyy_max)
print('xy max =', txy_max)

print('t1' , s1_max)
print('t3' , s3_max)
print('theta' , theta_max)

T_mle.to_csv('../results/T_best.csv')


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
