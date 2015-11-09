import numpy as np
import pandas as pd
import time
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import halfspace.scripts as hss
import sys
sys.path.append('/Users/itchy/research/stress/aux_scripts/')
sys.path.append('/Users/itchy/research/stress/nepal_2015/scripts/')
from tect_stress_functions import *
from stress_plots import rose, scatter_w_marginals


wnfs = pd.read_csv('../data/fault_data/wnfs_tris.csv')

wnfs['slip_m'] = 1.5

eqs = pd.read_csv('/Users/itchy/research/stress/nepal_2015/data/fault_data/eq_pts_stresses.csv')

eqs['slip_m'] *= 2

rup = pd.read_csv('/Users/itchy/research/stress/nepal_2015/data/fault_data/fault_pts_stresses.csv',
                  index_col=0)

rup = rup[rup.fault_name == 'MHT_rupture']

fault_df = pd.concat((wnfs, eqs, rup), ignore_index=True, axis=0)

# start tect stress inversion
t_start = time.time()

np.random.seed(69)

n_mcs = 100000

trials_per_loop = 1000

n_loops = int(n_mcs / trials_per_loop)

T_likes = {}

print('doing T inv')
for i in range(n_loops):
    
    first_iter = i * trials_per_loop
    
    T_likes[i] = do_stress_calcs(fault_df, n_trials=trials_per_loop, 
                                s1_range=(0,2), s3_range=(-1,1),
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
np.random.seed(1010)

rand_filter = np.random.uniform(size=n_mcs)

T_keeps = T_like_df[T_like_df.rel_likelihood > rand_filter]
T_keeps[['s1', 's3', 'theta']]  = T_keeps.apply(pandas_car_to_eigs, axis=1)

print('done in', int((t_end - t_start) / 60 ), 'm')

del T_likes

### 

dens_data = T_keeps[['s1', 's3', 'theta']].values.T
kde = gaussian_kde(dens_data)
density = kde(dens_data)

s1_max, s3_max, theta_max = dens_data[:,np.argmax(density)]

txx_max = scv.xx_stress_from_s1_s3_theta(s1_max, s3_max, theta_max)
tyy_max = scv.yy_stress_from_s1_s3_theta(s1_max, s3_max, theta_max)
txy_max = scv.xy_stress_from_s1_s3_theta(s1_max, s3_max, theta_max)


print('xx max =', txx_max)
print('yy max =', tyy_max)
print('xy max =', txy_max)

print('t1' , s1_max)
print('t3' , s3_max)
print('theta' , theta_max)

T_max_series = pd.Series({'txx':txx_max, 'tyy':tyy_max, 'txy':txy_max})
T_max_series.to_csv('../results/T_best.csv')


fdf = fault_df.copy(deep=True)

fdf['txx'] = -fdf.depth * 2700 * 9.81 * txx_max / 1e6
fdf['tyy'] = -fdf.depth * 2700 * 9.81 * tyy_max / 1e6
fdf['txy'] = -fdf.depth * 2700 * 9.81 * txy_max / 1e6

fdf['xx_stress'] += fdf.txx
fdf['yy_stress'] += fdf.tyy
fdf['xy_stress'] += fdf.txy

fdf = hss.resolve_stresses(fdf)
fdf['tau_rake'] = hsp.get_rake_from_shear_components(fdf.tau_ss, fdf.tau_dd)
fdf.to_csv('../data/fault_data/fault_df_inv.csv', index=False)


T_fig = plt.figure(figsize=(5, 12))

plt.subplot(311, polar=True)
rose(T_keeps.theta, bidirectional=True, bins=120, lw=0)
plt.gca().set_yticks([])
plt.gca().set_xticks(np.radians(np.arange(12) * 30))
plt.title('azimuth of T_max', fontsize=24)

plt.subplot(312)
T_keeps.s1.plot(kind='kde', lw=2)
plt.title('T_max posterior PDF', fontsize=24)
plt.xlabel('T_max, * rho g z', fontsize=18)
plt.gca().tick_params(labelsize=18)

plt.subplot(313)
T_keeps.s3.plot(kind='kde', lw=2)
plt.title('T_min posterior PDF', fontsize=24)
plt.xlabel('T_min, * rho g z', fontsize=18)
plt.gca().tick_params(labelsize=18)


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

