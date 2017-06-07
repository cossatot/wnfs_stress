import numpy as np
import pandas as pd
import halfspace.projections as hsp
import halfspace.stress_comps_vectorized as scv


rho = 2700
g = 9.81

def cat_t_priors(num_pts, n_trials, s1_range, s3_range, theta_range, first_iter):
    
    t_priors = sample_T_priors(n_trials, s1_range, s3_range, theta_range)
    run_ind = np.arange(n_trials) + first_iter
    t_priors = np.hstack(( t_priors, run_ind.reshape([n_trials, 1]) ))
    t_priors = np.repeat(t_priors, num_pts, axis=0)
    
    t_prior_df = pd.DataFrame(t_priors, columns=['txx', 'tyy', 'txy', 'iter'])
    
    return t_prior_df
    
    
def sample_T_priors(n_trials, s1_range, s3_range, theta_range):
    s1s = np.random.uniform(s1_range[0], s1_range[1], n_trials)
    s3s = np.random.uniform(s3_range[0], s3_range[1], n_trials) * s1s
    thetas = np.random.uniform(theta_range[0], theta_range[1], n_trials)
    
    xxs = scv.xx_stress_from_s1_s3_theta(s1s, s3s, thetas)
    yys = scv.yy_stress_from_s1_s3_theta(s1s, s3s, thetas)
    xys = scv.xy_stress_from_s1_s3_theta(s1s, s3s, thetas)
    
    del s1s, s3s, thetas  # save some RAM (important for large n_trials)
    
    xxs = xxs.reshape([n_trials, 1])
    yys = yys.reshape([n_trials, 1])
    xys = xys.reshape([n_trials, 1])

    t_priors = np.concatenate((xxs, yys, xys), axis=1)
    
    return t_priors 


def make_mc_df(in_df, n_trials=1, s1_range=(0,2), s3_range=(-1,1),
               theta_range=(0, np.pi), first_iter=0):
    
    num_pts = len(in_df.index)
    
    important_cols = ['strike', 'dip', 'rake', 'depth', 'slip_m', 'xx_stress',
                      'yy_stress', 'zz_stress', 'xy_stress', 'xz_stress',
                      'yz_stress']
    
    mc_df = pd.DataFrame( np.tile(in_df[important_cols].values, [n_trials, 1]),
                         columns=important_cols)
    
    t_prior_df = cat_t_priors(num_pts, n_trials, s1_range, s3_range,
                              theta_range, first_iter)
    
    mc_df = pd.concat((mc_df, t_prior_df), axis=1)
    
    del t_prior_df
    
    mc_df.rename(columns={'xx_stress':'mxx', 'yy_stress':'myy',
                          'zz_stress':'mzz', 'xy_stress':'mxy',
                          'xz_stress':'mxz', 'yz_stress':'myz'}, inplace=True)
    return mc_df


def get_total_stresses(mc_df, rho, g):
    
    mc_df['tau_s'] = scv.strike_shear(strike=mc_df.strike,
                                      dip=mc_df.dip, rho=rho, g=g,
                                      mxx=mc_df.mxx*1e6,
                                      myy=mc_df.myy*1e6,
                                      mzz=mc_df.mzz*1e6,
                                      mxy=mc_df.mxy*1e6,
                                      mxz=mc_df.mxz*1e6,
                                      myz=mc_df.myz*1e6,
                                      txx=mc_df.txx,
                                      tyy=mc_df.tyy,
                                      txy=mc_df.txy,
                                      depth=mc_df.depth*-1)
    
    mc_df['tau_d'] = scv.dip_shear(strike=mc_df.strike,
                                      dip=mc_df.dip, rho=rho, g=g,
                                      mxx=mc_df.mxx*1e6,
                                      myy=mc_df.myy*1e6,
                                      mzz=mc_df.mzz*1e6,
                                      mxy=mc_df.mxy*1e6,
                                      mxz=mc_df.mxz*1e6,
                                      myz=mc_df.myz*1e6,
                                      txx=mc_df.txx,
                                      tyy=mc_df.tyy,
                                      txy=mc_df.txy,
                                      depth=mc_df.depth*-1)
    
    mc_df['tau_rake'] = hsp.get_rake_from_shear_components(
                                                      strike_shear=mc_df.tau_s,
                                                      dip_shear=mc_df.tau_d)
    
    mc_df['rake_misfit_rad'] = np.radians(hsp.angle_difference(mc_df.rake,
                                                               mc_df.tau_rake,
                                                               return_abs=True))
    return mc_df


def get_litho_tect_stresses(mc_df, rho, g):
    
    mc_df['tau_s'] = scv.strike_shear(strike=mc_df.strike,
                                      dip=mc_df.dip, rho=rho, g=g,
                                      mxx=0.,
                                      myy=0.,
                                      mzz=0.,
                                      mxy=0.,
                                      mxz=0.,
                                      myz=0.,
                                      txx=mc_df.txx,
                                      tyy=mc_df.tyy,
                                      txy=mc_df.txy,
                                      depth=mc_df.depth*-1)
    
    mc_df['tau_d'] = scv.dip_shear(strike=mc_df.strike,
                                      dip=mc_df.dip, rho=rho, g=g,
                                      mxx=0.,
                                      myy=0.,
                                      mzz=0.,
                                      mxy=0.,
                                      mxz=0.,
                                      myz=0.,
                                      txx=mc_df.txx,
                                      tyy=mc_df.tyy,
                                      txy=mc_df.txy,
                                      depth=mc_df.depth*-1)
    
    mc_df['tau_rake'] = hsp.get_rake_from_shear_components(
                                                      strike_shear=mc_df.tau_s,
                                                      dip_shear=mc_df.tau_d)
    
    mc_df['rake_misfit_rad'] = np.radians(hsp.angle_difference(mc_df.rake,
                                                               mc_df.tau_rake,
                                                               return_abs=True))
    return mc_df


    
def do_stress_calcs(in_df, n_trials=1, s1_range=(0,3), s3_range=(-1,1),
                    theta_range=(0, np.pi), topo_stress=True, first_iter=0,
                    rho=2700, g=9.81, l_norm=1):
   
    mc_df = make_mc_df(in_df, n_trials=n_trials, s1_range=s1_range,
                       s3_range=s3_range, theta_range=theta_range,
                       first_iter=first_iter)
    
    if topo_stress == True:
        mc_df = get_total_stresses(mc_df, rho, g)
    else:
        mc_df = get_litho_tect_stresses(mc_df, rho, g)
    
    #calculate misfits 
    max_slip = in_df.slip_m.max()
    sum_weights = np.sum(in_df.slip_m) 
    mean_weights = np.sum(max_slip / in_df.slip_m)
    #rake_err = np.pi/9
    kappa = 8.529
    
    #mc_df['weighted_diff'] = mc_df.rake_misfit_rad / (mc_df.slip_m / max_slip)
    mc_df['weighted_diff'] = mc_df.rake_misfit_rad * mc_df.slip_m / max_slip

    mc_df['weighted_diff_sq'] = mc_df.weighted_diff**2
    
    iters = mc_df.groupby('iter')

    del mc_df
    
    
    if l_norm==1: 
        misfit = iters.weighted_diff.mean()
        likelihood = np.exp(kappa * np.cos(misfit) )
    
    elif l_norm==2: 
        misfit = np.sqrt(iters.weighted_diff_sq.sum()) / sum_weights**2
        likelihood = np.exp(kappa * np.cos(misfit) )
    
    like_df = pd.DataFrame(index=iters['iter'].mean(), columns=['txx', 'tyy',
                                                                'txy',
                                                                'misfit',
                                                                'likelihood'])
    like_df['txx'] = iters.txx.mean() 
    like_df['tyy'] = iters.tyy.mean()
    like_df['txy'] = iters.txy.mean() 
    like_df['misfit'] = misfit
    like_df['likelihood'] = likelihood
    
    return like_df


def weighted_misfit(misfits, weight_col):

    weights = weight_col / np.sum(weight_col)

    return np.mean( misfits * weights)


def cart_stresses_to_eigs(xx=0., yy=0., xy=0.):
    
    T = hsp.make_xy_stress_tensor(xx, yy, xy)
    
    vals, vecs = hsp.sorted_eigens(T)
    
    s1, s3 = vals[1], vals[0]
    
    max_x = vecs[0,1] * s1
    max_y = vecs[1,1] * s1
    
    theta = hsp.angle_to_azimuth( np.arctan2(max_y, max_x))
    
    return s1, s3, theta

def pandas_car_to_eigs(row):
    s1, s3, theta = cart_stresses_to_eigs(row['txx'], row['tyy'], row['txy'])
    
    return pd.Series([s1, s3, theta], index=['s1', 's3', 'theta'])



def calc_S(row):
    
    rgz = np.abs(row.depth * 2700 * 9.81)
    
    T = hsp.make_xyz_stress_tensor(row.txx, row.txy, row.tyy) * rgz
    
    M = hsp.make_xyz_stress_tensor(row.mxx, row.myy, row.mzz, row.mxy,
                                   row.mxz, row.myz) * 1e6
    
    L = hsp.make_xyz_stress_tensor(sig_xx=rgz, sig_yy=rgz, sig_zz=rgz)
    
    return T + M + L


def calc_sig3(row):
    
    S = calc_S(row)
    
    vals, vecs = hsp.sorted_eigens(S)
    
    sig_3 = vals[0] / 1e6
    
    return sig_3


def check_hydrofrac(S, phi):
    
    sig_3 = calc_sig3(S)
    
    pressure = 1/3. * hsp.first_tensor_invariant(S)
    
    # phi * pressure > sig_3 -> hydrofrac
    return phi * pressure > sig_3


def get_pressure(mxx=0., myy=0., mzz=0., txx=0., tyy=0., depth=0.,
                rho=2700., g=9.81):

    rgz = np.abs(depth * rho * g) / 1e6
    
    return ((mxx + (txx * rgz) + rgz) 
            + (myy + (tyy * rgz) + rgz)
            + (mzz +  rgz) )  / 3.
    
    
def check_mc_failure(T_df, fdf, phi, mu, tau_dyn_s=0., tau_dyn_d=0., tau_dyn=0.):
    strike_shear = scv.strike_shear()
    dip_shear = scv.dip_shear()
    
    total_shear = np.sqrt((strike_shear**2 + dip_shear**2)) + tau_dyn
    
    sig_n_eff = scv.eff_normal_stress()
    
    return mu < total_shear / sig_n_eff
    
    
def check_intact_failure(T_df, fdf, phi, mu=0.6, tau_dyn_s=0., tau_dyn_d=0., tau_dyn=0.):
    strike_shear = scv.strike_shear()
    dip_shear = scv.dip_shear()
    
    total_shear = np.sqrt((strike_shear**2 + dip_shear**2)) + tau_dyn
    
    sig_n_eff = scv.eff_normal_stress()
    
    return mu < total_shear / sig_n_eff


def sample_fault_property_priors(n_trials, mu_range, phi_range):
    if np.isscalar(mu_range):
        mu = np.ones((n_trials,1)) * mu_range
    elif len(mu_range) == 2:
        mu = np.random.uniform(mu_range[0], mu_range[1], (n_trials,1))
    else:
        mu = sample_mu() # doesn't work
        
    if np.isscalar(phi_range):
        phi = np.ones((n_trials,1)) * phi_range
    elif len(phi_range) == 2:
        phi = np.random.uniform(phi_range[0], phi_range[1], (n_trials,1))
    else:
        phi = sample_phi() # doesn't work
        
    fault_properties = np.concatenate([mu, phi], axis=1)
    
    return fault_properties


def cat_fault_prop_priors(num_pts, n_trials, mu_range, phi_range, first_iter):
    
    fault_prop_priors = sample_fault_property_priors(n_trials, mu_range, phi_range)
    
    #run_ind = np.arange(n_trials) + first_iter
    
    #fault_prop_priors = np.hstack([fault_prop_priors, run_ind.reshape((n_trials, 1))])
    
    fault_prop_priors = np.repeat(fault_prop_priors, num_pts, axis=0)
    
    fault_prop_df = pd.DataFrame(fault_prop_priors, columns=['mu', 'phi'])
    
    return fault_prop_df
    

def make_params_mc_df(in_df, T_df, n_trials=1,  
               mu_range=0.6, phi_range=0.3,
               first_iter=0):
    
    num_pts = len(in_df.index)
    
    important_cols = ['strike', 'dip', 'rake', 'depth', 'slip_m',
                      'xx_stress', 'yy_stress', 'zz_stress',
                      'xy_stress', 'xz_stress', 'yz_stress']
    
    mc_df = pd.DataFrame( np.tile(in_df[important_cols].values, [n_trials, 1]),
                         columns=important_cols)
    
    fault_prop_df = cat_fault_prop_priors(num_pts, n_trials, mu_range, phi_range,
                                         first_iter)
    t_prior_df = cat_t_priors_params_inv(num_pts, n_trials, T_df, first_iter)
    
    mc_df = pd.concat((mc_df, fault_prop_df, t_prior_df), axis=1)
    
    del fault_prop_df
    
    mc_df.rename(columns={'xx_stress':'mxx', 'yy_stress':'myy', 
                          'zz_stress':'mzz',
                          'xy_stress':'mxy', 'xz_stress':'mxz', 
                          'yz_stress':'myz'},
                 inplace=True)
    return mc_df
    
    
def cat_t_priors_params_inv(num_pts, n_trials, T_df, first_iter):
    
    t_priors = sample_T_param_priors(n_trials, T_df)
    run_ind = np.arange(n_trials) + first_iter
    t_priors = np.hstack(( t_priors, run_ind.reshape([n_trials, 1]) ))
    t_priors = np.repeat(t_priors, num_pts, axis=0)
    
    t_prior_df = pd.DataFrame(t_priors, columns=['txx', 'tyy', 'txy', 'iter'])
    
    return t_prior_df
    
    
def sample_T_param_priors(n_trials, T_df):
    
    T_samples = T_df.ix[np.random.choice(T_df.index.values, n_trials), 
                        ('txx', 'tyy', 'txy')]
    
    return T_samples.values


def get_total_stresses_params_inv(mc_df, rho, g):
    
    mc_df['tau_s'] = scv.strike_shear(strike=mc_df.strike,
                                      dip=mc_df.dip, rho=rho, g=g,
                                      mxx=mc_df.mxx*1e6,
                                      myy=mc_df.myy*1e6,
                                      mzz=mc_df.mzz*1e6,
                                      mxy=mc_df.mxy*1e6,
                                      mxz=mc_df.mxz*1e6,
                                      myz=mc_df.myz*1e6,
                                      txx=mc_df.txx,
                                      tyy=mc_df.tyy,
                                      txy=mc_df.txy,
                                      depth=mc_df.depth*-1) / 1e6
    
    mc_df['tau_d'] = scv.dip_shear(strike=mc_df.strike,
                                      dip=mc_df.dip, rho=rho, g=g,
                                      mxx=mc_df.mxx*1e6,
                                      myy=mc_df.myy*1e6,
                                      mzz=mc_df.mzz*1e6,
                                      mxy=mc_df.mxy*1e6,
                                      mxz=mc_df.mxz*1e6,
                                      myz=mc_df.myz*1e6,
                                      txx=mc_df.txx,
                                      tyy=mc_df.tyy,
                                      txy=mc_df.txy,
                                      depth=mc_df.depth*-1) / 1e6
    
    mc_df['pressure'] = get_pressure(mxx=mc_df.mxx, myy=mc_df.myy, mzz=mc_df.mzz,
                                     txx=mc_df.txx, tyy=mc_df.tyy, depth=mc_df.depth)
    
    mc_df['sig_n'] = scv.eff_normal_stress(strike=mc_df.strike,
                                      dip=mc_df.dip, rho=rho, g=g,
                                      mxx=mc_df.mxx*1e6,
                                      myy=mc_df.myy*1e6,
                                      mzz=mc_df.mzz*1e6,
                                      mxy=mc_df.mxy*1e6,
                                      mxz=mc_df.mxz*1e6,
                                      myz=mc_df.myz*1e6,
                                      txx=mc_df.txx,
                                      tyy=mc_df.tyy,
                                      txy=mc_df.txy,
                                      depth=mc_df.depth*-1,
                                      phi=0.) / 1e6
    
    mc_df['sigma_3'] = mc_df.apply(calc_sig3, axis=1)
    
    return mc_df


def intact_failure(row):
    S = calc_S(row)
    
    sig_n = hsp.normal_stress_on_optimal_plane(S)
    tau = hsp.shear_stress_on_optimal_plane(S)[0]
    
    return tau > 0.6 * (sig_n - (row.pressure * row.phi))


def hit_or_miss(mc_df, slip_err=0.1):
    
    slip_cond = ( (mc_df.mc_failure & (mc_df.slip_m > 0.05) ) |
                  (~ mc_df.mc_failure & (mc_df.slip_m < slip_err)))
    
    frac_cond = ~ mc_df.frac
    intact_fail_cond = ~ mc_df.intact_failure
    
    return (slip_cond & frac_cond & intact_fail_cond)


def do_fault_prop_calcs(mht_df, T_df, n_trials=100, mu_range=0.6, phi_range=0.3,
                       first_iter=0):
    
    mc_df = make_params_mc_df(mht_df, T_df, n_trials, mu_range, phi_range, 
                              first_iter)
    
    mc_df = get_total_stresses_params_inv(mc_df, rho, g)
    
    mc_df['mc_failure'] = (np.sqrt( (mc_df.tau_s**2 + mc_df.tau_d**2) ) > 
                       mc_df.mu * (mc_df.sig_n - (mc_df.pressure * mc_df.phi)))

    mc_df['intact_failure'] = mc_df.apply(intact_failure, axis=1)

    mc_df['frac'] = mc_df.phi * mc_df.pressure > mc_df.sigma_3
    
    mc_df['hit'] = hit_or_miss(mc_df)
    
    iters = mc_df.groupby('iter')
    
    like_df = pd.DataFrame(index = iters['iter'].mean(),
                          columns = ['txx', 'tyy', 'txy',
                                     'mu', 'phi', 'num_hits'])
    
    like_df.mu = iters.mu.mean()
    like_df.phi = iters.phi.mean()
    like_df.num_hits = iters.hit.sum()
    like_df.txx = iters.txx.mean()
    like_df.txy = iters.txy.mean()
    like_df.tyy = iters.tyy.mean()
    
    return like_df


def get_tau_rake_and_misfits(in_df, s1=0., s3=0., theta=0., rho=rho, g=g):

    '''
    This function returns a whole dataframe instead of just a rake
    misfit column.  Maybe change return val, maybe change name
    '''

    stress_df = in_df.copy(deep=True)


    stress_df['txx'] = scv.xx_stress_from_s1_s3_theta(s1, s3, theta)
    stress_df['tyy'] = scv.yy_stress_from_s1_s3_theta(s1, s3, theta)
    stress_df['txy'] = scv.xy_stress_from_s1_s3_theta(s1, s3, theta)


    stress_df.rename(columns={'xx_stress':'mxx', 'yy_stress':'myy',
                              'zz_stress':'mzz', 'xy_stress':'mxy',
                              'xz_stress':'mxz', 'yz_stress':'myz'}, 
                     inplace=True)

    stress_df = get_total_stresses(stress_df, rho=rho, g=g)
    stress_df['rake_misfit'] = np.degrees(stress_df.rake_misfit_rad)

    return stress_df[['tau_rake', 'rake_misfit']]


def trend_from_sd_rake(strike=0., dip=0., rake=0., angle='degrees', 
                       aki_richards=False):
    
    if angle == 'degrees':
        strike = np.radians(strike)
        dip = np.radians(dip)
        rake = np.radians(rake)
    
    trend = -1 * (np.arctan( np.cos(dip) * np.tan(rake)) ) + strike
    
    if aki_richards == True:
        
        if np.isscalar(trend):
            if rake > np.pi/2:
                trend -= np.pi
            elif rake < - np.pi/2:
                trend += np.pi
        else:
            trend[rake > np.pi/2] -= np.pi
            trend[rake < -np.pi/2] += np.pi

    if angle == 'degrees':
        trend = np.degrees(trend)# + 180
        trend = hsp.unwrap_angle(trend)

    return trend


def trend_to_cart(trend, length=1.):
    angle = hsp.azimuth_to_angle(trend)

    xx = length * np.cos(angle)
    yy = length * np.sin(angle)

    return xx, yy
