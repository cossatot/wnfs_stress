import numpy as np
import pandas as pd
import halfspace.projections as hsp

def triangle_area(coords):
    ''' returns area in sq km, takes coords in m'''
    ab = np.array(coords[0]) - np.array(coords[1])
    ac = np.array(coords[0]) - np.array(coords[2])
    
    return np.linalg.norm(np.cross(ab, ac)) / 2e6


def triangle_center(coords):
    return ((np.array(coords[0]) + np.array(coords[1]) + np.array(coords[2]))
            / 3).tolist()


def add_tri_centers_area(tri_features):
    for feat in tri_features:
        coords = feat['geometry']['coordinates'][0]
        feat['properties']['area_sq_km'] = triangle_area(coords)
    
        feat['properties']['center'] = triangle_center(coords)


def add_tri_z_values(tri_features, point_df, z_col):
    ''' Operates in place on geojson triangle feature arrays from QGIS'''
    for feat in tri_features:
        coords = feat['geometry']['coordinates'][0]
        zs = point_df.loc[(int(feat['properties']['POINTA']),
                           int(feat['properties']['POINTB']),
                           int(feat['properties']['POINTC']),
                           int(feat['properties']['POINTA'])),
                          z_col].values
        for i, coord in enumerate(coords):
            coord.append(float(zs[i]))


def add_strike_dip(tri_features):
    for feat in tri_features:
        coords = feat['geometry']['coordinates'][0]
    
        coo = [[co[0], co[1], co[2]] for co in coords[:-1]]
    
        s, d = hsp.strike_dip_from_3_xyz(coo[0], coo[1], coo[2])
    
        feat['properties']['strike'], feat['properties']['dip'] = s, d


def add_rake_from_trend(tri_features, trend, rake_err):
    for feat in tri_features:
        prop = feat['properties']
        prop['rake'] = hsp.rake_from_sd_trend(prop['strike'], prop['dip'],
                                              trend, aki_richards=True)
        prop['rake_err'] = rake_err


def feature_to_df(feature, tri_num=0):
    coords = feature['geometry']['coordinates'][0]

    df = pd.DataFrame(columns = ['x', 'y', 'z', 'tri'],
                      index = ['a', 'b', 'c', 'm'])
    df.index.name = 'vertex'

    df.loc[('a','b','c'), ('x', 'y', 'z')] = coords[:3]
    df.loc['m', ('x', 'y', 'z')] = feature['properties']['center']
    df['tri'] = tri_num

    df.set_index('tri', append=True, inplace=True)
    df = df.swaplevel(0,1)

    return df


def tri_dict_to_df(tri_dict, fault_name=None):
    tri_df = pd.concat((feature_to_df(feat, i) 
                        for i, feat in enumerate(tri_dict['features'])))

    if fault_name is not None:
        tri_df['fault_name'] = fault_name
    
    return tri_df
