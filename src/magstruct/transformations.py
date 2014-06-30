from numpy import sin, cos

__all__ = [
    'radec2xyz'
]

def radec2xyz(observation, center):
    t_obs, ra_obs, dec_obs = observation
    t_cen, ra_cen, dec_cen = center

    delta_ra = ra_obs - ra_cen
    
    x = -t_obs * cos(dec_obs) * sin(delta_ra)
    y = t_obs * sin(dec_obs) * cos(dec_cen) \
        - t_obs * cos(dec_obs) * sin(dec_cen) * cos(delta_ra)
    z = t_cen - t_obs * cos(dec_obs) * cos(dec_cen) * cos(delta_ra)

    return x, y, z

def rotate2plane(x, y, z, i, theta):
    x_ = x*cos(theta) + y*sin(theta)
    y_ = -x*sin(theta)*cos(i) + y*cos(theta)*cos(i) - z*sin(i)
    z_ = -x*sin(theta)*sin(i) + y*cos(theta)*sin(i) + z*cos(i)

    return x_, y_, z_
