# -*- coding: utf-8 -*-
"""
sphere_rotate.py : Module defining the spherical coordinate rotations.

function rotate_sphere(): rotation of spherical coordinates (ra, dec) about Y axis 
    which is normal to both Z {dec= 90 deg} and X {ra =0, dec =0}
    by angle alpha in [0, 180 deg].
class SpherePolarRotate: inherits matplotlib.transforms.Transform
    a polar spherical rotation by angle alpha along the zero-meridian
class SphereAziRotate:
    a spherical rotation about the pole by angle phi along the equator
class CastModulo: 
    transformation that casts values into a defined range (e.g. (-180,+180) deg)

Version:  v1.3 (2024-10-21)
Created on Thu May 23 18:42:02 2024

@author: jgazer (Jonathan Gasser, SwRI)
"""

import numpy as np
import matplotlib.transforms as trf
import matplotlib.path as pth

""" rotate_sphere(): 
    rotation of spherical coordinates (ra, dec) about Y axis 
    which is normal to both Z {dec= 90 deg} and X {ra =0, dec =0}
    by angle alpha in [0, 180 deg].
All input angles in deg.

Cartesian:
    x = [ cos ra cos dec ; -sin ra cos dec ; sin dec]
    x'= Ry(a)*x = [cosa 0 sina ; 0 1 0 ; -sina 0 cosa] * [x; y; z]
      = [ cos ra cos dec cosa + sin dec sina ;
         -sin ra cos dec ;
         -cos ra cos dec sina + sin dec cosa ]
Spherical:
    dec' = arsin[ sin dec cosa -cos ra cos dec sina ]
    ra'  = artan[ sin ra cos dec / ( cos ra cos dec cosa + sin ra sina ) ]
            + 180 * if( {denominator} <0 ) * -sgn( sin ra sin dec )
Version:  v1.0 (2024-03-15)
"""
def rotate_sphere(ra, dec, alpha, radians=False):
    '''
    TODO.
    '''
    
    ra = np.array(ra)
    dec = np.array(dec)
    ra_shape = ra.shape
    dec_shape = dec.shape
    
    if ra.ndim >1:
        ra = np.reshape(ra, (1,-1))[0]
    if dec.ndim >1:
        dec = np.reshape(dec, (1,-1))[0]
    
    if not radians:
        sina = np.sin( np.radians(alpha) )
        cosa = np.cos( np.radians(alpha) )
        sin_ra = np.sin( np.radians(ra) )
        cos_ra = np.cos( np.radians(ra) )
        sin_dec = np.sin( np.radians(dec) )
        cos_dec = np.cos( np.radians(dec) )
    
    sin_dec_out = sin_dec * cosa -cos_ra * cos_dec * sina
    # prevent machine precision level errors when calculate asin(1+epsilon) in case alfa=45°, dec=45°, ra=0°
    sin_dec_out[ np.where( np.abs(sin_dec_out -1) <1.e-15 ) ] = 1
    sin_dec_out[ np.where( np.abs(sin_dec_out +1) <1.e-15 ) ] = -1
    dec_out = np.arcsin( sin_dec_out )
    
    denom = cos_ra *cos_dec *cosa + sin_dec *sina
    denom[ np.where(denom ==0)] = 1e-99 #avoid zero division for ra= +-pi/2
    ra_out = np.arctan( sin_ra *cos_dec /denom )
    
    # deconvolution of angles from [-pi/2,pi/2] into full 4pi range
    for i in np.where( denom <0 )[0] :
         ra_out[i] = ra_out[i] +np.pi *( -1 +2* int(sin_ra[i] *sin_dec[i] <0) )
    
    if not radians:
        dec_out= np.degrees(dec_out)
        ra_out= np.degrees(ra_out)
        ra_out[ np.where( np.abs(ra_out)< 1e-15 )] =360
        ra_out = ra_out %360
    else:
        ra_out = ra_out %(2*np.pi)
    
    ra_out = np.reshape(ra_out, ra_shape)
    dec_out = np.reshape(dec_out, dec_shape)
    return ra_out, dec_out
# END OF FUNCTION


class SpherePolarRotate(trf.Transform):
    """
    Performs a rotation of spherical coordinates along the zero meridian (azi=0).
    """
    #operates in radians internally
    
    input_dims = 2
    output_dims = 2
    
    is_separable = False
    has_inverse = True
    
    def __init__(self, angle_deg=None, radians=False, angle_rad=None, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._in_radians = bool(radians)
        if angle_deg:
            self._angle = np.radians(angle_deg)
        elif angle_rad:
            self._angle = angle_rad
        else:
            self._angle = 0 # radians
            
        self._sina = np.sin(self._angle)
        self._cosa = np.cos(self._angle)

    
    def transform(self, values):
        ''' Apply this transformation on the given array of *values*.

        Parameters
        ----------
        values : array-like
            The input values as an array of length :attr:`input_dims` or
            shape (N, :attr:`input_dims`).

        Returns
        -------
        array
            The output values as an array of length :attr:`output_dims` or
            shape (N, :attr:`output_dims`), depending on the input.
        '''
        
        values = np.asanyarray(values)
        ndim = values.ndim
        values = values.reshape((-1, self.input_dims))
        
        # Convert input to radians
        if not self._in_radians:
            values = np.radians(values)
        
        sin_azi = np.sin( values[0:,0] )
        cos_azi = np.cos( values[0:,0] )
        sin_elv = np.sin( values[0:,1] )
        cos_elv = np.cos( values[0:,1] )
        
        # Calculate rotated elevation:
        sin_elv_out = sin_elv * self._cosa - cos_azi *cos_elv *self._sina
        # prevent machine precision level errors when calculate asin(1+epsilon) for angle=pi/4, elv=pi/4, azi=0
        sin_elv_out[ np.where( np.abs(sin_elv_out -1) <1.e-15 ) ] = 1
        sin_elv_out[ np.where( np.abs(sin_elv_out +1) <1.e-15 ) ] = -1
        elv_out = np.arcsin( sin_elv_out )
        
        # Calculate rotated azimuth:
        denom = cos_azi *cos_elv *self._cosa + sin_elv *self._sina
        denom[ np.where(denom ==0) ] = 1.e-99 # avoid zero division for azi = +-pi/2
        azi_out = np.arctan( sin_azi *cos_elv /denom )
        # Resolve degeneracy of angles into full 4pi range
        for i in np.where( denom <0 )[0] :
             azi_out[i] = azi_out[i] + np.pi *( -1 +2* int(sin_azi[i] *sin_elv[i] <0) )
        
        # Convert output from radians back to degrees
        if not self._in_radians:
            elv_out = np.degrees(elv_out)
            azi_out = np.degrees(azi_out)
            azi_out[ np.where( np.abs(azi_out) <1.e-15 )] = 360
            azi_out = (azi_out+180) %360 -180
        else:
            azi_out = (azi_out +np.pi) %(2*np.pi) -np.pi
        res = np.array([azi_out, elv_out]).transpose()
        
        # Convert the result back to the shape of the input values
        if ndim == 0:
            assert not np.ma.is_masked(res)
            return res[0, 0]
        if ndim == 1:
            return res.reshape(-1)
        elif ndim == 2:
            return res
        raise ValueError(
            "Input values must have shape (N, {dims}) or ({dims},)"
            .format( dims= self.input_dims ) )
    # END OF FUNCITON transform()
    
    def inverted(self):
        return SpherePolarRotate(angle_rad= -self._angle, radians= self._in_radians, **self.__kwargs)
    
    def transform_non_affine(self, values):
        return self.transform(values)
    
    def transform_path_nonaffine(self,path):
        return self.transform_path(path)
    
    def transform_path(self, path):
        nsteps = path._interpolation_steps
        vert1 = path.vertices[0]
        vert2 = path.vertices[-1]
        vertices_new = self._great_circle_points(vert1[0], vert1[1], vert2[0], vert2[1], num_points=nsteps)
        
        vertices_trf = self.transform(vertices_new)
        
        if not path.codes:
            return pth.Path(vertices_trf, path.codes )#, _interpolation_steps= nsteps)
        
        # lift pen whenever vertices neighbor across +180deg.
        codes_buf = []
        for k, v_1, v in zip(
                path.codes[1:],
                              vertices_trf[:-1],
                              vertices_trf[1:]):
            if ( v[0] %(2*np.pi) <np.pi ) != ( v_1[0] %(2*np.pi) <np.pi ) :
                codes_buf.append(path.MOVETO)
            else:
                codes_buf.append(k)
        return pth.Path(vertices_trf, codes_buf, _interpolation_steps= nsteps)


    def _great_circle_points(self, lon1, lat1, lon2, lat2, npts=100):
        """
        Compute points along the great circle route between two points on a sphere.
        Args:
            lon1 (double): Longitude of the starting point (in degrees)
            lat1 (double): Latitude of the starting point (in degrees)
            lon2 (double): Longitude of the ending point (in degrees) 
            lat2 (double): Latitude of the ending point (in degrees)
            npts (int): Number of points to compute along the great circle route.
        
        Returns:
            points (double arr): Numpy (npts x 2) Array of longitude and latitude coordinates of the points along the great circle route.
        """
        
        # Convert degrees to radians
        if self._in_radians:
            lon1, lat1, lon2, lat2 = np.radians([lon1, lat1, lon2, lat2])
    
        # Compute the angular distance *sigma1* between the two points
        delta_lon = lon2 - lon1
        sigma1 = np.arctan2(np.sin(delta_lon) * np.cos(lat2),
                            np.cos(lat1) * np.sin(lat2) - np.sin(lat1) * np.cos(lat2) * np.cos(delta_lon))
    
        # Compute points along the great circle route
        lon_points = np.linspace(lon1, lon2, npts)
        lat_points = np.arcsin(np.sin(lat1) * np.cos(lon_points - lon1) * np.cos(sigma1) +
                               np.cos(lat1) * np.sin(lon_points - lon1) * np.cos(sigma1) +
                               np.sin(sigma1) * np.sin(lon_points - lon1))
        # Convert radians to degrees
        if self._in_radians:
            lon_points = np.degrees(lon_points)
            lat_points = np.degrees(lat_points)
        
        return np.array([lon_points, lat_points]).T

## =======================================

class SphereAziRotate(trf.Transform):
    """
    Performs a rotation of spherical coordinates along the equator (lat=0).
    """
    
    # Operates in degrees internally: _angle in deg
    
    input_dims = 2
    output_dims = 2
    
    is_separable = True
    has_inverse = True
    
    def __init__(self, angle_deg=None, radians=False, angle_rad=None, *args, **kwargs):

        super().__init__(*args, **kwargs)
        self._in_radians = bool(radians)
        if angle_deg:
            self._angle_deg = angle_deg
        elif angle_rad:
            self._angle_deg = np.deg2rad(angle_rad)
        else:
            self._angle_deg = 0

    def transform(self, values):
        '''
        Apply this transformation on the given array of *values*.

        Parameters
        ----------
        values : array-like
            The input values as an array of length :attr:`input_dims` or
            shape (N, :attr:`input_dims`).

        Returns
        -------
        array
            The output values as an array of length :attr:`output_dims` or
            shape (N, :attr:`output_dims`), depending on the input.
        '''
        
        values = np.asanyarray(values)
        ndim = values.ndim
        values = values.reshape((-1, self.input_dims))

        azi = values[0:,0]
        
        # Convert input to degrees
        if self._in_radians:
            azi = np.degrees(azi)
        
        # Calculate rotated azimuth:
        azi_out = azi + self._angle_deg
        azi_out = (azi_out +180) %360 -180

        # Convert output from radians back to degrees
        if self._in_radians:
            azi_out = np.radians(azi_out)

        res = np.array([azi_out, values[0:,1]]).transpose()
        
        # Convert the result back to the shape of the input values
        if ndim == 0:
            assert not np.ma.is_masked(res)
            return res[0, 0]
        if ndim == 1:
            return res.reshape(-1)
        elif ndim == 2:
            return res
        raise ValueError(
            "Input values must have shape (N, {dims}) or ({dims},)"
            .format( dims= self.input_dims ) )
    ## END OF FUNCITON transform()
    
    def inverted(self):
        return SphereAziRotate(angle_deg= -self._angle_deg, radians= self._in_radians, **self.__kwargs)
    
    def transform_non_affine(self, values):
        return self.transform(values)
    
    '''transform_affine(), get_affine(), get_matrix(), transform_path_affine()
        are defined in super() as
        > def get_affine(self):
        > """Get the affine part of this transform."""
        > return IdentityTransform()'''
    
    def get_angle(self):
        if self._in_radians:
            return np.radians(self._angle_deg)
        return self._angle_deg
    
    #TODO define rotate_path(...)
    
    
class CastModulo(trf.Transform):
    """
    Cast all values into some number range by applying %modulo operation.
    """
    
    input_dims = 2
    output_dims = 2
    
    is_separable = True
    has_inverse = False
    
    def __init__(self, range_size, range_min, trf_dim=0, *args, **kwargs):

        super().__init__(*args, **kwargs)
        self._range = range_size
        self._min = range_min
        self._trf_dim = trf_dim

    
    def transform(self, values):
        '''
        Apply this transformation on the given array of *values*.

        Parameters
        ----------
        values : array-like
            The input values as an array of length :attr:`input_dims` or
            shape (N, :attr:`input_dims`).

        Returns
        -------
        array
            The output values as an array of length :attr:`output_dims` or
            shape (N, :attr:`output_dims`), depending on the input.
        '''
        
        values = np.asanyarray(values)
        ndim = values.ndim
        values = values.reshape((-1, self.input_dims))

        # Calculate modulation into target range
        trf_vals = values[0:,self._trf_dim]
        trf_vals = (trf_vals -self._min) %self._range +self._min
        values[0:,self._trf_dim] = trf_vals 
        
        # Convert the result back to the shape of the input values
        if ndim == 0:
            assert not np.ma.is_masked(values)
            return values[0, 0]
        if ndim == 1:
            return values.reshape(-1)
        elif ndim == 2:
            return values
        raise ValueError(
            "Input values must have shape (N, {dims}) or ({dims},)"
            .format(dims=self.input_dims))

    def transform_non_affine(self, values):
        return self.transform(values)

