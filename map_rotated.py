# -*- coding: utf-8 -*-
"""
Routine for plotting celestial maps on Mollweide projection with arbitrary orientation of the map center

Version: v1.7 (2025-10-20)
Created on Tue Apr 16 10:44:22 2024

@author: jgazer (Jonathan Gasser, SwRI)
"""

import numpy as np
import matplotlib.pyplot as plt

from matplotlib.path import Path
from matplotlib.patches import PathPatch
import matplotlib.transforms as trf
# from data_io import load_json, load_data
import sphere_rotate as spr


def plot_rotated_map(lons, lats, data, alf=5, phi=255, th=0,
                     orientation_kw=None, center_meridian=False, **kwargs_plot):
    """
    Plot *data matrix on Mollweide map with ecliptic coordinates bin edges *Lon,*Lat;
    with view rotated by angles *phi (azimuthal), *alf (polar) such that:
        for *th =0, ecl. coordinates lon=phi, lat=alf are at the map center.
        for *th !=0 the view is rotated horizontally again towards positive Lon after applied phi,alf.
    if center_meridian=True is set, *th is interpreted as the Lon of the meridian that goes through the center of the map view.
    
    Args:
        lons (array-like): m+1 array of longitude pixel edges
        lats (array-like): n+1 array of latitude pixel edges
        data (numpy.array): m x n matrix of image data values.
        alf (float): angle of polar rotation in [deg]
        phi (float): angle of init azimuthal rotation in [deg]
        th (float): angle of final azimuthal rotation in [deg]
        orientation_kw (str): keyword to set the map orientation. Overrides any passed arguments for alf,phi,th.
        center_meridian (bool): if True, the 'th' argument is interpreted as the Ecliptic Longitude coordinate that should go through the center of the map.
    
    Returns:
        fig (matplotlib.pyplot Figure): Figure object containing the plotted map
        ax (matplotlib.pyplot Axes): Axes object of the plotted map (Child of *fig*)
    """
    
    # Toggle predefined map orientations such as heliospheric nose, tail, ..
    if not orientation_kw is None:
        phi,alf,th, center_meridian = orientate_map(orientation_kw)

    phi = cast_into_360range(phi, -180)
    alf = cast_into_360range(alf, -180)
    th = cast_into_360range(th, -180)
    
    # with this block, the map center will lie on the zero meridian:
    # "delta_th = atan( tan(phi) * cos(alf) )"  # derivation: see end of file.
    if center_meridian:
        th_mod = (np.abs(phi-th) >90) and not (np.abs(phi-th) >270)
        tanphi_cosalf = np.tan(np.radians( phi-th )) *np.cos(np.radians(alf))
        
        th = np.degrees(np.arctan( tanphi_cosalf ))  # in [-90,90 deg]
        if th_mod: # expand to full 360
            th = th +180
        th = cast_into_360range(th, -180)

    # if Lon range does not cover a 360deg range, add a NAN column to data
    if np.nanmax(lons) -np.nanmin(lons) <360 :
        data = np.hstack(( data, np.nan*np.ones((data.shape[0], 1)) ))
        lons = np.hstack(( lons, np.nanmin(lons) +360 ))

    if alf !=0:
        # add 'alfa' row in nonrotated meshgrid & data
        aq= np.abs(90-np.abs(alf))
        if np.min(lats) < -aq :
            lats, data = insert_coordinate( -aq, lats, data=data)
        if aq < np.max(lats) :
            lats, data = insert_coordinate( aq, lats, data=data)
    
    # add 'phi' column in nonrotated meshgrid & data
    ins = cast_into_360range(phi, np.min(lons))
    lons, data = insert_coordinate( ins, lons, is_lon=True, data=data)
    oppo = cast_into_360range( phi+180, np.min(lons) )
    lons, data = insert_coordinate( oppo, lons, is_lon=True, data=data)
    lons= lons -ins
    
    # Get cut indices == where lon=n*180 but not at beginning or end
    lon_cut_idx = np.where( (lons%180 == 0) * (lons != lons[0]) * (lons != lons[-1]) )[0]

    lons,lats = np.meshgrid(lons, lats)
    # Rotate coordinates
    lons, lats = spr.rotate_sphere( np.array(lons), np.array(lats), alf )
    
    if len(lon_cut_idx) ==0:
        fig,ax,cbar = draw_mollweide_map([lons],[lats],[data], [alf,phi,th], **kwargs_plot)
        return fig,ax
        
    # Cut the map portions
    c= lon_cut_idx[0]
    if len(lon_cut_idx) >1:
        c2= lon_cut_idx[1]
        lons_plt = ( lons.copy()[0:,0:c+1], lons.copy()[0:,c:c2+1], lons.copy()[0:,c2:] )
        lats_plt = ( lats[0:,0:c+1], lats[0:,c:c2+1], lats[0:,c2:] )
        data_plt = ( data[0:,0:c], data[0:,c:c2], data[0:,c2:] )
    else:
        lons_plt = ( lons.copy()[0:,0:c+1], lons.copy()[0:,c:] )
        lats_plt = ( lats[0:,0:c+1], lats[0:,c:] )
        data_plt = ( data[0:,0:c], data[0:,c:] )
        if cast_into_360range(phi, -180) <0:
            lons_plt = lons_plt[ ::-1] #reverse order
            lats_plt = lats_plt[ ::-1]
            data_plt = data_plt[ ::-1]

    def replace_values(arr, old_val, new_val, precision2 =1e-20):
        ''' Replace all values near *old_val* in numpy.Array *arr* by *new_val*
        '''
        arr[ np.where( (arr -old_val )**2 < precision2 )] = new_val

    # replace marginal values mod_360 to prevent 'smearing over' the map
    if cast_into_360range(phi, -180) >0 :  # case phi >0
        replace_values( lons_plt[0], 0, 360)
        replace_values( lons_plt[1], 0, 0)
        if len(lons_plt) >2:
            replace_values( lons_plt[2], 0, 360)
    else: # case phi <=0
        replace_values( lons_plt[0], 360, 0)
        replace_values( lons_plt[1], 0, 360)
        if len(lons_plt) >2:
            replace_values( lons_plt[2], 360, 0)
    
    fig, ax, cbar = draw_mollweide_map(lons_plt, lats_plt, data_plt, [alf,phi,th], **kwargs_plot)
    return fig, ax
#END OF FUNCTION
    
##  =======  plotting   =================

def draw_mollweide_map( lons_plt, lats_plt, data_plt, angles,
                figsize= (12,8), axes=None,
                title=None, xlabel=None,
                ylabel=None, clabel=None,
                lon_ticks= np.linspace(-150,180,12),
                lat_ticks= np.linspace(-75,75,11),
                lon_ticklabel_offset= (3.5, 1.5), #deg
                lat_ticklabel_offset= (-6, -8), #deg
                title_font=50, label_font=30, tick_font=20, font_scale=1,
                title_color='k', xlabel_color='k', ylabel_color='k', clabel_color='k',
                xtick_color='k', ytick_color='k', ctick_color='k',
                face_color='lightgray', grid_color=[.35,.35,.35], grid_on=True, cbar_on=True,
                cmap='turbo', cmap_bad=[0,0,0,0], cmap_under=[0,0,0,1], cmap_over=[1,1,1,1],
                vmin=None, vmax=None, cnorm=False, cbar_kwargs={}, map_alpha=1):
    """
    Draws a map of data in rotated coordinates in Mollweide projection, including colorbar, grid, labels...
    
    Args:
        lons_plt (list): set of i rotated longitude bin-edge (mi+1,ni+1) matrices, one for each data array
        lats_plt (list): set of i rotated latitude bin-edge (mi+1,ni+1) matrices, one for each data array
        data_plt (list): set of i split-up (mi,ni) data matrices
        angles ([int,int,int]): [alf, phi, th], where alf= polar rotation angle, phi= azimuthal rotation angle, th= view-centering angle
    
    Optional Args:
        axes (matplotlib.pyplot.Axes): plot axes to draw into. If given, parent figure will be used as fig.
        lon_ticks, lat_ticks (array): list of Lon & Lat axis tick values in [deg]
        lon_ticklabel_offset, lat_ticklabel_offset: (dx,dy) offsets of Lon & Lat ticks in [deg]
        title_font, label_font, tick_font (float): font sizes
        font_scale (float): scale factor for all fonts (default: 1)
        title_color, xlabel_color, ylabel_color, clabel_color: colors of title text and axes texts
        xtick_color, ytick_color, ctick_color: colors of axes tick labels
        face_color: background color of Mollweide plot area
        grid_color: color of grid lines
        grid_on (bool): if True grid is drawn
        cbar_on (bool): if True add colorbar
        cmap (str|colormap): colormap to be used for plot: mpl.pyplot.cm.Colormap or keyword or other colormap object.
        cmap_bad, cmap_under, cmap_over (list): RGB+alpha values for missing data, under or over the colorbar range
        vmin, vmax (float): colorbar limits (default: data min & max)
        cbar_kwargs (dict): keywords & values for colorbar (see matplotlib.pyplot.colorbar)
        
    Returns:
        fig (matplotlib.pyplot.Figure): Figure object
        ax (matplotlib.pyplot.Axes): Axes object containing the drawn map
        cbar (matplotlib.Colorbar): Colorbar object of the plot
    """

    alf,phi,th = angles
    
    if font_scale !=1 and font_scale >0:
        title_font = font_scale * title_font
        label_font = font_scale * label_font
        tick_font = font_scale * tick_font
    
    if axes:
        fig = axes.get_figure()
        ax = axes
    else:
        fig = plt.figure(figsize= figsize, layout='constrained')
        ax = plt.subplot(111, projection='mollweide')
    ax.set_facecolor(face_color)
    
    # determine colorbar min and max:
    if cnorm:
        vmin = 0
        vmax = 1
    elif vmin is None or vmax is None:
        minVals = [np.nanmin( M ) for M in data_plt]
        maxVals = [np.nanmax( M ) for M in data_plt]
        if vmax is None:
            vmax = np.nanmax( maxVals )
            if np.isnan(vmax):
                vmax = 1
        if vmin is None:
            vmin = np.nanmin( minVals )
            if np.isnan(vmin):
                vmin = 0
            if vmin > 0 and vmax > 10*vmin:
                vmin = 0
    
    # colorbar keywords: location, orientation, fraction, shrink, aspect, ..
    cbar_kwargs_default = {'location':'bottom', 'shrink': 0.8}
    assert isinstance(cbar_kwargs, dict), " 'cbar_kwargs' should be a dictionary."
    cbar_kwargs = {**cbar_kwargs_default, **cbar_kwargs}
    
    # define the colormap usage:
    if cmap is None:
        cmap = plt.cm.turbo.copy()
    elif isinstance(cmap, str):
        cmap = plt.get_cmap(cmap)
    cmap.set_extremes( bad=cmap_bad, under=cmap_under, over=cmap_over) #rgba
    
    # draw the data patches:
    for i in range(len( lons_plt )):
        plt.pcolormesh( -np.radians( lons_plt[i] +th ), np.radians(lats_plt[i]),
                       data_plt[i], cmap=cmap, vmin=vmin, vmax=vmax, alpha=map_alpha)
        plt.pcolormesh( -np.radians( lons_plt[i]-360 +th ), np.radians(lats_plt[i]),
                       data_plt[i], cmap=cmap, vmin=vmin, vmax=vmax, alpha=map_alpha)
    cbar= plt.colorbar( ax.get_children()[0], **cbar_kwargs)
    
    # get default labels
    if np.any([ (s is None) for s in [title, xlabel, ylabel, clabel] ]):
        from data_io import load_json
        lbls = load_json()['default_labels']
        if title is None:
            title = lbls['title']
        if xlabel is None:
            xlabel = lbls['xlabel']
        if ylabel is None:
            ylabel = lbls['ylabel']
        if clabel is None:
            clabel = lbls['clabel']
    
    # make labels
    ax.set_title(title, fontsize= title_font, color=title_color)
    ax.set_xlabel(xlabel, fontsize= label_font, color=xlabel_color)
    ax.set_ylabel(ylabel, fontsize= label_font, color=ylabel_color)
    cbar.set_label( clabel, size= label_font, color=clabel_color )
    ax.set_xticks( np.radians( lon_ticks ))
    ax.set_yticks( np.radians( lat_ticks ))
    ax.tick_params( axis='x', labelsize= tick_font, labelcolor=xtick_color)
    ax.tick_params( axis='y', labelsize= tick_font, labelcolor=ytick_color)
    cbar.ax.tick_params( labelsize= tick_font, labelcolor=clabel_color)
    cbar.ax.xaxis.get_offset_text().set_fontsize(tick_font)

    # transform X axis labels positions
    polar_rot_tr = spr.SpherePolarRotate( angle_deg= alf, radians=True )
    cast_tr = spr.CastModulo(2*np.pi,-np.pi)
    lon_shift_tr = trf.Affine2D.from_values( 1,0,0,1, -np.radians(th),0)
    txt_xoff_tr = trf.Affine2D.from_values( -1,0,0,1, 
                    np.radians(phi + lon_ticklabel_offset[0]),
                    np.radians(lon_ticklabel_offset[1]) )
    xTickTransform = txt_xoff_tr +polar_rot_tr +lon_shift_tr +cast_tr +ax.transData
    for xt in ax.get_xticklabels():
        xt.set_transform( xTickTransform )
    
    # transform Y axis labels positions
    pi_off_tr = trf.Affine2D.from_values( 1,0,0,1, np.radians(phi),0)
    txt_yoff_tr = trf.Affine2D.from_values( 1,0,0,1, *lat_ticklabel_offset)
    ytick_transform = pi_off_tr +polar_rot_tr +lon_shift_tr +cast_tr +ax.transData +txt_yoff_tr
    for yt in ax.get_yticklabels():
        yt.set_transform( ytick_transform )
    
    if grid_on:
        # transformed grid parallels
        for elv in lat_ticks:
            add_grid_patch(ax, elv, False, alf, th=th, color=grid_color)
        # transformed grid meridians
        for az in lon_ticks:
            add_grid_patch(ax, az+phi, True, alf, th=th, color=grid_color)
    
    return fig, ax, cbar
# END OF 'MAIN' FUNCTION

# ===================================================================

def add_grid_patch(axes, coord, is_meridian, alf, th=0, color=[0.4,0.4,0.4]):
    """
    Draws one line of the grid to the map.
    
    Args: 
        axes (matplotlib.pyplot Axes): axes to draw into.
        coord (double): ecl.Lon or ecl.Lat coordinate value of the grid line
        is_meridian (bool): True= meridians, False= parallels
        alf (double): polar rotation angle of the drawn grid
        th (double): horizontal view-centering angle of the drawn grid
        color (str|list): grid color rgb or keyword
    """
    
    if is_meridian:
        lat0 = np.linspace(-90,90,91)
        lon0 = np.linspace(coord,coord,len(lat0))
    else:
        lon0 = np.linspace(0,360,181)
        lat0 = np.linspace(coord,coord,len(lon0))
    lons,lats = spr.rotate_sphere(lon0,lat0, alf)
    lons -= th
    lons = cast_into_360range(lons,-180)
    
    # split gridline at each "+180 --> -180" jump by inserting
    #   (-180,NAN,+180) in Lons, and (lat_end,NAN,lat_end) in Lats.
    jumps = np.where( np.abs(lons[:-1] - lons[1:]) >90 )[0] +1 #indices after jumps
    lon_end = 180*np.round( lons[ jumps ]/180 )
    lat_end = np.array( (lats[jumps] +lats[jumps-1] )/2 )
    
    ins_lons = np.array([-lon_end, [np.nan]*lon_end.size, lon_end ]).flatten()
    ins_lats = np.array([lat_end, [np.nan]*lat_end.size, lat_end ]).flatten()
    lons = np.insert(lons, np.array([jumps,jumps,jumps]).flatten(), ins_lons )
    lats = np.insert(lats, np.array([jumps,jumps,jumps]).flatten(), ins_lats )
    
    gridline_path = Path(np.radians(np.array([lons,lats]).T))
    gridpatch = PathPatch( gridline_path, transform=axes.transData, fill=False, color=color )
    axes.add_patch(gridpatch)
# END OF FUNCTION

def cast_into_360range(angle, lo_lim):
    ''' modulates angle into 360deg range {lo_lim, lo_lim+360 }'''
    if np.array(angle).size ==1 and angle == lo_lim+360:
        return angle
    return ( np.array(angle) -lo_lim )%360 +lo_lim
# END OF FUNCTION


def insert_coordinate(new_coord, lons_or_lats, is_lon=False, data=None):
    """
    Insert a new longitude or latitude value to given Lon or Lat array,
    and duplicate corresponding col or row in data array.
    
    Args:
        new_coord (float): new coordinate value to be inserted
        lons_or_lats (numpy.array): n+1 array of initial longitude or latitude coordinates, sorted in ascending order.
        is_lon (bool): if True, coordinates are treated as longitude. if False, coordinates are treated as latitude.
        data (numpy.array): m x n array of image data values. if None, coordinates are processed independent of data, and data are not processed. 
            This may lead to array size mismatch between coordinates and data.
    
    Returns:
        coords (numpy.array): n+2 array of lon or lat coordinates including the inserted coord value.
        data (numpy.array): array of image data values with pixel row or col duplicated. Array size: (m+1,n) if is_lon=True, (m,n+1) if is_lon=False.
    """
    
    # if new_coord outside range of existing coords: return unchanged
    if np.all( lons_or_lats > new_coord ) or np.all( lons_or_lats < new_coord ):
        return lons_or_lats, data
    
    has_data = not(data is None or len(data) == 0)
    if is_lon:
        lons = lons_or_lats
        if has_data:
            data = np.concatenate(( data[ 0:, np.where(lons < new_coord)[0] ],
                              data[ 0:, np.where(new_coord < lons)[0] -1 ] ), axis=1)
        coords = np.concatenate(( lons[ np.where(lons < new_coord)[0] ], [new_coord],
                              lons[ np.where(new_coord < lons)[0] ] ))
    else:
        lats = lons_or_lats
        if has_data:
            data = np.concatenate(( data[ np.where(lats < new_coord)[0] ,0:],
                                  data[ np.where(new_coord < lats)[0] -1 ,0:] ))
        coords = np.concatenate(( lats[ np.where(lats < new_coord)[0] ], [new_coord],
                              lats[ np.where(new_coord < lats)[0] ] ))
    return coords, data
# END OF FUNCTION

def orientate_map(kw):
    ''' Set predefined map orientation based on keywords.
    Args:
        kw (str): keyword in {ecl, nose, tail, ribbon, ribbon_c} or few aliases
    Returns: 
        phi, alf, th, center_meridian.
    '''
    from data_io import load_json
    jsn = load_json()
    
    kw_dict = jsn['orientation_keywords']
    assert kw in kw_dict, 'orientation_keyword not found: "'+str(kw)+'"'
    
    coord = jsn['Ecl_Coords']
    th = 0
    center_meridian = False
    if kw_dict[kw] == 0: # ecliptic orientation
        alf =0
        phi =0
    elif kw_dict[kw] == 1: # nose-centered
        alf = coord['nose_lat']
        phi = coord['nose_lon']
    elif kw_dict[kw] == 2: # tail-centered
        alf = coord['tail_lat']
        phi = coord['tail_lon']
    elif kw_dict[kw] == 3: # ribbon center as pole of view / ribbon as line
        alf = -coord['ribbon_lat']
        phi = coord['ribbon_lon']
        th = -120
        center_meridian = True
    elif kw_dict[kw] == 4: # ribbon-centered / ribbon as circle
        alf = coord['ribbon_lat']
        phi = coord['ribbon_lon']
    elif kw_dict[kw] == 5: # galactic
        # from : 
        alf = coord['NGP_lat']-90
        phi = coord['NGP_lon']
        th= coord['Gal_center_lon']
        center_meridian = True
    return phi, alf, th, center_meridian
#END OF FUNCTION

def tag_date(fig):
    ''' Put a date tag at the upper right corner of the figure.
    '''
    from datetime import datetime
    date_tag = datetime.now().strftime('%Y-%m-%d %H:%M')
    fig.text(.98, .999, date_tag, ha='right', va='top', fontsize=8)


## ======== Execute this script if called directly: ========

if __name__ == "__main__":
    from data_io import load_data
    lons, lats, data = load_data('ibex-hi-3_noSP-ram-flux_2009-2019.csv')
    data[ np.where(data<=0) ] =np.nan
    import cmcrameri.cm as cmc
    plot_rotated_map(lons,lats,data, orientation_kw='nose', font_scale=0.75, title='')
    plt.plot(0,0,'.k')
    plt.text(0.02,0.02,'Nose', fontsize=14, color='k')
    plt.show()
    input("Press Enter to close the image and continue...")
