import os
import sys

sys.path.append('/mnt/sda2/Common_dir/Tools/ToolBox_Diedrichsen/PythonToolbox')

import numpy as np
import nibabel as nib
import surfAnalysisPy as surf
import matplotlib.pyplot as plt

from surfplot import Plot
from nilearn import plotting


def get_surface_data(atlas_dir, hemi_side, neuromap = False):
    """
    params:
        atlas_dir (str): Path to the directory with the DiedrichsenLab surface atlas
    Return: 
        dictionnaries: {hemisphere Name: hemi} for pial, white surface and flat surface
    """
    pial_surf = {}
    white_surf = {}
    flatsurf = {}
    inflated = {}
    
    for i, (hemiName, hemi) in enumerate(hemi_side.items()):
        inflated[hemiName] = os.path.join(atlas_dir,'fs_LR.32k.' + hemi + '.inflated.surf.gii')
        flatsurf[hemiName] = os.path.join(atlas_dir,'fs_LR.32k.' + hemi + '.flat.surf.gii')
        pial_surf[hemiName] = os.path.join(atlas_dir, 'fs_LR.32k.' + hemi + '.pial.surf.gii')
        white_surf[hemiName] = os.path.join(atlas_dir, 'fs_LR.32k.' + hemi + '.white.surf.gii')

    return pial_surf, white_surf, flatsurf, inflated
    


def get_vertices(flatsurf, hemi_side):
    """
    params:

    return:

    """
    vertices = {}
    
    for i, (hemiName, hemi) in enumerate(hemi_side.items()):
        surface = nib.load(flatsurf[hemiName])
        vertices[hemiName] = surface.darrays[0].data
        
        print(f'Maximum Framing for data visualization in {hemiName} hemisphere: ',
        '\n\tLeft: ', np.nanmin(vertices[hemiName][:, 0]), 
        '\n\tRight: ', np.nanmax(vertices[hemiName][:, 0]), 
        '\n\tBottom: ', np.nanmin(vertices[hemiName][:, 1]),
        '\n\tTop: ', np.nanmax(vertices[hemiName][:, 1]),
        '\n')
    
    return vertices



def extract_vertices(vertices, data, framing):
    """
    params:

    return:

    """    
    vertex_in = (vertices[:,0]>=(framing[0])) & \
                (vertices[:,0]<=(framing[1])) & \
                (vertices[:,1]>=(framing[2])) & \
                (vertices[:,1]<=(framing[3]))
    
    array = np.zeros(data.shape)
    vertices_zoom = np.where(vertex_in == True)[0]
    array[vertices_zoom] = 1
    
    return array



def surface_surfplot_3D(array_mask, inflated, hemiName, hemi):
    """
    params:
    
    return:
    
    """
    
    if hemiName == 'left':
        p = Plot(surf_lh = inflated, views = 'lateral', size = (1200, 500))
    elif hemiName == 'right':
        p = Plot(surf_rh = inflated, views = 'lateral', size = (1200, 500))
        
    p.add_layer({hemiName: f'/mnt/sda2/Common_dir/Tools/ToolBox_Diedrichsen/PythonToolbox/surfAnalysisPy/standard_mesh/fs_{hemi}/fs_LR.32k.{hemi}.shape.gii'},
             cmap='binary_r', cbar=False)
    p.add_layer({hemiName: array_mask}, cbar = False, alpha = .5)
        
    return p   



def visualize_zoom(data_path, atlas_dir, overlay, hemisphere = 'both', frame = {'left': None, 'right': None}, **kwargs):
    """
    params:
        hemisphere: 'left', 'right' or 'both'
        frame: Left (x min), Right (x max), Bottom (y min), Top (y max)
        **kwargs: arguments to pass to the surfplot function from Diedrichsen lab (scale, color, etc.) or the surface 3D plot function
    return:
    """    
    if hemisphere == 'both':
        hemi_side = {'left': 'L', 'right': 'R'}
    elif hemisphere == 'left':
        hemi_side = {'left': 'L'}
    elif hemisphere == 'right':
        hemi_side = {'right': 'R'}

    pial_surf, white_surf, flat_surf, inflated = get_surface_data(atlas_dir, hemi_side)

    vertices = get_vertices(flat_surf, hemi_side)
    
    for i, (hemiName, hemi) in enumerate(hemi_side.items()):  
        if frame[hemiName]:
            framing = frame[hemiName]
        else:
            framing = None
            
        if type(data_path) == list:
            data = surf.map.vol_to_surf(data_path, pial_surf[hemiName], white_surf[hemiName])
        elif type(data_path) == str:
            data = surf.map.vol_to_surf([data_path], pial_surf[hemiName], white_surf[hemiName])
        
        if framing:
            surf.plot.plotmap(data, f'fs32k_{hemi}', overlay_type = overlay, 
                              frame = framing, **kwargs)
            plt.show()

            array_mask = extract_vertices(vertices[hemiName], data, framing)
            
            surf.plot.plotmap(array_mask, f'fs32k_{hemi}', overlay_type = 'label', alpha = 0.4)
            plt.show()
            
            plotting.plot_surf_roi(surf_mesh = inflated[hemiName], roi_map = array_mask,
                       bg_map = f'/mnt/sda2/Common_dir/Tools/ToolBox_Diedrichsen/PythonToolbox/surfAnalysisPy/standard_mesh/fs_{hemi}/fs_LR.32k.{hemi}.shape.gii',
                       hemi = hemiName, cmap='Purples', alpha = 1, vmax = 1.4, bg_on_data=True)
            plt.show()
        
         # Make use of surfplot package
         # Currently not usable because some issues when it comes to p.build()
            # p = surface_surfplot_3D(array_mask, inflated[hemiName], hemiName, hemi)
            # fig = p.build()    
            # fig.show(embed_nb = True)
            
        else:
            surf.plot.plotmap(data, f'fs32k_{hemi}', overlay_type = overlay, **kwargs)
            plt.show()  