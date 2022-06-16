# -*- coding: utf-8 -*-
"""
Created on Wed Jun  1 21:41:03 2022

@author: grace
"""
#%% LOAD IMAGES
def loadimages(full_path):
       
    from skimage import io
    
    # full_path=(img_path +'/'+filename) #+'.tif')
    
    img=io.imread(full_path, fastij=False)
    
    
    return img

#%% RUN SEGMENTATION
def run_seg(full_path, script_path):
    
    # import os
    
    # os.chdir(script_path)
    
    from Run_Cellpose import cellpose_driver
    
    cellpose_driver(full_path)
    #print('Minutes elapsed='+str(datetime.timedelta(seconds = elapsed)))
    
    return

#%% LOAD SEGMENTATION DATA
def seg_load(img_path, filename):
    import os
    import numpy as np
    
    pathtosegmentation=os.path.join(img_path);
    seg_file=(filename[0:-4]+'_seg.npy')
    dat=np.load(os.path.join(pathtosegmentation,seg_file),allow_pickle=True).item();
    masks=dat['masks'];
    img=dat['img']
    
    return masks, img

#%%  CLEAN SEGMENTATION
def clean_seg(masks):
    
    from skimage import segmentation
    from skimage import morphology
    
    # Remove objects touching the border of the image
    masksclean=segmentation.clear_border(masks)

    # Remove small objects
    masksclean=morphology.remove_small_objects(masksclean, min_size=200)
    
    return masksclean
    
    
    
    