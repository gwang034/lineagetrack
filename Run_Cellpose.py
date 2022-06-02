# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 11:45:29 2022

@author: grace
"""
#%%
def cellpose_driver(full_path):
    import numpy as np
    import time, os, sys
    from urllib.parse import urlparse
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    #%matplotlib inline
    mpl.rcParams['figure.dpi'] = 300
    from cellpose import utils, io
    import time
    
    # # I will download images from website
    # urls = ['http://www.cellpose.org/static/images/img02.png',
    #         'http://www.cellpose.org/static/images/img03.png',
    #         'http://www.cellpose.org/static/images/img05.png']
    # files = []
    # for url in urls:
    #     parts = urlparse(url)
    #     filename = os.path.basename(parts.path)
    #     if not os.path.exists(filename):
    #         sys.stderr.write('Downloading: "{}" to {}\n'.format(url, filename))
    #         utils.download_url_to_file(url, filename)
    #     files.append(filename)
    
    # REPLACE FILES WITH YOUR IMAGE PATHS
    # files = [r"C:\Users\grace\OneDrive\Documents\Warmflash Lab\Week 2 Max Proj Practice\MaxProj\P1_W1_1_120200911_20855 PM_maxproj.tif"]
    #print(files)
    # view 1 image
    t = time.time()
    img = io.imread(full_path)
    plt.figure(figsize=(2,2))
    plt.imshow(img)
    plt.axis('off')
    plt.show()
    

    #%%
    # RUN CELLPOSE
    
    from cellpose import models, io
    
    # DEFINE CELLPOSE MODEL
    # model_type='cyto' or model_type='nuclei'
    model = models.Cellpose(gpu=False, model_type='cyto2', )
    
    # define CHANNELS to run segementation on
    # grayscale=0, R=1, G=2, B=3
    # channels = [cytoplasm, nucleus]
    # if NUCLEUS channel does not exist, set the second channel to 0
    # channels = [0,0]
    # IF ALL YOUR IMAGES ARE THE SAME TYPE, you can give a list with 2 elements
    channels = [[0,0]] # IF YOU HAVE GRAYSCALE
    # channels = [2,3] # IF YOU HAVE G=cytoplasm and B=nucleus
    # channels = [2,1] # IF YOU HAVE G=cytoplasm and R=nucleus
    
    # or if you have different types of channels in each image
    # channels = [[2,3], [0,0], [0,0]]
    
    # if diameter is set to None, the size of the cells is estimated on a per image basis
    # you can set the average cell `diameter` in pixels yourself (recommended) 
    # diameter can be a list or a single number for all images
    
    # you can run all in a list e.g.
    # >>> imgs = [io.imread(filename) in for filename in files]
    # >>> masks, flows, styles, diams = model.eval(imgs, diameter=None, channels=channels)
    # >>> io.masks_flows_to_seg(imgs, masks, flows, diams, files, channels)
    # >>> io.save_to_png(imgs, masks, flows, files)
    
    full_path=[full_path]
    
    from skimage import io as skio
    
    # or in a loop
    for chan, filename in zip(channels, full_path):
        img = skio.imread(filename, img_num=1)
        masks, flows, styles, diams = model.eval(img, diameter=None, channels=chan)
    
        # save results so you can load in gui
        io.masks_flows_to_seg(img, masks, flows, diams, filename, chan)
    
        # save results as png
        io.save_to_png(img, masks, flows, filename)
        
        
    elapsed = time.time() - t
        #%%
        # DISPLAY RESULTS
    from cellpose import plot
    
    fig = plt.figure(figsize=(12,5))
    plot.show_segmentation(fig, img, masks, flows[0], channels=chan)
    plt.tight_layout()
    plt.show()
    
    return elapsed
