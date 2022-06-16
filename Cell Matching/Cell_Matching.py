# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 10:40:08 2022

@author: grace
"""

#%% CELL MATCHING: MAX PROJ IMAGES

# For the image to be aligned, a max projection image should be used.
# Run this code if images are not currently max projection images.


#============================** CHANGE VARIABLES **============================

img_loc=r"C:\Users\grace\OneDrive\Documents\Warmflash Lab\Video Analysis Project\Data\photo frames"

Nchannels=4

#==============================================================================

# initialize
import os
from skimage import io
import tifffile

img_name=[]

# make new folder for max proj images
namempfolder = os.path.join(img_loc,"MaxProj")  # Specifies name of folder
    
if not os.path.exists(namempfolder):    #check if directory does not exist

    os.mkdir(namempfolder)              # make new folder for maxproj imgs
        
    
# create list of image names
for images in os.listdir(img_loc):
    
    if (images.endswith(".tif")):       # check if the image ends with .tif
    
        img_name.append(images)         # adds to list of image names


# make max projection images
for f in range(len(img_name)):
    
    filename=img_name[f]                   # choose image of interest
    
    img=io.imread(img_loc+"/"+filename, fastij=False)
        
    NZstacks=int(len(img)/Nchannels)        # compute # of z stacks
    
    maxprojchan = []
    
    for channelnum in range(0,Nchannels):
        
        z=1
        
        idx = channelnum+(z-1)*Nchannels
        
        imaux = img[idx,:,:]
        
        maxprojchan = imaux
                
        for z in range(1,NZstacks):
            
            idx = channelnum+z*Nchannels;
            
            imaux = img[idx,:,:]
            
            maxprojchan = np.maximum(maxprojchan,imaux)
            
        maxprojchan[:,:]=maxprojchan
        
        plt.figure()
        
        plt.grid(False)
        
        plt.imshow(maxprojchan)
        
        os.chdir(namempfolder)
        
        final_file=(filename[0:-4]+'_maxproj')
        
        if channelnum==0:
            
            tifffile.imwrite(str(final_file+'.tif'), maxprojchan)
            
        else:
            
            with tifffile.TiffWriter(str(final_file+'.tif'), append=True) as tif2write:
                tif2write.save(maxprojchan)
    
#%% CELL MATCHING: IMAGE REGISTRATION

# Aligns the max proj image and the last frame of the movie

# align = the image to be aligned; max projection image
# reference = the image to be referenced; the last frame of the movie

# ===========================** CHANGE PATHS **================================

# location of photos to be aligned
align_loc=r"C:\Users\grace\OneDrive\Documents\Warmflash Lab\Video Analysis Project\Data\photo frames\MaxProj"

# location of reference photos
ref_loc=r"C:\Users\grace\OneDrive\Documents\Warmflash Lab\Video Analysis Project\Data\movie frames"

# create a new folder to save final photos
import os

save_loc = os.path.join(align_loc,"Aligned Images")  # Name of folder
    
if not os.path.exists(save_loc):    #check if directory does not exist

    os.mkdir(save_loc)              # make new folder for final imgs
    
#==============================================================================


# create list of reference image names
ref_imgs=[]

for images in os.listdir(ref_loc):
    
    if (images.endswith(".tif")):       # check if the image ends with .tif
    
        ref_imgs.append(images)         # adds to list of image names
        
        
# create list of images to be aligned
align_imgs=[]
for images in os.listdir(align_loc):
    
    if (images.endswith(".tif")):       # check if the image ends with .tif
    
        align_imgs.append(images)         # adds to list of image names


from Image_Registration import run_script

# complete image registration for all images

for f in range(len(align_imgs)):
    
    reference=ref_imgs[f]

    align=align_imgs[f]    

    run_script(align, reference, save_loc, align_loc, ref_loc)


# aligned images will be saved to the folder specified by save_loc

#%% CELL MATCHING: SEGMENT IMAGES AND CREATE IMAGE DATAFRAME

# Note: Must run in an environment with Cellpose installed.

# ====================== ** CHANGE PARAMETERS ** ==============================

Nchannels=4

script_path=r"C:\Users\grace\OneDrive\Documents\Warmflash Lab\Video Analysis Project\code\lineagetrack"

# location of photos to be segmented
img_path=r"C:\Users\grace\OneDrive\Documents\Warmflash Lab\Video Analysis Project\Data\photo frames\MaxProj"

# =============================================================================

from skimage import io
from skimage import measure
import os
import numpy as np
import pandas as pd


# change path to where you have Image_Processing saved
os.chdir(script_path)

# import all functions and packages from Image_Processing
from Cell_Segmenter import loadimages, run_seg, seg_load, clean_seg

# collect all image names
image_name=[]       # initialize array image_name

for images in os.listdir(img_path):
        
    if (images.endswith(".tif")):       # check if the image ends with .tif
        image_name.append(images)

# initialize list to store data from experimental photos
experiment_photos=list()

for f in range(len(image_name)):
    
    # chooses image to focus on
    filename=image_name[f]
    
    # full file path
    full_path=(img_path +'/'+filename)
    
    # loads image
    img = loadimages(full_path)
    
    # complete segmentation if it does not already exist
    if not os.path.exists(img_path +'/'+filename[0:-4]+'_seg.npy'):

        run_seg(full_path, script_path)
    
    # load segmentation data
    [masks, img]=seg_load(img_path, filename)
    
    # clean masks
    masksclean = clean_seg(masks)
    
    maxprojimage=np.zeros((img.shape[0], img.shape[1], Nchannels)) # initialize array
    
    photo_df=pd.DataFrame()
    
    from skimage import io
    
    for chan in range(Nchannels):
        # load each channel of the maxproj image into array
        maxprojimage[:,:, chan]=io.imread(full_path, img_num=chan)
        
        if chan==0:
            centroid = measure.regionprops_table(masksclean, maxprojimage[:,:, chan], 
                                                 properties=['centroid'])
            
            dapi_int = measure.regionprops_table(masksclean, maxprojimage[:,:, chan], 
                                              properties=['mean_intensity'])
            
            area = measure.regionprops_table(masksclean, maxprojimage[:,:, chan], 
                                              properties=['area'])
            
            # create dataframe
            photo_df=pd.DataFrame(centroid)
            
            photo_df=pd.concat([photo_df, pd.DataFrame(dapi_int)], axis=1)
            
            photo_df.rename(columns={"mean_intensity":"dapi"}, inplace=True)
            
            photo_df=pd.concat([photo_df, pd.DataFrame(area)], axis=1)

            
        if chan==1:
            
            h2b_int = measure.regionprops_table(masksclean, maxprojimage[:,:, chan], 
                                              properties=['mean_intensity'])
            
            photo_df=pd.concat([photo_df, pd.DataFrame(h2b_int)], axis=1)
            
            photo_df.rename(columns={"mean_intensity":"h2b"}, inplace=True)
            
        if chan==2:
            
            sox2_int = measure.regionprops_table(masksclean, maxprojimage[:,:, chan], 
                                              properties=['mean_intensity'])
            
            photo_df=pd.concat([photo_df, pd.DataFrame(sox2_int)], axis=1)
            
            photo_df.rename(columns={"mean_intensity":"sox2"}, inplace=True)
            
        if chan==3:
            bra_int = measure.regionprops_table(masksclean, maxprojimage[:,:, chan], 
                                              properties=['mean_intensity'])
            
            photo_df=pd.concat([photo_df, pd.DataFrame(bra_int)], axis=1)
            
            photo_df.rename(columns={"mean_intensity":"bra"}, inplace=True)
            
        
    experiment_photos.append(photo_df)


#%% CELL MATCHING: MATCH CENTROIDS

# import data from MATLAB code

# find all files ending in .csv in folder
match_path=r"C:\Users\grace\OneDrive\Documents\Warmflash Lab\Video Analysis Project\Data\image alignment data\matches"

matches = glob.glob(os.path.join(match_path, "*.csv"))


# initialize list to store all data frames
allmatches=list()

# loop over the list of csv files
for f in matches:
      
    # read the csv file
    df = pd.read_csv(f)
    
    allmatches.append(df)
    
# initialize list to store dataframes containing matched data
matched_data=list()
    
for g in range(len(allmatches)):
    
    # load dataframe containing matches
    match_df=allmatches[g]
    
    # load dataframe containing movie data
    mov_df=fincell_only[g]
    
    # load dataframe containing photo data
    phot_df=experiment_photos[g]
    
    # initialize new dataframe for final match data
    finmatch_df=pd.DataFrame()
    
    # add to dataframe
    for k in range(len(match_df)):
        
        # find the matching key
        match_num=match_df["mov_phot_match"][k]
        
        # remove cells matched to dummy particle
        if match_num==-1:
            
            match_df=match_df.drop([k])
            
            mov_df=mov_df.drop([k])
        
        # for cells with a match, add match data to dataframe
        else:
            finmatch_df=pd.concat([finmatch_df, phot_df.iloc[[match_num-1]].reset_index(drop="True")], ignore_index=True, sort=False)
    
    finmatch_df=pd.concat([match_df, finmatch_df, mov_df], axis=1)
    
    # choose columns of interest
    finmatch_df=finmatch_df[["Final_Cell_Number","mov_phot_match","mov_cent_0", 
                             "mov_cent_1", "centroid-0", "centroid-1", 
                             "dapi", "area", "h2b","sox2", "bra", 
                             "H2B_Intensity", "SOX2_Intensity"]]
    
    # rename columns for clarity
    finmatch_df.rename(columns={"centroid-0":"phot_cent_0", 
                                "centroid-1":"phot_cent_1", "dapi":"phot_dapi",
                                "area":"phot_area", "h2b":"phot_h2b", 
                                "sox2":"phot_sox2", "bra":"phot_bra",
                                "H2B_Intensity": "mov_h2b", 
                                "SOX2_Intensity": "mov_sox2"}, inplace=True)
    
    matched_data.append(finmatch_df)
    
#%% PLOT: CELL MATCHES

import skimage.io
import skimage.util

# load in the movie frames
frame_path = r"C:\Users\grace\OneDrive\Documents\Warmflash Lab\Video Analysis Project\Data\movie frames"

# collect all image names
movframes=[]       # initialize array image_name

for images in os.listdir(frame_path):
        
    if (images.endswith(".tif")):       # check if the image ends with .tif
        
        movframes.append(images)
 

for j in range(len(image_name)):
    
    # place still image and movie frame side by side
    photo=namempfolder+'/'+image_name[j]
    
    photo=io.imread(photo, img_num=1)
    
    movie=frame_path+'/'+movframes[j]
    
    movie=io.imread(movie, img_num=1)

    combine = skimage.util.montage([photo, movie], grid_shape=(1, 2))
       
    # plot connections between centroids
    # initialize plot
    plt.figure(dpi=500)
    
    plt.grid(False)
    
    plt.axis('off')
    
    # select data to plot
    data=matched_data[j]
    
    movx=np.array(data["mov_cent_0"])+1048
    
    movy=np.array(data["mov_cent_1"])
    
    photx=np.array(data["phot_cent_1"])
    
    photy=np.array(data["phot_cent_0"])
    
    for k in range(len(movx)):
        
        x_values=[photx[k], movx[k]]
        
        y_values=[photy[k], movy[k]]
        
        plt.plot(x_values, y_values, linewidth=0.5)
    
    io.imshow(combine)
    
    plt.figure(dpi=500)
    
    plt.grid(False)
    
    plt.axis('off')
    
    plt.scatter(photx, photy)
    
    plt.scatter(movx, movy)

    io.imshow(combine)
















