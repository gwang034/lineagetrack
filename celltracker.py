# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

#%% IMPORT PACKAGES

import numpy as np
import seaborn as sns
import pandas as pd
import os
import glob
import matplotlib.pyplot as plt
sns.set()


#%% CHANGE VARIABLES

# path to folder containing movies
path = r"C:\Users\grace\OneDrive\Documents\Warmflash Lab\Video Analysis Project\Data"

# number of experimental conditions
nexpcon = 4


#%% LOAD DATA

# find all files ending in .csv in folder
movie_names = glob.glob(os.path.join(path, "*.csv"))


# initialize list to store all data frames
allmovies=list()

nmov=0 

# loop over the list of csv files
for f in movie_names:
      
    # read the csv file
    df = pd.read_csv(f)
    
    allmovies.append(df)
    
    # count number of movies
    nmov= nmov+1 

# compute number of movies in one experimental condition
clus=int(nmov/nexpcon)
    
#%% FIND FINAL CELLS IN EACH MOVIE

# initialize list to store all final experimental data
experiments=list()

    
for g in range(nexpcon):
    
    # focus on experimental condition of interest
    allmovies_con=allmovies[g*clus:((g+1)*clus)]
    
    # initialize data frame to store all data from given experimental condition
    fincell_df=pd.DataFrame()
    
    # create new column to store final cell number
    fincell_df["Final_Cell_Number"]=0 
    
    # create new column to store movie number
    fincell_df["Movie_Number"]=0
    
    exp_mov_num=0    

    
    for f in range(len(allmovies_con)):
        
        # movie of interest
        movie=allmovies_con[f]
        
        # finds final frame
        finframe=max(movie["frame"])
        
        # finds cells in final frame
        fincells=movie[movie["frame"]==finframe]
        
        # remove cells with unknown lineage
        fincells=fincells[fincells["trackId"] != -1]
        
        
        # create data frame tracking lineage of each final cell
        for i in range(len(fincells)):
            
            # add data from last frame for each final cell to data frame
            fincell_df=pd.concat([fincell_df, pd.DataFrame(fincells.iloc[i]).transpose()])
            
            # reset index
            fincell_df=fincell_df.reset_index(drop="True")
            
            # initialize variables
            trackId=fincell_df["trackId"][len(fincell_df)-1].item()
            
            frame=finframe-1
            
            row=len(fincell_df)-1
            
            # trace back the lineage for each final cell
            while frame >= 0:
                
                # if current frame is not first appearance of cell
                if int(fincell_df["parentTrackId"][row]) == 0:
                    
                    # add data to data frame for cell
                    fincell_df=pd.concat([fincell_df, movie[(movie["frame"]==frame) & (movie["trackId"]==trackId)]])
                    
                    fincell_df=fincell_df.reset_index(drop="True")
                    
                    # store final cell number
                    fincell_df["Final_Cell_Number"][row]=i
                    
                    # store movie number
                    fincell_df["Movie_Number"][row]=exp_mov_num
                    
                # if current frame is first appearance of cell
                else:
                    
                    # add data to data frame for cell
                    trackId=fincell_df["parentTrackId"][row].item()
                    
                    fincell_df=pd.concat([fincell_df, movie[(movie["frame"]==frame) & (movie["trackId"]==trackId)]])
    
                    fincell_df=fincell_df.reset_index(drop="True")
                    
                    # store final cell number
                    fincell_df["Final_Cell_Number"][row]=i
                    
                    # store movie number
                    fincell_df["Movie_Number"][row]=exp_mov_num
                    
                frame=frame-1
                
                row=len(fincell_df)-1
    
        exp_mov_num=exp_mov_num+1
      
    # CLEAN DATA
    
    # keep columns of interest
    fincell_df=fincell_df[["frame", "trackId", "lineageId", "parentTrackId", 
                           "Mean_Intensity_0", "Mean_Intensity_1", 
                           "Object_Area_0", "Final_Cell_Number",
                           "Movie_Number"]]
    
    # rename "frame" column "hours"
    fincell_df.rename(columns={"frame":"hours"}, inplace=True)
    
    # convert frames to hours c 
    fincell_df["hours"]=15*fincell_df["hours"]/60
    
    # rename H2B and SOX2 Columns
    fincell_df.rename(columns={"Mean_Intensity_0":"H2B_Intensity", 
                               "Mean_Intensity_1":"SOX2_Intensity"},
                      inplace=True)
    
    # remove nan values
    for j in range(len(fincell_df)):
        
        if pd.isna(fincell_df["Final_Cell_Number"][j])==True:
            
            fincell_df["Final_Cell_Number"][j]=fincell_df["Final_Cell_Number"][j-1]
            
        if pd.isna(fincell_df["Movie_Number"][j])==True:
            
            fincell_df["Movie_Number"][j]=fincell_df["Movie_Number"][j-1]
                
    # add data frame for experimental condition to experiments list
    experiments.append(fincell_df)
    
    
#%% MEAN SOX2 AND H2B INTENSITIES

# initialize empty list to store dataframes
experiments_means=list()

# iterate through experimental conditions
for j in range(nexpcon):
    
    # initialize new dataframe to store means for each exp. con.
    fincellmeans_df=pd.DataFrame()
    
    # iterate through movies in each exp. con.
    for k in range(clus):
                
        for l in range(finframe):
            
            row=0
            
            # mean intensity data for all cells in a frame
            mean_data=experiments[j][(experiments[j]["Movie_Number"]==k) & 
                                 (experiments[j]["hours"]==l/4)].mean()
            
            # keep columns of interest
            mean_data=mean_data[["hours", "H2B_Intensity", "SOX2_Intensity", 
                                 "Movie_Number"]]
            
            # append row to dataframe
            fincellmeans_df=pd.concat([fincellmeans_df, 
                                       pd.DataFrame(mean_data).transpose()])
            
            # reset index
            fincellmeans_df=fincellmeans_df.reset_index(drop="True")

        
    # rename H2B and SOX2 Columns
    fincellmeans_df.rename(columns={"H2B_Intensity":"H2B_Mean", 
                                "SOX2_Intensity":"SOX2_Mean"},
                      inplace=True)
    
    # add dataframe to list
    experiments_means.append(fincellmeans_df)
    
#%% PLOTS: SOX2 INTENSITIES

# create figure with 4 subplots for each experimental condition
fig, axes = plt.subplots(2,2, figsize=(15,10))

# set style
sns.set_context("paper", font_scale=2)

# create line plots for each final cell in each movie
sns.lineplot(x="hours", y="SOX2_Intensity", hue="Final_Cell_Number", 
             data=experiments[0], palette="tab10", ax=axes[0,0], legend=False)

sns.lineplot(x="hours", y="SOX2_Intensity", hue="Final_Cell_Number", 
             data=experiments[1], palette="tab10", ax=axes[0,1], legend=False)

sns.lineplot(x="hours", y="SOX2_Intensity", hue="Final_Cell_Number", 
             data=experiments[2], palette="tab10", ax=axes[1,0], legend=False)

sns.lineplot(x="hours", y="SOX2_Intensity", hue="Final_Cell_Number", 
             data=experiments[3], palette="tab10", ax=axes[1,1], legend=False)

# plot the mean SOX2 expression in all cells at each time point
sns.lineplot(x="hours", y="SOX2_Intensity", data=experiments[0], color="Black", 
             err_style="bars", err_kws={'capsize':3}, ax=axes[0,0])

sns.lineplot(x="hours", y="SOX2_Intensity", data=experiments[1], color="Black", 
             err_style="bars", err_kws={'capsize':3}, ax=axes[0,1])

sns.lineplot(x="hours", y="SOX2_Intensity", data=experiments[2], color="Black", 
             err_style="bars", err_kws={'capsize':3}, ax=axes[1,0])

sns.lineplot(x="hours", y="SOX2_Intensity", data=experiments[3], color="Black", 
             err_style="bars", err_kws={'capsize':3}, ax=axes[1,1])

# set titles for plots
axes[0,0].set_title('mTeSR 0-48', fontweight="bold")

axes[0,1].set_title('BMP 10ng/ml 0-30, Noggin 30-48', fontweight="bold")

axes[1,0].set_title('BMP 50ng/ml 0-30, Noggin 30-48', fontweight="bold")

axes[1,1].set_title('BMP 50ng/ml 0-48', fontweight="bold")

# set same y axis for all
plt.setp(axes, ylim=(550,750))

# fix formatting
fig.tight_layout()


#%% PLOTS: H2B INTENISTIES

# create figure with 4 subplots for each experimental condition
fig, axes = plt.subplots(2,2, figsize=(15,10))

# set style
sns.set_context("paper", font_scale=2)

# create line plots for each final cell in each movie
sns.lineplot(x="hours", y="H2B_Intensity", hue="Final_Cell_Number", 
             data=experiments[0], palette="tab10", ax=axes[0,0], legend=False)

sns.lineplot(x="hours", y="H2B_Intensity", hue="Final_Cell_Number", 
             data=experiments[1], palette="tab10", ax=axes[0,1], legend=False)

sns.lineplot(x="hours", y="H2B_Intensity", hue="Final_Cell_Number", 
             data=experiments[2], palette="tab10", ax=axes[1,0], legend=False)

sns.lineplot(x="hours", y="H2B_Intensity", hue="Final_Cell_Number", 
             data=experiments[3], palette="tab10", ax=axes[1,1], legend=False)

# set titles for plots
axes[0,0].set_title('mTeSR 0-48', fontweight="bold")

axes[0,1].set_title('BMP 10ng/ml 0-30, Noggin 30-48', fontweight="bold")

axes[1,0].set_title('BMP 50ng/ml 0-30, Noggin 30-48', fontweight="bold")

axes[1,1].set_title('BMP 50ng/ml 0-48', fontweight="bold")

# set same y axis for all
plt.setp(axes, ylim=(500,1000))

# fix formatting
fig.tight_layout()


#%% PLOTS: SOX2 EXPRESSION IN A LINEAGE

# Exp. Conditions of interest: 
    # BMP 10ng/ml 0-30, Noggin 30-48 (expcon = 1)
    # BMP 50ng/ml 0-30, Noggin 30-48 (expcon = 2)
    
expcon=1

orig_cells=pd.unique(experiments[expcon]["lineageId"])

for f in range(len(orig_cells)):
    
    # chooses original cell of interest
    orig_cell=orig_cells[f]
    
    # focuses only on cells in the lineage of the original cell
    lineage=experiments[expcon][experiments[expcon]["lineageId"]==orig_cell]
    
    # plot SOX2 expression of cells in lineage
    plt.figure(dpi=500)
    
    sns.set_context("paper", font_scale=1.5)
    
    sns.lineplot(x="hours", y="SOX2_Intensity", data=lineage, hue="trackId",
                 palette="tab10", legend=False)
    
    # plot vertical line at hour 0 and hour 30
    plt.axvline(x = 30, color = 'k', ls='--')
    
    plt.axvline(x = 0, color = 'k', ls='--')
    
    # set title
    plt.suptitle('SOX2 Intensities within Cell Lineage', fontsize=18)
    
    titles=[0, "BMP 10ng/ml 0-30, Noggin 30-48", 
            "BMP 50ng/ml 0-30, Noggin 30-48"]
    
    plt.title(titles[expcon], fontsize=12)

    # set same y axis for all
    plt.ylim(550,700)
    
    # fix formatting
    fig.tight_layout()
               
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
        
        plt.imshow(maxprojchan)
        
        os.chdir(namempfolder)
        
        final_file=(filename[0:-4]+'_maxproj')
        
        if channelnum==0:
            
            tifffile.imwrite(str(final_file+'.tif'), maxprojchan)
            
        else:
            
            with tifffile.TiffWriter(str(final_file+'.tif'), append=True) as tif2write:
                tif2write.save(maxprojchan)
    
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
            
            # create dataframe
            photo_df=pd.DataFrame(centroid)
            
            photo_df=pd.concat([photo_df, pd.DataFrame(dapi_int)], axis=1)
            
            photo_df.rename(columns={"mean_intensity":"dapi"}, inplace=True)
            
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













