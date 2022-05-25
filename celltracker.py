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
    fincell_df=fincell_df[["frame", "parentTrackId", "Mean_Intensity_0",
                    "Mean_Intensity_1", "Object_Area_0", "Final_Cell_Number",
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
            mean_data=mean_data[["hours", "H2B_Intensity",
                            "SOX2_Intensity", "Movie_Number"]]
            
            # 
            fincellmeans_df=pd.concat([fincellmeans_df, 
                                       pd.DataFrame(mean_data).transpose()])
        
    # rename H2B and SOX2 Columns
    fincell_df.rename(columns={"H2B_Intensity":"H2B_Mean", 
                               "SOX2_Intensity":"SOX2_Mean"},
                      inplace=True)
    
    experiments_means.append(fincellmeans_df)
    
#%% PLOTS: SOX2 EXPRESSION

# create figure with 4 subplots for each experimental condition
fig, axes = plt.subplots(2,2, figsize=(15,10))

# set style
sns.set_context("paper", font_scale=2)

# create line plots for each final cell in each movie
sns.lineplot(x="hours", y="SOX2_Intensity", hue="Final_Cell_Number", 
             data=experiments[0], ax=axes[0,0], legend=False)

sns.lineplot(x="hours", y="SOX2_Intensity", hue="Final_Cell_Number", 
             data=experiments[1], ax=axes[0,1], legend=False)

sns.lineplot(x="hours", y="SOX2_Intensity", hue="Final_Cell_Number", 
             data=experiments[2], ax=axes[1,0], legend=False)

sns.lineplot(x="hours", y="SOX2_Intensity", hue="Final_Cell_Number", 
             data=experiments[3], ax=axes[1,1], legend=False)

# plot the mean SOX2 expression in all cells at each time point
sns.lineplot(x=experiments[0]["hours"], 
             y=experiments[0][experiments[0]["hours"]==experiments[0]["hours"]]["Mean_Intensity_1"].mean(),
             ax=axes[0,0])


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
















