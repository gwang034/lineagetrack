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
path = r"C:\Users\grace\OneDrive\Documents\Warmflash Lab\Video Analysis Project\Data\Full Data\CSV Files"

# number of experimental conditions
nexpcon = 4

# number of expected movies (may be different from real number)
nmov = 64

# compute number of movies in one experimental condition
clus=int(nmov/nexpcon)


#%% LOAD DATA

# find all files ending in .csv in folder
movie_names = glob.glob(os.path.join(path, "*.csv"))


# initialize list to store all data frames
allmovies=list()

# loop over the list of csv files
for f in movie_names:
      
    # read the csv file
    df = pd.read_csv(f)
    
    allmovies.append(df)


    
#%% DATA FRAME: LINEAGES FROM FINAL CELLS

# initialize list to store all final experimental data
experiments=list()

# make new folder for dataframes
namefindffolder = os.path.join(path,"Final Dataframes")  # Specifies name of folder
    
if not os.path.exists(namefindffolder): #check if directory does not exist

    os.mkdir(namefindffolder)           # make new folder for final dataframes
        
    
for g in range(nexpcon):
    
    # initialize list to store all movies of the same experimental condition
    allmovies_con=list()
    
    # create an array of all possible movie numbers
    nmov_ar = np.array(range(nmov))
    
    nmov_st=list()
    
    # change all to two digit number strings
    for j in range(len(nmov_ar)):
        
        nmov_st.append(str(nmov_ar[j]).zfill(2))
    
    # focus on experimental condition of interest
    mov_nums=nmov_st[g*clus//2:((g+1)*clus//2)]

    mov_nums=np.append(mov_nums, nmov_st[nmov-(clus//2*(g+1)):nmov-(clus//2*(g))])
    
    # store all movies of one experimental condition in a dataframe
    for h in range(len(mov_nums)):
        
        for i in range(len(movie_names)):
            
            if movie_names[i][-16:-14]==mov_nums[h]:
                
                allmovies_con.append(allmovies[i])
                
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
                           "Movie_Number", "Object_Center_0", 
                           "Object_Center_1", "Object_Area_0"]]
    
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
            
            fincell_df.loc[j].at["Final_Cell_Number"]=fincell_df["Final_Cell_Number"][j-1]
            
        if pd.isna(fincell_df["Movie_Number"][j])==True:
            
            fincell_df.loc[j].at["Movie_Number"]=fincell_df["Movie_Number"][j-1]
    
    print(g)
    
    # save dataframe
    fincell_df.to_csv(namefindffolder+'/'+'condition_'+str(g)+'.csv')
    
    # add data frame for experimental condition to experiments list
    experiments.append(fincell_df)
    
#%% DATA FRAME: FINAL CELLS ONLY

# initialize list to contain all dataframes
fincell_only=list()

for f in range(len(experiments)):
    
    exp=experiments[f]
    
    # only keep cells in the final frame
    new_df=exp[exp["hours"]==50.50]
    
    new_df=new_df.reset_index(drop="True")
    
    fincell_only.append(new_df)
    
    
    
    
#%% DATA FRAME: MEAN SOX2 AND H2B INTENSITIES

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
               















