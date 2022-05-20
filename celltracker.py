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
from matplotlib import pyplot as plt
sns.set()


#%% CHANGE VARIABLES

# path to folder containing movies
path = r"C:\Users\grace\OneDrive\Documents\Warmflash Lab\Video Analysis Project\Data"

# number of experimental conditions
nexpcon = 42
#dsdfasfs

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

    
#%% FIND FINAL CELLS IN EACH MOVIE

# initialize list to store all final experimental data
experiments=list()

for f in range(len(allmovies)):
    
    # movie of interest
    movie=allmovies[f]
    
    # finds final frame
    finframe=max(movie["frame"])
    
    # finds cells in final frame
    fincells=movie[movie["frame"]==finframe]
    
    # remove cells with unknown lineage
    fincells=fincells[fincells["trackId"] != -1]
    
    # initialize matrix containing final cellular data from the movie
    finmovie_matrix=np.zeros([finframe+1,5,len(fincells)])
    
    # create data frame tracking lineage of each final cell
    for i in range(len(fincells)):
        
        # new data frame with final cell as first row
        fincell_df=pd.DataFrame(fincells.iloc[i]).transpose()
        
        # reset index
        fincell_df=fincell_df.reset_index(drop="True")
        
        # initialize variables
        trackId=fincell_df["trackId"].item()
        
        frame=finframe-1
        
        row=0
        
        # trace back the lineage for each final cell
        while frame >= 0:
            
            # if current frame is not first appearance of cell
            if int(fincell_df["parentTrackId"][row]) == 0:
                
                # add data to data frame for cell
                fincell_df=pd.concat([fincell_df, movie[(movie["frame"]==frame) & (movie["trackId"]==trackId)]])
                                
                fincell_df=fincell_df.reset_index(drop="True")


            # if current frame is first appearance of cell
            else:
                
                # add data to data frame for cell
                trackId=fincell_df["parentTrackId"][row].item()
                
                fincell_df=pd.concat([fincell_df, movie[(movie["frame"]==frame) & (movie["trackId"]==trackId)]])

                fincell_df=fincell_df.reset_index(drop="True")

                
                
            frame=frame-1
            
            row=row+1
        
        fincell_df=fincell_df[["frame", "parentTrackId", "Mean_Intensity_0",
                              "Mean_Intensity_1", "Object_Area_0"]]
        
        finmovie_matrix[:,:,i]=np.array(fincell_df)

    experiments.append(finmovie_matrix)
    
    
    
    
#%% PLOTS: SOX2 EXPRESSION

fig, axes = plt.subplots(2,2, dpi=1000)

axes[0].set_title('mTeSR 0-48')

axes[1].set_title('BMP 10ng/ml 0-30, Noggin 30-48')

axes[2].set_title('BMP 50ng/ml 0-30, Noggin 30-48')

axes[3].set_title('BMP 50ng/ml 0-48')


















