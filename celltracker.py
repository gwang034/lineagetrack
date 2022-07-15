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
        
# create an array of all possible movie numbers
nmov_ar = np.array(range(nmov))

nmov_st=list()

# change all to two digit number strings
for j in range(len(nmov_ar)):
    
    nmov_st.append(str(nmov_ar[j]).zfill(2))

    
for g in range(nexpcon):
    
    if not os.path.exists(namefindffolder+'/'+'condition_'+str(g)+'.csv'):
        # initialize list to store all movies of the same experimental condition
        allmovies_con=list()
        
        # focus on experimental condition of interest
        mov_nums=nmov_st[g*clus//2:((g+1)*clus//2)]
    
        mov_nums=np.append(mov_nums, nmov_st[nmov-(clus//2*(g+1)):nmov-(clus//2*(g))])
        
        # store all movies of one experimental condition in a dataframe
        for h in range(len(mov_nums)):
            
            for i in range(len(movie_names)):
                
                if movie_names[i][-16:-14]==mov_nums[h]:
                    
                    # add column storing movie number
                    allmovies[i]["Movie_Number"]=mov_nums[h]
                    
                    allmovies_con.append(allmovies[i])
             
        # initialize data frame to store all data from given experimental condition
        fincell_df=pd.DataFrame()
        
        # create new column to store final cell number
        fincell_df["Final_Cell_Number"]=0 
        
        
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
                trackId=fincell_df["trackId"][len(fincell_df)-1]
                
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
                        fincell_df.loc[row, "Final_Cell_Number"]=i
                        
                    # if current frame is first appearance of cell
                    else:
                        
                        # add data to data frame for cell
                        trackId=fincell_df["parentTrackId"][row]
                        
                        fincell_df=pd.concat([fincell_df, movie[(movie["frame"]==frame) & (movie["trackId"]==trackId)]])
        
                        fincell_df=fincell_df.reset_index(drop="True")
                        
                        # store final cell number
                        fincell_df.loc[row, "Final_Cell_Number"]=i
                        
                        
                    frame=frame-1
                    
                    row=len(fincell_df)-1
                    
                    # print(str(g)+'_'+str(row))
        
          
        # CLEAN DATA
        
        # keep columns of interest
        fincell_df=fincell_df[["frame", "trackId", "lineageId", "parentTrackId", 
                               "Mean_Intensity_0", "Mean_Intensity_1", 
                               "Mean_Intensity_2", "Object_Area_0", 
                               "Final_Cell_Number", "Movie_Number", 
                               "Object_Center_0", "Object_Center_1", 
                               "Object_Area_0"]]
        
        # rename "frame" column "hours"
        fincell_df.rename(columns={"frame":"hours"}, inplace=True)
        
        # convert frames to hours c 
        fincell_df["hours"]=15*fincell_df["hours"]/60
        
        # rename H2B and SOX2 Columns
        fincell_df.rename(columns={"Mean_Intensity_0":"H2B_Intensity", 
                                   "Mean_Intensity_1":"SOX2_Intensity",
                                   "Mean_Intensity_2":"BRA_Intensity"},
                          inplace=True)
        
        # remove nan values
        for j in range(len(fincell_df)):
            
            if pd.isna(fincell_df["Final_Cell_Number"][j])==True:
                
                fincell_df.at[j, "Final_Cell_Number"]=fincell_df["Final_Cell_Number"][j-1]        
            
                
        # save dataframe
        fincell_df.to_csv(namefindffolder+'/'+'condition_'+str(g)+'.csv')
        
        print(g)
        
        # add data frame for experimental condition to experiments list
        experiments.append(fincell_df)

    else:
        
        # find all files ending in .csv in folder
        df_names = glob.glob(os.path.join(namefindffolder, "*.csv"))
        
        # read the csv file
        fincell_df = pd.read_csv(df_names[g])
            
        experiments.append(fincell_df)
    
#%% DATA FRAME: FINAL CELLS ONLY

# initialize list to contain all dataframes
fincell_only=list()

for f in range(len(experiments)):
    
    exp=experiments[f]
    
    # only keep cells in the final frame
    new_df=exp[exp["hours"]==50.75]
    
    new_df=new_df.reset_index(drop="True")
    
    fincell_only.append(new_df)
    
    
#%% DATA FRAME: STILL IMAGE DATA

# initialize list to contain all dataframes
still_data=list()

for f in range(len(experiments)):
    
    exp=experiments[f]
    
    # only keep cells in the final frame
    still_df=exp[exp["hours"]==51.00]
        
    still_df.rename(columns={"Unnamed: 0":"Original_Index"}, inplace=True)
    
    still_df.reset_index(inplace=True)
    
    still_data.append(still_df)    
    
    
#%% DATA FRAME: CELL FATES PT. 2

# initialize DataFrames
bra_pos=pd.DataFrame()

cdx2_pos=pd.DataFrame()

sox2_pos=pd.DataFrame()


for f in range(len(experiments)):
    
    exp=experiments[f]
    
    # add column storing experimental condition
    exp["exp_con"]=f
    
    # remove cells with SOX2 < 550
    exp=exp[exp["SOX2_Intensity"]>550]
    
    # find movie numbers in condition
    mov_in_con=np.array(pd.unique(exp["Movie_Number"]))
    
    for g in range(len(mov_in_con)):
        
        # focus on specific movie number
        mov=exp[exp["Movie_Number"]==mov_in_con[g]]
        
        # Identify lineages for each cell fate
        bra_pos_cells=np.array(pd.unique(mov[(mov["BRA_Intensity"]>850) & (mov["hours"]==51)]["Final_Cell_Number"]))
    
        cdx2_pos_cells=np.array(pd.unique(mov[(mov["BRA_Intensity"]<700) & (mov["SOX2_Intensity"]<570) & (mov["hours"]==51)]["Final_Cell_Number"]))
    
        sox2_pos_cells=np.array(pd.unique(mov[(mov["BRA_Intensity"]<800) & (mov["SOX2_Intensity"]>600) & (mov["hours"]==51)]["Final_Cell_Number"]))
        
        # store lineages
        for h in range(len(bra_pos_cells)):
            # store lineage of BRA positive cells
            bra_pos=pd.concat([bra_pos, mov[mov["Final_Cell_Number"]==bra_pos_cells[h]]])
            
        for h in range(len(cdx2_pos_cells)):
            # store lineage of CDX2 positive cells
            cdx2_pos=pd.concat([cdx2_pos, mov[mov["Final_Cell_Number"]==cdx2_pos_cells[h]]])
            
        for h in range(len(sox2_pos_cells)):
            # store lineage of SOX2 positive cells
            sox2_pos=pd.concat([sox2_pos, mov[mov["Final_Cell_Number"]==sox2_pos_cells[h]]])
        
bra_pos["fate"]="BRA"

cdx2_pos["fate"]="CDX2"

sox2_pos["fate"]="SOX2"

bra_pos.reset_index(drop=True, inplace=True)

cdx2_pos.reset_index(drop=True, inplace=True)

sox2_pos.reset_index(drop=True, inplace=True)

all_fates=pd.concat([bra_pos, sox2_pos, cdx2_pos])

all_fates.reset_index(drop=True, inplace=True)

#%% FULL DATA PLOTS: SOX2 INTENSITIES

plt.figure(dpi=500)

# set style
sns.set_context("paper")

# plot the mean SOX2 expression in all cells at each time point
sns.lineplot(x="hours", y="SOX2_Intensity", data=experiments[0], color="Blue", 
             label="mTeSR 0-48")

sns.lineplot(x="hours", y="SOX2_Intensity", data=experiments[1], color="Green",
             label="BMP 10ng/ml 0-30, Noggin 30-48")

sns.lineplot(x="hours", y="SOX2_Intensity", data=experiments[2], color="Red",
             label="BMP 50ng/ml 0-30, Noggin 30-48")

sns.lineplot(x="hours", y="SOX2_Intensity", data=experiments[3], color="Orange",
             label="BMP 50ng/ml 0-48")

plt.title("Mean SOX2 Expression")

plt.ylim(550,690)

plt.legend(loc="lower left")

#%% FULL DATA PLOTS: H2B INTENSITIES

plt.figure(dpi=500)

# set style
sns.set_context("paper")

# plot the mean SOX2 expression in all cells at each time point
sns.lineplot(x="hours", y="H2B_Intensity", data=experiments[0], color="Blue", 
             label="mTeSR 0-48")

sns.lineplot(x="hours", y="H2B_Intensity", data=experiments[1], color="Green",
             label="BMP 10ng/ml 0-30, Noggin 30-48")

sns.lineplot(x="hours", y="H2B_Intensity", data=experiments[2], color="Red",
             label="BMP 50ng/ml 0-30, Noggin 30-48")

sns.lineplot(x="hours", y="H2B_Intensity", data=experiments[3], color="Orange",
             label="BMP 50ng/ml 0-48")

plt.title("Mean H2B Expression")



#%% FULL DATA PLOTS: SCATTER PLOT - FINAL CELL SOX2 AND BRA

plt.figure(dpi=500)

# set style
sns.set_context("paper")

sns.scatterplot(x="SOX2_Intensity", y="BRA_Intensity", data=still_data[0],
                color="Blue", label="mTeSR 0-48")

sns.scatterplot(x="SOX2_Intensity", y="BRA_Intensity", data=still_data[1],
                color="Green", label="BMP 10ng/ml 0-30, Noggin 30-48")

sns.scatterplot(x="SOX2_Intensity", y="BRA_Intensity", data=still_data[2],
                color="Red", label="BMP 50ng/ml 0-30, Noggin 30-48")

sns.scatterplot(x="SOX2_Intensity", y="BRA_Intensity", data=still_data[3],
                color="Orange", label="BMP 50ng/ml 0-48")

plt.title("SOX2 and BRA Expression in Final Frame Cells")

plt.ylim(550,1300)

# Create individual subplots

# set style
sns.set_context("paper", font_scale=2)


fig, axes = plt.subplots(2,2, figsize=(15,10))

sns.scatterplot(x="SOX2_Intensity", y="BRA_Intensity", data=still_data[0],
                color="Blue", legend = False, ax=axes[0,0])

sns.scatterplot(x="SOX2_Intensity", y="BRA_Intensity", data=still_data[1],
                color="Green", legend = False, ax=axes[0,1])

sns.scatterplot(x="SOX2_Intensity", y="BRA_Intensity", data=still_data[2],
                color="Red", legend = False, ax=axes[1,0])

sns.scatterplot(x="SOX2_Intensity", y="BRA_Intensity", data=still_data[3],
                color="Orange", legend = False, ax=axes[1,1])

# set titles for plots
axes[0,0].set_title('mTeSR 0-48', fontweight="bold")

axes[0,1].set_title('BMP 10ng/ml 0-30, Noggin 30-48', fontweight="bold")

axes[1,0].set_title('BMP 50ng/ml 0-30, Noggin 30-48', fontweight="bold")

axes[1,1].set_title('BMP 50ng/ml 0-48', fontweight="bold")

# fix formatting
fig.tight_layout()


#%% FULL DATA PLOTS: SOX2 WITHIN A LINEAGE COLOR CODED WITH BRA

# Plots SOX2 expression throughout a lineage, with final cells color coded
# based on their BRA expression

for expcon in range(0, nexpcon):

    # find all movies in experimental condition
    mov_in_con=pd.unique(experiments[expcon]["Movie_Number"])
        
    for j in range(len(mov_in_con)):
        
        # select movie number
        movie=experiments[expcon][experiments[expcon]["Movie_Number"]==mov_in_con[j]]
        
        # find original cells
        orig_cells=pd.unique(movie["lineageId"])
        
        for f in range(len(orig_cells)):
            
            # chooses original cell of interest
            orig_cell=orig_cells[f]
            
            # focuses only on cells in the lineage of the original cell
            lineage=movie[movie["lineageId"]==orig_cell]
            
            # saves trackId of final cells
            fin_ids=np.array(pd.unique(lineage[lineage["hours"]==51]["trackId"]))
            
            fin_lin=pd.DataFrame()
            
            # saves final cell data from previous frames
            for k in range(len(fin_ids)):
                add=lineage[lineage["trackId"]==fin_ids[k]]
                
                fin_lin=pd.concat([fin_lin, add])
                        
            fin_lin.reset_index(inplace=True)
            
            # updates BRA expression for final cells for seaborn "hue" separation
            for l in range(len(fin_lin)):
                
                if fin_lin.at[l,"BRA_Intensity"]==0:
                    
                    fin_lin.at[l, "BRA_Intensity"]=fin_lin["BRA_Intensity"][l-1]
                
            # plot SOX2 expression of cells in lineage
            plt.figure(dpi=500)
            
            sns.set_context("paper", font_scale=1.5)
            
            sns.lineplot(x="hours", y="SOX2_Intensity", data=lineage, hue="trackId",
                         palette=sns.dark_palette("Blue", n_colors=len(pd.unique(lineage["trackId"])), reverse=True), legend=False)

            
            sns.lineplot(x="hours", y="SOX2_Intensity", data=fin_lin,
                            hue="BRA_Intensity", palette=sns.dark_palette("Red", n_colors=len(pd.unique(fin_lin["Final_Cell_Number"]))))
            
            # plot vertical line at hour 0 and hour 30
            plt.axvline(x = 30, color = 'k', ls='--')
            
            plt.axvline(x = 0, color = 'k', ls='--')
            
            # set title
            plt.suptitle('SOX2 Intensities within Cell Lineage', fontsize=18)
            
            titles=['mTeSR 0-48', "BMP 10ng/ml 0-30, Noggin 30-48", 
                    "BMP 50ng/ml 0-30, Noggin 30-48", 'BMP 50ng/ml 0-48']
            
            plt.title(titles[expcon], fontsize=12)
        
            # set same y axis for all
            plt.ylim(550,750)
            
            # fix formatting
            plt.legend(title="BRA Intensity", bbox_to_anchor =(1.02, 1), loc="upper left")
            
            plt.tight_layout()
            
            plt.savefig(r"C:\Users\grace\OneDrive\Documents\Warmflash Lab\Video Analysis Project\code\figures\SOX2 Lineage with BRA Sorting\Condition "+str(expcon)+'/sox2_c'+str(expcon)+'_m'+str(j)+'_'+str(f)+'.png')
    

#%% FULL DATA PLOTS: MAX AND MIN BRA AND SOX2 LINEAGES

# cycle through experimental conditions
for j in range(len(experiments)):
    
    # sort from greatest to least BRA Intensity
    last_cells=still_data[j].sort_values(by=["BRA_Intensity"], inplace=True)
    
    exp=experiments[j]
    
    last_cells.reset_index(inplace="True")
    
    # initialize variables
    i=0
    
    row_last=0; row_top=len(last_cells)-1
    
    last_lineage=list(); top_lineage=list()
    
    last_num=[]; top_num=[]
    
    last_sox_num=[]; top_sox_num=[]
    
    while i<5:
        
        last=last_cells["Original_Index"][row_last]
        
        sox=last_cells["SOX2_Intensity"][row_last]
        
        # find index of cell in larger dataframe
        idx=exp.iloc[last]
        
        # find data for that cell
        least_cell=pd.DataFrame(idx).transpose()
        
        # focuses only on cells in the lineage of the given cell
        lineage=exp[(exp["lineageId"]==least_cell["lineageId"].item()) & (exp["Movie_Number"]==least_cell["Movie_Number"].item()) & (exp["Final_Cell_Number"]==least_cell["Final_Cell_Number"].item())]
        
        if min(lineage["hours"])==0:
            
            last_lineage.append(lineage)
            
            i=i+1
            
            row_last=row_last+1
            
            last_num=np.append(last_num, last_cells["BRA_Intensity"][row_last].item())
        
            last_sox_num=np.append(last_sox_num, sox)
        
        else:
            
            row_last=row_last+1
            
    
    i=0
    
    while i<5:
        
        top=last_cells["Original_Index"][row_top]
        
        sox=last_cells["SOX2_Intensity"][row_top]
        
        # find index of cell in larger dataframe
        idx=exp.iloc[top]
        
        # find data for that cell
        top_cell=pd.DataFrame(idx).transpose()
        
        # focuses only on cells in the lineage of the given cell
        lineage=exp[(exp["lineageId"]==top_cell["lineageId"].item()) & (exp["Movie_Number"]==top_cell["Movie_Number"].item()) & (exp["Final_Cell_Number"]==top_cell["Final_Cell_Number"].item())]
        
        if min(lineage["hours"])==0:
            
            top_lineage.append(lineage)
            
            i=i+1
            
            row_top=row_top-1
            
            top_num=np.append(top_num, last_cells["BRA_Intensity"][row_top].item())

            top_sox_num=np.append(top_sox_num, sox)

        else:
            
            row_top=row_top-1
            
    
    
    for k in range(0,5):
        
        # plot SOX2 expression of cells in lineage
        fig, ax = plt.subplots(dpi=500)
                
        sns.set_context("paper", font_scale=1.5)
        
        sns.lineplot(x="hours", y="SOX2_Intensity", data=last_lineage[k], hue="trackId",
                     palette=sns.dark_palette("Blue", n_colors=len(pd.unique(last_lineage[k]["trackId"])), reverse=True), legend=False)
        
        # plot vertical line at hour 0 and hour 30
        plt.axvline(x = 30, color = 'k', ls='--')
        
        plt.axvline(x = 0, color = 'k', ls='--')
        
        # set title
        plt.suptitle('SOX2 Intensities within Lineage of Cells with Least and Greatest Final BRA Intensities', fontsize=14)
        
        titles=['mTeSR 0-48', "BMP 10ng/ml 0-30, Noggin 30-48", 
                "BMP 50ng/ml 0-30, Noggin 30-48",
                'BMP 50ng/ml 0-48']
        
        plt.title(titles[j], fontsize=12)
            
        # set same y axis for all
        # plt.ylim(550,750)
        
        # fix formatting
        plt.tight_layout()
        
        
        
        # Plot lineage for cells with greatest BRA expression
        
        
        ## UNCOMMENT TO CREATE SEPARATE PLOTS FOR LEAST & GREATEST
        
        # # plot SOX2 expression of cells in lineage
        # plt.figure(dpi=500)
        
        # sns.set_context("paper", font_scale=1.5)
        
        # # plot vertical line at hour 0 and hour 30
        # plt.axvline(x = 30, color = 'k', ls='--')
        
        # plt.axvline(x = 0, color = 'k', ls='--')
        
        # # set title
        # plt.suptitle('SOX2 Intensities within Lineage of Cells with Greatest Final BRA Intensities', fontsize=14)
        
        # titles=['mTeSR 0-48', "BMP 10ng/ml 0-30, Noggin 30-48", 
        #         "BMP 50ng/ml 0-30, Noggin 30-48",
        #         'BMP 50ng/ml 0-48']
        
        # plt.title(titles[j], fontsize=12)
        
        ## END OF UNCOMMENT SECTION
        
        sns.lineplot(x="hours", y="SOX2_Intensity", data=top_lineage[k], hue="trackId",
                     palette=sns.dark_palette("Red", n_colors=len(pd.unique(top_lineage[k]["trackId"])), reverse=True), legend=False)
        
        plt.plot([51], [top_sox_num[k]], 'lightcoral', marker='o', markersize=5)
        
        plt.plot([51], [last_sox_num[k]], 'dodgerblue', marker='o', markersize=5)
        
        # create legend (COMMENT OUT IF PLOTTING SEPARATELY)
        from matplotlib.lines import Line2D
        
        custom_lines = [Line2D([0], [0], color="Red", lw=4),
                        Line2D([0], [0], color="Blue", lw=4)]
        
        
        ax.legend(custom_lines, ["BRA="+str(top_num[k]),
                                 "BRA="+str(last_num[k])])
        
    
        # set same y axis for all
        # plt.ylim(550,750)
        
        # fix formatting
        plt.tight_layout()
        
        plt.savefig(r"C:\Users\grace\OneDrive\Documents\Warmflash Lab\Video Analysis Project\code\figures\SOX2 Min & Max BRA\Condition "+str(j)+'/bra_c'+str(j)+'_'+str(k)+'.png')
        
#%% FULL DATA PLOTS: TRAJECTORIES GREATER THAN 850

for j in range(len(experiments)):
    
    exp=experiments[j]
    
    pos=pd.DataFrame()
    
    # save the data for BRA positive cells from the last frame
    for i in range(len(still_data[j])):
        
        if still_data[j]["BRA_Intensity"][i]>850:
            
            pos=pd.concat([pos, pd.DataFrame(still_data[j].iloc[i]).transpose()])
            
            pos=pos.reset_index(drop=True)
    
    for k in range(len(pos)):
        
        cell=int(pos.loc[k, "Original_Index"])
        
        # find index of cell in larger dataframe
        idx=exp.iloc[cell]
        
        # find data for that cell
        cell_data=pd.DataFrame(idx).transpose()
        
        # focuses only on cells in the lineage of the given cell
        lineage=exp[(exp["lineageId"]==cell_data["lineageId"].item()) & (exp["Movie_Number"]==cell_data["Movie_Number"].item()) & (exp["Final_Cell_Number"]==cell_data["Final_Cell_Number"].item())]
        
        plt.figure(dpi=500)
        
        sns.lineplot(x="hours", y="SOX2_Intensity", data=lineage, hue="trackId",
                     palette=sns.dark_palette("Red", n_colors=len(pd.unique(lineage["trackId"])), reverse=True), legend=False)
        
        # plot vertical line at hour 0 and hour 30
        plt.axvline(x = 30, color = 'k', ls='--')
        
        plt.axvline(x = 0, color = 'k', ls='--')
        
        
        plt.suptitle("SOX2 Intensity of BRA Positive Cell", fontsize=14)
        
        titles=['mTeSR 0-48', "BMP 10ng/ml 0-30, Noggin 30-48", 
                "BMP 50ng/ml 0-30, Noggin 30-48",
                'BMP 50ng/ml 0-48']
        
        plt.title(titles[j]+'; BRA='+str(pos.loc[k, "BRA_Intensity"]), fontsize=12)
        
        plt.savefig(r"C:\Users\grace\OneDrive\Documents\Warmflash Lab\Video Analysis Project\code\figures\BRA Positive\Condition "+str(j)+'/bra_c'+str(j)+'_'+str(k)+'.png')

#%% FULL DATA PLOT: EACH CELL FATES IN EACH EXPERIMENTAL CONDITION

# run DATA FRAME: CELL FATES first
    
for f in range(len(experiments)):
    
    bra_pos_plot=bra_pos[bra_pos["exp_con"]==f]
    
    # Plot BRA Positive Cells first
    if bra_pos_plot.empty:
        
        continue
    
    else:
        plt.figure(dpi=500)
        
        y_ops=[(550,775), (525,750), (525, 750), (550, 775)]
        
        plt.ylim(y_ops[f])
        
        
        hue=bra_pos_plot[["trackId", "Movie_Number"]].apply(
            lambda row: f"{row.trackId}, {row.Movie_Number}", axis=1)
        
        sns.lineplot(x="hours", y="SOX2_Intensity", data=bra_pos_plot, hue=hue,
                     palette=sns.dark_palette("Red", n_colors=len(pd.unique(hue)), reverse=True), legend=False)
        
        sns.lineplot(x="hours", y="SOX2_Intensity", data=bra_pos_plot, color="paleturquoise",
                     err_style="bars", err_kws={'capsize':3, 'ecolor':"deepskyblue"}, legend=False)
        
        
        # plot vertical line at hour 0 and hour 30
        plt.axvline(x = 30, color = 'k', ls='--')
        
        plt.axvline(x = 0, color = 'k', ls='--')
        
        
        plt.suptitle("SOX2 Intensity of BRA Positive Cells", fontsize=14)
        
        titles=['mTeSR 0-48', "BMP 10ng/ml 0-30, Noggin 30-48", 
                "BMP 50ng/ml 0-30, Noggin 30-48",
                'BMP 50ng/ml 0-48']
        
        plt.title(titles[f], fontsize=12)
        
for f in range(len(experiments)):
    
    cdx2_pos_plot=cdx2_pos[cdx2_pos["exp_con"]==f]
    
    # Plot CDX2 Positive Cells
    if cdx2_pos_plot.empty:
        
        continue
    
    else:
        plt.figure(dpi=500)
        
        y_ops=[(550,775), (525,750), (525, 750), (550, 775)]
        
        plt.ylim(y_ops[f])
        
        hue=cdx2_pos_plot[["trackId", "Movie_Number"]].apply(
            lambda row: f"{row.trackId}, {row.Movie_Number}", axis=1)
        
        sns.lineplot(x="hours", y="SOX2_Intensity", data=cdx2_pos_plot, hue=hue,
                     palette=sns.dark_palette("Red", n_colors=len(pd.unique(hue)), reverse=True), legend=False)
        
        sns.lineplot(x="hours", y="SOX2_Intensity", data=cdx2_pos_plot, color="paleturquoise",
                     err_style="bars", err_kws={'capsize':3, 'ecolor':"deepskyblue"}, legend=False)
        
        
        # plot vertical line at hour 0 and hour 30
        plt.axvline(x = 30, color = 'k', ls='--')
        
        plt.axvline(x = 0, color = 'k', ls='--')
        
        
        plt.suptitle("SOX2 Intensity of CDX2 Positive Cells", fontsize=14)
        
        titles=['mTeSR 0-48', "BMP 10ng/ml 0-30, Noggin 30-48", 
                "BMP 50ng/ml 0-30, Noggin 30-48",
                'BMP 50ng/ml 0-48']
        
        plt.title(titles[f], fontsize=12)
        
        
for f in range(len(experiments)):
    
    sox2_pos_plot=sox2_pos[sox2_pos["exp_con"]==f]
    
    # Plot SOX2 Positive Cells last
    if sox2_pos_plot.empty:
        
        continue
    
    else:
        plt.figure(dpi=500)
        
        y_ops=[(550,775), (525,750), (525, 750), (550, 775)]
        
        plt.ylim(y_ops[f])
        
        hue=sox2_pos_plot[["trackId", "Movie_Number"]].apply(
            lambda row: f"{row.trackId}, {row.Movie_Number}", axis=1)
        
        sns.lineplot(x="hours", y="SOX2_Intensity", data=sox2_pos_plot, hue=hue,
                     palette=sns.dark_palette("Red", n_colors=len(pd.unique(hue)), reverse=True), legend=False)
        
        sns.lineplot(x="hours", y="SOX2_Intensity", data=sox2_pos_plot, color="paleturquoise",
                     err_style="bars", err_kws={'capsize':3, 'ecolor':"deepskyblue"}, legend=False)
        
        
        # plot vertical line at hour 0 and hour 30
        plt.axvline(x = 30, color = 'k', ls='--')
        
        plt.axvline(x = 0, color = 'k', ls='--')
        
        
        plt.suptitle("SOX2 Intensity of SOX2 Positive Cells", fontsize=14)
        
        titles=['mTeSR 0-48', "BMP 10ng/ml 0-30, Noggin 30-48", 
                "BMP 50ng/ml 0-30, Noggin 30-48",
                'BMP 50ng/ml 0-48']
        
        plt.title(titles[f], fontsize=12)
    
#%% FULL DATA PLOT: ALL CELL FATES IN EACH EXPERIMENTAL CONDITION

for f in range(len(experiments)):
    
    # plot cell fates in experimental condition
    plt.figure(dpi=500)
    
    plt.suptitle("SOX2 Intensity of Different Cell Fates", fontsize=14)
    
    titles=['mTeSR 0-48', "BMP 10ng/ml 0-30, Noggin 30-48", 
            "BMP 50ng/ml 0-30, Noggin 30-48",
            'BMP 50ng/ml 0-48']
    
    plt.title(titles[f], fontsize=12)
    
    if bra_pos[bra_pos["exp_con"]==f].empty==False:
       
        sns.lineplot(x="hours", y="SOX2_Intensity", data=bra_pos[bra_pos["exp_con"]==f], color="indigo",
                 err_style="bars", err_kws={'capsize':3, 'ecolor':"indigo"}, label="BRA Positive")
    
    if sox2_pos[sox2_pos["exp_con"]==f].empty==False:    
        sns.lineplot(x="hours", y="SOX2_Intensity", data=sox2_pos[sox2_pos["exp_con"]==f], color="r",
                 err_style="bars", err_kws={'capsize':3, 'ecolor':"r"}, label="SOX2 Positive")
    
    if cdx2_pos[cdx2_pos["exp_con"]==f].empty==False:
        sns.lineplot(x="hours", y="SOX2_Intensity", data=cdx2_pos[cdx2_pos["exp_con"]==f], color="b",
                 err_style="bars", err_kws={'capsize':3, 'ecolor':"b"}, label="CDX2 Positive")
    

#%% FULL DATA PLOT: HISTOGRAM OF SOX2 FOR EACH CELL FATE

dat=[bra_pos, sox2_pos, cdx2_pos]

for f in range(len(dat)):

    plt.figure(dpi=500)
    
    titles=["BRA", "SOX2", "CDX2"]
    
    plt.suptitle("SOX2 Intensity of "+titles[f]+" Positive Cells", fontsize=14)
    
    colors=["indigo", "r", "b"]
    
    dat_plot=dat[f][dat[f]["hours"]==51]
    
    sns.histplot(x="SOX2_Intensity", data=dat_plot, color=colors[f])


plt.figure(dpi=500)

plt.suptitle("SOX2 Intensity of Different Cell Fates", fontsize=14)

sns.histplot(x="SOX2_Intensity", data=all_fates[all_fates["hours"]==51], palette=colors, hue="fate")
 
#%% PLOT CENTROIDS OF CELLS WITH SOX2 EXPRESSION LESS THAN 550

# Necessary for checking that these cells are being rightfully removed

for f in range(len(experiments)):
    
    exp=still_data[f]
    
    mov_nums=np.unique()
    
    exp=exp[exp["SOX2_Intensity"]<550]
    
    mov_nums=np.array(pd.unique(exp["Movie_Number"]))
    
    for g in range(len(mov_nums)):
        
        mov=mov_nums[g]
        
        plt.imshow(r"")
        
#%% SMALL DATA 4 SUBPLOTS: SOX2 EXPRESSION FOR FINAL CELLS IN A MOVIE

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


#%% SMALL DATA PLOTS: H2B INTENISTIES

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


#%% SMALL DATA PLOTS: SOX2 EXPRESSION IN A LINEAGE

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
               















