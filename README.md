# lineagetrack
Code for tracking lineages of human embryonic stem cells for the Warmflash Lab @ Rice University.

celltracker processes .csv files from cellular movies analyzed in Ilastik and produces figures pertaining to SOX2 and BRA expression. celltracker should be used with data taken from movies analyzed in Ilastik with a still image of the cells appended as the last frame. 

Within the Cell Matching folder, Cell_Matcher performs cell matching between a still image and the last frame of the cellular movie through one of three methods: image registration with code given in Image_Registration, matching centroids following segmentation using Cellpose (performed by scripts Cell_Segmenter and Run_Cellpose), or by accepting data processed by code written by the Warmflash lab given by the Warmflash Lab (https://github.com/warmflasha/simpleTracking). 
