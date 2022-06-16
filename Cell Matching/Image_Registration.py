# -*- coding: utf-8 -*-
"""
Created on Mon May 30 10:43:15 2022

@author: grace
"""

# Code modified from code provided by GeeksforGeeks
# source: https://www.geeksforgeeks.org/image-registration-using-opencv-python/

# Align max projection images of cells at the final time point in the 
# experiment with the last frame of the movie in order to collect data from
# more fluorescense channels.

def run_script(align, reference, save_loc, align_loc, ref_loc):

    import cv2
    import numpy as np
    import os
    from skimage import io
    import matplotlib as plt
    
    
    # Open the image files.
    
    # image to be aligned: second channel of the max proj image
    os.chdir(align_loc)
    img1=io.imread(align, img_num=1)
    img1=np.uint8(img1)
    
    # reference image: first channel of the last frame of the movie
    os.chdir(ref_loc)
    img2=io.imread(reference, img_num=0)
    img2=np.uint8(img2)
    
    height, width = img2.shape
    
    # Create ORB detector with 5000 features.
    orb_detector = cv2.ORB_create(1000)
    
    # Find keypoints and descriptors.
    # The first arg is the image, second arg is the mask
    # (which is not required in this case).
    kp1, d1 = orb_detector.detectAndCompute(img1, None)
    kp2, d2 = orb_detector.detectAndCompute(img2, None)
    
    # plot key points
    imgplot = cv2.drawKeypoints(img1, kp1, None, color=(0,255,0), flags=0)
    plt.figure(), io.imshow(imgplot)
    imgplot = cv2.drawKeypoints(img2, kp2, None, color=(0,255,0), flags=0)
    plt.figure(), io.imshow(imgplot)
    
    # Match features between the two images.
    # We create a Brute Force matcher with
    # Hamming distance as measurement mode.
    matcher = cv2.BFMatcher(cv2.NORM_HAMMING, crossCheck = True)
    
    # Match the two sets of descriptors.
    matches = matcher.match(d1, d2)
    matches = list(matches)
    
    # Sort matches on the basis of their Hamming distance.
    matches.sort(key = lambda x: x.distance)
    
    # Take the top 90 % matches forward.
    matches = matches[:int(len(matches)*0.1)]
    no_of_matches = len(matches)
    
    # Define empty matrices of shape no_of_matches * 2.
    p1 = np.zeros((no_of_matches, 2))
    p2 = np.zeros((no_of_matches, 2))
    
    for i in range(len(matches)):
        p1[i, :] = kp1[matches[i].queryIdx].pt
        p2[i, :] = kp2[matches[i].trainIdx].pt
    
    # Find the homography matrix.
    homography, mask = cv2.findHomography(p1, p2, cv2.RANSAC)
    
    # Use this matrix to transform the
    # colored image wrt the reference image.
    transformed_img = cv2.warpPerspective(img1,
    					homography, (width, height))
    
    os.chdir(save_loc)
    
    # Save the output.
    cv2.imwrite(align+'.tif', transformed_img)
    
    plt.figure(), io.imshow(transformed_img)
    
    plt.figure(dpi=1000)    
    result = cv2.drawMatches(img1, kp1, img2, kp2, matches, img1, flags = 2)
    plt.imshow(result)
    plt.show()
