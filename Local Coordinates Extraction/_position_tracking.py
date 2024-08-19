#!/usr/bin/env python
# coding: utf-8

# In[126]:


import csv
from skimage import io
import numpy as np
from pathlib import Path


# In[1]:


def count_samples(path):
    files = Path(path).glob('*.tif')
    
    nfiles = 0
        
    for file in files:
        nfiles += 1
        
    return nfiles


# In[123]:


# returns a list of time steps
def get_pos_data(path, t_step=0.5):
    files = Path(path).glob('*.tif')
    
    pos_data = []
    for file in files:

        file_id = int(str(file).split('-')[-1].split('.')[0]) # get integer ID of file (by splitting file name)

        img = io.imread(str(file)) # creates a 3D array (nframes, x-coord, y-coord) with RGB values

        nframes = img.shape[0] # total number of frames within the sample

        coor = np.where(img==255) 
                # coor[0] = frame numbers (numbered 0 to n_max)
                # coor[1] = x-coordinates of white pixel
                # coor[2] = y-coordinates of white pixel

        # create a new entry [file_id, time_point, x, y]    
        for ii in range(nframes):
            if(ii in coor[0]):
                jj = np.where(coor[0]== ii)[0][0]
                t = ii * t_step
                x = coor[1][jj]
                y = coor[2][jj]
            else:
                t = ii * t_step
                x = -1
                y = -1
            pos_data.append([file_id, t, x, y])
    return pos_data

# returns an array of standard deviations; number of entries is equal to the number of block sizes
def get_mean_stdev(path, block_size=60):
    files = Path(path).glob('*.tif')
    data=[] # initialize - [file name, max stdev]
    stdevs=[] # dummy variable to be changed
    for file in files: # iterates depending on the number of nodes

        img = io.imread(str(file)) # creates a 3D array (nframes, x-coord, y-coord) with RGB values
        
        nframes = img.shape[0] # total number of frames within the sample
        
        # determines number of blocks by dividing total frames by block size
        nblocks = int(np.ceil(nframes/block_size)) 
        
        coor = np.where(img==255) 
            # coor[0] = frame numbers (numbered 0 to n_max)
            # coor[1] = x-coordinates of white pixel
            # coor[2] = y-coordinates of white pixel
        
        coor_blocks = []
        for parameter in coor:
            # split each parameter of the coor array
            # coor_blocks = [ array(split nframes), array(split x), array(split y)]
            coor_blocks.append(np.array_split(parameter[0:nframes], nblocks))
        
        for ii in range(nblocks):
            x_block = coor_blocks[1][ii]
            y_block = coor_blocks[2][ii]
            stdevs.append(np.sqrt(np.std(x_block)**2+ np.std(y_block)**2))

    return stdevs