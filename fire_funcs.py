# -*- coding: utf-8 -*-
"""
Created on Sun Apr 17 22:14:30 2022

Functions for fire simulations and detection

@author: hills
"""

import numpy as np
from hotstuff import totalfinder

# initialise global variables
temp_bg = 32.123
temp_grid = np.full((100,100),temp_bg)
bmass_grid = np.random.uniform(0,1,(100,100))
check_grid = np.zeros((100,100))
N_frames = input('Number of frames for simulation? ')

def temp_curve(temp_0, tpeak, burn_time):
    """ Returns the temperature evolution curve for all time, needs params to control peak temp and how long it burns for """
    time = np.linspace(0,1000,N_frames)

    B = burn_time
    A = - B**2 / 4 / (tpeak - temp_0)

    y = A*time**2 + B*time + temp_0
    y[y<temp_bg] = temp_bg

    return y

# note for future: improve the fire simulation algo

def temp_evolution(time, temp_0, tpeak, burn_time):
    """ Increments timestep from time = 0 by time """
    temp_ls = temp_curve(temp_0, tpeak, burn_time)
#     print(temp_ls[time])
    return temp_ls[time]

# at every timestep, check the array for temp values. make it a binomial process with chance to spread proportional to T of a location

def spread_fire(temp_self, k, coords):
    """ Given a temp_self, determines probability of spreading and spreads the fire to other grids """
    prob_spread = k * min((temp_self - temp_bg)/600,1)
    x_coord, y_coord = coords
    for i in [-1,0,1]:
        for j in [-1,0,1]:
            if i==0 and j==0: continue
            if x_coord+i<0 or y_coord+j< 0: continue
            if np.random.rand() < prob_spread:
                try:
                    if check_grid[x_coord+i, y_coord+j]==0:
                        temp_grid[x_coord+i, y_coord+j] = np.random.normal(450,50)
                        check_grid[x_coord+i, y_coord+j]=1
                except IndexError: pass

    return temp_grid

def start_sim():
  """ Starts simulation """
  temp_grid[50,50] = 500 # seed with fire, temperature = 500 celsius

  frames = [] # make frames for the gif
  firefrac_frames = [] # list of fire frac for each frame

  for _ in range(N_frames):
      old_temp_grid = np.copy(temp_grid)
      for coord, item in np.ndenumerate(old_temp_grid):
          k = bmass_grid[coord]
          temp_grid = spread_fire(item, k, coord)
          if check_grid[coord]!=0:
              temp_grid[coord] = temp_evolution(int(check_grid[coord]-1), 500, 1100, 8)
              check_grid[coord]+=1
      firefrac_frames.append(np.sum(check_grid!=0)/np.size(check_grid))
      frames.append(temp_grid.copy()+273) # convert to Kelvins

  return frames, firefrac_frames

# for the detection part
def modis4_detect(T4, T11):
    """ If the temperature detected by satellite exceeds 310K. """
    if T11>360: return True
    elif T4>310 and T4-T11>10: return True
    else: return False

def make_detected_temp(frames):
  """ Find the temperature detected by the satellite """
  detected_list = []
  for temp_arr in frames:
      detected_list.append(totalfinder(temp_arr))
  return np.array(detected_list)