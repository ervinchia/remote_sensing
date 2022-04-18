# -*- coding: utf-8 -*-
"""
Created on Sun Apr 17 22:16:01 2022

Running simulations and plotting

@author: hills
"""
import os
os.chdir(__file__[:-23]) # put the things in the same directory to make sure no weirdness with the importing
from fire_funcs import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

### make simulation
print('Starting simulation...')
frames, firefrac_frames = start_sim()
print('Simulation complete')
detected_list = make_detected_temp(frames)

### Start the fire and plot
fig,ax = plt.subplots(figsize=(8,8))
im = plt.imshow(frames[0], cmap = 'coolwarm', vmin = 30, vmax = 1200)
fig.colorbar(im, ax= ax, fraction=.045)
def animate(i):
    plt.title(f'$\\Delta t$ = {i}, mean temp = {round(np.mean(frames[i]),5)}K\nFire fraction = {firefrac_frames[i]}')
    im.set_array(frames[i])
    return [im]
anim = FuncAnimation(fig, animate, frames=N_frames, repeat=False)
if input("Save fire_spread.gif?(Y/N)\n")=="Y":
  file_name = input('Name? ')
  anim.save(file_name+'.gif', writer=PillowWriter(fps=10), dpi=200)

### Detections plots
fire_detected = [[1000,1000],[1000,1000]] # need this to draw the colors

### Pure saturation
fire_list = np.transpose([detected_list[:,0]>331,detected_list[:,1]>400])
### initialise animation
detector_fig = plt.figure()
im_detector = plt.imshow(fire_detected, cmap = 'coolwarm', vmin = 30, vmax = 1200)
plt.vlines(0.5,-0.5,1.5,color='white',linewidth=2)
plt.text(-0.25,1.4,'4 microns',fontsize=16,c='white')
plt.text(.75,1.4,'11 microns',fontsize=16,c='white')
def animate_saturator(i):
    plt.title(f'$\\Delta t$ = {i}, Does it saturate?\nFire fraction = {firefrac_frames[i]}')
    im_detector.set_array((-1+2*fire_list[i,:])*fire_detected)
    return [im_detector]
### start animation
anim_saturator = FuncAnimation(detector_fig, animate_saturator, frames=N_frames, repeat=False)
if input("Save saturation detected gif?(Y/N)\n")=="Y":
  file_name = input('Name? ')
  anim_saturator.save(file_name+'.gif', writer=PillowWriter(fps=10), dpi=200)
### 10 percent and 50 percent deprecated since they don't make sense
##### 10 percent
##fire_list_10 = np.transpose([detected_list[:,0]>275,detected_list[:,1]>40])
##### initialise animation
##detector_fig_10 = plt.figure()
##im_detector = plt.imshow(fire_detected, cmap = 'coolwarm', vmin = 30, vmax = 1200)
##plt.vlines(0.5,-0.5,1.5,color='white',linewidth=2)
##plt.text(-0.25,1.4,'4 microns',fontsize=16,c='white')
##plt.text(.75,1.4,'11 microns',fontsize=16,c='white')
##def animate_saturator(i):
##    plt.title(f'$\\Delta t$ = {i}, Does it saturate?\nFire fraction = {firefrac_frames[i]}')
##    im_detector.set_array((-1+2*fire_list_10[i,:])*fire_detected)
##    return [im_detector]
##### start animation
##anim_10_percent = FuncAnimation(detector_fig_10, animate_saturator, frames=N_frames, repeat=False)
##if input("Save 10% saturation detected gif?(Y/N)\n")=="Y":
##  file_name = input('Name? ')
##  anim_10_percent.save(file_name+'.gif', writer=PillowWriter(fps=10), dpi=200)
##
##### 50 percent
##fire_list_50 = np.transpose([detected_list[:,0]>40,detected_list[:,1]>70])
##### initialise animation
##detector_fig_50 = plt.figure()
##im_detector = plt.imshow(fire_detected, cmap = 'coolwarm', vmin = 30, vmax = 1200)
##plt.vlines(0.5,-0.5,1.5,color='white',linewidth=2)
##plt.text(-0.25,1.4,'4 microns',fontsize=16,c='white')
##plt.text(.75,1.4,'11 microns',fontsize=16,c='white')
##def animate_saturator(i):
##    plt.title(f'$\\Delta t$ = {i}, Does it saturate?\nFire fraction = {firefrac_frames[i]}')
##    im_detector.set_array((-1+2*fire_list_50[i,:])*fire_detected)
##    return [im_detector]
##### start animation
##anim_50_percent = FuncAnimation(detector_fig_50, animate_saturator, frames=N_frames, repeat=False)
##if input("Save 50% saturation detected gif?(Y/N)\n")=="Y":
##  file_name = input('Name? ')
##  anim_50_percent.save(file_name+'.gif', writer=PillowWriter(fps=10), dpi=200)

### MODIS algo (making the figure was a bit weird here, different from the saturation detection since only one channel)
fire_list_modis = np.vectorize(modis4_detect)(detected_list[:,0], detected_list[:,1]).astype(int)
fire_detected = np.array([[1,1],[1,1]])
### initialise animation
detector_fig = plt.figure()
plt.axis('off')
im_detector = plt.imshow(fire_list_modis[0]*fire_detected, cmap = 'coolwarm',vmin=-1,vmax=1.5)
def animate_saturator(i):
    plt.title(f'$\\Delta t$ = {i}, Does it detect?\nFire fraction = {firefrac_frames[i]}')
    im_detector.set_array(-1+2*fire_list_modis[i]*fire_detected)
    return [im_detector]
### start animation
anim_modis = FuncAnimation(detector_fig, animate_saturator, frames=N_frames, repeat=False)
if input("Save MODIS detected gif?(Y/N)\n")=="Y":
  file_name = input('Name? ')
  anim_modis.save('fire_detected_modis.gif', writer=PillowWriter(fps=10), dpi=200)

input('Done, press enter to exit')
