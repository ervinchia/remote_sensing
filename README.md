# remote_sensing
Simulation of satellite detection of unresolved subpixel fire for PC4262 Remote Sensing.

# files
1. `hotstuff.py`

  A helper `python` script that when imported, initializes common physical constants and imports multiple functions that compute, transfer and integrate radiances, designed for rapid computation as matrices via the `numpy` module. Written with grid simulation in mind.
  
 
2. `fire_modelling.ipynb`

  An interactive Jupyter notebook, contains everything from `fire_funcs.py` and `fire_detection_plots.py`.
  
  
3. `fire_funcs.py`

  All helper functions to do the simulations are found here.
  
  
4. `fire_detection_plots.py`

  Run this script for a example of the simulation. Outputs: 1) GIF of fire spread, 2) detection at saturation temperature, 3) detection at 10% saturation temperature, 4) detection at 50% saturation temperature, and 5) detection using MODIS v4 detection algorithm. Uses `matplotlib.animation` module to make the GIFs.
  
#

# Additional links

* https://code.earthengine.google.com/53303233dc0d27a44ed6cee64288a7ba?hideCode=true

Google Earth Engine script to obtain land surface temperature and vegetation data from MODIS.
