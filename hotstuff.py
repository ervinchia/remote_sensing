from array import array
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
h = 6.626E-34
c = 299792458
k = 1.38E-23
sigma = 5.67E-8
matplotlib.use('TKAgg')

def renorm(Temp:float, lda:float) -> float:
    '''
    Takes temperature `T` and wavelength in microns `lda`, returns the thermalized ratio `x`.
    '''
    return((h*c)/(k*Temp*lda))

def longwave(x:float)-> float:
    '''
    takes a `renorm` value `x` and returns the longwave approximation.
    '''
    f = ((15*x**3)/np.pi**4)*((1/3) - (x/8) + (x**2)/60)
    return(f)

def shortwave(x:float)-> float:
    '''
    takes a `renorm` value `x` and returns the shortwave term approximation.
    '''
    g = 15/(np.pi**4) * (np.exp(-1*x)*(x**3+3*x**2+6*x+6) + (np.exp(-2*x)/8)*(4*x**3+6*x**2+6*x+3)) 
    return(g)

def planck_estimate(x_1:float, x_2:float, temp:float, propagate:bool = True, attenuate:int = 4, xflip:float = 2.25, legacy:bool = True)-> float:
    '''
    Uses `shortwave` and `longwave` approximation function to yield the estimated blackbody radiance over `x_2` and `x_1` at temperature `temp`.
    `attenuate` takes either `4` or `11` microns, which attenuates the radiance at that wavelength.
    '''
    if legacy == False:
        if type(x_1) == float:
            if max(x_1, x_2) <= xflip:
                m = sigma*(temp**4)*(longwave(x_2) - longwave(x_1))
            else:
                m = sigma*(temp**4)*(shortwave(x_2) - shortwave(x_1))
        else:
            m=np.where(x_2 <= xflip, sigma*(temp**4)*(longwave(x_1)-longwave(x_2)), sigma*(temp**4)*(shortwave(x_2)-shortwave(x_1)))
    else:
        m= sigma*(temp**4)*(shortwave(x_2) - shortwave(x_1))
    if attenuate == 4:
        m = 0.98*m
    elif attenuate == 11:
        m = 0.99*m
    return(m)

def radiatedpower(T, wl ,bw = 1E-7, emis:float = 1):
    '''
    Accepts temperature `T` and wavelength `wl`, applies the `planck_estimate` at the appropriate shortwave/longwave limit. Also accepts different bandwidths `bw`, defaults to `1E-7` and emissivity `emis`, default to 1. Emissivity is assumed to be constant over a short `bw`.'''
    if wl == 11E-6:
        atn = 11
    else:
        atn = 4
    x = emis*planck_estimate(renorm(Temp = T, lda = wl-bw), renorm(Temp =T, lda = wl+bw), temp = T, attenuate=atn)
    return x

def brighttemp(lda, power, bw = 1E-7, err_check = None):
    '''
    Returns the brightness temperature based on the shortwave approximation. bandwidth `bw` is onesided to be consistent with other functions; the function will double it when computing spectral radiance. `err_check` accepts `float` or `int` as a input temperature for checking. If given, the function returns a tuple of (`T_b`, `x`-`err_check`). 
    '''
    x = (h*c/k/lda)/(np.log((bw*2*2*np.pi*h*c**2/lda**5/power)+1))
    if isinstance((err_check),(float, int)):
        x = (x, x - err_check)
    return x

def gridpowerfinder(grid:np.ndarray, wl:float, bw = 1E-7)-> np.ndarray:
    '''
    accepts an arbitrary grid and applies the `planck_estimate` function element-wise at the wavelength `wl`. Bandwidth `bw` is onesided to be consistent with other functions; the function will double it when computing spectral radiance. The computed element-wise grid power is normalized to grid size with equal weightage, reflecting an isotropic pixel detector on the satellite. Returns a `radiancegrid` of the same dimensions of input `grid`. 
    '''
    gridc = grid.copy()
    gridsize = np.size(gridc)
    if wl == 11E-6:
        atn = 11
    else:
        atn = 4
    radiancegrid = planck_estimate(renorm(gridc, lda=wl-bw), renorm(gridc, lda=wl+bw), gridc, attenuate = atn)/gridsize
    return(radiancegrid)

def gridbrightness(rgrid: np.ndarray, wl:float, bandwidth = 1E-7)-> list[float]:
    '''
    accepts a normalized radiance grid from `gridpowerfinder` and applies the `brighttemp` function element-wise at the wavelength `wl`. Bandwidth `bandwidth` (default `1E-7`) is onesided to be consistent with other functions; the function will double it when computing temperature. The individual radiances in each element of the grid is summed and used as the total power emitted by the grid for brightness temperature computation.
    '''
    rgridc = rgrid.copy()
    x = brighttemp(lda = wl, power = rgridc.sum(axis = None),bw = bandwidth, err_check= None)
    return(x)

def totalfinder(tempgrid:np.ndarray, bandwidth = 1E-7, wllist:list = [4E-6, 11E-6])-> list[float]:
    '''
    Accepts the temperature grid and parses it through `gridpowerfinder` and `gridbrightness` together at wavelengths `wllist` that defaults to 4 and 11 microns. Returns a list of brightness temperatures at the wavelengths in `wllist` that is detected by the satellite.'''
    result = []
    tempgridc = tempgrid.copy()
    for i in wllist:
        result.append(gridbrightness(gridpowerfinder(tempgridc, i),i))
    return(result)

# %%
# spare code to get apparent brightness with fire temp, recreating prof's plot
plt.rcParams.update({'font.size': 16})

firefraction = np.logspace(-4,0,40) 
firetemp = 650
bgtemp = 305.273
averaged_power_4 = firefraction*radiatedpower(firetemp, wl = 4E-6) + (1-firefraction)*radiatedpower(bgtemp, wl = 4E-6)
averaged_power_11 = firefraction*radiatedpower(firetemp, wl = 11E-6) + (1-firefraction)*radiatedpower(bgtemp, wl = 11E-6)
fig, ax = plt.subplots(figsize = (10,7))
ax.plot(firefraction, brighttemp(4E-6, averaged_power_4), label = '4 μm')
ax.plot(firefraction, brighttemp(11E-6, averaged_power_11), label = '11 μm')
ax.set(xlabel = 'fire fraction', ylabel ='brightness temperature', title = f'Apparent brightness temperature at TOA, Fire temperature = {firetemp}K')
ax.legend()
ax.set_xscale('log')
ax.grid(alpha = 0.8)
plt.tight_layout()
plt.savefig('firefraction', dpi = 300)
plt.show()

# %%
# more spare code, thermal discrepancies plot
# templist = np.linspace(300,1500)
# v = 1E-7
# plist4 = planck_estimate(renorm(templist, 4E-6-v), renorm(templist,4E-6+v), temp= templist, xflip= 1.75, legacy = True)
# plist11 = planck_estimate(renorm(templist, 11E-6 -v), renorm(templist,11E-6 +v), temp= templist, xflip =1.75, legacy = True)
# bt4 = brighttemp(4E-6, plist4)
# bt11 = brighttemp(11E-6, plist11)
# fig, ax = plt.subplots(figsize = (10,8))
# ax.plot(templist, (bt4-templist), label = '4 μm')
# ax.plot(templist, (bt11-templist), label = '11 μm')
# ax.set(xlabel = 'Absolute Temperature, K', ylabel = 'Absolute error', title = r'Shortwave approximation thermal discrepancies')
# ax.grid(alpha =0.8)
# ax.legend()
# plt.tight_layout()
# plt.savefig('sw_abs_error.png', dpi = 300)
# plt.show()

# %%
