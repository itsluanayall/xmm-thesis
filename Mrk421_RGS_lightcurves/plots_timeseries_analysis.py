from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import glob
from array import array

MJDREF = 50814.0

def plot(x, y, title, xlabel, ylabel, output_folder, dy, dx=[]):
    """
    """
    fig = plt.figure(figsize=(60,20))
    ax = fig.add_subplot(1, 1, 1)
    plt.autoscale()
    #ratio = 20
    #ax.set_aspect(1.0/ax.get_data_ratio()*ratio)
    #plt.axis([51330, 59040, 0, 60])

    #Drop NaN values by making a numpy mask
    mask_nan = np.invert(np.isnan(y)) 
    x = x[mask_nan]
    y = y[mask_nan]
    dy = dy[mask_nan]

    if len(dx)==0:
        
        plt.errorbar(x,y, yerr=dy, color='black', marker='.', ecolor='gray', linestyle='')
    
    else:

        dx = dx[mask_nan]
        plt.errorbar(x,y, yerr=dy, xerr=dx, color='black', marker='.', markersize=0.01, ecolor='gray', linestyle='')
    
    plt.grid(True)
    plt.title(title, fontsize=20)
    plt.xlabel(xlabel, fontsize=10)
    plt.ylabel(ylabel, fontsize=10)
    #plt.xlim(min(x),53000)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.margins(0)

    # Separate in years
    year = 2001
    for mjd in [51910, 52275, 52640, 53005, 53371, 	53736, 54101, 54466,54832, 55197, 55562, 55927, 56293, 56658, 57023, 57388, 57754, 58119, 58484]: 
        plt.text(x=mjd+70, y=57, s=year)
        year+=1
    plt.vlines([51910, 52275, 52640, 53005, 53371, 	53736, 54101, 54466,54832, 55197, 55562, 55927, 56293, 56658, 57023, 57388, 57754, 58119, 58484], min(y), max(y), colors='r', linestyles='solid')
 
    #Save figure 
    plt.savefig(f'{output_folder}/{title}.png')
    plt.show()

if __name__ == "__main__":
    """
    # Go where the result log fits file is placed
    path_log = "/media/xmmsas/thesis/Markarian421/Products"
    os.chdir(path_log)

    hdul = fits.open('obs_table.fits')

    x = hdul[1].data['MJD_avg_time']    
    y = hdul[1].data['RGS_Rate[count/s]']
    dy = hdul[1].data['Stdev_rate']
    
    title = 'Long-term X-ray Variability Lightcurve'
    plot(x,y, dy=dy, title=title, xlabel='MJD', ylabel='Mean Rate [ct/s]' )


    x = y
    dx = dy
    y = hdul[1].data['F_var']
    dy = hdul[1].data['F_var_sigma']
    mask_fvar = np.invert(np.equal(y,-1.0))
    x = x[mask_fvar]
    y = y[mask_fvar]
    y = y*100
    dy = dy[mask_fvar]
    dy = dy*100
    dx = dx[mask_fvar]

    title = 'Fractional Variability vs Rate'
    plot(x, y, dx=dx, dy=dy, title=title, xlabel='Mean rate [ct/s]', ylabel='Fractional Variability [%]')
    """

    #total light curve
    target_dir = "/media/xmmsas/thesis/Markarian421"
    os.chdir(target_dir)
    total_lightcurve_rates = []
    total_lightcurve_errates = []
    total_lightcurve_times = []

    for directory in os.listdir(target_dir):
        if directory.startswith('0'):
            os.chdir(f"{target_dir}/{directory}/rgs")
            
            for filename in glob.glob('*_RGS_rates.ds'):
                with fits.open(filename) as hdul: 
                    x = hdul[1].data['TIME']
                    y = hdul[1].data['RATE']               
                    yerr = hdul[1].data['ERROR']

                #Drop NaN values by making a numpy mask
                mask_nan = np.invert(np.isnan(y)) 
                x = x[mask_nan]
                y = y[mask_nan]
                yerr = yerr[mask_nan]
                
                total_lightcurve_rates.extend(y)
                total_lightcurve_errates.extend(yerr)
                total_lightcurve_times.extend(x)
        
    total_lightcurve_rates = np.asarray(total_lightcurve_rates)
    total_lightcurve_errates = np.asarray(total_lightcurve_errates)
    total_lightcurve_times = np.asarray(total_lightcurve_times)
    
    # Conversion of times (from MET to MJD)
    total_lightcurve_times_mjd = MJDREF + (total_lightcurve_times/86400.0)
    print(total_lightcurve_times_mjd)
    plot(total_lightcurve_times_mjd, total_lightcurve_rates, 
        title="Historical lightcurve evolution Mrk421", xlabel='MJD', 
        ylabel='Rate [ct/s]', dy=total_lightcurve_errates, output_folder=f"{target_dir}/Products")
                

