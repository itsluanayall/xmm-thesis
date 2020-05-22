from astropy.io import fits
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import glob
from array import array

import pandas as pd
from matplotlib.patches import Rectangle
target_dir = "/media/xmmsas/thesis/Markarian421"
MJDREF = 50814.0

def plot(x, y, title, xlabel, ylabel, output_folder, dy, dx=[]):
    """
    """
    fig = plt.figure(figsize=(60,20))
    ax = fig.add_subplot(1, 1, 1)

    #Drop NaN values by making a numpy mask
    mask_nan = np.invert(np.isnan(y)) 
    x = x[mask_nan]
    y = y[mask_nan]
    dy = dy[mask_nan]

    if len(dx)==0:
        
        plt.errorbar(x,y, yerr=dy, color='black', marker='.', ecolor='gray', linestyle='-')
    
    else:

        dx = dx[mask_nan]
        plt.errorbar(x,y, yerr=dy, xerr=dx, color='black', marker='.', ecolor='gray', linestyle='-')
    
    plt.grid(True)
    plt.title(title, fontsize=20)
    plt.xlabel(xlabel, fontsize=10)
    plt.ylabel(ylabel, fontsize=10)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.margins(0)

    #Save figure 
    plt.savefig(f'{output_folder}/{title}.png')
    plt.show()

def plot_total_lc(x, y, title, xlabel, ylabel, output_folder, dy, dx=[]):
    """
    """


    fig = plt.figure(figsize=(80,10))
    ax = fig.add_subplot(1, 1, 1)
    #plt.autoscale()
    #ratio = 0.2
    #ax.set_aspect(1.0/ax.get_data_ratio()*ratio)
    #plt.axis([51330, 59040, 0, 60])
    #Drop NaN values by making a numpy mask
    mask_nan = np.invert(np.isnan(y)) 
    x = x[mask_nan]
    y = y[mask_nan]
    dy = dy[mask_nan]

    if len(dx)==0:
        
        plt.errorbar(x,y, yerr=dy, color='r', marker='.', ecolor='gray', linestyle='',)
    
    else:

        dx = dx[mask_nan]
        plt.errorbar(x,y, yerr=dy, xerr=dx, color='black', marker='.', markersize=0.05, ecolor='gray', linestyle='-')
    plt.yscale('log')
    plt.grid(True)
    plt.title(title, fontsize=20)
    plt.xlabel(xlabel, fontsize=10)
    plt.ylabel(ylabel, fontsize=10)
    plt.xlim(51840,51870)
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
    # LONG TERM VARIABILITY PLOT
    path_log = "/media/xmmsas/thesis/Markarian421/Products"
    os.chdir(path_log)

    hdul = Table.read('obs_table.fits', hdu=1)    
    data = hdul.to_pandas()
    data = data.sort_values(by=['MJD_avg_time'])

    title = 'Long-term X-ray Variability Lightcurve'
    plot(data['MJD_avg_time'].values ,data['RGS_Rate'].values, dy=data['Stdev_rate'].values, title=title, xlabel='MJD', ylabel='Mean Rate [ct/s]', output_folder=f"{target_dir}/Products" )

    

    # FRACTIONAL VARIABILITY PLOTS
    hdul = Table.read('obs_table.fits', hdu=1)    
    data = hdul.to_pandas()
    data = data[data['F_var']!=-1.000]
    data['F_var'] = data['F_var'].apply(lambda x: x*100)
    data['F_var_sigma'] = data['F_var_sigma'].apply(lambda x: x*100)
    data = data.sort_values(by=['RGS_Rate'])

    title = 'Fractional Variability vs Rate'
    plot(data['RGS_Rate'].values, data['F_var'].values, dx=data['F_var_sigma'].values, dy=data['Stdev_rate'].values, title=title, output_folder=f"{target_dir}/Products", xlabel='Mean rate [ct/s]', ylabel='Fractional Variability [%]')
    
    data = data.sort_values(by=['MJD_avg_time'])

    title = 'Fractional Variability vs time'
    plot(data['MJD_avg_time'].values, data['F_var'].values, dx=data['F_var_sigma'].values, dy=data['Stdev_rate'].values, title=title, output_folder=f"{target_dir}/Products", xlabel='MJD [d]', ylabel='Fractional Variability [%]')
    
    """
    
    #TOTAL LIGHT CURVE
    os.chdir(target_dir)
    total_lightcurve_rates = []
    total_lightcurve_errates = []
    total_lightcurve_times = []

    for directory in os.listdir(target_dir):
        if directory.startswith('0'):
            os.chdir(f"{target_dir}/{directory}/rgs")
            
            for filename in glob.glob('*_RGS_rates.ds'):
                hdul = Table.read(filename, hdu=1)    
                data = hdul.to_pandas()
                data.dropna()

                #with fits.open(filename) as hdul: 
                #    x = hdul[1].data['TIME']
                #    y = hdul[1].data['RATE']               
                #    yerr = hdul[1].data['ERROR']

                #Drop NaN values by making a numpy mask
                #mask_nan = np.invert(np.isnan(y)) 
                #x = x[mask_nan]
                #y = y[mask_nan]
                #yerr = yerr[mask_nan]
                
                total_lightcurve_rates.extend(data['RATE'].values)
                total_lightcurve_errates.extend(data['ERROR'].values)
                total_lightcurve_times.extend(data['TIME'].values)
        
    total_lightcurve_rates = np.asarray(total_lightcurve_rates)
    total_lightcurve_errates = np.asarray(total_lightcurve_errates)
    total_lightcurve_times = np.asarray(total_lightcurve_times)
    
    # Conversion of times (from MET to MJD)
    total_lightcurve_times_mjd = MJDREF + (total_lightcurve_times/86400.0)
    
    data_lc = pd.DataFrame({"RATE":total_lightcurve_rates, "MJD":total_lightcurve_times_mjd, "ERROR":total_lightcurve_errates})
    data_lc = data_lc.sort_values(by=['MJD'])
    #sns.catplot(x='MJD', y='RATE',  data=data_lc)
    #plt.show()
    #plt.savefig(f"{target_dir}/Products/seaborn.png")
    """
    plot_total_lc(total_lightcurve_times_mjd, total_lightcurve_rates, 
        title="Historical lightcurve evolution Mrk421", xlabel='MJD', 
        ylabel='Rate [ct/s]', dy=total_lightcurve_errates, output_folder=f"{target_dir}/Products")
                
    """
    plot_total_lc(data_lc['MJD'].values, data_lc['RATE'].values, dy=data_lc['ERROR'].values,title="Historical lightcurve evolution Mrk421", xlabel='MJD', 
        ylabel='Rate [ct/s]',  output_folder=f"{target_dir}/Products")

