"""
Welcome to plot_timeseries_analysis.py! This script allows to generate some interesting plots concerning 
timeseries, lightcurves and variability quantities. Here follow the possible combinations:


 --fvar  
 --fvar --correlation
 --lightcurve total
 --lightcurve mean
 
"""

from astropy.io import fits
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import glob
from array import array
from brokenaxes import brokenaxes
from config import CONFIG
from tools import *
from scipy.stats import linregress
from argparse import ArgumentParser
import pandas as pd
from matplotlib.patches import Rectangle

target_dir = CONFIG['target_dir']
target_name = CONFIG['target_name']
timescale_fvar = CONFIG['timescale_fvar']
MJDREF = CONFIG['MJDREF']

parser = ArgumentParser(description=__doc__)

parser.add_argument('--fvar', action='store_true',
                    help='make plots of Fractional variability vs rate and vs time.')
parser.add_argument("--correlation", action="store_true", 
                    help="add correlation line to fvar plot")
parser.add_argument("--lightcurve", action="store_true", 
                    help="make lightcurve plot")
parser.add_argument("--total", action="store_true", 
                    help="make historical lightcurve plot")    
parser.add_argument("--mean", action="store_true", 
                    help="make mean lightcurve plot")                   
parser.add_argument("--xs", action="store_true", 
                    help="make excess variance plot")     
parser.add_argument("--brokenaxis", action="store_true", 
                    help="make lightcurve with broken axis (only for markarian 421)")    
args = parser.parse_args()

def plot(x, y, title, xlabel, ylabel, output_folder, dy, dx=[]):
    """
    Makes a plot with the given arguments. One-liner.
    """

    fig = plt.figure(figsize=(10,5))
    ax = fig.add_subplot(1, 1, 1)

    if len(dx)==0:    
        plt.errorbar(x,y, yerr=dy, color='black', ecolor='gray', linestyle='', fmt='.', markersize=1, elinewidth=1, capsize=2, capthick=1)
    
    else:
        plt.errorbar(x,y, yerr=dy, xerr=dx, color='black', ecolor='gray', fmt='.', linestyle='', markersize=1, elinewidth=1, capsize=2, capthick=1)
    
    plt.grid(True)
    plt.title(title, fontsize=20)
    plt.xlabel(xlabel, fontsize=10)
    plt.ylabel(ylabel, fontsize=10)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    
    if args.fvar and args.correlation:

        #Correlation coefficient
        slope, intercept, r, p, stderr = linregress(x,y)
        line = f'Correlation: {r:.2f}'
        plt.plot(x, intercept + slope * x, label=line, color='red')
        plt.legend()

def plot_total_lc(x, y, title, xlabel, ylabel, output_folder, dy, dx=[]):
    """
    """
    plt.errorbar(x,y, yerr=dy, color='b', marker='.', ecolor='gray', linestyle='', markersize=0.5)  
    plt.show()

    fig, axs = plt.subplots(3,11, sharey=True, figsize=(30,12), gridspec_kw = {'wspace':0.2, 'hspace':0.2})
    fig.suptitle(title, fontsize=20, y=0.92)
    #Drop NaN values by making a numpy mask
    mask_nan = np.invert(np.isnan(y)) 
    x = x[mask_nan]
    y = y[mask_nan]
    dy = dy[mask_nan]

    # Subplot 1
    axs[0,0].errorbar(x,y, yerr=dy, color='b', marker='.', ecolor='gray', linestyle='', markersize=0.5)  
    axs[0,0].grid(True)
    axs[0,0].set_xlim(51683,51697)
    axs[0,0].annotate('Year 2000', xy=(0, 1), xycoords='axes fraction', fontsize=8, xytext=(5, -5), textcoords='offset points',
                horizontalalignment='left', verticalalignment='top', color='r')
    axs[0,0].margins(0)

    # Subplot 2
    axs[0,1].errorbar(x,y, yerr=dy, color='b', marker='.', ecolor='gray', linestyle='', markersize=0.5)
    axs[0,1].grid(True)
    axs[0,1].set_xlim(51850,51863)
    axs[0,1].margins(0)

    #Subplot 3
    axs[0,2].errorbar(x,y, yerr=dy, color='b', marker='.', ecolor='gray', linestyle='', markersize=0.5)
    axs[0,2].grid(True)
    axs[0,2].set_xlim(52035,52040)
    axs[0,2].annotate('Year 2001', xy=(0, 1), xycoords='axes fraction', fontsize=8, xytext=(5,-5), textcoords='offset points',
                horizontalalignment='left', verticalalignment='top', color='r')
    axs[0,2].margins(0)

    
    #Subplot 4
    axs[0,3].errorbar(x,y, yerr=dy, color='b', marker='.', ecolor='gray', linestyle='', markersize=0.5)
    axs[0,3].grid(True)
    axs[0,3].set_xlim(52393,52400)
    axs[0,3].margins(0) 
    axs[0,3].get_xaxis().get_major_formatter().set_useOffset(False)
    axs[0,3].annotate('Year 2002', xy=(0, 1), xycoords='axes fraction', fontsize=8, xytext=(5,-5), textcoords='offset points',
                horizontalalignment='left', verticalalignment='top', color='r')


    #Subplot 5
    axs[0,4].errorbar(x,y, yerr=dy, color='b', marker='.', ecolor='gray', linestyle='', markersize=0.5)
    axs[0,4].grid(True)
    axs[0,4].set_xlim(52581,52612)
    axs[0,4].margins(0)

    #Subplot 6
    axs[0,5].errorbar(x,y, yerr=dy, color='b', marker='.', ecolor='gray', linestyle='', markersize=0.5)
    axs[0,5].grid(True)
    axs[0,5].set_xlim(52789,52800)
    axs[0,5].margins(0)
    axs[0,5].annotate('Year 2003', xy=(0, 1), xycoords='axes fraction', fontsize=8, xytext=(5,-5), textcoords='offset points',
                horizontalalignment='left', verticalalignment='top', color='r')
    axs[0,5].get_xaxis().get_major_formatter().set_useOffset(False)
    
    #Subplot 7
    axs[0,6].errorbar(x,y, yerr=dy, color='b', marker='.', ecolor='gray', linestyle='', markersize=0.5)
    axs[0,6].grid(True)
    axs[0,6].set_xlim(52956,52985)
    axs[0,6].margins(0)

    
    #Subplot 8
    axs[0,7].errorbar(x,y, yerr=dy, color='b', marker='.', ecolor='gray', linestyle='', markersize=0.5)
    axs[0,7].grid(True)
    axs[0,7].annotate('Year 2004', xy=(0, 1), xycoords='axes fraction', fontsize=8, xytext=(5,-5), textcoords='offset points',
                horizontalalignment='left', verticalalignment='top', color='r')
    axs[0,7].set_xlim(53125,53136)
    axs[0,7].margins(0)

    #Subplot 9
    axs[0,8].errorbar(x,y, yerr=dy, color='b', marker='.', ecolor='gray', linestyle='', markersize=0.5)
    axs[0,8].grid(True)
    axs[0,8].annotate('Year 2005', xy=(0, 1), xycoords='axes fraction', fontsize=8, xytext=(5,-5), textcoords='offset points',
                horizontalalignment='left', verticalalignment='top', color='r')
    axs[0,8].set_xlim(53678,53690)
    axs[0,8].margins(0)

    #Subplot 10
    axs[0,9].errorbar(x,y, yerr=dy, color='b', marker='.', ecolor='gray', linestyle='', markersize=0.5)
    axs[0,9].grid(True)
    axs[0,9].annotate('Year 2006', xy=(0, 1), xycoords='axes fraction', fontsize=8, xytext=(5,-5), textcoords='offset points',
                horizontalalignment='left', verticalalignment='top', color='r')
    axs[0,9].set_xlim(53854,53884)
    axs[0,9].margins(0)

    #Subplot 11
    axs[0,10].errorbar(x,y, yerr=dy, color='b', marker='.', ecolor='gray', linestyle='', markersize=0.5)
    axs[0,10].grid(True)
    axs[0,10].set_xlim(54068, 54080)
    axs[0,10].margins(0)
    

    #Subplot 12
    axs[1,0].errorbar(x,y, yerr=dy, color='b', marker='.', ecolor='gray', linestyle='', markersize=0.5)
    axs[1,0].grid(True)
    axs[1,0].annotate('Year 2007', xy=(0, 1), xycoords='axes fraction', fontsize=8, xytext=(5,-5), textcoords='offset points',
                horizontalalignment='left', verticalalignment='top', color='r')
    axs[1,0].set_xlim(54220, 54235)
    axs[1,0].set_ylabel(ylabel, fontsize=14)
    axs[1,0].margins(0)


    #Subplot 13
    axs[1,1].errorbar(x,y, yerr=dy, color='b', marker='.', ecolor='gray', linestyle='', markersize=0.5)
    axs[1,1].grid(True)
    axs[1,1].set_xlim(54418,54430)
    axs[1,1].margins(0)


    #Subplot 14
    axs[1,2].errorbar(x,y, yerr=dy, color='b', marker='.', ecolor='gray', linestyle='', markersize=0.5)
    axs[1,2].grid(True)
    axs[1,2].annotate('Year 2008', xy=(0, 1), xycoords='axes fraction', fontsize=8, xytext=(5,-5), textcoords='offset points',
                horizontalalignment='left', verticalalignment='top', color='r')
    axs[1,2].set_xlim(54592,54618)
    axs[1,2].margins(0)

    #Subplot 15
    axs[1,3].errorbar(x,y, yerr=dy, color='b', marker='.', ecolor='gray', linestyle='', markersize=0.5)
    axs[1,3].grid(True)
    axs[1,3].set_xlim(54785,54800)
    axs[1,3].margins(0)

    #Subplot 16
    axs[1,4].errorbar(x,y, yerr=dy, color='b', marker='.', ecolor='gray', linestyle='', markersize=0.5)
    axs[1,4].grid(True)
    axs[1,4].annotate('Year 2009', xy=(0, 1), xycoords='axes fraction', fontsize=8, xytext=(5,-5), textcoords='offset points',
                horizontalalignment='left', verticalalignment='top', color='r')
    axs[1,4].set_xlim(54965, 54980)
    axs[1,4].margins(0)

    #Subplot 17
    axs[1,5].errorbar(x,y, yerr=dy, color='b', marker='.', ecolor='gray', linestyle='', markersize=0.5)
    axs[1,5].grid(True)
    axs[1,5].set_xlim(55145, 55160)
    axs[1,5].margins(0)

    #Subplot 18
    axs[1,6].errorbar(x,y, yerr=dy, color='b', marker='.', ecolor='gray', linestyle='', markersize=0.5)
    axs[1,6].grid(True)
    axs[1,6].annotate('Year 2010', xy=(0, 1), xycoords='axes fraction', fontsize=8, xytext=(5,-5), textcoords='offset points',
                horizontalalignment='left', verticalalignment='top', color='r')
    axs[1,6].set_xlim(55310, 55325)
    axs[1,6].margins(0)

    #Subplot 19
    axs[1,7].errorbar(x,y, yerr=dy, color='b', marker='.', ecolor='gray', linestyle='', markersize=0.5)
    axs[1,7].grid(True)
    axs[1,7].set_xlim(55505, 55520)
    axs[1,7].margins(0)

    #Subplot 20
    axs[1,8].errorbar(x,y, yerr=dy, color='b', marker='.', ecolor='gray', linestyle='', markersize=0.5)
    axs[1,8].grid(True)
    axs[1,8].annotate('Year 2011', xy=(0, 1), xycoords='axes fraction', fontsize=8, xytext=(5,-5), textcoords='offset points',
                horizontalalignment='left', verticalalignment='top', color='r')
    axs[1,8].set_xlim(55694, 55708)
    axs[1,8].get_xaxis().get_major_formatter().set_useOffset(False)
    axs[1,8].margins(0)

    #Subplot 21
    axs[1,9].errorbar(x,y, yerr=dy, color='b', marker='.', ecolor='gray', linestyle='', markersize=0.5)
    axs[1,9].grid(True)
    axs[1,9].set_xlim(55873,55901)
    axs[1,9].margins(0)

    #Subplot 22
    axs[1,10].errorbar(x,y, yerr=dy, color='b', marker='.', ecolor='gray', linestyle='', markersize=0.5)
    axs[1,10].grid(True)
    axs[1,10].annotate('Year 2014', xy=(0, 1), xycoords='axes fraction', fontsize=8, xytext=(5,-5), textcoords='offset points',
                horizontalalignment='left', verticalalignment='top', color='r')
    axs[1,10].set_xlim(56770, 56785)
    axs[1,10].margins(0)

    #Subplot 23
    axs[2,0].errorbar(x,y, yerr=dy, color='b', marker='.', ecolor='gray', linestyle='', markersize=0.5)
    axs[2,0].grid(True)
    axs[2,0].annotate('Year 2015', xy=(0, 1), xycoords='axes fraction', fontsize=8, xytext=(5,-5), textcoords='offset points',
                horizontalalignment='left', verticalalignment='top', color='r')
    axs[2,0].set_xlim(57175, 57190)
    axs[2,0].margins(0)
    

    #Subplot 24
    axs[2,1].errorbar(x,y, yerr=dy, color='b', marker='.', ecolor='gray', linestyle='', markersize=0.5)
    axs[2,1].grid(True)
    axs[2,1].set_xlim(57328, 57340)
    axs[2,1].margins(0)

    #Subplot 25
    axs[2,2].errorbar(x,y, yerr=dy, color='b', marker='.', ecolor='gray', linestyle='', markersize=0.5)
    axs[2,2].grid(True)
    axs[2,2].set_xlim(57360, 57375)
    axs[2,2].margins(0)

    #Subplot 26
    axs[2,3].errorbar(x,y, yerr=dy, color='b', marker='.', ecolor='gray', linestyle='', markersize=0.5)
    axs[2,3].grid(True)
    axs[2,3].annotate('Year 2016', xy=(0, 1), xycoords='axes fraction', fontsize=8, xytext=(5,-5), textcoords='offset points',
                horizontalalignment='left', verticalalignment='top', color='r')
    axs[2,3].set_xlim(57513, 57535)
    axs[2,3].margins(0)

    #Subplot 27
    axs[2,4].errorbar(x,y, yerr=dy, color='b', marker='.', ecolor='gray', linestyle='', markersize=0.5)
    axs[2,4].grid(True)
    axs[2,4].set_xlim(57692, 57707)
    axs[2,4].get_xaxis().get_major_formatter().set_useOffset(False)
    axs[2,4].margins(0)

    #Subplot 28
    axs[2,5].errorbar(x,y, yerr=dy, color='b', marker='.', ecolor='gray', linestyle='', markersize=0.5)
    axs[2,5].grid(True)
    axs[2,5].annotate('Year 2017', xy=(0, 1), xycoords='axes fraction', fontsize=8, xytext=(5,-5), textcoords='offset points',
                horizontalalignment='left', verticalalignment='top', color='r') 
    axs[2,5].set_xlim(57870, 57885)
    axs[2,5].margins(0)

    #Subplot 29
    axs[2,6].errorbar(x,y, yerr=dy, color='b', marker='.', ecolor='gray', linestyle='', markersize=0.5)
    axs[2,6].grid(True)
    axs[2,6].set_xlim(58069, 58080)
    axs[2,6].margins(0)


    #Subplot 30
    axs[2,7].errorbar(x,y, yerr=dy, color='b', marker='.', ecolor='gray', linestyle='', markersize=0.5)
    axs[2,7].grid(True)
    axs[2,7].annotate('Year 2018', xy=(0, 1), xycoords='axes fraction', fontsize=8, xytext=(5,-5), textcoords='offset points',
                horizontalalignment='left', verticalalignment='top', color='r')
    axs[2,7].set_xlim(58235, 58250)
    axs[2,7].margins(0)

    #Subplot 31
    axs[2,8].errorbar(x,y, yerr=dy, color='b', marker='.', ecolor='gray', linestyle='', markersize=0.5)
    axs[2,8].grid(True)
    axs[2,8].set_xlim(58438, 58449)
    axs[2,8].get_xaxis().get_major_formatter().set_useOffset(False)
    axs[2,8].margins(0)

    #Subplot 32
    axs[2,9].errorbar(x,y, yerr=dy, color='b', marker='.', ecolor='gray', linestyle='', markersize=0.5)
    axs[2,9].grid(True)
    axs[2,9].annotate('Year 2019', xy=(0, 1), xycoords='axes fraction', fontsize=8, xytext=(5,-5), textcoords='offset points',
                horizontalalignment='left', verticalalignment='top', color='r')
    axs[2,9].set_xlim(58603, 58612)
    axs[2,9].margins(0)

    #Subplot 33
    axs[2,10].errorbar(x,y, yerr=dy, color='b', marker='.', ecolor='gray', linestyle='', markersize=0.5)
    axs[2,10].grid(True)
    axs[2,10].set_xlabel(xlabel, fontsize=14)
    axs[2,10].set_xlim(58812, 58822)
    axs[2,10].margins(0)

    #Save figure 
    plt.savefig(f'{output_folder}/{title}.png')
    plt.show()

if __name__ == "__main__":

    #Create Products directory
    if not os.path.isdir(f'{target_dir}/Products/Plots_timeseries'):
        os.makedirs(f'{target_dir}/Products/Plots_timeseries')
    
    #----------------Mean lightcurve plot------------#
    if args.lightcurve and args.mean:
        
        # Open table containing lightcurve information
        hdul = Table.read(os.path.join(target_dir, 'Products', 'RGS_Lightcurves', 'obs_table.fits'), hdu=1)    
        
        # Put table in pandas dataframe 
        data = hdul.to_pandas()
        data.dropna()
        data = data.sort_values(by=['MJD_avg_time'])
        time_err = data['Duration_Obs'].values/86400
       
        #Make plot
        title = 'Long-term X-ray Variability Lightcurve'
        plot(x=data['MJD_avg_time'].values, y=data['RGS_Rate'].values, dx=time_err, dy=data['RGS_erate'].values, title=title, xlabel='MJD', ylabel='Mean Rate [ct/s]', output_folder=os.path.join(target_dir,'Products','Plots_timeseries'))
        plt.savefig(os.path.join(target_dir,'Products','Plots_timeseries', 'mean_lightcurve.png'))

    #--------------Fvar plots------------------------#
    if args.fvar:
        
        # Open table containing lightcurve and variability information
        hdul = Table.read(os.path.join(target_dir,'Products','RGS_Lightcurves','obs_table.fits'), hdu=1)    
        data = hdul.to_pandas()
        data = data.dropna()
        #data = data[data['ObsId'] != 150498701]
        data = data[data['F_var']!=-1.000]
        data['F_var'] = data['F_var'].apply(lambda x: x*100)  #for % units
        data['F_var_sigma'] = data['F_var_sigma'].apply(lambda x: x*100)
        data = data.sort_values(by=['RGS_Rate'])
        
        #Make fvar vs rate plot
        title = f'Fractional Variability vs Rate'
        print("# datapoints fvar =", len(data))
        plot(data['RGS_Rate'].values, data['F_var'].values, dx=data['RGS_erate'].values, dy=data['F_var_sigma'].values, title=title, output_folder=os.path.join(target_dir,'Products','Plots_timeseries'), xlabel='Mean rate [ct/s]', ylabel='Fractional Variability [%]')
        plt.savefig(os.path.join(target_dir,'Products','Plots_timeseries', 'fvar_vs_rate.png'))


        #Make fvar vs time plot
        data = data.sort_values(by=['MJD_avg_time'])
        title = 'Fractional Variability vs time'
        plot(data['MJD_avg_time'].values, data['F_var'].values, dy=data['F_var_sigma'].values, title=title, output_folder=f"{target_dir}/Products/Plots_timeseries", xlabel='MJD [d]', ylabel='Fractional Variability [%]')
        plt.savefig(os.path.join(target_dir,'Products','Plots_timeseries', 'fvar_vs_time'))

    #--------------Excess variance plots-------------#
    if args.xs:

        #Open table containing variability information
        hdul = Table.read(os.path.join(target_dir,'Products', 'RGS_Lightcurves','obs_table.fits'), hdu=1)    
        data = hdul.to_pandas()
        data = data.dropna()
        #data = data[data['ObsId'] != 150498701]
        data = data[data['F_var']!=-1.000]
        data = data.sort_values(by=['MJD_avg_time'])
        print("# datapoints excess variance =", len(data))
        
        #Square root of xs 
        xs_linear = np.sqrt(data['Excess_Variance'].values)
        xs_linear_err = 0.5*xs_linear * (data['xs_sigma'].values/data['Excess_Variance'].values)
        
        #Make xs vs time plot
        plot(x=data['MJD_avg_time'].values, y=xs_linear, dy=xs_linear_err, title='linear xs vs time', output_folder=f"{target_dir}/Products/Plots_timeseries", xlabel='MJD', ylabel='$\sqrt{\sigma_{XS}^2}$')
        plt.savefig(os.path.join(target_dir,'Products','Plots_timeseries', 'linear_xs_vs_time.png'))
        
        #Make xs vs rate plot
        plot(x=data['RGS_Rate'].values, y=xs_linear, dx=data['RGS_erate'].values, dy=xs_linear_err, title='linear xs vs rate', output_folder=f"{target_dir}/Products/Plots_timeseries", xlabel='Rate[ct/s]', ylabel='$\sqrt{\sigma_{XS}^2}$')
        plt.savefig(os.path.join(target_dir,'Products','Plots_timeseries', 'linear_xs_vs_rate.png'))

    #-------------Historic RGS lightcurve------------#
    if args.lightcurve and args.total and target_name=="Markarian 421":

        total_lightcurve_rates = []
        total_lightcurve_errates = []
        total_lightcurve_times = []

        for directory in os.listdir(target_dir):
            if directory.startswith('0'):  #then it's a directory of an observation
                os.chdir(os.path.join(target_dir, directory, 'rgs'))

                for filename in glob.glob('*_RGS_rates.ds'):
                    x, y, yerr, fracexp, y_bg, yerr_bg = mask_fracexp15(filename)
                    total_lightcurve_rates.extend(y)
                    total_lightcurve_errates.extend(yerr)
                    total_lightcurve_times.extend(x)
            
        total_lightcurve_rates = np.asarray(total_lightcurve_rates)
        total_lightcurve_errates = np.asarray(total_lightcurve_errates)
        total_lightcurve_times = np.asarray(total_lightcurve_times)
        
        # Conversion of times (from MET to MJD)
        total_lightcurve_times_mjd = MJDREF + (total_lightcurve_times/86400.0)
        
        # Make lightcurve pandas dataframe 
        data_lc = pd.DataFrame({"RATE":total_lightcurve_rates, "MJD":total_lightcurve_times_mjd, "ERROR":total_lightcurve_errates})
        data_lc = data_lc.sort_values(by=['MJD'])
        data_lc.to_csv(os.path.join(target_dir, 'Products', 'Plots_timeseries', 'data_lc.csv'))

        # Separate data into years
        year_array = []
        for mjd in data_lc['MJD'].values:
            if mjd<51910:
                year_array.append(int(2000))
            elif mjd<52275:
                year_array.append(int(2001))
            elif mjd<52640:
                year_array.append(int(2002))
            elif mjd<53005:
                year_array.append(int(2003))
            elif mjd<53371:
                year_array.append(int(2004))
            elif mjd<53736:
                year_array.append(int(2005))
            elif mjd<54101:
                year_array.append(int(2006))
            elif mjd<54466:
                year_array.append(int(2007))
            elif mjd<54832:
                year_array.append(int(2008))
            elif mjd<55197:
                year_array.append(int(2009))
            elif mjd<55562:
                year_array.append(int(2010))
            elif mjd<55927:
                year_array.append(int(2011))
            elif mjd<56293:
                year_array.append(int(2012))
            elif mjd<56658:
                year_array.append(int(2013))
            elif mjd<57023:
                year_array.append(int(2014))
            elif mjd<57388:
                year_array.append(int(2015))
            elif mjd<57754:
                year_array.append(int(2016))
            elif mjd<58119:
                year_array.append(int(2017))
            elif mjd<58484:
                year_array.append(int(2018))
            elif mjd<58849:
                year_array.append(int(2019))

        # Add 'YEAR' column to dataframe
        data_lc['YEAR'] = year_array  
        data_lc = data_lc.reset_index(drop=True)
        data_lc = data_lc.reset_index()
        #sample_data = data_lc[0:792]
        #sample_data.to_csv(f'{target_dir}/Products/data_lc_mrk421.csv', index=False)    

        #Broken axis plot - define limits of observations (only implemented for Marakarian 421!)
        if args.brokenaxis:

            lims = ((51689,51690), (51849,51851), (51861, 51864), (52037,52038),
            (52398, 52400),(52582, 52584),(52592, 52594), (52609, 52612),
            (52791,52793), (52797, 52798),(52957, 52959), (52983,52985),
            (53131, 53132), (53681, 53685), (53854,53856), (53883,53884), 
            (54074, 54075), (54228, 54231), (54423, 54426), (54593, 54594),
            (54617, 54618), (54792, 54794), (54976, 54977), (55151, 55153),
            (55319, 55320), (55510, 55516), (55698, 55703), (55874, 55875),
            (55893, 55895), (55900, 55901), (56776, 56781), (57179, 57187.5),
            (57334, 57337), (57364, 57365), (57370, 57371), (57514, 57515), 
            (57518, 57522.5), (57527.5, 57530.5), (57533.5, 57534.5),
            (57695, 57698.5), (57705, 57706), (57876.5, 57878), (57880.5, 57881.5),
            (58076, 58077.5), (58239.5, 58241), (58247.5, 58249), (58443, 58444.5),
            (58609, 58610), (58816, 58817))

            fig_brkax = plt.figure(figsize=(60,8))
            bax = brokenaxes(xlims=lims, wspace=0.2, tilt=90, diag_color='red', d=0.0015)
            bax.errorbar(data_lc['MJD'].values, data_lc['RATE'].values, data_lc['ERROR'].values, linestyle='', markersize=0.05, marker='.')
            my_ticks = [51689,51690, 51849,51851, 51861, 51864, 52037,52038,
                                    52398, 52400,52582, 52584,52592, 52594, 52609, 52612,
                                    52791,52793, 52797, 52798,52957, 52959, 52983,52985,
                                    53131, 53132, 53681, 53685, 53854,53856, 53883,53884, 
                                    54074, 54075, 54228, 54231, 54423, 54426, 54593, 54594,
                                    54617, 54618, 54792, 54794, 54976, 54977, 55151, 55153,
                                    55319, 55320, 55510, 55516, 55698, 55703, 55874, 55875,
                                    55893, 55895, 55900, 55901, 56776, 56781, 57179, 57187.5,
                                    57334, 57337, 57364, 57365, 57370, 57371, 57514, 57515, 
                                    57518, 57522.5, 57527.5, 57530.5, 57533.5, 57534.5,
                                    57695, 57698.5, 57705, 57706, 57876.5, 57878, 57880.5, 57881.5,
                                    58076, 58077.5, 58239.5, 58241, 58247.5, 58249, 58443, 58444.5,
                                    58609, 58610, 58816, 58817]
            bax.set_xlabel('MJD', labelpad=60, ha="right", fontsize=20)
            bax.set_ylabel('Rate [ct/s]', fontsize=20)
            bax.tick_params(rotation=60)
            bax.ticklabel_format(useOffset=False, style='plain')
            bax.grid(axis='both', which='major')
            bax.grid(axis='both', which='minor', alpha=0.4)
            
            year = 2001
            year_endpoints = []
            for i in range(1, len(data_lc)):
                if data_lc['YEAR'][i] != data_lc['YEAR'][i-1]:
                    year_endpoints.append(data_lc['MJD'][i])
                    
            bax.vlines(year_endpoints, 0, 60, colors='r', linestyles='solid')
            bax.axs[0].text(x=data_lc['MJD'][0], y=58, s=str(2000), c='r')
            bax.text(x=year_endpoints[0]+0.10, y=58, s=str(2001), c='r')
            bax.text(x=year_endpoints[1]+0.10, y=58, s=str(2002), c='r')
            bax.axs[8].text(x=year_endpoints[2]+0.10, y=58, s=str(2003), c='r')
            bax.text(x=year_endpoints[3]+0.10, y=58, s=str(2004), c='r')
            bax.text(x=year_endpoints[4]+0.10, y=58, s=str(2005), c='r')
            bax.axs[14].text(x=year_endpoints[5]+0.10, y=58, s=str(2006), c='r')
            bax.text(x=year_endpoints[6]+0.10, y=58, s=str(2007), c='r')
            bax.axs[19].text(x=year_endpoints[7]+0.10, y=58, s=str(2008), c='r')
            bax.text(x=year_endpoints[8]+0.10, y=58, s=str(2009), c='r')
            bax.text(x=year_endpoints[9]+0.10, y=58, s=str(2010), c='r')
            bax.text(x=year_endpoints[10]+0.10, y=58, s=str(2011), c='r')
            bax.text(x=year_endpoints[11]+0.10, y=58, s=str(2014), c='r')
            bax.text(x=year_endpoints[12]+0.10, y=58, s=str(2015), c='r')
            bax.axs[35].text(x=year_endpoints[13]+0.10, y=58, s=str(2016), c='r')
            bax.axs[41].text(x=year_endpoints[14]+0.10, y=58, s=str(2017), c='r')
            bax.axs[44].text(x=year_endpoints[15]+0.10, y=58, s=str(2018), c='r')
            bax.text(x=year_endpoints[16]+0.10, y=58, s=str(2019), c='r')

            bax.margins(0)
            plt.suptitle('Historical Lightcurve Mrk421', fontsize=25)
            plt.savefig(f"{target_dir}/Products/Plots_timeseries/lc_broken_axis.png")

        # Make distribution of observations 
        from collections import Counter
        import seaborn as sns
        cnt = Counter(year_array)
        labels = ['2000', '2001', '2002', '2003', '2004', '2005', '2006', '2007',
            '2008', '2009', '2010', '2011',  '2014', '2015', '2016',
            '2017', '2018', '2019']

        fig = plt.figure(figsize=(15,15))
        ax = fig.add_subplot(111)
        ax.bar(labels, cnt.values())
        ax.tick_params( rotation = 60)
        ax.set_ylabel('# datapoints') 
        plt.savefig(os.path.join(target_dir, 'Products', 'Plots_timeseries', 'distrib_data.png'))
        
        #Plot rate distribution 
        fig_hist_rate = plt.figure(figsize=(6,6))
        counts, bins, _ = plt.hist(data_lc['RATE'].values, bins=10,color='skyblue' ,edgecolor='black', linewidth=1.2, histtype='stepfilled')
        bin_centers = 0.5 * (bins[:-1] + bins[1:])
        
        #Errors for bins
        bin_prob = counts/len(data_lc['RATE'].values)
        bin_up = (bin_prob + 1.96*np.sqrt(bin_prob*(1-bin_prob)/len(data_lc['RATE'].values)) ) * len(data_lc['RATE'].values)
        bin_low = (bin_prob - 1.96*np.sqrt(bin_prob*(1-bin_prob)/len(data_lc['RATE'].values)) ) * len(data_lc['RATE'].values)
        bin_top = bin_up-counts
        bin_bot = counts - bin_low
        plt.errorbar(bin_centers, counts, yerr=(bin_bot, bin_top), fmt='.', color='black')

        # Lognormal fit
        from scipy.optimize import curve_fit
        from scipy.stats import lognorm

        def lognorm_func(x, mu, sigma):
            return np.exp(-(np.log(x) - mu)**2 / (2 * sigma**2))/ (x * sigma * np.sqrt(2 * np.pi))

        # Restore data from histogram: counts multiplied bin centers
        restored = [[d]*int(counts[n]) for n,d in enumerate((bins[1:]+bins[:-1])/2)]
        
        # Flatten the result
        restored = [item for sublist in restored for item in sublist]

        print(lognorm.fit(restored, floc=0))

        dist = lognorm(*lognorm.fit(restored, floc=0))
        x = np.arange(1,60)
        y = dist.pdf(x)

        # The pdf is normalized, so we need to scale it to match the histogram
        y = y/y.max()
        y = y*counts.max()

        # Chi2 lognormal 
        y2 = dist.pdf(bin_centers) 
        y2 = y2/y2.max()
        y2 = y2*counts.max()
        chisq =(((counts-y2)/bin_top)**2).sum()
        ndof = len(counts) - 3
        print('Chisquare/ndof lognormal = %f/%d' % (chisq, ndof))

        # Gauss fit
        def Gauss(x, a, x0, sigma):
            return a * np.exp(-(x - x0)**2 / (2 * sigma**2))

        popt,pcov = curve_fit(Gauss, bin_centers, counts, p0=[max(counts), 23, 1])

        # Chi2 gauss
        chisq_gauss =(((counts-Gauss(bin_centers, popt[0], popt[1], popt[2]))/bin_top)**2).sum()
        ndof_gauss = len(counts) - 3
        print('Chisquare/ndof gauss = %f/%d' % (chisq_gauss, ndof_gauss))

        plt.plot(x,y,'r',linewidth=2, label='Lognormal Fit $\chi^2_{red}$ = '+ f'{chisq/ndof:.2f}')
        plt.plot(x,Gauss(x, popt[0], popt[1], popt[2]), 'g',linewidth=2, label='Gauss Fit $\chi^2_{red}$ = '+  f'{chisq_gauss/ndof_gauss:.2f}')
        plt.xlim(0,60)
       
        plt.xlabel("Rate [ct/s]")
        plt.ylabel('# datapoints')
        plt.xticks()
        plt.legend()
        plt.savefig(os.path.join(target_dir, 'Products', 'Plots_timeseries', 'rate_histogram.png'))
        
        #Plot compact lightcurve (x=index, y=rate, hue=year)
        g = sns.FacetGrid(data=data_lc, hue='YEAR', height=8, aspect=2)
        g.map(plt.errorbar, 'index', 'RATE', 'ERROR', fmt='.', ecolor='gray', elinewidth=1, capsize=2, capthick=1)
        g.fig.suptitle("Mrk421 total data")
        plt.grid(True)
        g.add_legend()
        plt.savefig(os.path.join(target_dir, 'Products', 'Plots_timeseries', "compact_lightcurve.png"))
        plt.show()
        

    if args.lightcurve and args.total:
        #TOTAL LIGHT CURVE
        os.chdir(target_dir)
        total_lightcurve_rates = []
        total_lightcurve_errates = []
        total_lightcurve_times = []
        observations = []
        for directory in os.listdir(target_dir):
            if directory.startswith('0'):
                os.chdir(f"{target_dir}/{directory}/rgs")
                
                for filename in glob.glob('*_RGS_rates.ds'):
                    x, y, yerr, fracexp, y_bg, yerr_bg = mask_fracexp15(filename)
                    total_lightcurve_rates.extend(y[:-1])
                    total_lightcurve_errates.extend(yerr[:-1])
                    total_lightcurve_times.extend(x[:-1])
                    observations.extend([int(directory) for i in range(len(x)-1)])
            
        total_lightcurve_rates = np.asarray(total_lightcurve_rates)
        total_lightcurve_errates = np.asarray(total_lightcurve_errates)
        total_lightcurve_times = np.asarray(total_lightcurve_times)
        observations = np.asarray(observations)
      
        # Conversion of times (from MET to MJD)
        total_lightcurve_times_mjd = MJDREF + (total_lightcurve_times/86400.0)

        data_lc = pd.DataFrame({"RATE":total_lightcurve_rates, "TIME":total_lightcurve_times, "ERROR":total_lightcurve_errates, "MJD": total_lightcurve_times_mjd, "OBSERVATION":observations})
        data_lc = data_lc.sort_values(by=['TIME'])
        
        #Add year column to dataframe
        year_array = []
        for mjd in data_lc['MJD'].values:
            if mjd<51910:
                year_array.append(int(2000))
            elif mjd<52275:
                year_array.append(int(2001))
            elif mjd<52640:
                year_array.append(int(2002))
            elif mjd<53005:
                year_array.append(int(2003))
            elif mjd<53371:
                year_array.append(int(2004))
            elif mjd<53736:
                year_array.append(int(2005))
            elif mjd<54101:
                year_array.append(int(2006))
            elif mjd<54466:
                year_array.append(int(2007))
            elif mjd<54832:
                year_array.append(int(2008))
            elif mjd<55197:
                year_array.append(int(2009))
            elif mjd<55562:
                year_array.append(int(2010))
            elif mjd<55927:
                year_array.append(int(2011))
            elif mjd<56293:
                year_array.append(int(2012))
            elif mjd<56658:
                year_array.append(int(2013))
            elif mjd<57023:
                year_array.append(int(2014))
            elif mjd<57388:
                year_array.append(int(2015))
            elif mjd<57754:
                year_array.append(int(2016))
            elif mjd<58119:
                year_array.append(int(2017))
            elif mjd<58484:
                year_array.append(int(2018))
            elif mjd<58849:
                year_array.append(int(2019))

        data_lc['YEAR'] = year_array  
        data_lc = data_lc.reset_index(drop=True)
        data_lc = data_lc.reset_index()
        
        #indexes where to place the xticks in the plot
        year_endpoints = []
        for i in range(1, len(data_lc)):
            if data_lc['YEAR'][i] != data_lc['YEAR'][i-1]:
                year_endpoints.append(data_lc['index'][i])
        labels = np.linspace(start=2000,stop=2019,num=20,dtype=int)
        labels = np.delete(labels, [12,13])  #delete year 2012 and 2013
        year_endpoints = np.array(year_endpoints)
        year_endpoints = np.insert(year_endpoints, 0, values=0)

        fig = plt.figure(figsize=(15,5))
        plt.errorbar('index', 'RATE', 'ERROR', data=data_lc, ecolor='black', linestyle='', color='black')
        #Add vertical lines separating years
        plt.vlines(year_endpoints, 0, 20, colors='r', linestyles='solid')
        plt.grid(axis='y')
        plt.xlabel('Year', labelpad=0.1, fontsize=15)
        plt.ylabel('Rate [ct/s]', labelpad=0.1, fontsize=15)
        plt.margins(0)
        plt.tight_layout()
        plt.xticks(ticks=year_endpoints, labels=labels, rotation=50, fontsize=11)
        plt.savefig(os.path.join(target_dir, "Products", "Plots_timeseries", "compact_LC_vlines.png"))