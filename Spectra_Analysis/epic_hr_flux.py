import numpy as np
import os
from config import CONFIG
from astropy.io import fits
from astropy.table import Table
from tools import *
import pandas as pd
import glob
import random
from scipy.stats import linregress
from scipy.optimize import curve_fit


if __name__ == "__main__":

    target_dir = CONFIG['target_dir']

    #Open epic_obs_table.fits to get the hardness ratio values
    hdul_time = Table.read(os.path.join(target_dir, 'Products', "EPIC_Lightcurves", "EPIC_obs_table.fits"), hdu=1)
    data_time = hdul_time.to_pandas()
    data_time = data_time.sort_values(by='ObsId')

    #Open epic_spectra.fits to get the flux values  (0.2-10)
    os.chdir(os.path.join(target_dir, 'Products', 'EPIC_Spectra'))

    hdul_spec = Table.read(os.path.join(target_dir, 'Products', "EPIC_Spectra", "EPIC_spectra_table.fits"), hdu=1)
    data_spec = hdul_spec.to_pandas()
    data_spec_logpar = data_spec[data_spec['model']=='TBabs*zlogpar']
    data_spec_logpar = data_spec_logpar.sort_values(by='obsid')

    #-------------------------------Plot HR vs flux-----------------------------------
    figure = plt.figure(figsize=(6,5))
    plt.errorbar(x=data_spec_logpar['flux'].values, xerr=(data_spec_logpar['flux'].values-data_spec_logpar['flux_low'].values, data_spec_logpar['flux_up'].values-data_spec_logpar['flux'].values),
                 y=data_time['hr'].values,  color='black', linestyle='', marker='.', markersize=5, ecolor='gray', elinewidth=1, capsize=2, capthick=1)
    
    plt.xlabel('Flux (0.2 - 10  keV )[erg/cm2/s]')
    plt.ylabel('HR')
    plt.savefig(os.path.join(target_dir, 'Products', 'EPIC_Spectra', 'HR_vs_flux.png'))
    plt.close()

    #---------------------------------HR vs rate (soft + hard)--------------------------
    os.chdir(os.path.join(target_dir, 'Products', 'EPIC_Lightcurves'))
    figure2 = plt.figure(figsize=(6,6))
    df_total = pd.DataFrame()
    
    for filename in glob.glob('*.{}'.format("csv")):
        
        observation = filename[0:10]
        rgb = '#%06X' % random.randint(0, 0xFFFFFF)  #create random color
        df = pd.read_csv(filename) #read csv file of single observation
        df_total = df_total.append(df)
        plt.errorbar(data=df, x='rate', y='hr', yerr='hr_err', xerr='erate', linestyle='', ecolor=rgb, color=rgb, label=observation)

    #linear function between x and y
    def fitfunc (x,a,b):
        return a*x + b
    
    # Initial values
    initial_values =(1.,-1)
    pars, covm = curve_fit(fitfunc, df_total['rate'].values, df_total['hr'].values, initial_values, df_total['hr_err'].values) 

    a0, b0 = pars    #parameter of fit
    da, db = np.sqrt(covm.diagonal())   #and its error (from covariance matrix)

    # Print fit results
    print('a = %f +- %f' % (a0, da))
    print('b = %f +- %f' % (b0, db))

    #chi2
    chisq =(((df_total['hr'].values-fitfunc(df_total['rate'].values, a0,b0) )/df_total['hr_err'].values)**2).sum()
    ndof = len(df_total['rate'].values) - 1
    print('Chisquare/ndof = %f/%d' % (chisq, ndof))

    #Correlation coefficient
    slope, intercept, r, p, stderr = linregress(df_total['rate'].values, df_total['hr'].values)
    line = f'Correlation= {r:.2f}'


    func_grid = np.linspace(200, 650, 100)
    #plt.plot(func_grid, fitfunc(func_grid, a0, b0), color = 'blue')
    plt.plot(df_total['rate'].values, intercept + slope * df_total['rate'].values, label=line, color='red')
    plt.xlabel('Total Rate [ct/s]')
    plt.grid()
    plt.ylabel('HR: (H-S)/(H+S)')
    plt.legend(ncol=2)
    plt.savefig(os.path.join(target_dir, 'Products', 'EPIC_Lightcurves', 'hr_vs_rate.png'))
    plt.close()
   