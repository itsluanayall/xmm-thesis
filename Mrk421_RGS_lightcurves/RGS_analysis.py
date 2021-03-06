"""
Welcome to the RGS Analysis Software for XMM-Newton! 

At the current version, running this script will process RGS data of the XMM-Newton observatory 
and plot the lightcurves and spectra.
This software was made to analyse the source (target) Markarian421. The data has been previously 
downloaded in the target directory and each observation is labelled with an ID (obsid). Each
observation is a directory inside the target directory.
Inside each observation directory, there should be the 'odf' and 'rgs' directories, where the 
ODFs and the RGS products will be stored respectively.\n

Because of the large number of observations to be processed, this software is driven by
the object-oriented paradigm, so each observation is an instance of the class 'Observation', and 
all the operations made on the data are encapsulated in the methods of this class. 
Note that it will take a few hours to finish executing if you're running the code for the first time.\n


The python code shouldn't be modified by the user, but in order to make it run on you computer,
make sure to set all the user-dependent variables in the config.json file of the package.
In particular:
"sas_dir" is the directory where SAS is installed;
"ccf_dir" is the directory where the CCF files are stored;
"target_name" is the name ID of the source to analyse;
"target_dir" is the directory of the source to analyse, e.g. 'Markarian421', that contains the directories of the observation(s);
"target_RA" right ascension coordinate of the chosen source;
"target_DEC" declination coordinate of the chosen source;
"target_REDSHIFT" redshift of the chosen source.

"""

import logging
import os
from observation import Observation
from tools import run_command, setupSAS, constant
from config import CONFIG
from astropy.table import Table
import matplotlib.pyplot as plt
from astropy.io import ascii
import pandas as pd
import glob
import numpy as np
import random
from scipy.optimize import curve_fit
logging.basicConfig(level=logging.INFO)

#Introduction message for the user
print(__doc__)

if __name__ == "__main__":

    #Set configuration variables of the config.json file. Don't forget to change them according to your needs!
    version = CONFIG['version']
    sas_dir = CONFIG['sas_dir']
    ccf_dir = CONFIG['ccf_dir']
    mjdref = CONFIG['MJDREF']
    target_name = CONFIG['target_name']
    target_dir = CONFIG['target_dir']
    target_REDSHIFT = CONFIG['target_REDSHIFT']
    timescale_fvar = CONFIG['timescale_fvar']
    logging.info(f'{target_name} Analysis - version:{version}')
    setupSAS(sas_dir=sas_dir, ccf_dir=ccf_dir)
    logging.info(f'Timescale chosen for fractional variability: {timescale_fvar} ks')

    #Remove the .tar files from which we extracted the data
    for directory in os.listdir(target_dir):
        if directory.endswith('.tar.gz'):
            os.remove(os.path.join(target_dir, directory))
    
    #Create Products directory
    if not os.path.isdir(os.path.join(target_dir, 'Products')):
        os.makedirs(os.path.join(target_dir, 'Products'))
        os.makedirs(os.path.join(target_dir, 'Products', 'RGS_Lightcurves'))
        os.makedirs(os.path.join(target_dir, 'Products', 'Plots_timeseries'))
        os.makedirs(os.path.join(target_dir, 'Products', 'RGS_Spectra'))
        os.makedirs(os.path.join(target_dir, 'Products', 'Backgrnd_LC'))

    #Loop analysis for each observation

    mrk421_observation_list = []
    obs_table = Table(names=('ObsId', 'RevolutionId', 'ExposureID', 'Start', 
                            'End', 'Duration_Obs', 'RGS_Rate',
                            'RGS_erate', 'MJD_avg_time',
                            'F_var', 'F_var_sigma', 'Excess_Variance', 
                             'xs_sigma',   'Norm_excess_variance', 'nxs_sigma', 'VA', 'VA_sigma'), 
                    dtype=('i', 'i', 'U9', 'U30',
                            'U30', 'd', 'd', 
                            'd', 'd',
                            'd', 'd', 'd',
                            'd', 'd', 'd', 'd', 'd'))

    counter = 0
    mrk421_problematic_obs = ['0658802001', '0411082701']
    pks_problematic_obs = ['0411780834', '0411780501']
    duration_lc_ks = []
    
    for obsid in os.listdir(target_dir):
        
        if obsid.startswith('0'):   #All observation folders start with 0
            print('----------------')
            if obsid in mrk421_problematic_obs or obsid in pks_problematic_obs:
                print(f'The observation {obsid} is tricky. Please analyse individually. Moving forward to next observation.')
                print('----------------')
                continue

            obs = Observation(obsid=obsid, target_dir=target_dir)   #instance of the observation
            mrk421_observation_list.append(obs)
            
            #Process each observation
            obs.cifbuild()
            obs.odfingest()
            obs.rgsproc()
            obs.create_pairs_exposures()
            if obs.skipobs == True:
                print(f'The observation {obsid} is tricky. Please analyse individually. Moving forward to next observation.')
                print('----------------')
                continue
            obs.bkg_lightcurve()
            obs.check_flaring_particle_bkgr()
            obs.rgslccorr()
            obs.lightcurve(mjdref=mjdref)
            obs.fracvartest(screen=True, instrument='rgs')
            obs.vaughan_panel(N=30, M=20, timescale=timescale_fvar, timebinsize=25)
            obs.divide_spectrum()
            obs.xspec_divided_spectra_average(target_REDSHIFT)
            obs.xspec_divided_spectra(target_REDSHIFT)

            #Save attributes of observable into a table
            if len(obs.rgsrate)==0:
                obs_table.add_row((str(obs.obsid), str(obs.revolution),  None , str(obs.starttime), str(obs.endtime), str(int(obs.duration)),
                                None, None, None, None,
                                None , None, None,
                                None, None, None, None))
                

            else:
    
                obs_table.add_row((str(obs.obsid), str(obs.revolution),  f"{obs.expoid[0][0]}+{obs.expoid[0][1]}" , str(obs.starttime), str(obs.endtime), str(int(obs.duration)),
                                    "{:.3f}".format(obs.rgsrate[0]), "{:.3f}".format(obs.stdev[0]), "{:.3f}".format(obs.longterm_lc_times[0]), "{:.3f}".format(obs.fracvardict[0].get('Fractional Variability')), "{:.3f}".format(obs.fracvardict[0].get('Fractional Variability Error')),
                                    "{:.3f}".format(obs.fracvardict[0].get('Excess variance')), "{:.3f}".format(obs.fracvardict[0].get('Excess variance error')),
                                    "{:.5f}".format(obs.fracvardict[0].get('Normalized excess variance')), "{:.5f}".format(obs.fracvardict[0].get('Normalized excess variance error')),
                                    "{:.3f}".format(obs.fracvardict[0].get('Variability Amplitude')), "{:.3f}".format(obs.fracvardict[0].get('Variability amplitude error'))))
                duration_lc_ks.append(obs.duration_lc_ks[0])

                if len(obs.rgsrate)>1:
                    for i in range(1, len(obs.rgsrate)):
                        obs_table.add_row((str(obs.obsid), str(obs.revolution),  f"{obs.expoid[i][0]}+{obs.expoid[i][1]}" , str(obs.starttime), str(obs.endtime), str(int(obs.duration)),
                                            "{:.3f}".format(obs.rgsrate[i]), "{:.3f}".format(obs.stdev[i]), "{:.3f}".format(obs.longterm_lc_times[i]), "{:.3f}".format(obs.fracvardict[i].get('Fractional Variability')), "{:.3f}".format(obs.fracvardict[i].get('Fractional Variability Error')),
                                            "{:.3f}".format(obs.fracvardict[i].get('Excess variance')), "{:.3f}".format(obs.fracvardict[i].get('Excess variance error')),
                                            "{:.5f}".format(obs.fracvardict[i].get('Normalized excess variance')), "{:.5f}".format(obs.fracvardict[i].get('Normalized excess variance error')),
                                            "{:.3f}".format(obs.fracvardict[i].get('Variability Amplitude')), "{:.3f}".format(obs.fracvardict[i].get('Variability amplitude error'))))
                        duration_lc_ks.append(obs.duration_lc_ks[i])

            #Keep track of number of observations that have been processed so far
            counter += 1
            logging.info(f'Processed {counter} observations!')
            print('----------------')
       
    #Show and save in csv format the table with all the attributes of the observations
    print(obs_table)
    obs_table['Duration_Obs'].unit = 's'
    obs_table['RGS_Rate'].unit = 'ct/s'
    obs_table['RGS_erate'].unit = 'ct/s'
    obs_table['MJD_avg_time'].unit = 'd'
    obs_table['Excess_Variance'].unit = 'ct2/s2'
    obs_table['xs_sigma'].unit = 'ct2/s2'
    obs_table.write(output=f'{target_dir}/Products/RGS_Lightcurves/obs_table.fits', format='fits', overwrite=True)

    # Duration_lc_ks distribution
    hist = plt.figure(figsize=(10,10))
    plt.hist(duration_lc_ks, bins=40)
    plt.savefig(f'{target_dir}/Products/RGS_Lightcurves/distribution_expos_duration_ks.png')
    plt.close()
    
    #Combined long observations xs vs rate
    os.chdir(os.path.join(target_dir, 'Products', 'Plots_timeseries'))
    fig_xs_rate, axs  = plt.subplots(2, 1, figsize=(8, 8), sharex=True, gridspec_kw={'hspace':0.1})   
    
    rate_fit = []
    fvar_fit = []
    fvar_err_fit = []
    for filename in glob.glob('*.{}'.format("csv")):
        if filename!='data_lc.csv':
            df_xs_rate = pd.read_csv(filename) #read csv file of single observation
            rgb = '#%06X' % random.randint(0, 0xFFFFFF)  #create random color
            try:
                rate_fit.append(df_xs_rate['rate'].values)
                fvar_fit.append(df_xs_rate['fvar'].values)
                fvar_err_fit.append(df_xs_rate['fvar_err'].values)
                obsid_xs = df_xs_rate['observation'].values[0]

                initial_values =[np.mean(df_xs_rate['fvar'].values)]
                pars, covm = curve_fit(constant, df_xs_rate['rate'].values, df_xs_rate['fvar'].values, initial_values,df_xs_rate['fvar_err'].values) 

                q0 = pars    #parameter of fit
                dq = np.sqrt(covm.diagonal())   #and its error (from covariance matrix)

                # Print fit results
                print(f'Fit constant for obs {obsid_xs} Fvar vs rate:')
                print('q = %f +- %f' % (q0, dq))
                
                #chi2
                chisq =(((df_xs_rate['fvar'].values-constant(df_xs_rate['rate'].values, q0) )/df_xs_rate['fvar_err'].values)**2).sum()
                ndof = len(df_xs_rate['rate'].values) - 1
                print('Chisquare/ndof = %f/%d' % (chisq, ndof))


                axs[0].errorbar(data=df_xs_rate, x='rate', y='xs', yerr='xs_err', xerr='erate', fmt='.', markersize=10, ecolor='gray', elinewidth=1, capsize=2, capthick=1, color=rgb,label=obsid_xs)
                axs[1].errorbar(data=df_xs_rate, x='rate', y='fvar', yerr='fvar_err', xerr='erate', fmt='.', markersize=10, ecolor='gray', elinewidth=1, capsize=2, capthick=1, color=rgb, label=obsid_xs)
            
            except IndexError as e:
                print(f'The file {filename} is likely empty. Please check the .csv files in the Products directory.')

    axs[0].legend(title='Observation ID', fancybox=True)
    axs[0].set_ylabel('$<\sigma_{XS}^2>$')
    axs[0].grid()
    axs[1].set_xlabel('Rate [ct/s]')
    axs[1].set_ylabel('$<F_{var}>$')
    axs[1].grid(True)
    '''
    # Initial values
    fvar_fit = np.concatenate(fvar_fit).ravel()
    fvar_err_fit = np.concatenate(fvar_err_fit).ravel()
    rate_fit = np.concatenate(rate_fit).ravel()

    initial_values =[np.mean(fvar_fit)]
    pars, covm = curve_fit(constant, rate_fit, fvar_fit, initial_values,fvar_err_fit) 

    q0 = pars    #parameter of fit
    dq = np.sqrt(covm.diagonal())   #and its error (from covariance matrix)

    # Print fit results
    print('Fit constant for Fvar vs rate:')
    print('q = %f +- %f' % (q0, dq))
    
    #chi2
    chisq =(((fvar_fit-constant(rate_fit, q0) )/fvar_err_fit)**2).sum()
    ndof = len(rate_fit) - 1
    print('Chisquare/ndof = %f/%d' % (chisq, ndof))

    axs[1].hlines(q0, min(rate_fit), max(rate_fit), color='red', label=f"Constant fit: {q0[0]:.3f} +- {dq[0]:.3f} \n $\chi^2$/ndof = {chisq:.2f} / {ndof}")
    '''

    plt.savefig(os.path.join(target_dir, 'Products', 'Plots_timeseries','xs_rate_combined.png'))
    plt.close()
    '''
    #For a single observation
    #sample_obs = ['0099280101', '0153951201', '0670920301', '0810860701'] #for spectra
    #sample_obs = ['0153951201']
    sample_obs = ['0136541001', '0158971201', '0810860201', '0411080301', '0560980101', '0791781401', '0810860701', '0791782001', '0791780201'] #for vaughan panels
    #sample_obs = ['0791782001']
    #sample_obs = ['0136540701']
    
    for observation in sample_obs:
            
        obs = Observation(obsid=observation, target_dir=target_dir)   
        
        #Process each observation
        obs.cifbuild()
        obs.odfingest()
        obs.rgsproc()
        obs.create_pairs_exposures()
        obs.bkg_lightcurve()
        obs.check_flaring_particle_bkgr()
        obs.rgslccorr()
        obs.lightcurve(mjdref=mjdref, use_grace=use_grace)
        obs.fracvartest(screen=True)
        obs.vaughan_panel(N=15, M=15, timescale=60, timebinsize=25)
        obs.divide_spectrum()
        obs.xspec_divided_spectra_average(target_REDSHIFT)
        obs.xspec_divided_spectra(target_REDSHIFT)
    
    os.chdir(os.path.join(target_dir, 'Products', 'Plots_timeseries'))

    fig_xs_rate, axs  = plt.subplots(2, 1, figsize=(8, 8), sharex=True, gridspec_kw={'hspace':0.1})   
    
    for filename in glob.glob('*.{}'.format("csv")):
        if filename!='data_lc.csv':
            df_xs_rate = pd.read_csv(filename) #read csv file of single observation
            rgb = '#%06X' % random.randint(0, 0xFFFFFF)  #create random color
            axs[0].errorbar(data=df_xs_rate, x='rate', y='xs', yerr='xs_err', xerr='erate', fmt='.', markersize=10, ecolor='gray', elinewidth=1, capsize=2, capthick=1, color=rgb, label=df_xs_rate['observation'][0])
            axs[1].errorbar(data=df_xs_rate, x='rate', y='fvar', yerr='fvar_err', xerr='erate', fmt='.', markersize=10, ecolor='gray', elinewidth=1, capsize=2, capthick=1, color=rgb, label=df_xs_rate['observation'][0])
        
    axs[0].legend(title='Observation ID', fancybox=True)
    axs[0].set_ylabel('$<\sigma_{XS}^2>$')
    axs[0].grid()
    axs[1].set_xlabel('Rate [ct/s]')
    axs[1].set_ylabel('$<F_{var}>$')
    axs[1].grid(True)
    plt.savefig(f'{target_dir}/Products/Plots_timeseries/xs_rate_combined.png')
    '''