"""
Welcome to the Mrk421 Analysis Software! 

At the current version, running this script will process RGS data of the XMM-Newton observatory 
and plot the lightcurves.
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
"target_dir" is the directory of the source to analyse, e.g. 'Markarian421', that contains the directories of the observation(s);
"use_grace" is a boolean that you can set to True is you wish to save the lightcurves using the Xmgrace plotting software.

"""

import logging
import os
from observation import Observation
from tools import run_command, setupSAS
from config import CONFIG
from astropy.table import Table
import matplotlib.pyplot as plt
from astropy.io import ascii

logging.basicConfig(level=logging.INFO)

#Introduction message for the user
print(__doc__)

if __name__ == "__main__":

    #Set configuration variables of the config.json file. Don't forget to change them according to your needs!
    version = CONFIG['version']
    sas_dir = CONFIG['sas_dir']
    ccf_dir = CONFIG['ccf_dir']
    mjdref = CONFIG['MJDREF']
    target_dir = CONFIG['target_dir']
    use_grace = CONFIG['use_grace']
    logging.info(f'MRK421 Analysis - version:{version}')
    setupSAS(sas_dir=sas_dir, ccf_dir=ccf_dir)

    #Remove the .tar files from which we extracted the data
    for directory in os.listdir(target_dir):
        if directory.endswith('.tar.gz'):
            os.remove(os.path.join(target_dir, directory))
    
    #Create Products directory
    if not os.path.isdir(f'{target_dir}/Products'):
        os.makedirs(f'{target_dir}/Products')
        os.makedirs(f'{target_dir}/Products/RGS_Lightcurves')
    
    #Loop analysis for each observation

    mrk421_observation_list = []
    obs_table = Table(names=('ObsId', 'RevolutionId', 'ExposureID', 'Start', 
                            'End', 'Duration_Obs', 'RGS_Rate',
                            'Stdev_rate', 'MJD_avg_time',
                            'F_var', 'F_var_sigma', 'Excess_Variance', 
                             'xs_sigma',   'Norm_excess_variance', 'nxs_sigma', 'VA', 'VA_sigma'), 
                    dtype=('i', 'i', 'U9', 'U30',
                            'U30', 'd', 'd', 
                            'd', 'd',
                            'd', 'd', 'd',
                            'd', 'd', 'd', 'd', 'd'))

    counter = 0
    mrk421_problematic_obs = ['']
    duration_lc_ks = []
    for obsid in os.listdir(target_dir):
        
        if obsid.startswith('0'):   #All observation folders start with 0
            print('----------------')
            if obsid in mrk421_problematic_obs:
                print('This observation is tricky. Please analyse individually. Moving forward to next observation.')
                print('----------------')
                continue
            obs = Observation(obsid=obsid, target_dir=target_dir)   #instance of the observation
            mrk421_observation_list.append(obs)
            
            #Process each observation
            obs.cifbuild()
            obs.odfingest()
            obs.rgsproc()
            obs.rgslccorr()
            obs.lightcurve(mjdref=mjdref, use_grace=use_grace)
            obs.fracvartest(screen=True, timescale=9)
            #obs.bkg_lightcurve()

            #Save attributes of observable into a table
            if len(obs.rgsrate)==0:
                obs_table.add_row((str(obs.obsid), str(obs.revolution),  None , str(obs.starttime), str(obs.endtime), str(int(obs.duration)),
                                None, None, None, None,
                                None , None, None,
                                None, None, None, None))
                

            else:
    
                obs_table.add_row((str(obs.obsid), str(obs.revolution),  f"{obs.expoid[0][0]}+{obs.expoid[0][1]}" , str(obs.starttime), str(obs.endtime), str(int(obs.duration)),
                                    obs.rgsrate[0], obs.stdev[0], obs.longterm_lc_times[0], obs.fracvardict[0].get('Fractional Variability'), obs.fracvardict[0].get('Fractional Variability Error'),
                                    obs.fracvardict[0].get('Excess variance'),obs.fracvardict[0].get('Excess variance error'),
                                    obs.fracvardict[0].get('Normalized excess variance'), obs.fracvardict[0].get('Normalized excess variance error'),
                                    obs.fracvardict[0].get('Variability Amplitude'), obs.fracvardict[0].get('Variability amplitude error')))
                duration_lc_ks.append(obs.duration_lc_ks[0])

                if len(obs.rgsrate)>1:
                    for i in range(1, len(obs.rgsrate)):
                        obs_table.add_row((str(obs.obsid), str(obs.revolution),  f"{obs.expoid[i][0]}+{obs.expoid[i][1]}" , str(obs.starttime), str(obs.endtime), str(int(obs.duration)),
                                            obs.rgsrate[i], obs.stdev[i], obs.longterm_lc_times[i], obs.fracvardict[i].get('Fractional Variability'), obs.fracvardict[i].get('Fractional Variability Error'),
                                            obs.fracvardict[i].get('Excess variance'),obs.fracvardict[i].get('Excess variance error'),
                                            obs.fracvardict[i].get('Normalized excess variance'), obs.fracvardict[i].get('Normalized excess variance error'),
                                            obs.fracvardict[i].get('Variability Amplitude'), obs.fracvardict[i].get('Variability amplitude error')))
                        duration_lc_ks.append(obs.duration_lc_ks[i])

            #Keep track of number of observations that have been processed so far
            counter += 1
            logging.info(f'Processed {counter} observations!')
            print('----------------')
            
    #Show and save in csv format the table with all the attributes of the observations
    print(obs_table)
    obs_table['Duration_Obs'].unit = 's'
    obs_table['RGS_Rate'].unit = 'ct/s'
    obs_table['Stdev_rate'].unit = 'ct/s'
    obs_table['MJD_avg_time'].unit = 'd'
    obs_table['Excess_Variance'].unit = 'ct2/s2'
    obs_table['xs_sigma'].unit = 'ct2/s2'
    obs_table.write(output=f'{target_dir}/Products/obs_table.fits', format='fits', overwrite=True)

    # Duration_lc_ks distribution
    plt.hist(duration_lc_ks, bins=40)
    plt.show()
    plt.savefig(f'{target_dir}/Products/RGS_Lightcurves/distribution_expos_duration_ks.png')
    
    '''
    #For a single observation
    obs = Observation(obsid='0658802001', target_dir=target_dir)   #instance of the observation
    
    
    #Process each observation
    obs.cifbuild()
    obs.odfingest()
    obs.rgsproc()
    obs.rgslccorr()
    
    obs.lightcurve(use_grace=use_grace)
    obs.fracvartest(screen=True, netlightcurve=True)
    '''