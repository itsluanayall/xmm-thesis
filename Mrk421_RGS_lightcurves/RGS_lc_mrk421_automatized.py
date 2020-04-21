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
from astropy.io import ascii

logging.basicConfig(level=logging.INFO)

#Introduction message for the user
print(__doc__)

if __name__ == "__main__":

    #Set configuration variables of the config.json file. Don't forget to change them according to your needs!
    version = CONFIG['version']
    sas_dir = CONFIG['sas_dir']
    ccf_dir = CONFIG['ccf_dir']
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
    #nobs = len(os.listdir(target_dir))
    mrk421_observation_list = []
    obs_table = Table(names=('ObsId', 'RevolutionId', 'Start', 'End', 'Duration[s]', 'RGS_Rate[count/s]', 'Discarded_Exposures', 'Fractional_Variability', 'Excess_Variance', 'Norm_excess_variance'), 
                    dtype=('object', 'object', 'object', 'object', 'object', 'object', 'object', 'object', 'object', 'object'))
    counter = 0

    for obsid in os.listdir(target_dir):
        
        if obsid.startswith('0'):   #All observation folders start with 0
            print('----------------')
            obs = Observation(obsid=obsid, target_dir=target_dir)   #instance of the observation
            mrk421_observation_list.append(obs)
            
            #Process each observation
            obs.cifbuild()
            obs.odfingest()
            obs.rgsproc()
            obs.rgslccorr()
            obs.lightcurve(use_grace=use_grace)
            obs.fracvartest(screen=True)
            #obs.bkg_lightcurve()

            #Save attributes of observable into a table
            try:
                rgsrates = ""
                fracvar_str = ""
                xs_str = ""
                nxs_str = ""
                for el in obs.rgsrate:
                    rgsrates += f"{el},"
                rgsrates = rgsrates[:-1]

                for el in obs.fracvardict:
                    str_var = f"{el.get('Fractional Variability'):.4f} +- {el.get('Fractional Variability Error'):.4f} "
                    str_xs = f"{el.get('Excess variance'):.4f}"
                    str_nxs = f"{el.get('Normalized excess variance'):.4f} +- {el.get('Normalized excess variance error'):.4f}"

                    fracvar_str += f"{str_var},"
                    xs_str += f"{str_xs},"
                    nxs_str += f"{str_nxs},"

                fracvar_str = fracvar_str[:-1]
                xs_str = xs_str[:-1]
                nxs_str = nxs_str[:-1]
                discarded_expos_str = ', '.join(obs.discarded_expos)

            except TypeError as e:
                logging.error(e)
                fracvar_str = 'None'
                discarded_expos_str = 'None'

            #Write observation info into w csv table row
            obs_table.add_row((str(obs.obsid), str(obs.revolution), str(obs.starttime), str(obs.endtime), str(int(obs.duration)),
                             rgsrates,  discarded_expos_str, fracvar_str, xs_str, nxs_str ))

            #Keep track of number of observations that have been processed so far
            counter += 1
            logging.info(f'Processed {counter} observations!')
            print('----------------')
            
    #Show and save in csv format the table with all the attributes of the observations
    print(obs_table)
    ascii.write(table=obs_table, output=f'{target_dir}/Products/obs_table.csv', format='csv', overwrite=True)
            
    '''
    #For a single observation
    obs = Observation(obsid='0162960101', target_dir=target_dir)   #instance of the observation
    
    
    #Process each observation
    obs.cifbuild()
    obs.odfingest()
    obs.rgsproc()
    obs.rgslccorr()
    
    obs.lightcurve(use_grace=use_grace)
    obs.fracvartest(screen=True, netlightcurve=True)
    '''