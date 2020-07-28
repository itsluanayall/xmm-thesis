import logging
import os
from observation import Observation
from tools import run_command, setupSAS
from config import CONFIG
from astropy.table import Table, vstack
import matplotlib.pyplot as plt
from astropy.io import ascii
import pandas as pd
import seaborn as sns
import glob
import numpy as np
import random
import sys
logging.basicConfig(level=logging.INFO)

if __name__ == "__main__":

    #Set configuration variables of the config.json file. Don't forget to change them according to your needs!
    version = CONFIG['version']
    sas_dir = CONFIG['sas_dir']
    ccf_dir = CONFIG['ccf_dir']
    mjdref = CONFIG['MJDREF']
    target_dir = CONFIG['target_dir']
    target_REDSHIFT = CONFIG['target_REDSHIFT']
    use_grace = CONFIG['use_grace']
    timescale_fvar = CONFIG['timescale_fvar']
    setupSAS(sas_dir=sas_dir, ccf_dir=ccf_dir)
    logging.info(f'Timescale chosen for fractional variability: {timescale_fvar} ks')

    epic_observations = ['0670920301', '0670920401', '0670920501', '0302180101', '0150498701', '0502030101']


    #Make final table
    EPIC_obs_table = Table(names=('ObsId', 'RevolutionId', 'ExposureID', 'Start', 
                            'End', 'Duration_Obs', 'EPIC_rate_soft', 'EPIC_erate_soft',
                            'EPIC_rate_hard', 'EPIC_erate_hard', 'hr',
                            'fvar_soft', 'efvar_soft', 'xs_soft', 'exs_soft', 
                            'nxs_soft', 'enxs_soft', 'VA_soft', 'eVA_soft',
                            'fvar_hard', 'efvar_hard', 'xs_hard', 'exs_hard', 
                            'nxs_hard', 'enxs_hard', 'VA_hard', 'eVA_hard'), 
                    dtype=('i', 'i', 'U9', 'U30',
                            'U30', 'd', 'd', 'd',
                            'd', 'd', 'd',
                            'd', 'd', 'd', 'd',
                            'd', 'd', 'd', 'd',
                            'd', 'd', 'd', 'd',
                            'd', 'd', 'd', 'd') 
                        )
    i = 0
    for obsid in epic_observations:
        obs = Observation(obsid=obsid, target_dir=target_dir)
        print('-----------------------------------------------------------')

        #Process each observation
        obs.cifbuild()
        obs.odfingest()
        obs.epproc()
        obs.filter_epic(pileup=True)
        obs.epiclccorr(pileup=True) #always correct for pile-up
        obs.epic_lightcurve()
        obs.fracvartest(instrument='epic')
        obs.pn_spectrum(pileup=True)
        obs.pn_xspec(target_REDSHIFT)
        
        #Save attributes of observation into the EPIC_table (to edit)
        EPIC_obs_table.add_row((str(obs.obsid), str(obs.revolution),  f"{obs.epic_expid}" , str(obs.starttime), 
                                str(obs.endtime), str(int(obs.duration)), obs.epicrate[0], obs.epic_erate[0], obs.mean_hr,
                                obs.epicrate[1], obs.epic_erate[1],
                                obs.fracvardict[0].get('Fractional Variability'), obs.fracvardict[0].get('Fractional Variability Error'),
                                obs.fracvardict[0].get('Excess variance'), obs.fracvardict[0].get('Excess variance error'),
                                obs.fracvardict[0].get('Normalized excess variance'), obs.fracvardict[0].get('Normalized excess variance error'),
                                obs.fracvardict[0].get('Variability Amplitude'), obs.fracvardict[0].get('Variability amplitude error'),
                                obs.fracvardict[1].get('Fractional Variability'), obs.fracvardict[1].get('Fractional Variability Error'),
                                obs.fracvardict[1].get('Excess variance'), obs.fracvardict[1].get('Excess variance error'),
                                obs.fracvardict[1].get('Normalized excess variance'), obs.fracvardict[1].get('Normalized excess variance error'),
                                obs.fracvardict[1].get('Variability Amplitude'), obs.fracvardict[1].get('Variability amplitude error')))
                            
        i+=1

    EPIC_obs_table.write(output=f'{target_dir}/Products/EPIC_Lightcurves/EPIC_obs_table.fits', format='fits', overwrite=True)
    print(EPIC_obs_table)

    # Combine all spectra tables of individual observation into a single table
    os.chdir(os.path.join(target_dir, 'Products', 'EPIC_Spectra'))
    
    #Initialize final spectral table and its path
    spectra_table = Table()
    spectra_table_path = os.path.join(target_dir, 'Products', 'EPIC_Spectra','EPIC_spectra_table.fits')
    
    for obs in epic_observations:
        #Read individual spectral table
        table = Table.read(os.path.join(target_dir, 'Products', 'EPIC_Spectra',f'{obs}_table.fits'), format='fits')

        #Append to final table
        spectra_table = vstack([spectra_table, table])
        spectra_table.write(spectra_table_path, format='fits', overwrite=True)
        print(f'Appended observation {obs}!')


