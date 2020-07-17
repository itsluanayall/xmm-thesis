import logging
import os
from observation import Observation
from tools import run_command, setupSAS
from config import CONFIG
from astropy.table import Table
import matplotlib.pyplot as plt
from astropy.io import ascii
import pandas as pd
import seaborn as sns
import glob
import numpy as np
import random
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
    epic_observations = ['0670920301']

    #Make final table
    EPIC_obs_table = Table(names=('ObsId', 'RevolutionId', 'ExposureID', 'Start', 
                            'End', 'Duration_Obs', 'EPIC_rate_soft', 'EPIC_erate_soft',
                            'EPIC_rate_hard', 'EPIC_erate_hard', 'HR',
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
    for obsid in epic_observations:
        obs = Observation(obsid=obsid, target_dir=target_dir)
        print('-----------------------------------------------------------')

        #Process each observation
        obs.cifbuild()
        obs.odfingest()
        obs.epproc()
        obs.filter_epic(pileup=True)
        obs.epiclccorr()
        obs.epic_lightcurve()
        obs.fracvartest(instrument='epic')
        obs.pn_spectrum()
        obs.pn_xspec(target_REDSHIFT)
        
        #Save attributes of observation into the EPIC_table (to edit)
        EPIC_obs_table.add_row((str(obs.obsid), str(obs.revolution),  f"{obs.epic_expid}" , str(obs.starttime), 
                                str(obs.endtime), str(int(obs.duration)), obs.epicrate[0], obs.epic_erate[0],
                                obs.epicrate[1], obs.epic_erate[1], obs.hardness_ratio,
                                obs.fracvardict[0].get('Fractional Variability'), obs.fracvardict[0].get('Fractional Variability Error'),
                                obs.fracvardict[0].get('Excess variance'), obs.fracvardict[0].get('Excess variance error'),
                                obs.fracvardict[0].get('Normalized excess variance'), obs.fracvardict[0].get('Normalized excess variance error'),
                                obs.fracvardict[0].get('Variability Amplitude'), obs.fracvardict[0].get('Variability amplitude error'),
                                obs.fracvardict[1].get('Fractional Variability'), obs.fracvardict[1].get('Fractional Variability Error'),
                                obs.fracvardict[1].get('Excess variance'), obs.fracvardict[1].get('Excess variance error'),
                                obs.fracvardict[1].get('Normalized excess variance'), obs.fracvardict[1].get('Normalized excess variance error'),
                                obs.fracvardict[1].get('Variability Amplitude'), obs.fracvardict[1].get('Variability amplitude error')))
                            
        print(EPIC_obs_table)