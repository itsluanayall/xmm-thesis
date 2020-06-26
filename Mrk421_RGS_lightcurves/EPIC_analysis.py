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

    for obsid in epic_observations:
        obs = Observation(obsid=obsid, target_dir=target_dir)

        #Process each observation
        obs.cifbuild()
        obs.odfingest()
        obs.emproc()
        obs.epproc()
        obs.filter_epic()