"""
Welcome to ftest.py! This scripts allows you to perform the F-test on the results of the fits using the 
XSpec. Please append "--epic" if you want to perform the test on EPIC data.
"""

import logging
import os
import glob
import xspec
from astropy.table import Table
from config import CONFIG
import numpy as np
import sys
import matplotlib.pyplot as plt
from argparse import ArgumentParser

parser = ArgumentParser(description=__doc__)
parser.add_argument('--epic', action='store_true',
                    help='perform ftest on EPIC data')

if __name__ == "__main__":

    target_dir = CONFIG['target_dir'] 
    products_dir = os.path.join(target_dir, "Products")
    if not os.path.isdir(os.path.join(target_dir, "Products", "Plots_spectra")):
        os.makedirs(os.path.join(target_dir, "Products", "Plots_spectra"))

    #Read spectral data 
    if EPIC_ftest:
        hdul_spec = Table.read(os.path.join(products_dir, "EPIC_Spectra", "EPIC_spectra_table.fits"), hdu=1)
        data_spec = hdul_spec.to_pandas()

    else:
        hdul_spec = Table.read(os.path.join(products_dir, "RGS_Spectra", "spectra_table.fits"), hdu=1)
        data_spec = hdul_spec.to_pandas()
        data_spec = data_spec[data_spec['tbinid'] == 0]   #Use only average spectra

    #Make ftest array
    ftest_array = []
    ftest_dict = {'logpar': 0, 'powerlaw': 0}
    ftest_logpar_obs = []
    ftest_powerlaw_obs = []
    i = 0
    while i<len(data_spec):
        if 'TBabs*zlogp' in data_spec['model'].values[i]:
            chi2_lp = data_spec['chi2'].values[i]
            chi2_pw = data_spec['chi2'].values[i+1]
            dof_lp = data_spec['dof'].values[i]
            dof_pw = data_spec['dof'].values[i+1]
        elif 'TBabs*zpowe' in data_spec['model'].values[i]:
            chi2_lp = data_spec['chi2'].values[i+1]
            chi2_pw = data_spec['chi2'].values[i]
            dof_lp = data_spec['dof'].values[i+1]
            dof_pw = data_spec['dof'].values[i]

        print(f"Observation {data_spec['obsid'].values[i]}, exposures {data_spec['exposures_id'].values[i]}, logparabola vs powerlaw fstatistic" )
        print('chi2 lp:', chi2_lp, 'dof lp:', dof_lp )
        print('chi2 pw:', chi2_pw, 'dof pw:', dof_pw)

        fstatistic = xspec.Fit.ftest(chi2_lp, int(dof_lp), chi2_pw, int(dof_pw))
        if (fstatistic<0.1 and fstatistic!= -999.):
            ftest_dict['logpar'] += 1
            ftest_logpar_obs.append(data_spec['obsid'].values[i])
        else:
            ftest_dict['powerlaw'] += 1
            ftest_powerlaw_obs.append(data_spec['obsid'].values[i])

        ftest_array.append(fstatistic)
        ftest_array.append(fstatistic)
        i+=2

    #Print on screen results
    print('Final results of ftest:', ftest_dict)
    print('Observations which favour logparabola model:', ftest_logpar_obs)
    print('Observations which favour powerlaw model:', ftest_powerlaw_obs)

    data_spec['ftest'] = np.array(ftest_array)
    final_table = Table.from_pandas(data_spec)
    if EPIC_ftest:
        final_table.write(output=os.path.join(target_dir, 'Products', 'EPIC_Spectra', 'spectra_table_average.fits'), format='fits', overwrite=True)
    else:
        final_table.write(output=os.path.join(target_dir, 'Products', 'RGS_Spectra', 'spectra_table_average.fits'), format='fits', overwrite=True)

    #histogram
    data_spec_clean = data_spec.dropna()
    data_spec_LP = data_spec_clean[(data_spec_clean['ftest']<0.1) & (data_spec_clean['ftest']!=-999.)]  #threshold comes from Bhutta et al. 2018
    data_spec_PL = data_spec_clean[(data_spec_clean['ftest']>=0.1) | (data_spec_clean['ftest']==-999.)]
    
    
    x = np.arange(2)
    height = np.array([len(data_spec_LP), len(data_spec_PL)])
    print(height)
    sys.exit()
    fig = plt.figure()
    
    plt.bar(x=x, height=height)
    plt.xticks(x, ('LogParabola', 'Powerlaw'))
    plt.savefig(os.path.join(target_dir, "Products", "Plots_spectra", "ftest_barchart.png"))