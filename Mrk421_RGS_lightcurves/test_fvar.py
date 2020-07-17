import os
from tools import *
from config import CONFIG
from astropy.table import Table
import matplotlib.pyplot as plt
from astropy.io import ascii, fits
import pandas as pd
import glob
import numpy as np

target_dir = CONFIG['target_dir']
#os.chdir(os.path.join(target_dir, 'Products', 'RGS_Lightcurves'))
mrk421_problematic_obs = ['0658802001', '0411082701']
obs_list = []
exposure_list = []
mean_error_list = []     
fvar_list = []  
efvar_list = []
for obsid in os.listdir(target_dir):
        
    if obsid.startswith('0'):   #All observation folders start with 0
        if obsid in mrk421_problematic_obs:
            print(f'This observation ({obsid}) is tricky. Please analyse individually. Moving forward to next observation.')
            continue
        
        os.chdir(os.path.join(target_dir, obsid, 'rgs'))
        rate_files = glob.glob("*RGS_rates.ds")

        for rate_file in rate_files:

            with fits.open(rate_file) as hdul:
                obs_list.append(obsid)
                exposure_list.append(rate_file[11:18])
                mean_error = hdul['RATE'].data['ERROR'][~np.isnan(hdul['RATE'].data['ERROR'])].mean()
                mean_error_list.append(mean_error)
                
                fvar, fvar_err = fractional_variability(hdul['RATE'].data['RATE'], hdul['RATE'].data['ERROR'], hdul['RATE'].data['BACKV'], hdul['RATE'].data['BACKE'], netlightcurve=True)
                fvar_list.append(fvar)
                efvar_list.append(fvar_err)
df_final = pd.DataFrame({'obsid': obs_list, 'exposures': exposure_list, 'mean_rate_error': mean_error_list, 'fvar': fvar_list, 'efvar': efvar_list})
df_error_sorted = df_final.sort_values(by=['mean_rate_error'])
df_final =df_final.reset_index()
print(df_error_sorted.to_string())


high_error_fvar = [df_final[df_final['obsid']=='0510610101']['fvar'].values[0], df_final[df_final['obsid']=='0162960101']['fvar'].values[0],
                    df_final[df_final['obsid']=='0670920301']['fvar'].values[0], df_final[df_final['obsid']=='0656380801']['fvar'].values[0],
                    df_final[df_final['obsid']=='0136540601']['fvar'].values[0], df_final[df_final['obsid']=='0510610101']['fvar'].values[0] ]

low_error_fvar = [df_final[df_final['obsid']=='0658800601']['fvar'].values[0], df_final[df_final['obsid']=='0153950801']['fvar'].values[0],
                    df_final[df_final['obsid']=='0658800301']['fvar'].values[0], df_final[df_final['obsid']=='0658800101']['fvar'].values[0],
                    df_final[df_final['obsid']=='0153950601']['fvar'].values[0], df_final[df_final['obsid']=='0791780101']['fvar'].values[0]]


high_error_efvar =  [df_final[df_final['obsid']=='0510610101']['efvar'].values[0], df_final[df_final['obsid']=='0162960101']['efvar'].values[0],
                    df_final[df_final['obsid']=='0670920301']['efvar'].values[0], df_final[df_final['obsid']=='0656380801']['efvar'].values[0],
                    df_final[df_final['obsid']=='0136540601']['efvar'].values[0], df_final[df_final['obsid']=='0510610101']['efvar'].values[0] ]

low_error_efvar = [df_final[df_final['obsid']=='0658800601']['efvar'].values[0], df_final[df_final['obsid']=='0153950801']['efvar'].values[0],
                    df_final[df_final['obsid']=='0658800301']['efvar'].values[0], df_final[df_final['obsid']=='0658800101']['efvar'].values[0],
                    df_final[df_final['obsid']=='0153950601']['efvar'].values[0], df_final[df_final['obsid']=='0791780101']['efvar'].values[0]]



x=[1,2,3,4,5, 6]
fig, ax = plt.subplots(2,1, figsize=(8,11), gridspec_kw={'wspace':7})
ax[0].errorbar(y=high_error_fvar, x=x, yerr=high_error_efvar, label='High error on lightcurve', linestyle='', fmt='.')
ax[0].errorbar(y=low_error_fvar, x=x, yerr=low_error_efvar , label='Low error on lightcurve',  linestyle='', fmt='.')
xticks_labels = ['0510610101 \n vs\n0658800601', '0162960101\n vs \n0153950801', '0670920301\n vs \n 0658800301', '0656380801\n vs \n 0658800101', '0136540601 \n vs \n 0153950601', '0510610101 \n vs\n  0791780101']
ax[0].set_xticks(x)
ax[0].set_xticklabels(xticks_labels, rotation=60, fontsize=7)
ax[0].legend()
ax[0].xaxis.grid() # vertical lines
ax[0].set_ylabel('Fvar')


variable_lc_fvar = [df_final[df_final['obsid']=='0153950601']['fvar'].values[0], df_final[df_final['obsid']=='0658800101']['fvar'].values[0],
                    df_final[df_final['obsid']=='0136540601']['fvar'].values[0], df_final[df_final['obsid']=='0658801501']['fvar'].values[0]]
flat_lc_fvar = [df_final[df_final['obsid']=='0658800601']['fvar'].values[0], df_final[df_final['obsid']=='0791780101']['fvar'].values[0],
                    df_final[df_final['obsid']=='0510610101']['fvar'].values[0], df_final[df_final['obsid']=='0153951301']['fvar'].values[0]]

variable_lc_efvar = [df_final[df_final['obsid']=='0153950601']['efvar'].values[0], df_final[df_final['obsid']=='0658800101']['efvar'].values[0],
                    df_final[df_final['obsid']=='0136540601']['efvar'].values[0], df_final[df_final['obsid']=='0658801501']['efvar'].values[0]]
flat_lc_efvar = [df_final[df_final['obsid']=='0658800601']['efvar'].values[0], df_final[df_final['obsid']=='0791780101']['efvar'].values[0],
                    df_final[df_final['obsid']=='0510610101']['efvar'].values[0], df_final[df_final['obsid']=='0153951301']['efvar'].values[0]]



x1 = [1,2,3,4]

ax[1].errorbar(y=variable_lc_fvar, x=x1, yerr=variable_lc_efvar, label='Variable lightcurve', linestyle='', fmt='.')
ax[1].errorbar(y=flat_lc_fvar, x=x1, yerr=flat_lc_efvar, label='Flat lightcurve', linestyle='', fmt='.')
ax[1].set_xticks(x1)
xticks_labels1 = ['0153950601 \n vs\n 0658800601', '0658800101\n vs \n 0791780101', '0136540601\n vs \n 0510610101', '0658801501\n vs \n 0153951301']
ax[1].set_xticklabels(xticks_labels1, rotation=60, fontsize=7)
ax[1].legend()
ax[1].set_ylabel('Fvar')
ax[1].xaxis.grid() # vertical lines

plt.tight_layout()
plt.savefig(os.path.join(target_dir, "Products", "test_fvar"))
