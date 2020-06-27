import os
import matplotlib.pyplot as plt
from config import CONFIG
from astropy.table import Table
import pandas as pd

missing_xsa_obs = [658802001]

if __name__ == "__main__":

    target_dir = CONFIG['target_dir'] 
    products_dir = os.path.join(target_dir, "Products")
    if not os.path.isdir(os.path.join(target_dir, "Products", "Plots_spectra")):
        os.makedirs(os.path.join(target_dir, "Products", "Plots_spectra"))

    #Read Fvar table and rates
    hdul_lc = Table.read(os.path.join(products_dir, "RGS_Lightcurves", "obs_table.fits"), hdu=1)   
    data_lc = hdul_lc.to_pandas()
    data_lc = data_lc.dropna(subset=['RGS_Rate'])   #do not consider NaN Rates
    print(len(data_lc))

    #Read Flux table for logpar and powerlaw
    hdul_spec = Table.read(os.path.join(products_dir, "RGS_Spectra", "spectra_table.fits"), hdu=1)
    data_spec = hdul_spec.to_pandas()
    data_spec = data_spec[data_spec['tbinid'] == 0]   #Use only average spectra
    data_spec_clean = data_spec[data_spec['obsid'] != 658802001]
    data_spec_clean = data_spec[data_spec['obsid'] != 411082701]
    data_spec_zlogp = data_spec_clean[data_spec_clean['model']=='constant*TBabs*zlogp']
    data_spec_zpowe = data_spec_clean[data_spec_clean['model']=='constant*TBabs*zpowe']
    print(len(data_spec_zlogp))

    #Plot panel
    fig, axs =plt.subplots(2, 1, figsize=(20,10), gridspec_kw={'hspace':0.4})
    #fig.suptitle(f'Correlation Fractional Variability and Flux', fontsize=15)

    # make error bars from upper and lower limits
    flux_top_logpar = data_spec_zlogp['flux_up'].values - data_spec_zlogp['flux'].values
    flux_bot_logpar = data_spec_zlogp['flux'].values - data_spec_zlogp['flux_low'].values 
    flux_top_pow = data_spec_zpowe['flux_up'].values - data_spec_zpowe['flux'].values
    flux_bot_pow = data_spec_zpowe['flux'].values - data_spec_zpowe['flux_low'].values 

    #Create dataframes for plot
    fvar_flux_logpar = pd.DataFrame({'flux': data_spec_zlogp['flux'].values, 'fvar': data_lc['F_var'].values, 'fvar_err': data_lc['F_var_sigma'].values, 'xerr_bot': flux_bot_logpar, 'xerr_top':flux_top_logpar})
    fvar_flux_pow = pd.DataFrame({'flux': data_spec_zpowe['flux'].values, 'fvar': data_lc['F_var'].values, 'fvar_err': data_lc['F_var_sigma'].values, 'xerr_bot': flux_bot_pow, 'xerr_top':flux_top_pow})
    fvar_flux_logpar = fvar_flux_logpar[fvar_flux_logpar['fvar']!=-1.]
    fvar_flux_pow = fvar_flux_pow[fvar_flux_pow['fvar']!= -1.]

    axs[0].errorbar(data=fvar_flux_logpar, x='flux', y='fvar', yerr='fvar_err', xerr=(fvar_flux_logpar['xerr_bot'].values, fvar_flux_logpar['xerr_top'].values), color='b', linestyle='', ecolor='cornflowerblue', fmt='.', markersize=3, elinewidth=1, capsize=2, capthick=1, label='zlogpar')
    axs[0].errorbar(data=fvar_flux_pow, x='flux', y='fvar', yerr='fvar_err', xerr=(fvar_flux_pow['xerr_bot'].values, fvar_flux_pow['xerr_top'].values), color='red', ecolor='lightcoral', fmt='.', markersize=3, elinewidth=1, capsize=2, capthick=1, linestyle='', label='zpowerlaw')

    axs[0].grid(True)
    axs[0].set_ylabel('$F_{var}$', fontsize=15)
    axs[0].set_xlabel('Flux [cm$^{-2}$ erg s$^{-1}$]', fontsize=15)
    axs[0].legend()

    data_lc = data_lc[data_lc['F_var']!=-1.]
    axs[1].errorbar(data_lc['RGS_Rate'].values, data_lc['F_var'].values, yerr=data_lc['F_var_sigma'].values, color='black', linestyle='', ecolor='gray', fmt='.', markersize=3, elinewidth=1, capsize=2, capthick=1)
    axs[1].grid(True)
    axs[1].set_ylabel('$F_{var}$', fontsize=15)
    axs[1].set_xlabel('Rate [ct/s]', fontsize=15)

    #plt.tight_layout()
    plt.savefig(os.path.join(products_dir, "Plots_spectra", "fvar_flux_correlations.png"))
