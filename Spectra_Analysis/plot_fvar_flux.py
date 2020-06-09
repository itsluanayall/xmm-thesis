import os
import matplotlib.pyplot as plt
from config import CONFIG
from astropy.table import Table


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

    #Read Flux table for logpar and powerlaw
    hdul_spec = Table.read(os.path.join(products_dir, "RGS_Spectra", "spectra_table.fits"), hdu=1)
    data_spec = hdul_spec.to_pandas()
    data_spec = data_spec[data_spec['tbinid'] == 0]   #Use only average spectra
    data_spec_clean = data_spec[data_spec['obsid'] != 658802001]

    data_spec_zlogp = data_spec_clean[data_spec_clean['model']=='constant*TBabs*zlogp']
    data_spec_zpowe = data_spec_clean[data_spec_clean['model']=='constant*TBabs*zpowe']

    #Plot panel
    fig, axs =plt.subplots(2, 1, figsize=(20,10), gridspec_kw={'hspace':0.4})
    #fig.suptitle(f'Correlation Fractional Variability and Flux', fontsize=15)

    flux_top_logpar = data_spec_zlogp['flux_up'].values - data_spec_zlogp['flux'].values
    flux_bot_logpar = data_spec_zlogp['flux'].values - data_spec_zlogp['flux_low'].values 
    flux_top_pow = data_spec_zpowe['flux_up'].values - data_spec_zpowe['flux'].values
    flux_bot_pow = data_spec_zpowe['flux'].values - data_spec_zpowe['flux_low'].values 

    axs[0].errorbar(data_spec_zlogp['flux'].values, data_lc['F_var'].values, yerr=data_lc['F_var_sigma'].values, xerr=(flux_bot_logpar, flux_top_logpar), color='b', linestyle='', marker='.', markersize=0.5, ecolor='b', label='zlogpar')
    axs[0].errorbar(data_spec_zpowe['flux'].values, data_lc['F_var'].values, yerr=data_lc['F_var_sigma'].values, xerr=(flux_bot_pow, flux_top_pow), color='red', linestyle='', marker='.', markersize=0.5, ecolor='red', label='zpowerlaw')
    axs[0].grid(True)
    axs[0].set_ylabel('$F_{var}$', fontsize=15)
    axs[0].set_xlabel('Flux [cm$^{-2}$ erg s$^{-1}$]', fontsize=15)
    axs[0].legend()

    axs[1].errorbar(data_lc['RGS_Rate'].values, data_lc['F_var'].values, yerr=data_lc['F_var_sigma'].values, color='black', linestyle='', marker='.', ecolor='gray')
    axs[1].grid(True)
    axs[1].set_ylabel('$F_{var}$', fontsize=15)
    axs[1].set_xlabel('Rate [ct/s]', fontsize=15)

    #plt.tight_layout()
    plt.savefig(os.path.join(products_dir, "Plots_spectra", "fvar_flux_correlations.png"))
