import os
import matplotlib.pyplot as plt
from config import CONFIG
from astropy.table import Table
import pandas as pd
import seaborn as sns


missing_xsa_obs = [658802001]

if __name__ == "__main__":

    # Go to Plots_spectra directory
    target_dir = CONFIG['target_dir'] 
    products_dir = os.path.join(target_dir, "Products")
    if not os.path.isdir(os.path.join(target_dir, "Products", "Plots_spectra")):
        os.makedirs(os.path.join(target_dir, "Products", "Plots_spectra"))

    #Read Fvar table and rates
    hdul_lc = Table.read(os.path.join(products_dir, "RGS_Lightcurves", "obs_table.fits"), hdu=1)   
    data_lc = hdul_lc.to_pandas()
    data_lc = data_lc.dropna(subset=['RGS_Rate'])   #do not consider NaN Rates


    #Read spectra table for logpar and powerlaw parameters
    hdul_spec = Table.read(os.path.join(products_dir, "RGS_Spectra", "spectra_table.fits"), hdu=1)
    data_spec = hdul_spec.to_pandas()
    data_spec = data_spec[data_spec['tbinid'] == 0]   #Use only average spectra
    data_spec_clean = data_spec[data_spec['obsid'] != 658802001]

    data_spec_zlogp = data_spec_clean[data_spec_clean['model']=='constant*TBabs*zlogp']
    data_spec_zpowe = data_spec_clean[data_spec_clean['model']=='constant*TBabs*zpowe']

    df_plot_zlogp = pd.DataFrame({"alpha": data_spec_zlogp['phoindex'].values, "alpha_top": data_spec_zlogp['phoindex_up'].values - data_spec_zlogp['phoindex'].values,
                                  "alpha_bot": data_spec_zlogp['phoindex'].values - data_spec_zlogp['phoindex_low'].values, 
                                  "fvar": data_lc['F_var'].values, "fvar_err": data_lc['F_var_sigma'].values,
                                  "obsid": data_lc['ObsId'].values})
    df_plot_zlogp = df_plot_zlogp[df_plot_zlogp['fvar'] != -1.]

    df_plot_zpowe = pd.DataFrame({"phoindex": data_spec_zpowe['phoindex'].values, "phoindex_top": data_spec_zpowe['phoindex_up'].values - data_spec_zpowe['phoindex'].values,
                                  "phoindex_bot": data_spec_zpowe['phoindex'].values - data_spec_zpowe['phoindex_low'].values, 
                                  "fvar": data_lc['F_var'].values, "fvar_err": data_lc['F_var_sigma'].values,
                                  "obsid": data_lc['ObsId'].values})
    df_plot_zpowe = df_plot_zpowe[df_plot_zpowe['fvar'] != -1.]

    fig, axs =plt.subplots(1, 1, figsize=(20,10), gridspec_kw={'hspace':0.4})
    axs.errorbar(df_plot_zlogp['alpha'].values, df_plot_zlogp['fvar'].values, yerr=df_plot_zlogp['fvar_err'].values,
                 xerr = (df_plot_zlogp['alpha_bot'].values, df_plot_zlogp['alpha_top'].values), color='b', linestyle='', marker='.', markersize=0.5, ecolor='b', label='zlogpar')
    axs.errorbar(df_plot_zpowe['phoindex'].values, df_plot_zpowe['fvar'].values, yerr=df_plot_zpowe['fvar_err'].values,
                 xerr = (df_plot_zpowe['phoindex_bot'].values, df_plot_zpowe['phoindex_top'].values), color='red', linestyle='', marker='.', markersize=0.5, ecolor='red', label='zpowerlaw')

    axs.grid(True)
    axs.set_ylabel('$F_{var}$', fontsize=15)
    axs.set_xlabel('phoindex', fontsize=15)
    axs.legend()
    
    #Linear fit?
    
    #sns.lmplot(x="alpha", y="fvar", data=df_plot_zlogp, lowess=True)

    plt.show()
    '''
    #Plot panel
    fig, axs =plt.subplots(1, 1, figsize=(20,10), gridspec_kw={'hspace':0.4})
    #fig.suptitle(f'Correlation Fractional Variability and Flux', fontsize=15)

    phoindex_top_logpar = data_spec_zlogp['phoindex_up'].values - data_spec_zlogp['phoindex'].values
    phoindex_bot_logpar = data_spec_zlogp['phoindex'].values - data_spec_zlogp['phoindex_low'].values 
    phoindex_top_pow = data_spec_zpowe['phoindex_up'].values - data_spec_zpowe['phoindex'].values
    phoindex_bot_pow = data_spec_zpowe['phoindex'].values - data_spec_zpowe['phoindex_low'].values 

    axs.errorbar(data_spec_zlogp['phoindex'].values, data_lc['F_var'].values, yerr=data_lc['F_var_sigma'].values, xerr=(phoindex_bot_logpar, phoindex_top_logpar), color='b', linestyle='', marker='.', markersize=0.5, ecolor='b', label='zlogpar')
    axs.errorbar(data_spec_zpowe['phoindex'].values, data_lc['F_var'].values, yerr=data_lc['F_var_sigma'].values, xerr=(phoindex_bot_pow, phoindex_top_pow), color='red', linestyle='', marker='.', markersize=0.5, ecolor='red', label='zpowerlaw')
    axs.grid(True)
    axs.set_ylabel('$F_{var}$', fontsize=15)
    axs.set_xlabel('Flux [cm$^{-2}$ erg s$^{-1}$]', fontsize=15)
    axs.legend()


    #plt.tight_layout()
    plt.savefig(os.path.join(products_dir, "Plots_spectra", "fvar_phoindex_correlations.png"))
    '''