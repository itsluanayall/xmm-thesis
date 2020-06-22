import os
import matplotlib.pyplot as plt
from config import CONFIG
from astropy.table import Table
import pandas as pd
import seaborn as sns
from argparse import ArgumentParser
import matplotlib.lines as mlines
from matplotlib.patches import Patch
parser = ArgumentParser()
parser.add_argument("--phoindex", action="store_true", 
                    help="make phoindex vs fvar plot")
parser.add_argument('--beta', action='store_true',
                    help='make beta vs fvar plot')
args = parser.parse_args()


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

    df_plot_zlogp = pd.DataFrame({"fvar": data_lc['F_var'].values, "fvar_err": data_lc['F_var_sigma'].values,
                                  "obsid": data_lc['ObsId'].values})
    df_plot_zlogp = df_plot_zlogp[df_plot_zlogp['fvar'] != -1.]

    df_plot_zpowe = pd.DataFrame("fvar": data_lc['F_var'].values, "fvar_err": data_lc['F_var_sigma'].values,
                                  "obsid": data_lc['ObsId'].values})
    df_plot_zpowe = df_plot_zpowe[df_plot_zpowe['fvar'] != -1.]


    if args.phoindex: #user wants phoindex vs fvar
        df_plot_zlogp = df_plot_zlogp.assign('alpha' = data_spec_zlogp['phoindex'].values, 'alpha_top' = data_spec_zlogp['phoindex_up'].values - data_spec_zlogp['phoindex'].values,
                        'alpha_bot' = data_spec_zlogp['phoindex'].values - data_spec_zlogp['phoindex_low'].values)
        df_plot_zpowe = df_plot_zpowe.assign("phoindex" = data_spec_zpowe['phoindex'].values, "phoindex_top"= data_spec_zpowe['phoindex_up'].values - data_spec_zpowe['phoindex'].values,
                                  "phoindex_bot"= data_spec_zpowe['phoindex'].values - data_spec_zpowe['phoindex_low'].values)
        
        # Distinguish between fvar high (+) and fvar low (x)
        df_plot_zlogp_high = df_plot_zlogp[df_plot_zlogp['fvar']>3*df_plot_zlogp['fvar_err']]
        df_plot_zlogp_low =  df_plot_zlogp[df_plot_zlogp['fvar']<=3*df_plot_zlogp['fvar_err']]

        df_plot_zpowe_high = df_plot_zpowe[df_plot_zpowe['fvar']>3*df_plot_zpowe['fvar_err']]
        df_plot_zpowe_low =  df_plot_zpowe[df_plot_zpowe['fvar']<=3*df_plot_zpowe['fvar_err']]
        
        # Plot
        fig, axs =plt.subplots(1, 1, figsize=(20,10), gridspec_kw={'hspace':0.4})
        axs.errorbar(df_plot_zlogp_high['alpha'].values, df_plot_zlogp_high['fvar'].values, yerr=df_plot_zlogp_high['fvar_err'].values,
                    xerr = (df_plot_zlogp_high['alpha_bot'].values, df_plot_zlogp_high['alpha_top'].values), color='b', linestyle='', marker='+', markersize=0.5, ecolor='b')
        axs.errorbar(df_plot_zlogp_low['alpha'].values, df_plot_zlogp_low['fvar'].values, yerr=df_plot_zlogp_low['fvar_err'].values,
                    xerr = (df_plot_zlogp_low['alpha_bot'].values, df_plot_zlogp_low['alpha_top'].values), color='b', linestyle='', marker='x', markersize=0.5, ecolor='b')
        
        
        axs.errorbar(df_plot_zpowe_high['phoindex'].values, df_plot_zpowe_high['fvar'].values, yerr=df_plot_zpowe_high['fvar_err'].values,
                    xerr = (df_plot_zpowe_high['phoindex_bot'].values, df_plot_zpowe_high['phoindex_top'].values), color='red', linestyle='', marker='+', markersize=0.5, ecolor='red')
        axs.errorbar(df_plot_zpowe_low['phoindex'].values, df_plot_zpowe_low['fvar'].values, yerr=df_plot_zpowe_low['fvar_err'].values,
                    xerr = (df_plot_zpowe_low['phoindex_bot'].values, df_plot_zpowe_low['phoindex_top'].values), color='red', linestyle='', marker='x', markersize=0.5, ecolor='red')

        axs.grid(True)
        axs.set_ylabel('$F_{var}$', fontsize=15)
        axs.set_xlabel('phoindex', fontsize=15)
        
        #Make legend 
        cross = mlines.Line2D([], [], color='black', marker='x', linestyle='None',
                          markersize=10, label='Low Fvar')
        plus = mlines.Line2D([], [], color='black', marker='+', linestyle='None',
                          markersize=10, label='High Fvar')
        red_patch =  Patch(facecolor='red', edgecolor='black', label='zpowerlaw')
        blue_patch =  Patch(facecolor='b', edgecolor='black', label='zlogpar')
        axs.legend(handles=[plus, cross, red_patch, blue_patch])
        
        #Linear fit?
        
        #sns.lmplot(x="alpha", y="fvar", data=df_plot_zlogp, lowess=True)
        plt.show()

    if args.beta: #user wants beta vs fvar
        df_plot_zlogp = df_plot_zlogp.assign('beta' = data_spec_zlogp['beta'].values, 'beta_top' = data_spec_zlogp['beta_up'].values - data_spec_zlogp['beta'].values,
                        'beta_bot' = data_spec_zlogp['beta'].values - data_spec_zlogp['beta_low'].values)
        figure = plt.figure(figsize=(10,10))
        plt.errorbar(df_plot_zlogp['beta'].values, df_plot_zlogp['fvar'].values, yerr=df_plot_zlogp['fvar_err'].values,
                    xerr = (df_plot_zlogp['beta_bot'].values, df_plot_zlogp['beta_top'].values), color='b', linestyle='', marker='.', markersize=0.5, ecolor='b', label='zlogpar')
        plt.grid()
        plt.xlabel('beta', fontsize=15)
        plt.ylabel('$F_{var}$', fontsize=15)
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