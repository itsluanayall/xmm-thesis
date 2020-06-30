"""
Welcome to plot_parameters.py! This scripts allows you to make some interesting plot regarding the spectral
parameters (phoindex, beta...) of blazars. Possible options are:

 --fvar --phoindex (--beta) --average, 
 --rate --phoindex (--beta) --bins, 
 --phoindex --beta, 
 --panel --bins --logpar (--powerlaw) 
"""

import os
import glob
import matplotlib.pyplot as plt
from config import CONFIG
from astropy.table import Table
import pandas as pd
import seaborn as sns
from argparse import ArgumentParser
import matplotlib.lines as mlines
from matplotlib.patches import Patch
from tools import *

parser = ArgumentParser(description=__doc__)
parser.add_argument("--phoindex", action="store_true", 
                    help="make phoindex  plot")
parser.add_argument('--beta', action='store_true',
                    help='make beta plot')
parser.add_argument('--fvar', action='store_true',
                    help='make fvar plot')
parser.add_argument('--rate', action='store_true',
                    help='make rate plot')
parser.add_argument('--flux', action='store_true',
                    help='make flux plot')
parser.add_argument('--lumin', action='store_true',
                    help='make lumionosity plot')
parser.add_argument('--distribution', action='store_true',
                    help='plot the distributions of low and high state parameters')
parser.add_argument('--average', action='store_true',
                    help='use only average spectra')
parser.add_argument('--bins', action='store_true',
                    help='use spectra of 1000s each')
parser.add_argument('--panel', action='store_true',
                    help='make panel with lightcurve and parameters')
parser.add_argument('--logpar', action='store_true',
                    help='use logparabola model to make the panel')
parser.add_argument('--powerlaw', action='store_true',
                    help='use powerlaw model to make the panel')
parser.add_argument('--state', action='store_true',
                    help='make plot differentiating between low and high state')
args = parser.parse_args()


missing_xsa_obs = [658802001]

if __name__ == "__main__":

    # Go to Plots_spectra directory
    target_dir = CONFIG['target_dir'] 
    MJDREF = 50814.0
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

    if args.average:
        data_spec = data_spec[data_spec['tbinid'] == 0]   #Use only average spectra
    elif args.bins:
        data_spec = data_spec[data_spec['tbinid'] != 0]   #Use divided bins

    #Clean data
    data_spec_clean = data_spec[data_spec['obsid'] != 658802001]

    # Make dataframe for the two different models
    data_spec_zlogp = data_spec_clean[data_spec_clean['model']=='constant*TBabs*zlogp']
    data_spec_zpowe = data_spec_clean[data_spec_clean['model']=='constant*TBabs*zpowe']


    if args.fvar:

        if args.phoindex: #user wants phoindex vs fvar
            df_plot_zlogp = pd.DataFrame({"fvar": data_lc['F_var'].values, "fvar_err": data_lc['F_var_sigma'].values, 'alpha': data_spec_zlogp['phoindex'].values, 
                                          'alpha_top': data_spec_zlogp['phoindex_up'].values - data_spec_zlogp['phoindex'].values,
                                          'alpha_bot': data_spec_zlogp['phoindex'].values - data_spec_zlogp['phoindex_low'].values,
                                          "obsid": data_lc['ObsId'].values})
            df_plot_zlogp = df_plot_zlogp[df_plot_zlogp['fvar'] != -1.]
            
            
            df_plot_zpowe = pd.DataFrame({"fvar": data_lc['F_var'].values, "fvar_err": data_lc['F_var_sigma'].values,
                                         "obsid": data_lc['ObsId'].values, "phoindex": data_spec_zpowe['phoindex'].values,
                                         "phoindex_top": data_spec_zpowe['phoindex_up'].values - data_spec_zpowe['phoindex'].values,
                                         "phoindex_bot": data_spec_zpowe['phoindex'].values - data_spec_zpowe['phoindex_low'].values})
            df_plot_zpowe = df_plot_zpowe[df_plot_zpowe['fvar'] != -1.]

            # Distinguish between fvar high (+) and fvar low (x)
            df_plot_zlogp_high = df_plot_zlogp[df_plot_zlogp['fvar']>3*df_plot_zlogp['fvar_err']]
            df_plot_zlogp_low =  df_plot_zlogp[df_plot_zlogp['fvar']<=3*df_plot_zlogp['fvar_err']]

            df_plot_zpowe_high = df_plot_zpowe[df_plot_zpowe['fvar']>3*df_plot_zpowe['fvar_err']]
            df_plot_zpowe_low =  df_plot_zpowe[df_plot_zpowe['fvar']<=3*df_plot_zpowe['fvar_err']]
            
            # Plot
            fig, axs =plt.subplots(1, 1, figsize=(5,5), gridspec_kw={'hspace':0.4})
            axs.errorbar(df_plot_zlogp_high['alpha'].values, df_plot_zlogp_high['fvar'].values, yerr=df_plot_zlogp_high['fvar_err'].values,
                        xerr = (df_plot_zlogp_high['alpha_bot'].values, df_plot_zlogp_high['alpha_top'].values), color='b', linestyle='', marker='+', markersize=10, ecolor='cornflowerblue', elinewidth=1, capsize=2, capthick=1)
            axs.errorbar(df_plot_zlogp_low['alpha'].values, df_plot_zlogp_low['fvar'].values, yerr=df_plot_zlogp_low['fvar_err'].values,
                        xerr = (df_plot_zlogp_low['alpha_bot'].values, df_plot_zlogp_low['alpha_top'].values), color='b', linestyle='', marker='x', markersize=10, ecolor='cornflowerblue', elinewidth=1, capsize=2, capthick=1)
            
            
            axs.errorbar(df_plot_zpowe_high['phoindex'].values, df_plot_zpowe_high['fvar'].values, yerr=df_plot_zpowe_high['fvar_err'].values,
                        xerr = (df_plot_zpowe_high['phoindex_bot'].values, df_plot_zpowe_high['phoindex_top'].values), color='red', linestyle='', marker='+', markersize=10, ecolor='lightcoral', elinewidth=1, capsize=2, capthick=1)
            axs.errorbar(df_plot_zpowe_low['phoindex'].values, df_plot_zpowe_low['fvar'].values, yerr=df_plot_zpowe_low['fvar_err'].values,
                        xerr = (df_plot_zpowe_low['phoindex_bot'].values, df_plot_zpowe_low['phoindex_top'].values), color='red', linestyle='', marker='x', markersize=10, ecolor='lightcoral', elinewidth=1, capsize=2, capthick=1)

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
            plt.savefig(os.path.join(target_dir, "Products", "Plots_spectra", "fvar_phoindex_correlations.png"))

        if args.beta: #user wants beta vs fvar
            df_plot_zlogp = pd.DataFrame({"fvar": data_lc['F_var'].values, "fvar_err": data_lc['F_var_sigma'].values, 'beta': data_spec_zlogp['beta'].values, 
                                          'beta_top': data_spec_zlogp['beta_up'].values - data_spec_zlogp['beta'].values,
                                          'beta_bot': data_spec_zlogp['beta'].values - data_spec_zlogp['beta_low'].values,
                                          "obsid": data_lc['ObsId'].values})
            df_plot_zlogp = df_plot_zlogp[df_plot_zlogp['fvar'] != -1.]

            figure = plt.figure(figsize=(5,5))
            plt.errorbar(df_plot_zlogp['beta'].values, df_plot_zlogp['fvar'].values, yerr=df_plot_zlogp['fvar_err'].values,
                        xerr = (df_plot_zlogp['beta_bot'].values, df_plot_zlogp['beta_top'].values), color='b', linestyle='', marker='.', markersize=5, ecolor='cornflowerblue', elinewidth=1, capsize=2, capthick=1, label='zlogpar')
            plt.grid()
            plt.xlabel('beta', fontsize=15)
            plt.ylabel('$F_{var}$', fontsize=15)
            plt.savefig(os.path.join(target_dir, "Products", "Plots_spectra", "fvar_beta_correlations.png"))


    if args.rate:

        if args.phoindex: #user wants phoindex vs rate
            data_spec_zlogp = data_spec_zlogp[data_spec_zlogp['phoindex_up']!=0.]
            df_plot_zlogp = pd.DataFrame({"rate": data_spec_zlogp['rate'].values, "erate": data_spec_zlogp['erate'].values, 'alpha': data_spec_zlogp['phoindex'].values, 
                                          'alpha_top': data_spec_zlogp['phoindex_up'].values - data_spec_zlogp['phoindex'].values,
                                          'alpha_bot': data_spec_zlogp['phoindex'].values - data_spec_zlogp['phoindex_low'].values,
                                          "obsid": data_spec_zlogp['obsid'].values})
            
            data_spec_zpowe = data_spec_zpowe[data_spec_zpowe['phoindex_up']!=0.]
            df_plot_zpowe = pd.DataFrame({"rate": data_spec_zpowe['rate'].values, "erate": data_spec_zpowe['erate'].values,
                                          "phoindex": data_spec_zpowe['phoindex'].values,
                                         "phoindex_top": data_spec_zpowe['phoindex_up'].values - data_spec_zpowe['phoindex'].values,
                                         "phoindex_bot": data_spec_zpowe['phoindex'].values - data_spec_zpowe['phoindex_low'].values,
                                         "obsid": data_spec_zpowe['obsid'].values})
            
            figure = plt.figure(figsize=(10,5))
            plt.errorbar(y=df_plot_zlogp['alpha'].values, x=df_plot_zlogp['rate'].values, xerr=df_plot_zlogp['erate'].values,
                        yerr = (df_plot_zlogp['alpha_bot'].values, df_plot_zlogp['alpha_top'].values), fmt='.', markersize=5, ecolor='cornflowerblue', elinewidth=1, capsize=2, capthick=1, color='b', label='zlogpar')
            plt.errorbar(y=df_plot_zpowe['phoindex'].values, x=df_plot_zpowe['rate'].values, xerr=df_plot_zpowe['erate'].values,
                        yerr = (df_plot_zpowe['phoindex_bot'].values, df_plot_zpowe['phoindex_top'].values), fmt='.', color='red', markersize=5, elinewidth=1, capsize=2, capthick=1, ecolor='lightcoral', label='zpowerlaw')

            plt.grid()
            plt.xlabel('Rate [ct/s]', fontsize=15)
            plt.ylabel('phoindex', fontsize=15)
            plt.legend()
            plt.savefig(os.path.join(target_dir, "Products", "Plots_spectra", "phoindex_vs_rate.png"))


        if args.beta:
            df_plot_zlogp = pd.DataFrame({"rate": data_spec_zlogp['rate'].values, "erate": data_spec_zlogp['erate'].values, 'beta': data_spec_zlogp['beta'].values, 
                                          'beta_top': data_spec_zlogp['beta_up'].values - data_spec_zlogp['beta'].values,
                                          'beta_bot': data_spec_zlogp['beta'].values - data_spec_zlogp['beta_low'].values,
                                          "obsid": data_spec_zlogp['obsid'].values})
            
            figure = plt.figure(figsize=(10,5))
            plt.errorbar(y=df_plot_zlogp['beta'].values, x=df_plot_zlogp['rate'].values, xerr=df_plot_zlogp['erate'].values,
                        yerr = (df_plot_zlogp['beta_bot'].values, df_plot_zlogp['beta_top'].values), fmt='.', markersize=5, ecolor='cornflowerblue', elinewidth=1, capsize=2, capthick=1, color='b', label='zlogpar')
            plt.grid()
            plt.xlabel('Rate [ct/s]', fontsize=15)
            plt.ylabel('beta', fontsize=15)
            plt.savefig(os.path.join(target_dir, "Products", "Plots_spectra", "beta_vs_rate.png"))
    

    if args.flux:

        if args.phoindex: #user wants phoindex vs flux

            #Logpar dataframe
            data_spec_zlogp = data_spec_zlogp[data_spec_zlogp['flux_low']!=0.]
            data_spec_zlogp = data_spec_zlogp[data_spec_zlogp['phoindex_up']!=0.]
            df_plot_zlogp = pd.DataFrame({"rate": data_spec_zlogp['rate'].values, "erate": data_spec_zlogp['erate'].values, 
                                          "flux": data_spec_zlogp['flux'].values, "flux_top": data_spec_zlogp['flux_up'].values - data_spec_zlogp['flux'],
                                          'flux_bot': data_spec_zlogp['flux'].values - data_spec_zlogp['flux_low'],
                                          "flux_up": data_spec_zlogp['flux_up'].values, "flux_low":  data_spec_zlogp['flux_low'].values,
                                          'alpha': data_spec_zlogp['phoindex'].values, 
                                          'alpha_top': data_spec_zlogp['phoindex_up'].values - data_spec_zlogp['phoindex'].values,
                                          'alpha_bot': data_spec_zlogp['phoindex'].values - data_spec_zlogp['phoindex_low'].values,
                                          "obsid": data_spec_zlogp['obsid'].values})
            
            #Powerlaw dataframe
            data_spec_zpowe = data_spec_zpowe[data_spec_zpowe['phoindex_up']!=0.]
            data_spec_zpowe = data_spec_zpowe[data_spec_zpowe['flux_low']!=0]
            df_plot_zpowe = pd.DataFrame({"rate": data_spec_zpowe['rate'].values, "erate": data_spec_zpowe['erate'].values,
                                          "flux": data_spec_zpowe['flux'].values, "flux_top": data_spec_zpowe['flux_up'].values - data_spec_zpowe['flux'],
                                          "flux_bot": data_spec_zpowe['flux'].values - data_spec_zpowe['flux_low'],
                                          "flux_up": data_spec_zpowe['flux_up'].values, "flux_low": data_spec_zpowe['flux_low'].values,
                                          "phoindex": data_spec_zpowe['phoindex'].values,
                                          "phoindex_top": data_spec_zpowe['phoindex_up'].values - data_spec_zpowe['phoindex'].values,
                                          "phoindex_bot": data_spec_zpowe['phoindex'].values - data_spec_zpowe['phoindex_low'].values,
                                          "obsid": data_spec_zpowe['obsid'].values})

            if args.state:  #separate in high and low state
                df_plot_zlogp_high = df_plot_zlogp[df_plot_zlogp['rate']>=20]
                df_plot_zlogp_low = df_plot_zlogp[df_plot_zlogp['rate']<20]
                df_plot_zpowe_high = df_plot_zpowe[df_plot_zpowe['rate']>=20]
                df_plot_zpowe_low = df_plot_zpowe[df_plot_zpowe['rate']<20]               

            #Plot phoindex vs flux for zlogpar and powerlaw
            if args.state:
                figure, axs = plt.subplots(2,1, figsize=(10,5), sharex=True, gridspec_kw={'hspace':0.3})
                axs[0].errorbar(y=df_plot_zlogp_high['alpha'].values, x=df_plot_zlogp_high['flux'].values, xerr=(df_plot_zlogp_high['flux_bot'].values, df_plot_zlogp_high['flux_top'].values),
                            yerr = (df_plot_zlogp_high['alpha_bot'].values, df_plot_zlogp_high['alpha_top'].values), fmt='.', markersize=3, ecolor='peachpuff', elinewidth=1, capsize=2, capthick=1, color='orange', label='High state')
                axs[0].errorbar(y=df_plot_zlogp_low['alpha'].values, x=df_plot_zlogp_low['flux'].values, xerr=(df_plot_zlogp_low['flux_bot'].values, df_plot_zlogp_low['flux_top'].values),
                            yerr = (df_plot_zlogp_low['alpha_bot'].values, df_plot_zlogp_low['alpha_top'].values),fmt='.', markersize=3, ecolor='mediumseagreen', elinewidth=1, capsize=2, capthick=1, color='g', label='Low state')
                axs[1].errorbar(y=df_plot_zpowe_high['phoindex'].values, x=df_plot_zpowe_high['flux'].values, xerr=(df_plot_zpowe_high['flux_bot'].values, df_plot_zpowe_high['flux_top'].values),
                            yerr = (df_plot_zpowe_high['phoindex_bot'].values, df_plot_zpowe_high['phoindex_top'].values), fmt='.', color='orange', markersize=3, elinewidth=1, capsize=2, capthick=1, ecolor='peachpuff', label='High state')
                axs[1].errorbar(y=df_plot_zpowe_low['phoindex'].values, x=df_plot_zpowe_low['flux'].values, xerr=(df_plot_zpowe_low['flux_bot'].values, df_plot_zpowe_low['flux_top'].values),
                            yerr = (df_plot_zpowe_low['phoindex_bot'].values, df_plot_zpowe_low['phoindex_top'].values), fmt='.', color='g', markersize=3, elinewidth=1, capsize=2, capthick=1, ecolor='mediumseagreen', label='Low state')
                
                axs[1].set_xlabel('Flux [cm$^{-2}$ erg s$^{-1}$]', fontsize=15)
                axs[1].set_ylabel('phoindex', fontsize=15)
                axs[0].set_ylabel('alpha', fontsize=15)
                axs[0].set_title('Logparabola')
                axs[1].set_title('Powerlaw')
                axs[0].legend(loc="upper right")
                plt.savefig(os.path.join(target_dir, "Products", "Plots_spectra", "phoindex_vs_flux_state.png"))

            else:  #plot without distinguishing betweeen low and high state
                figure = plt.figure(figsize=(10,5))
                plt.errorbar(y=df_plot_zlogp['alpha'].values, x=df_plot_zlogp['flux'].values, xerr=(df_plot_zlogp['flux_bot'].values, df_plot_zlogp['flux_top'].values),
                            yerr = (df_plot_zlogp['alpha_bot'].values, df_plot_zlogp['alpha_top'].values), fmt='.', markersize=5, ecolor='cornflowerblue', elinewidth=1, capsize=2, capthick=1, color='b', label='zlogpar')
                plt.errorbar(y=df_plot_zpowe['phoindex'].values, x=df_plot_zpowe['flux'].values, xerr=(df_plot_zpowe['flux_bot'].values, df_plot_zpowe['flux_top'].values),
                            yerr = (df_plot_zpowe['phoindex_bot'].values, df_plot_zpowe['phoindex_top'].values), fmt='.', color='red', markersize=5, elinewidth=1, capsize=2, capthick=1, ecolor='lightcoral', label='zpowerlaw')
                plt.grid()
                plt.xlabel('Flux [cm$^{-2}$ erg s$^{-1}$]', fontsize=15)
                plt.ylabel('phoindex', fontsize=15)
                plt.legend()
                plt.savefig(os.path.join(target_dir, "Products", "Plots_spectra", "phoindex_vs_flux.png"))


        if args.beta:
            data_spec_zlogp = data_spec_zlogp[data_spec_zlogp['flux_low']!=0.]
            df_plot_zlogp = pd.DataFrame({"rate": data_spec_zlogp['rate'].values, "erate": data_spec_zlogp['erate'].values,
                                          "flux": data_spec_zlogp['flux'].values, "flux_top": data_spec_zlogp['flux_up'].values - data_spec_zlogp['flux'].values, 
                                          'flux_bot': data_spec_zlogp['flux'].values - data_spec_zlogp['flux_low'].values,
                                          'beta': data_spec_zlogp['beta'].values, 
                                          'beta_top': data_spec_zlogp['beta_up'].values - data_spec_zlogp['beta'].values,
                                          'beta_bot': data_spec_zlogp['beta'].values - data_spec_zlogp['beta_low'].values,
                                          "obsid": data_spec_zlogp['obsid'].values})
            if args.state:  #separate in high and low state
                df_plot_zlogp_high = df_plot_zlogp[df_plot_zlogp['rate']>=20]
                df_plot_zlogp_low = df_plot_zlogp[df_plot_zlogp['rate']<20]           

            figure = plt.figure(figsize=(10,5))
            if args.state:
                
                plt.errorbar(y=df_plot_zlogp_high['beta'].values, x=df_plot_zlogp_high['flux'].values, xerr=(df_plot_zlogp_high['flux_bot'].values, df_plot_zlogp_high['flux_top'].values),
                            yerr = (df_plot_zlogp_high['beta_bot'].values, df_plot_zlogp_high['beta_top'].values), fmt='.', markersize=3, ecolor='peachpuff', elinewidth=1, capsize=2, capthick=1, color='orange', label='High state')
                plt.errorbar(y=df_plot_zlogp_low['beta'].values, x=df_plot_zlogp_low['flux'].values, xerr=(df_plot_zlogp_low['flux_bot'].values, df_plot_zlogp_low['flux_top'].values),
                            yerr = (df_plot_zlogp_low['beta_bot'].values, df_plot_zlogp_low['beta_top'].values),fmt='.', markersize=3, ecolor='mediumseagreen', elinewidth=1, capsize=2, capthick=1, color='g', label='Low state')
                plt.legend()

            else:
                plt.errorbar(y=df_plot_zlogp['beta'].values, x=df_plot_zlogp['flux'].values, xerr=(df_plot_zlogp['flux_bot'].values, df_plot_zlogp['flux_top'].values),
                            yerr = (df_plot_zlogp['beta_bot'].values, df_plot_zlogp['beta_top'].values), fmt='.', markersize=5, ecolor='cornflowerblue', elinewidth=1, capsize=2, capthick=1, color='b', label='zlogpar')
            plt.grid()
            plt.xlabel('Flux [cm$^{-2}$ erg s$^{-1}$]', fontsize=15)
            plt.ylabel('beta', fontsize=15)

            #Save figure
            if args.state:
                plt.savefig(os.path.join(target_dir, "Products", "Plots_spectra", "beta_vs_flux_state.png"))
            else:
                plt.savefig(os.path.join(target_dir, "Products", "Plots_spectra", "beta_vs_flux.png"))
    

    if not args.fvar and not args.rate:
        if args.phoindex and args.beta:
            data_spec_zlogp = data_spec_zlogp[data_spec_zlogp['phoindex_up']!=0.]
            df_plot_zlogp = pd.DataFrame({'alpha': data_spec_zlogp['phoindex'].values, 
                                            'alpha_top': data_spec_zlogp['phoindex_up'].values - data_spec_zlogp['phoindex'].values,
                                            'alpha_bot': data_spec_zlogp['phoindex'].values - data_spec_zlogp['phoindex_low'].values,
                                            'beta': data_spec_zlogp['beta'].values, 
                                            'beta_top': data_spec_zlogp['beta_up'].values - data_spec_zlogp['beta'].values,
                                            'beta_bot': data_spec_zlogp['beta'].values - data_spec_zlogp['beta_low'].values,
                                            "obsid": data_spec_zlogp['obsid'].values})
            
            data_spec_zpowe = data_spec_zpowe[data_spec_zpowe['phoindex_up']!=0.]

            figure = plt.figure(figsize=(10,5))
            plt.errorbar(y=df_plot_zlogp['alpha'].values, x=df_plot_zlogp['beta'].values, yerr=(df_plot_zlogp['alpha_bot'].values, df_plot_zlogp['alpha_top'].values),
                        xerr = (df_plot_zlogp['beta_bot'].values, df_plot_zlogp['beta_top'].values), fmt='.', markersize=5, ecolor='gray', elinewidth=1, capsize=2, capthick=1, color='black')
            plt.grid()
            plt.xlabel('beta', fontsize=15)
            plt.ylabel('alpha', fontsize=15)
            plt.savefig(os.path.join(target_dir, "Products", "Plots_spectra", "alpha_vs_beta.png"))


    if args.panel and args.bins:

        #TOTAL LIGHT CURVE
        os.chdir(target_dir)
        total_lightcurve_rates = []
        total_lightcurve_errates = []
        total_lightcurve_times = []
        observations = []
        for directory in os.listdir(target_dir):
            if directory.startswith('0'):
                os.chdir(f"{target_dir}/{directory}/rgs")
                
                for filename in glob.glob('*_RGS_rates.ds'):
                    x, y, yerr, fracexp, y_bg, yerr_bg = mask_fracexp15(filename)
                    total_lightcurve_rates.extend(y[:-1])
                    total_lightcurve_errates.extend(yerr[:-1])
                    total_lightcurve_times.extend(x[:-1])
                    observations.extend([int(directory) for i in range(len(x)-1)])
            
        total_lightcurve_rates = np.asarray(total_lightcurve_rates)
        total_lightcurve_errates = np.asarray(total_lightcurve_errates)
        total_lightcurve_times = np.asarray(total_lightcurve_times)
        observations = np.asarray(observations)
      
        # Conversion of times (from MET to MJD)
        total_lightcurve_times_mjd = MJDREF + (total_lightcurve_times/86400.0)

        data_lc = pd.DataFrame({"RATE":total_lightcurve_rates, "TIME":total_lightcurve_times, "ERROR":total_lightcurve_errates, "MJD": total_lightcurve_times_mjd, "OBSERVATION":observations})
        data_lc = data_lc.sort_values(by=['TIME'])
        
        #Add year column to dataframe
        year_array = []
        for mjd in data_lc['MJD'].values:
            if mjd<51910:
                year_array.append(int(2000))
            elif mjd<52275:
                year_array.append(int(2001))
            elif mjd<52640:
                year_array.append(int(2002))
            elif mjd<53005:
                year_array.append(int(2003))
            elif mjd<53371:
                year_array.append(int(2004))
            elif mjd<53736:
                year_array.append(int(2005))
            elif mjd<54101:
                year_array.append(int(2006))
            elif mjd<54466:
                year_array.append(int(2007))
            elif mjd<54832:
                year_array.append(int(2008))
            elif mjd<55197:
                year_array.append(int(2009))
            elif mjd<55562:
                year_array.append(int(2010))
            elif mjd<55927:
                year_array.append(int(2011))
            elif mjd<56293:
                year_array.append(int(2012))
            elif mjd<56658:
                year_array.append(int(2013))
            elif mjd<57023:
                year_array.append(int(2014))
            elif mjd<57388:
                year_array.append(int(2015))
            elif mjd<57754:
                year_array.append(int(2016))
            elif mjd<58119:
                year_array.append(int(2017))
            elif mjd<58484:
                year_array.append(int(2018))
            elif mjd<58849:
                year_array.append(int(2019))

        data_lc['YEAR'] = year_array  
        data_lc = data_lc.reset_index(drop=True)
        data_lc = data_lc.reset_index()
        
        #indexes where to place the xticks in the plot
        year_endpoints = []
        for i in range(1, len(data_lc)):
            if data_lc['YEAR'][i] != data_lc['YEAR'][i-1]:
                year_endpoints.append(data_lc['index'][i])
        labels = np.linspace(start=2000,stop=2019,num=20,dtype=int)
        labels = np.delete(labels, [12,13])  #delete year 2012 and 2013
        year_endpoints = np.array(year_endpoints)
        year_endpoints = np.insert(year_endpoints, 0, values=0)
        
        if args.logpar:
            data_spec_zlogp = data_spec_zlogp[data_spec_zlogp['phoindex_up']!=0]
            df_plot_zlogp = pd.DataFrame({'tbinstart': data_spec_zlogp['tbinstart'].values, 'tbinstop': data_spec_zlogp['tbinstop'].values , 'phoindex': data_spec_zlogp['phoindex'].values,
                                            "phoindex_top": data_spec_zlogp['phoindex_up'].values - data_spec_zlogp['phoindex'].values,
                                            "phoindex_bot": data_spec_zlogp['phoindex'].values - data_spec_zlogp['phoindex_low'].values, 
                                            "beta": data_spec_zlogp['beta'].values, "beta_top": data_spec_zlogp['beta_up'].values - data_spec_zlogp['beta'].values,
                                            "beta_bot": data_spec_zlogp['beta'].values - data_spec_zlogp['beta_low'].values,
                                            "nH": data_spec_zlogp['nH'].values, "nH_up": data_spec_zlogp['nH_up'].values, "nH_low": data_spec_zlogp['nH_low'].values,
                                            "nH_top": data_spec_zlogp['nH_up'].values - data_spec_zlogp['nH'].values,
                                            "nH_bot": data_spec_zlogp['nH'].values - data_spec_zlogp['nH_low'].values,
                                            "obsid": data_spec_zlogp['obsid'].values})
            #Order dataframe and reset index
            df_plot_zlogp = df_plot_zlogp.sort_values(by=['tbinstart'])
            df_plot_zlogp = df_plot_zlogp[df_plot_zlogp['nH_up']!=0]
            df_plot_zlogp = df_plot_zlogp.reset_index(drop=True)
            df_plot_zlogp = df_plot_zlogp.reset_index()

            #Upper limits on nH
            df_plot_nH_pegged = df_plot_zlogp[df_plot_zlogp['nH_low']<=1e-4]
            df_plot_zlogp = df_plot_zlogp[df_plot_zlogp['nH_low']>1e-4]

            #Plot panel            
            fig_logp, axs_logp = plt.subplots(4, 1, figsize=(16,13), sharex=True, gridspec_kw={'hspace':0, 'wspace':0})
    
            axs_logp[0].errorbar('index', 'RATE', 'ERROR', data=data_lc, ecolor='black', linestyle='', color='black')
            axs_logp[1].errorbar('index', 'phoindex', yerr=(df_plot_zlogp['phoindex_bot'].values, df_plot_zlogp['phoindex_top'].values),
                                data=df_plot_zlogp, ecolor='black', linestyle='', color='black', capthick=1, elinewidth=1)
            axs_logp[2].errorbar('index', 'beta', yerr=(df_plot_zlogp['beta_bot'].values, df_plot_zlogp['beta_top'].values),
                                data=df_plot_zlogp, ecolor='black', linestyle='', color='black', capthick=1, elinewidth=1)
            axs_logp[3].errorbar('index', 'nH', yerr=(df_plot_zlogp['nH_bot'].values, df_plot_zlogp['nH_top'].values), data=df_plot_zlogp,
                                 ecolor='black', linestyle='', color='black', capthick=1, elinewidth=1)
            axs_logp[3].errorbar('index', 'nH_up', yerr='nH_top', data=df_plot_nH_pegged, uplims=True, linestyle='', capthick=1, elinewidth=1, ecolor='black')

            #Add vertical lines separating years
            axs_logp[0].vlines(year_endpoints, 0, 60, colors='r', linestyles='solid')
            axs_logp[1].vlines(year_endpoints, 1, 3.6, colors='r', linestyles='solid')
            axs_logp[2].vlines(year_endpoints, -0.7, 2, colors='r', linestyles='solid')
            axs_logp[3].vlines(year_endpoints, -0.015, 0.145, colors='r', linestyles='solid')
            
            #Bellurie
            axs_logp[0].grid()
            axs_logp[1].grid()
            axs_logp[2].grid()
            axs_logp[3].grid()
            axs_logp[0].set_ylabel("Rate [ct/s]")
            axs_logp[1].set_ylabel("phoindex")
            axs_logp[2].set_ylabel("beta")
            axs_logp[3].set_ylabel("nH [$10^{22}$ g/cm$^2$]")
            axs_logp[3].set_xlabel("Year")
            axs_logp[0].margins(0)
            axs_logp[1].margins(0)
            axs_logp[2].margins(0)
            axs_logp[3].margins(0)
            plt.xticks(ticks=year_endpoints, labels=labels, rotation=60)
            plt.savefig(os.path.join(target_dir, "Products", "Plots_spectra", "panel_logpar.png"))
            plt.show()
            

        if args.powerlaw:

            data_spec_zpowe = data_spec_zpowe[data_spec_zpowe['phoindex_up']!=0]
            df_plot_powerlaw = pd.DataFrame({'tbinstart': data_spec_zpowe['tbinstart'].values, 'tbinstop': data_spec_zpowe['tbinstop'].values , 'phoindex': data_spec_zpowe['phoindex'].values,
                                            "phoindex_top": data_spec_zpowe['phoindex_up'].values - data_spec_zpowe['phoindex'].values,
                                            "phoindex_bot": data_spec_zpowe['phoindex'].values - data_spec_zpowe['phoindex_low'].values, 
                                            "nH": data_spec_zpowe['nH'].values, "nH_up": data_spec_zpowe['nH_up'].values, "nH_low": data_spec_zpowe['nH_low'].values,
                                            "nH_top": data_spec_zpowe['nH_up'].values - data_spec_zpowe['nH'].values,
                                            "nH_bot": data_spec_zpowe['nH'].values - data_spec_zpowe['nH_low'].values,
                                            "obsid": data_spec_zpowe['obsid'].values})
            
            #Order dataframe and reset index
            df_plot_powerlaw = df_plot_powerlaw.sort_values(by=['tbinstart'])
            df_plot_powerlaw = df_plot_powerlaw[df_plot_powerlaw['nH_up']!=0]
            df_plot_powerlaw = df_plot_powerlaw.reset_index(drop=True)
            df_plot_powerlaw = df_plot_powerlaw.reset_index()

            #Upper limits on nH
            df_plot_nH_pegged = df_plot_powerlaw[df_plot_powerlaw['nH_low']<=1e-4]
            df_plot_powerlaw = df_plot_powerlaw[df_plot_powerlaw['nH_low']>1e-4]
            
            #Plot panel            
            fig_pw, axs_pw = plt.subplots(3, 1, figsize=(16,12), sharex=True, gridspec_kw={'hspace':0, 'wspace':0})
    
            axs_pw[0].errorbar('index', 'RATE', 'ERROR', data=data_lc, ecolor='black', linestyle='', color='black')
            axs_pw[1].errorbar('index', 'phoindex', yerr=(df_plot_powerlaw['phoindex_bot'].values, df_plot_powerlaw['phoindex_top'].values),
                                data=df_plot_powerlaw, ecolor='black', linestyle='', color='black', capthick=1, elinewidth=1)
            axs_pw[2].errorbar('index', 'nH', yerr=(df_plot_powerlaw['nH_bot'].values, df_plot_powerlaw['nH_top'].values), data=df_plot_powerlaw,
                                 ecolor='black', linestyle='', color='black', capthick=1, elinewidth=1)
            axs_pw[2].errorbar('index', 'nH_up', yerr='nH_top', data=df_plot_nH_pegged, uplims=True, linestyle='', capthick=1, elinewidth=1, ecolor='black')

            #Add vertical lines separating years
            axs_pw[0].vlines(year_endpoints, 0, 60, colors='r', linestyles='solid')
            axs_pw[1].vlines(year_endpoints, 1.5, 3.7, colors='r', linestyles='solid')
            axs_pw[2].vlines(year_endpoints, -0.015, 0.145, colors='r', linestyles='solid')
            
            #Bellurie
            axs_pw[0].grid()
            axs_pw[1].grid()
            axs_pw[2].grid()
            axs_pw[0].set_ylabel("Rate [ct/s]")
            axs_pw[1].set_ylabel("phoindex")
            axs_pw[2].set_ylabel("nH [$10^{22}$ g/cm$^2$]")
            axs_pw[2].set_xlabel("Year")
            axs_pw[0].margins(0)
            axs_pw[1].margins(0)
            axs_pw[2].margins(0)
            plt.xticks(ticks=year_endpoints, labels=labels, rotation=60)
            plt.savefig(os.path.join(target_dir, "Products", "Plots_spectra", "panel_powerlaw.png"))
            plt.show()
    

    if args.distribution:

        if args.logpar:

            # Separate data into low and high state (rate >20)
            data_spec_zlogp = data_spec_zlogp[data_spec_zlogp['phoindex_up']!=0]
            data_spec_zlogp_high = data_spec_zlogp[data_spec_zlogp['rate']>=20]
            data_spec_zlogp_low = data_spec_zlogp[data_spec_zlogp['rate']<20]

            # Plot distribution for alpha and beta 
            fig, axs = plt.subplots(2, 2, figsize=(5,7), sharex=True, gridspec_kw={'hspace':0.1})
            axs[0,0].hist(data_spec_zlogp_high['phoindex'].values, bins=15, color='red')
            axs[1,0].hist(data_spec_zlogp_low['phoindex'].values, bins=15, color='blue')
            axs[0,1].hist(data_spec_zlogp_high['beta'].values, bins=15, color='red')
            axs[1,1].hist(data_spec_zlogp_low['beta'].values, bins=15, color='blue')
            red_patch =  Patch(facecolor='red', edgecolor='black', label='high state')
            blue_patch =  Patch(facecolor='b', edgecolor='black', label='low state')
            axs[0,1].legend(handles=[red_patch, blue_patch], loc='upper right',  fancybox=True, shadow=True)
            axs[1,0].set_xlabel('phoindex logparabola')
            axs[1,1].set_xlabel('beta logparabola')
            plt.savefig(os.path.join(target_dir, "Products", "Plots_spectra", "distribution_logpar.png"))

        if args.powerlaw:

            # Separate data into low and high state (rate >20)
            data_spec_zpowe = data_spec_zpowe[data_spec_zpowe['phoindex_up']!=0]
            data_spec_zpowe_high = data_spec_zpowe[data_spec_zpowe['rate']>=20]
            data_spec_zpowe_low = data_spec_zpowe[data_spec_zpowe['rate']<20]

            #Plot distribution only for photon index
            fig, axs = plt.subplots(2, 1, figsize=(5,7), sharex=True, gridspec_kw={'hspace':0.1})
            y_high, x_high, _ = axs[0].hist(data_spec_zpowe_high['phoindex'].values, bins=15, color='red')
            y_low, x_low, _ = axs[1].hist(data_spec_zpowe_low['phoindex'].values, bins=15, color='blue')
            red_patch =  Patch(facecolor='red', edgecolor='black', label='high state')
            blue_patch =  Patch(facecolor='b', edgecolor='black', label='low state')
            axs[0].legend(handles=[red_patch, blue_patch])
            axs[0].set_xlabel('phoindex powerlaw')
            plt.savefig(os.path.join(target_dir, "Products", "Plots_spectra", "distribution_powerlaw.png"))

