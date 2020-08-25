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
import random
from scipy.optimize import curve_fit
from scipy.stats import linregress
from matplotlib.offsetbox import AnchoredText
from scipy.stats import ks_2samp, chisquare


parser = ArgumentParser(description=__doc__)
parser.add_argument('--epic', action='store_true',
                    help='use EPIC results')
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
parser.add_argument('--hysteresis', type=int,
                    help='make phoindex vs rate plot of a given observation ID to search for hysteresis cycles')
parser.add_argument('--vaughan', action='store_true',
                    help='use only long observations (the ones used to make vaughan variability panels)')

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
    if args.epic:
        hdul_lc = Table.read(os.path.join(products_dir, "EPIC_Lightcurves", "EPIC_obs_table.fits"), hdu=1)   
        data_lc = hdul_lc.to_pandas()
    else:
        hdul_lc = Table.read(os.path.join(products_dir, "RGS_Lightcurves", "obs_table.fits"), hdu=1)   
        data_lc = hdul_lc.to_pandas()
        data_lc = data_lc.dropna(subset=['RGS_Rate'])   #do not consider NaN Rates

    #Read spectra table for logpar and powerlaw parameters
    if args.epic:
        hdul_spec = Table.read(os.path.join(products_dir, "EPIC_Spectra", "EPIC_spectra_table.fits"), hdu=1)
        data_spec = hdul_spec.to_pandas()
    else:
        hdul_spec = Table.read(os.path.join(products_dir, "RGS_Spectra", "spectra_table.fits"), hdu=1)
        data_spec = hdul_spec.to_pandas()
        hdul_spec_avg = Table.read(os.path.join(products_dir, "RGS_Spectra", "spectra_table_average.fits"), hdu=1)
        data_spec_avg = hdul_spec_avg.to_pandas()

        if args.average:
            hdul_spec = Table.read(os.path.join(products_dir, "RGS_Spectra", "spectra_table_average.fits"), hdu=1)
            data_spec = hdul_spec.to_pandas()  #Use only average spectra
        elif args.bins:
            data_spec = data_spec[data_spec['tbinid'] != 0]   #Use divided bins
        


    #Clean data
    data_spec_clean = data_spec[data_spec['obsid'] != 658802001]

    # Make dataframe for the two different models
    if args.epic:
        data_spec_zlogp = data_spec_clean[data_spec_clean['model']=='TBabs*zlogpar']
        data_spec_zpowe = data_spec_clean[data_spec_clean['model']=='TBabs*zpowerlw']
    else:
        data_spec_zlogp = data_spec_clean[data_spec_clean['model']=='constant*TBabs*zlogp']
        data_spec_zpowe = data_spec_clean[data_spec_clean['model']=='constant*TBabs*zpowe']


    if args.fvar and not args.epic:

        if args.phoindex: #user wants phoindex vs fvar
            
            df_plot_zlogp = pd.DataFrame({"fvar": data_lc['F_var'].values, "fvar_err": data_lc['F_var_sigma'].values, 'alpha': data_spec_zlogp['phoindex'].values, 
                                          "xs": data_lc['Excess_Variance'].values, "rate": data_spec_zlogp['rate'].values,
                                          'alpha_top': data_spec_zlogp['phoindex_up'].values - data_spec_zlogp['phoindex'].values,
                                          'alpha_bot': data_spec_zlogp['phoindex'].values - data_spec_zlogp['phoindex_low'].values,
                                          "ftest": data_spec_zlogp['ftest'].values, "obsid": data_lc['ObsId'].values})
            df_plot_zlogp = df_plot_zlogp[df_plot_zlogp['xs'] >0]
            
            
            df_plot_zpowe = pd.DataFrame({"fvar": data_lc['F_var'].values, "fvar_err": data_lc['F_var_sigma'].values,
                                         "xs": data_lc['Excess_Variance'].values, "rate": data_spec_zpowe['rate'].values,
                                         "obsid": data_lc['ObsId'].values, "phoindex": data_spec_zpowe['phoindex'].values,
                                         "phoindex_top": data_spec_zpowe['phoindex_up'].values - data_spec_zpowe['phoindex'].values,
                                         "phoindex_bot": data_spec_zpowe['phoindex'].values - data_spec_zpowe['phoindex_low'].values,
                                         "ftest": data_spec_zpowe['ftest'].values})
            df_plot_zpowe = df_plot_zpowe[df_plot_zpowe['xs'] >0]

            if args.vaughan:
                vaughan_obs = ['0136541001', '0158971201', '0810860201', '0411080301', '0560980101', '0791781401', '0810860701', '0791782001']
                df_plot_zpowe = df_plot_zpowe[df_plot_zpowe.obsid.isin(vaughan_obs)] #filter on obsid
                df_plot_zlogp = df_plot_zlogp[df_plot_zlogp.obsid.isin(vaughan_obs)] #filter on obsid
            
            #Distinguish model favoring based on Ftest column
            df_plot_zlogp_better = df_plot_zlogp[(df_plot_zlogp['ftest']<0.1) & (df_plot_zlogp['ftest']!=-999.)]
            df_plot_zpowe_better = df_plot_zpowe[(df_plot_zpowe['ftest']>=0.1) | (df_plot_zpowe['ftest']==-999.)]

            # Plot
            fig, axs =plt.subplots(1, 2, figsize=(7,5), sharey=False, gridspec_kw={'wspace':0.3})
        
            axs[0].errorbar(x=df_plot_zlogp_better['fvar'].values, y=df_plot_zlogp_better['alpha'].values, xerr=df_plot_zlogp_better['fvar_err'].values,
                        yerr = (df_plot_zlogp_better['alpha_bot'].values, df_plot_zlogp_better['alpha_top'].values), color='b', fmt='.', markersize=5, ecolor='cornflowerblue', elinewidth=1, capsize=2, capthick=1, label='logpar best model')
   
            axs[1].errorbar(y=df_plot_zpowe_better['phoindex'].values, x=df_plot_zpowe_better['fvar'].values, xerr=df_plot_zpowe_better['fvar_err'].values,
                        yerr = (df_plot_zpowe_better['phoindex_bot'].values, df_plot_zpowe_better['phoindex_top'].values), color='red', fmt='.', markersize=5, ecolor='lightcoral', elinewidth=1, capsize=2, capthick=1, label='powerlaw best model')

            axs[0].grid(True)
            axs[1].grid(True)
            axs[0].set_title('Logparabola')
            axs[1].set_title('Powerlaw')
            axs[0].set_xlabel('$F_{var}$', fontsize=14)
            axs[0].set_ylabel(r'$ \alpha$', fontsize=14)
            axs[1].set_ylabel('$\Gamma$', fontsize=14)
            axs[1].set_xlabel('$F_{var}$', fontsize=14)

            plt.savefig(os.path.join(target_dir, "Products", "Plots_spectra", "fvar_phoindex_correlations.png"))

        if args.beta: #user wants beta vs fvar
            
            df_plot_zlogp = pd.DataFrame({"fvar": data_lc['F_var'].values, "fvar_err": data_lc['F_var_sigma'].values,
                                          "xs": data_lc['Excess_Variance'].values,
                                          'beta': data_spec_zlogp['beta'].values, 
                                          'beta_top': data_spec_zlogp['beta_up'].values - data_spec_zlogp['beta'].values,
                                          'beta_bot': data_spec_zlogp['beta'].values - data_spec_zlogp['beta_low'].values,
                                          "ftest": data_spec_zlogp['ftest'].values, "obsid": data_lc['ObsId'].values})
            df_plot_zlogp = df_plot_zlogp[df_plot_zlogp['xs']>0]

            # Distinguish between fvar high (+) and fvar low (x)
            df_plot_zlogp_better = df_plot_zlogp[(df_plot_zlogp['ftest']<0.1) & (df_plot_zlogp['ftest']!=-999.)]

            #Plot
            figure = plt.figure(figsize=(6,5))
            plt.errorbar(y=df_plot_zlogp_better['beta'].values,x= df_plot_zlogp_better['fvar'].values, xerr=df_plot_zlogp_better['fvar_err'].values,
                        yerr = (df_plot_zlogp_better['beta_bot'].values, df_plot_zlogp_better['beta_top'].values), color='b', linestyle='', marker='.', markersize=5, ecolor='cornflowerblue', elinewidth=1, capsize=2, capthick=1, label='logpar best model')
            
            plt.grid()
            plt.ylabel(r'$\beta$', fontsize=14)
            plt.xlabel('$F_{var}$', fontsize=14)
            plt.title('Logparabola', fontsize=15)
            plt.tight_layout()
            plt.savefig(os.path.join(target_dir, "Products", "Plots_spectra", "fvar_beta_correlations.png"))


    if args.fvar and args.epic:

        
        df_plot_zlogp = pd.DataFrame({"fvar_soft": data_lc['fvar_soft'].values, "efvar_soft": data_lc['efvar_soft'].values, 
                                        "fvar_hard": data_lc['fvar_hard'].values, "efvar_hard": data_lc['efvar_hard'].values,
                                        'alpha': data_spec_zlogp['phoindex'].values, 
                                        'alpha_top': data_spec_zlogp['phoindex_up'].values - data_spec_zlogp['phoindex'].values,
                                        'alpha_bot': data_spec_zlogp['phoindex'].values - data_spec_zlogp['phoindex_low'].values,
                                        'beta': data_spec_zlogp['beta'].values, 
                                        'beta_top': data_spec_zlogp['beta_up'].values - data_spec_zlogp['beta'].values,
                                        'beta_bot': data_spec_zlogp['beta'].values - data_spec_zlogp['beta_low'].values,
                                        "obsid": data_lc['ObsId'].values})
        
        
        df_plot_zpowe = pd.DataFrame({"fvar_soft": data_lc['fvar_soft'].values, "efvar_soft": data_lc['efvar_soft'].values, 
                                        "fvar_hard": data_lc['fvar_hard'].values, "efvar_hard": data_lc['efvar_hard'].values,
                                        "obsid": data_lc['ObsId'].values, "phoindex": data_spec_zpowe['phoindex'].values,
                                        "phoindex_top": data_spec_zpowe['phoindex_up'].values - data_spec_zpowe['phoindex'].values,
                                        "phoindex_bot": data_spec_zpowe['phoindex'].values - data_spec_zpowe['phoindex_low'].values})

        # Plot
        fig, axs =plt.subplots(1, 3, figsize=(7,6), sharey=True, gridspec_kw={'wspace':0.2})
        axs[0].errorbar(df_plot_zlogp['alpha'].values, df_plot_zlogp['fvar_soft'].values, yerr=df_plot_zlogp['efvar_soft'].values,
                    xerr = (df_plot_zlogp['alpha_bot'].values, df_plot_zlogp['alpha_top'].values), color='orange', fmt='.', markersize=5, ecolor='peachpuff', elinewidth=1, capsize=2, capthick=1, label='soft LC')
        axs[0].errorbar(df_plot_zlogp['alpha'].values, df_plot_zlogp['fvar_hard'].values, yerr=df_plot_zlogp['efvar_hard'].values,
                    xerr = (df_plot_zlogp['alpha_bot'].values, df_plot_zlogp['alpha_top'].values), color='g', fmt='.', markersize=5, ecolor='yellowgreen', elinewidth=1, capsize=2, capthick=1, label='hard LC')
    
        axs[1].errorbar(df_plot_zpowe['phoindex'].values, df_plot_zpowe['fvar_soft'].values, yerr=df_plot_zpowe['efvar_soft'].values,
                    xerr = (df_plot_zpowe['phoindex_bot'].values, df_plot_zpowe['phoindex_top'].values), color='orange', fmt='.', markersize=5, ecolor='peachpuff', elinewidth=1, capsize=2, capthick=1, label='soft LC')
        axs[1].errorbar(df_plot_zpowe['phoindex'].values, df_plot_zpowe['fvar_hard'].values, yerr=df_plot_zpowe['efvar_hard'].values,
                    xerr = (df_plot_zpowe['phoindex_bot'].values, df_plot_zpowe['phoindex_top'].values), color='g', fmt='.', markersize=5, ecolor='yellowgreen', elinewidth=1, capsize=2, capthick=1, label='hard LC')
            
        axs[2].errorbar(df_plot_zlogp['beta'].values, df_plot_zlogp['fvar_soft'].values, yerr=df_plot_zlogp['efvar_soft'].values,
                    xerr = (df_plot_zlogp['beta_bot'].values, df_plot_zlogp['beta_top'].values), color='orange', fmt='.', markersize=5, ecolor='peachpuff', elinewidth=1, capsize=2, capthick=1, label='soft LC')
        axs[2].errorbar(df_plot_zlogp['beta'].values, df_plot_zlogp['fvar_hard'].values, yerr=df_plot_zlogp['efvar_hard'].values,
                    xerr = (df_plot_zlogp['beta_bot'].values, df_plot_zlogp['beta_top'].values), color='g', fmt='.', markersize=5, ecolor='yellowgreen', elinewidth=1, capsize=2, capthick=1, label='hard LC')
        

        axs[0].grid(True)
        axs[1].grid(True)
        axs[2].grid(True)
        axs[0].set_title('Logparabola')
        axs[1].set_title('Powerlaw')
        axs[2].set_title('Logparabola')
        axs[0].set_ylabel('$F_{var}$', fontsize=10)
        axs[0].set_xlabel('alpha', fontsize=10)
        axs[1].set_xlabel('phoindex', fontsize=10)
        axs[2].set_xlabel('beta', fontsize=10)
        axs[2].legend(loc='upper left')
        plt.savefig(os.path.join(target_dir, "Products", "EPIC_Spectra", "fvar_param_correlations.png"))


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
    

    if args.hysteresis:
        
        # Make the dataframe for the plot of logpar
        data_spec_zlogp = data_spec_zlogp[data_spec_zlogp['phoindex_up']!=0.]   #we want to study details of light curve
        df_plot_zlogp = pd.DataFrame({"time":data_spec_zlogp['tbinstart'].values, "rate": data_spec_zlogp['rate'].values, "erate": data_spec_zlogp['erate'].values, 'alpha': data_spec_zlogp['phoindex'].values, 
                                        'alpha_top': data_spec_zlogp['phoindex_up'].values - data_spec_zlogp['phoindex'].values,
                                        'alpha_bot': data_spec_zlogp['phoindex'].values - data_spec_zlogp['phoindex_low'].values,
                                        "obsid": data_spec_zlogp['obsid'].values})
        
        # Take only data from observation specified by command line
        mjdref = 50814.0
        df_plot_zlogp = df_plot_zlogp[df_plot_zlogp['obsid']==args.hysteresis]
        df_plot_zlogp['alpha'] = df_plot_zlogp['alpha'].apply(lambda x : -x)
        df_plot_zlogp['time'] = df_plot_zlogp['time'].apply(lambda x : mjdref+ (x/86400.0))

        #Same for powerlaw
        data_spec_zpowe = data_spec_zpowe[data_spec_zpowe['phoindex_up']!=0.]
        df_plot_zpowe = pd.DataFrame({"time": data_spec_zpowe['tbinstart'].values,"rate": data_spec_zpowe['rate'].values, "erate": data_spec_zpowe['erate'].values,
                                        "phoindex": data_spec_zpowe['phoindex'].values,
                                        "phoindex_top": data_spec_zpowe['phoindex_up'].values - data_spec_zpowe['phoindex'].values,
                                        "phoindex_bot": data_spec_zpowe['phoindex'].values - data_spec_zpowe['phoindex_low'].values,
                                        "obsid": data_spec_zpowe['obsid'].values})
        
        df_plot_zpowe = df_plot_zpowe[df_plot_zpowe['obsid']==args.hysteresis]
        df_plot_zpowe['phoindex'] = df_plot_zpowe['phoindex'].apply(lambda x : -x)
        df_plot_zpowe['time'] = df_plot_zpowe['time'].apply(lambda x : mjdref+ (x/86400.0))

        if args.hysteresis==136540101:
            loop1 = df_plot_zpowe['phoindex'].values[10:14]
            x_loop1 = df_plot_zpowe['rate'].values[10:14]

            loop2 = df_plot_zpowe['phoindex'].values[15:20]
            x_loop2 = df_plot_zpowe['rate'].values[15:20]        

            loop3 = df_plot_zpowe['phoindex'].values[21:26]
            x_loop3 = df_plot_zpowe['rate'].values[21:26]

            loop4 = df_plot_zpowe['phoindex'].values[27:31]
            x_loop4 = df_plot_zpowe['rate'].values[27:31]

            fig = plt.figure(figsize=(8,8))
            axs1 = plt.subplot(313)
            axs1.errorbar(y=df_plot_zpowe['rate'].values[0:9], x=df_plot_zpowe['time'].values[0:9], yerr=df_plot_zpowe['erate'].values[0:9],
                            fmt='.',color='black', ecolor='gray', marker='.', markersize=5, elinewidth=1, capsize=2, capthick=1)
            axs1.errorbar(y=df_plot_zpowe['rate'].values[10:14], x=df_plot_zpowe['time'].values[10:14], yerr=df_plot_zpowe['erate'].values[10:14],
                            fmt='.',color='red', ecolor='gray', marker='o', markersize=5, elinewidth=1, capsize=2, capthick=1, label='Epoch 1')
            axs1.errorbar(y=df_plot_zpowe['rate'].values[15:20], x=df_plot_zpowe['time'].values[15:20], yerr=df_plot_zpowe['erate'].values[15:20],
                            fmt='.',color='blue', ecolor='gray', marker='v', markersize=5, elinewidth=1, capsize=2, capthick=1, label='Epoch 2')
            axs1.errorbar(y=df_plot_zpowe['rate'].values[21:26], x=df_plot_zpowe['time'].values[21:26], yerr=df_plot_zpowe['erate'].values[21:26],
                            fmt='.',color='g', ecolor='gray', marker='^', markersize=5, elinewidth=1, capsize=2, capthick=1, label='Epoch 3')
            axs1.errorbar(y=df_plot_zpowe['rate'].values[27:31], x=df_plot_zpowe['time'].values[27:31], yerr=df_plot_zpowe['erate'].values[27:31],
                            fmt='.',color='orange', ecolor='gray', marker='s', markersize=5, elinewidth=1, capsize=2, capthick=1, label='Epoch 4')
            axs1.errorbar(y=df_plot_zpowe['rate'].values[32:-1], x=df_plot_zpowe['time'].values[32:-1], yerr=df_plot_zpowe['erate'].values[32:-1],
                            fmt='.',color='black', ecolor='gray', marker='.', markersize=5, elinewidth=1, capsize=2, capthick=1,)
            axs1.set_xlabel('Time [MJD]', fontsize=12)
            axs1.set_ylabel('Rate [ct/s]', fontsize=12)
            axs1.ticklabel_format(useOffset=False)
            axs1.legend(ncol=4)

            axs2 = plt.subplot(321)
            axs2.set_title('Epoch 1', fontsize=11)
            axs2.set_ylabel('$\Gamma$ (powerlaw)', fontsize=12)
            axs2.errorbar(y=loop1, x=x_loop1, fmt='.',color='red', markersize=5, elinewidth=1, capsize=2, capthick=1, ecolor='gray', label='Epoch 1')
            axs2.quiver(x_loop1[:-1], loop1[:-1], x_loop1[1:]-x_loop1[:-1], loop1[1:]-loop1[:-1], scale_units='xy', angles='xy', scale=1, width=0.006, headwidth=10)

            axs3 = plt.subplot(322)
            axs3.set_title('Epoch 2', fontsize=11)
            axs3.errorbar(y=loop2, x=x_loop2, fmt='.',color='blue', markersize=5, elinewidth=1, capsize=2, capthick=1, ecolor='gray', label='Epoch 2')
            axs3.quiver(x_loop2[:-1], loop2[:-1], x_loop2[1:]-x_loop2[:-1], loop2[1:]-loop2[:-1], scale_units='xy', angles='xy', scale=1, width=0.006, headwidth=10)
        
            axs4 = plt.subplot(323)
            axs4.set_title('Epoch 3', fontsize=11)
            axs4.set_ylabel('$\Gamma$ (powerlaw)', fontsize=12)
            axs4.errorbar(y=loop3, x=x_loop3, fmt='.',color='g', markersize=5, elinewidth=1, capsize=2, capthick=1, ecolor='gray', label='Epoch 3')
            axs4.quiver(x_loop3[:-1], loop3[:-1], x_loop3[1:]-x_loop3[:-1], loop3[1:]-loop3[:-1], scale_units='xy', angles='xy', scale=1, width=0.006, headwidth=10)
            
            axs5 = plt.subplot(324)
            axs5.set_title('Epoch 4', fontsize=11)
            axs5.set_xlabel('Rate [ct/s]', fontsize=11)
            axs5.errorbar(y=loop4, x=x_loop4,fmt='.',color='orange', markersize=5, elinewidth=1, capsize=2, capthick=1, ecolor='gray', label='Epoch 3')
            axs5.quiver(x_loop4[:-1], loop4[:-1], x_loop4[1:]-x_loop4[:-1], loop4[1:]-loop4[:-1], scale_units='xy', angles='xy', scale=1, width=0.006, headwidth=10)
            
            plt.tight_layout()
            plt.savefig(os.path.join(target_dir, "Products", "Plots_spectra", f"hysteresis_{args.hysteresis}.png"))


        if args.hysteresis==150498701:
            loop1 = df_plot_zpowe['phoindex'].values[1:7]
            x_loop1 = df_plot_zpowe['rate'].values[1:7]

            loop2 = df_plot_zpowe['phoindex'].values[7:16]
            x_loop2 = df_plot_zpowe['rate'].values[7:16]

            loop3 = df_plot_zpowe['phoindex'].values[23:28]
            x_loop3 = df_plot_zpowe['rate'].values[23:28]

            loop4 = df_plot_zpowe['phoindex'].values[33:40]
            x_loop4 = df_plot_zpowe['rate'].values[33:40]

            loop5 = df_plot_zpowe['phoindex'].values[43:47]
            x_loop5 = df_plot_zpowe['rate'].values[43:47]


            fig = plt.figure(figsize=(11,13))
            axs1 = plt.subplot(331)
            axs1.set_title('Epoch 1', fontsize=12)
            axs1.set_ylabel('$\Gamma$ (powerlaw)', fontsize=14)
            axs1.errorbar(y=loop1, x=x_loop1, fmt='.',color='red', markersize=5, elinewidth=1, capsize=2, capthick=1, ecolor='gray', label='Epoch 1')
            axs1.quiver(x_loop1[:-1], loop1[:-1], x_loop1[1:]-x_loop1[:-1], loop1[1:]-loop1[:-1], scale_units='xy', angles='xy', scale=1, width=0.006, headwidth=10)
            
            axs2 = plt.subplot(332)
            axs2.set_title('Epoch 2', fontsize=12)
            axs2.errorbar(y=loop2, x=x_loop2, fmt='.',color='blue', markersize=5, elinewidth=1, capsize=2, capthick=1, ecolor='gray', label='Epoch 2')
            axs2.quiver(x_loop2[:-1], loop2[:-1], x_loop2[1:]-x_loop2[:-1], loop2[1:]-loop2[:-1], scale_units='xy', angles='xy', scale=1, width=0.006, headwidth=10)

            axs3 = plt.subplot(333)
            axs3.set_title('Epoch 3', fontsize=12)
            axs3.errorbar(y=loop3, x=x_loop3, fmt='.',color='g', markersize=5, elinewidth=1, capsize=2, capthick=1, ecolor='gray', label='Epoch 3')
            axs3.quiver(x_loop3[:-1], loop3[:-1], x_loop3[1:]-x_loop3[:-1], loop3[1:]-loop3[:-1], scale_units='xy', angles='xy', scale=1, width=0.006, headwidth=10)

            axs4 = plt.subplot(323)
            axs4.set_title('Epoch 4', fontsize=12)
            axs4.set_ylabel('$\Gamma$ (powerlaw)', fontsize=14)
            axs4.errorbar(y=loop4, x=x_loop4, fmt='.',color='orange', markersize=5, elinewidth=1, capsize=2, capthick=1, ecolor='gray', label='Epoch 3')
            axs4.quiver(x_loop4[:-1], loop4[:-1], x_loop4[1:]-x_loop4[:-1], loop4[1:]-loop4[:-1], scale_units='xy', angles='xy', scale=1, width=0.006, headwidth=10)

            axs5 = plt.subplot(324)
            axs5.set_title('Epoch 5', fontsize=12)
            axs5.set_xlabel('Rate [ct/s]', fontsize=14)
            axs5.errorbar(y=loop5, x=x_loop5, fmt='.',color='purple', markersize=5, elinewidth=1, capsize=2, capthick=1, ecolor='gray', label='Epoch 3')
            axs5.quiver(x_loop5[:-1], loop5[:-1], x_loop5[1:]-x_loop5[:-1], loop5[1:]-loop5[:-1], scale_units='xy', angles='xy', scale=1, width=0.006, headwidth=10)

            axs6 = plt.subplot(313)
            axs6.errorbar(y=df_plot_zpowe['rate'].values[1:7], x=df_plot_zpowe['time'].values[1:7], yerr=df_plot_zpowe['erate'].values[1:7],
                            fmt='.',color='r', ecolor='gray', marker='d', markersize=7, elinewidth=1, capsize=2, capthick=1, label='Epoch 1')
            axs6.errorbar(y=df_plot_zpowe['rate'].values[7:16], x=df_plot_zpowe['time'].values[7:16], yerr=df_plot_zpowe['erate'].values[7:16],
                            fmt='.',color='blue', ecolor='gray', marker='o', markersize=7, elinewidth=1, capsize=2, capthick=1, label='Epoch 2')
            axs6.errorbar(y=df_plot_zpowe['rate'].values[23:28], x=df_plot_zpowe['time'].values[23:28], yerr=df_plot_zpowe['erate'].values[23:28],
                            fmt='.',color='g', ecolor='gray', marker='v', markersize=7, elinewidth=1, capsize=2, capthick=1, label='Epoch 3')
            axs6.errorbar(y=df_plot_zpowe['rate'].values[33:40], x=df_plot_zpowe['time'].values[33:40], yerr=df_plot_zpowe['erate'].values[33:40],
                            fmt='.',color='orange', ecolor='gray', marker='^', markersize=7, elinewidth=1, capsize=2, capthick=1, label='Epoch 4')
            axs6.errorbar(y=df_plot_zpowe['rate'].values[43:47], x=df_plot_zpowe['time'].values[43:47], yerr=df_plot_zpowe['erate'].values[43:47],
                            fmt='.',color='purple', ecolor='gray', marker='s', markersize=7, elinewidth=1, capsize=2, capthick=1, label='Epoch 5')
            axs6.errorbar(y=df_plot_zpowe['rate'].values[16:23], x=df_plot_zpowe['time'].values[16:23], yerr=df_plot_zpowe['erate'].values[16:23],
                            fmt='.',color='black', ecolor='gray', marker='.', markersize=5, elinewidth=1, capsize=2, capthick=1)
            axs6.errorbar(y=df_plot_zpowe['rate'].values[40:43], x=df_plot_zpowe['time'].values[40:43], yerr=df_plot_zpowe['erate'].values[40:43],
                            fmt='.',color='black', ecolor='gray', marker='.', markersize=5, elinewidth=1, capsize=2, capthick=1)
            axs6.errorbar(y=df_plot_zpowe['rate'].values[28:33], x=df_plot_zpowe['time'].values[28:33], yerr=df_plot_zpowe['erate'].values[28:33],
                            fmt='.',color='black', ecolor='gray', marker='.', markersize=5, elinewidth=1, capsize=2, capthick=1)
            axs6.set_xlabel('Time [MJD]', fontsize=14)
            axs6.set_ylabel('Rate [ct/s]', fontsize=14)
            axs6.ticklabel_format(useOffset=False)
            axs6.legend(ncol=2)
            plt.savefig(os.path.join(target_dir, "Products", "Plots_spectra", f"hysteresis_{args.hysteresis}.png"))


        if args.hysteresis==99280601:
            loop1 = df_plot_zpowe['phoindex'].values[0:4]
            x_loop1 = df_plot_zpowe['rate'].values[0:4]

            loop2 = df_plot_zpowe['phoindex'].values[5:12]
            x_loop2 = df_plot_zpowe['rate'].values[5:12]

            loop3 = df_plot_zpowe['phoindex'].values[12:23]
            x_loop3 = df_plot_zpowe['rate'].values[12:23]

            axs1 = plt.subplot(231)
            axs1.set_title('Epoch 1', fontsize=10)
            axs1.set_ylabel('$\Gamma$ (powerlaw)', fontsize=12)
            axs1.errorbar(y=loop1, x=x_loop1, fmt='.',color='red', markersize=5, elinewidth=1, capsize=2, capthick=1, ecolor='gray', label='Epoch 1')
            axs1.quiver(x_loop1[:-1], loop1[:-1], x_loop1[1:]-x_loop1[:-1], loop1[1:]-loop1[:-1], scale_units='xy', angles='xy', scale=1, width=0.006, headwidth=10)
            
            axs2 = plt.subplot(232)
            axs2.set_title('Epoch 2', fontsize=10)
            axs2.errorbar(y=loop2, x=x_loop2, fmt='.',color='blue', markersize=5, elinewidth=1, capsize=2, capthick=1, ecolor='gray', label='Epoch 1')
            axs2.quiver(x_loop2[:-1], loop2[:-1], x_loop2[1:]-x_loop2[:-1], loop2[1:]-loop2[:-1], scale_units='xy', angles='xy', scale=1, width=0.006, headwidth=10)
             
            axs3 = plt.subplot(233)
            axs3.set_title('Epoch 3', fontsize=10)
            axs3.set_xlabel('Rate [ct/s]', fontsize=10)
            axs3.errorbar(y=loop3, x=x_loop3, fmt='.',color='g', markersize=5, elinewidth=1, capsize=2, capthick=1, ecolor='gray', label='Epoch 1')
            axs3.quiver(x_loop3[:-1], loop3[:-1], x_loop3[1:]-x_loop3[:-1], loop3[1:]-loop3[:-1], scale_units='xy', angles='xy', scale=1, width=0.006, headwidth=10)
                       
            axs4 = plt.subplot(212)
            axs4.errorbar(y=df_plot_zpowe['rate'].values[0:4], x=df_plot_zpowe['time'].values[0:4], yerr=df_plot_zpowe['erate'].values[0:4],
                            fmt='.',color='r', ecolor='gray', marker='^', markersize=4, elinewidth=1, capsize=2, capthick=1, label='Epoch 1')
            axs4.errorbar(y=df_plot_zpowe['rate'].values[5:12], x=df_plot_zpowe['time'].values[5:12], yerr=df_plot_zpowe['erate'].values[5:12],
                            fmt='.',color='blue', ecolor='gray', marker='o', markersize=4, elinewidth=1, capsize=2, capthick=1, label='Epoch 2')
            axs4.errorbar(y=df_plot_zpowe['rate'].values[12:23], x=df_plot_zpowe['time'].values[12:23], yerr=df_plot_zpowe['erate'].values[12:23],
                            fmt='.',color='g', ecolor='gray', marker='s', markersize=4, elinewidth=1, capsize=2, capthick=1, label='Epoch 3')           
            axs4.errorbar(y=df_plot_zpowe['rate'].values[4], x=df_plot_zpowe['time'].values[4], yerr=df_plot_zpowe['erate'].values[4],
                            fmt='.',color='black', ecolor='gray', marker='.', markersize=5, elinewidth=1, capsize=2, capthick=1)           
            plt.legend()
            axs4.set_xlabel('Time [MJD]', fontsize=12)
            axs4.set_ylabel('Rate [ct/s]', fontsize=12)
            axs4.ticklabel_format(useOffset=False)
            plt.tight_layout()
            plt.savefig(os.path.join(target_dir, "Products", "Plots_spectra", f"hysteresis_{args.hysteresis}.png"))


        if args.hysteresis==136540801 :
            #Make figure
            figure, axs = plt.subplots(2,1, figsize=(8,6), gridspec_kw={'hspace':0.3})
            
            #Separate data into 3 parts: stable, ascending, descending

            ascending = []
            descending = []
            stable = []
            
            i = 0
            while i<len(df_plot_zpowe['rate']):

                if i==0:
                    if df_plot_zpowe['rate'].values[i]<=df_plot_zpowe['rate'].values[i+1]:
                        ascending.append(df_plot_zpowe['rate'].values[i])
                        i+=1
                    elif  df_plot_zpowe['rate'].values[i]>df_plot_zpowe['rate'].values[i+1]:
                        descending.append(df_plot_zpowe['rate'].values[i])
                        i+=1
                    continue

                if df_plot_zpowe['rate'].values[i]<df_plot_zpowe['rate'].values[i-1]:
                    descending.append(df_plot_zpowe['rate'].values[i])
                    i+=1
                elif df_plot_zpowe['rate'].values[i]>df_plot_zpowe['rate'].values[i-1] :
                    ascending.append(df_plot_zpowe['rate'].values[i])
                    i+=1          
                else:
                    stable.append(df_plot_zpowe['rate'].values[i])
                    i+=1
            
            df_plot_zpowe1 = df_plot_zpowe[df_plot_zpowe['rate'].isin(ascending)]
            df_plot_zpowe2 = df_plot_zpowe[df_plot_zpowe['rate'].isin(stable)]
            df_plot_zpowe3 = df_plot_zpowe[df_plot_zpowe['rate'].isin(descending)]

            df_plot_zlogp1 = df_plot_zlogp[df_plot_zlogp['rate'].isin(ascending)]
            df_plot_zlogp2 = df_plot_zlogp[df_plot_zlogp['rate'].isin(stable)]
            df_plot_zlogp3 = df_plot_zlogp[df_plot_zlogp['rate'].isin(descending)]
            
            
            #Subplot 1: lightcurve
            axs[0].errorbar(y=df_plot_zpowe1['rate'].values, x=df_plot_zpowe1['time'].values, yerr=df_plot_zpowe1['erate'].values,
                            fmt='.',color='red', ecolor='lightcoral', markersize=1.5, elinewidth=1, capsize=2, capthick=1)
            axs[0].errorbar(y=df_plot_zpowe3['rate'].values, x=df_plot_zpowe3['time'].values, yerr=df_plot_zpowe3['erate'].values,
                            fmt='.',color='b', ecolor='cornflowerblue', markersize=1.5, elinewidth=1, capsize=2, capthick=1)
            
            #Subplot 2 photon index vs rate (powerlaw)
            axs[1].errorbar(y=df_plot_zpowe1['phoindex'].values, x=df_plot_zpowe1['rate'].values, xerr=df_plot_zpowe1['erate'].values,
                        yerr = (df_plot_zpowe1['phoindex_bot'].values, df_plot_zpowe1['phoindex_top'].values), fmt='.',color='red', markersize=7, elinewidth=1, capsize=2, capthick=1, ecolor='lightcoral', label='zpowerlaw')
            axs[1].errorbar(y=df_plot_zpowe2['phoindex'].values, x=df_plot_zpowe2['rate'].values, xerr=df_plot_zpowe2['erate'].values,
                        yerr = (df_plot_zpowe2['phoindex_bot'].values, df_plot_zpowe2['phoindex_top'].values), fmt='.',color='g', markersize=7, elinewidth=1, capsize=2, capthick=1, ecolor='yellowgreen', label='zpowerlaw')
            axs[1].errorbar(y=df_plot_zpowe3['phoindex'].values, x=df_plot_zpowe3['rate'].values, xerr=df_plot_zpowe3['erate'].values,
                        yerr = (df_plot_zpowe3['phoindex_bot'].values, df_plot_zpowe3['phoindex_top'].values), fmt='.',color='b', markersize=7, elinewidth=1, capsize=2, capthick=1, ecolor='cornflowerblue', label='zpowerlaw')
            '''
            #Subplot 3: alpha vs rate (logpar)      
            axs[2].errorbar(y=df_plot_zlogp1['alpha'].values, x=df_plot_zlogp1['rate'].values, xerr=df_plot_zlogp1['erate'].values,
                        yerr = (df_plot_zlogp1['alpha_bot'].values, df_plot_zlogp1['alpha_top'].values), fmt='.', linestyle='',markersize=7, ecolor='lightcoral', elinewidth=1, capsize=2, capthick=1, color='red', label='zlogpar')
            axs[2].errorbar(y=df_plot_zlogp2['alpha'].values, x=df_plot_zlogp2['rate'].values, xerr=df_plot_zlogp2['erate'].values,
                        yerr = (df_plot_zlogp2['alpha_bot'].values, df_plot_zlogp2['alpha_top'].values), fmt='.', linestyle='',markersize=7, ecolor='yellowgreen', elinewidth=1, capsize=2, capthick=1, color='green', label='zlogpar')
            axs[2].errorbar(y=df_plot_zlogp3['alpha'].values, x=df_plot_zlogp3['rate'].values, xerr=df_plot_zlogp3['erate'].values,
                        yerr = (df_plot_zlogp3['alpha_bot'].values, df_plot_zlogp3['alpha_top'].values), fmt='.', linestyle='',markersize=7, ecolor='cornflowerblue', elinewidth=1, capsize=2, capthick=1, color='b', label='zlogpar')
            '''
            #Linestyle arrow 
            axs[1].quiver(df_plot_zpowe['rate'].values[:-1], df_plot_zpowe['phoindex'].values[:-1], df_plot_zpowe['rate'].values[1:]-df_plot_zpowe['rate'].values[:-1], df_plot_zpowe['phoindex'].values[1:]-df_plot_zpowe['phoindex'].values[:-1], scale_units='xy', angles='xy', scale=1, width=0.003, headwidth=8)
            #axs[2].quiver(df_plot_zlogp['rate'].values[:-1], df_plot_zlogp['alpha'].values[:-1], df_plot_zlogp['rate'].values[1:]-df_plot_zlogp['rate'].values[:-1], df_plot_zlogp['alpha'].values[1:]-df_plot_zlogp['alpha'].values[:-1], scale_units='xy', angles='xy', scale=1, width=0.003, headwidth=8)

            #Labels and stuff
            axs[0].grid()
            axs[0].set_xlabel('Time [MJD]')
            axs[0].set_ylabel('Rate [ct/s]')
            axs[1].grid()
            axs[1].set_xlabel('Rate [ct/s]')
            axs[1].set_ylabel('$\Gamma$ (powerlaw)')
            axs[2].grid()
            #axs[2].set_xlabel('Rate [ct/s]')
            #axs[2].set_ylabel(r'$\alpha$ (logparabola)')
            axs[0].ticklabel_format(useOffset=False)

            plt.savefig(os.path.join(target_dir, "Products", "Plots_spectra", f"hysteresis_{args.hysteresis}.png"))
            

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
                df_plot_zlogp = df_plot_zlogp.sort_values(by=['rate'])
                plt.errorbar(y=df_plot_zlogp['alpha'].values, x=df_plot_zlogp['flux'].values, xerr=(df_plot_zlogp['flux_bot'].values, df_plot_zlogp['flux_top'].values),
                            yerr = (df_plot_zlogp['alpha_bot'].values, df_plot_zlogp['alpha_top'].values), fmt='.', linestyle='', markersize=5, ecolor='cornflowerblue', elinewidth=1, capsize=2, capthick=1, color='b', label='zlogpar')
                plt.errorbar(y=df_plot_zpowe['phoindex'].values, x=df_plot_zpowe['flux'].values, xerr=(df_plot_zpowe['flux_bot'].values, df_plot_zpowe['flux_top'].values),
                            yerr = (df_plot_zpowe['phoindex_bot'].values, df_plot_zpowe['phoindex_top'].values), fmt='.', color='red', markersize=5, elinewidth=1, capsize=2, capthick=1, ecolor='lightcoral', label='zpowerlaw')
                plt.grid()
                plt.xlabel('Flux [cm$^{-2}$ erg s$^{-1}$]', fontsize=15)
                plt.ylabel('photon index', fontsize=15)
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

        if args.phoindex and args.beta:  #user wants alpha vs beta plot (only possible for logpar)
            data_spec_zlogp = data_spec_zlogp[data_spec_zlogp['phoindex_up']!=0.]
            df_plot_zlogp = pd.DataFrame({"rate": data_spec_zlogp['rate'].values, "erate": data_spec_zlogp['erate'].values, 
                                            'alpha': data_spec_zlogp['phoindex'].values, 
                                            'alpha_top': data_spec_zlogp['phoindex_up'].values - data_spec_zlogp['phoindex'].values,
                                            'alpha_bot': data_spec_zlogp['phoindex'].values - data_spec_zlogp['phoindex_low'].values,
                                            'beta': data_spec_zlogp['beta'].values, 
                                            'beta_top': data_spec_zlogp['beta_up'].values - data_spec_zlogp['beta'].values,
                                            'beta_bot': data_spec_zlogp['beta'].values - data_spec_zlogp['beta_low'].values,
                                            "obsid": data_spec_zlogp['obsid'].values})
            
            if args.state:  #separate in high and low state
                df_plot_zlogp_high = df_plot_zlogp[df_plot_zlogp['rate']>=20]
                df_plot_zlogp_low = df_plot_zlogp[df_plot_zlogp['rate']<20]  

            figure = plt.figure(figsize=(10,5))

            if args.state:
                
                plt.errorbar(y=df_plot_zlogp_high['alpha'].values, x=df_plot_zlogp_high['beta'].values, xerr=(df_plot_zlogp_high['beta_bot'].values, df_plot_zlogp_high['beta_top'].values),
                            yerr = (df_plot_zlogp_high['alpha_bot'].values, df_plot_zlogp_high['alpha_top'].values), fmt='.', markersize=3, ecolor='peachpuff', elinewidth=1, capsize=2, capthick=1, color='orange', label='High state')
                plt.errorbar(y=df_plot_zlogp_low['alpha'].values, x=df_plot_zlogp_low['beta'].values, xerr=(df_plot_zlogp_low['beta_bot'].values, df_plot_zlogp_low['beta_top'].values),
                            yerr = (df_plot_zlogp_low['alpha_bot'].values, df_plot_zlogp_low['alpha_top'].values),fmt='.', markersize=3, ecolor='mediumseagreen', elinewidth=1, capsize=2, capthick=1, color='g', label='Low state')
                plt.legend()

            else:
                plt.errorbar(y=df_plot_zlogp['alpha'].values, x=df_plot_zlogp['beta'].values, yerr=(df_plot_zlogp['alpha_bot'].values, df_plot_zlogp['alpha_top'].values),
                            xerr = (df_plot_zlogp['beta_bot'].values, df_plot_zlogp['beta_top'].values), fmt='.', markersize=5, ecolor='gray', elinewidth=1, capsize=2, capthick=1, color='black')
            
            
            #Correlation coefficient
            slope, intercept, r, p, stderr = linregress(df_plot_zlogp['beta'].values, df_plot_zlogp['alpha'].values)
            line = f'Correlation: {r:.2f}'
            plt.plot(df_plot_zlogp['beta'].values, intercept + slope * df_plot_zlogp['beta'].values, label=line, color='red', zorder=3)
            
            plt.grid()
            plt.legend()
            plt.xlabel(r'$\beta$', fontsize=15)
            plt.ylabel(r'$\alpha$', fontsize=15)
            plt.tight_layout()

            if args.state:
                plt.savefig(os.path.join(target_dir, "Products", "Plots_spectra", "alpha_vs_beta_state.png"))
            else:
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
                                            "phoindex_top68": data_spec_zlogp['phoindex_up68'].values - data_spec_zlogp['phoindex'].values,
                                            "phoindex_bot68": data_spec_zlogp['phoindex'].values - data_spec_zlogp['phoindex_low68'].values, 
                                            "beta_top68": data_spec_zlogp['beta_up68'].values - data_spec_zlogp['beta'].values,
                                            "beta_bot68": data_spec_zlogp['beta'].values - data_spec_zlogp['beta_low68'].values,
                                            "nH_up68": data_spec_zlogp['nH_up68'].values, "nH_low68": data_spec_zlogp['nH_low68'].values,
                                            "nH_top68": data_spec_zlogp['nH_up68'].values - data_spec_zlogp['nH'].values,
                                            "nH_bot68": data_spec_zlogp['nH'].values - data_spec_zlogp['nH_low68'].values,
                                            "obsid": data_spec_zlogp['obsid'].values})
            #Order dataframe and reset index
            df_plot_zlogp = df_plot_zlogp.sort_values(by=['tbinstart'])
            df_plot_zlogp = df_plot_zlogp[df_plot_zlogp['nH_up']!=0]
            df_plot_zlogp = df_plot_zlogp.reset_index(drop=True)
            df_plot_zlogp = df_plot_zlogp.reset_index()

            #Upper limits on nH
            df_plot_nH_pegged = df_plot_zlogp[df_plot_zlogp['nH_low']<=1e-4]
            df_plot_zlogp = df_plot_zlogp[df_plot_zlogp['nH_low']>1e-4]
            print(f"Mean and standard deviation of NH evolution: {np.mean(df_plot_zlogp['nH'].values)}, {np.std(df_plot_zlogp['nH'].values, ddof=1)}")

            fig = plt.figure(figsize=(13,6))
            plt.errorbar('index', 'RATE', 'ERROR', data=data_lc, ecolor='black', linestyle='', color='black')
            #Add vertical lines separating years
            plt.vlines(year_endpoints, 0, 60, colors='r', linestyles='solid')
            plt.grid(axis='y')
            plt.xlabel('Year', labelpad=0.1, fontsize=15)
            plt.ylabel('Rate [ct/s]', labelpad=0.1, fontsize=15)
            plt.margins(0)
            plt.tight_layout()
            plt.xticks(ticks=year_endpoints, labels=labels, rotation=50, fontsize=11)
            plt.savefig(os.path.join(target_dir, "Products", "Plots_spectra", "compact_LC_vlines.png"))

            #Plot panel            
            fig_logp, axs_logp = plt.subplots(4, 1, figsize=(13,11), sharex=True, gridspec_kw={'hspace':0, 'wspace':0})


                # Total lightcurve
            axs_logp[0].errorbar('index', 'RATE', 'ERROR', data=data_lc, ecolor='blue', linestyle='', color='b')
              
                #alpha vs time  
                    #98% confidence intervals
            axs_logp[1].errorbar('index', 'phoindex', yerr=(df_plot_zlogp['phoindex_bot'].values, df_plot_zlogp['phoindex_top'].values),
                                data=df_plot_zlogp, ecolor='lightgray', linestyle='', marker='.', markersize=3, color='b', capthick=1, elinewidth=1)
                    #68% confidence intervals 
            axs_logp[1].errorbar('index', 'phoindex', yerr=(df_plot_zlogp['phoindex_bot68'].values, df_plot_zlogp['phoindex_top68'].values),
                                data=df_plot_zlogp, ecolor='cornflowerblue', linestyle='', color='b', capthick=1, elinewidth=1)
                #beta vs time
                    #90 confidence intervals   
            axs_logp[2].errorbar('index', 'beta', yerr=(df_plot_zlogp['beta_bot'].values, df_plot_zlogp['beta_top'].values),
                                data=df_plot_zlogp, ecolor='lightgray', linestyle='', marker='.', markersize=3, color='b', capthick=1, elinewidth=1)
                    #68% confidence intervals
            axs_logp[2].errorbar('index', 'beta', yerr=(df_plot_zlogp['beta_bot68'].values, df_plot_zlogp['beta_top68'].values),
                                data=df_plot_zlogp, ecolor='cornflowerblue', linestyle='', color='b', capthick=1, elinewidth=1)
                #nH vs time
                    #%90
            axs_logp[3].errorbar('index', 'nH', yerr=(df_plot_zlogp['nH_bot'].values, df_plot_zlogp['nH_top'].values), data=df_plot_zlogp,
                                 ecolor='lightgray', linestyle='', marker='.', markersize=3, color='b', capthick=1, elinewidth=1)
                    #68%
            axs_logp[3].errorbar('index', 'nH', yerr=(df_plot_zlogp['nH_bot68'].values, df_plot_zlogp['nH_top68'].values), data=df_plot_zlogp,
                                 ecolor='cornflowerblue', linestyle='', color='b', capthick=1, elinewidth=1)
                    #upper limits 90% confidence
            axs_logp[3].errorbar('index', 'nH_up', data=df_plot_nH_pegged, uplims=True, linestyle='', marker='v', markersize=3, capthick=1, elinewidth=1, color='b', markeredgewidth=0.3, markeredgecolor='white')

            '''
            #Locate vaughan panel observations
            obs_vaughuan = [136541001, 158971201, 810860201, 411080301, 560980101, 791781401, 810860701, 791782001]
            for obs in obs_vaughuan:
                rgb = '#%06X' % random.randint(0, 0xFFFFFF)  #create random color
                df_obs = data_lc[data_lc['OBSERVATION']==obs]
                axs_logp[0].errorbar('index', 'RATE', 'ERROR', data=df_obs, ecolor='black', marker='.', markersize='5',linestyle='', color=rgb, label=obs)
                axs_logp[0].annotate(obs, xy=(df_obs['index'].values[30],df_obs['RATE'].values[30]), xytext=(df_obs['index'].values[0]+20,  df_obs['RATE'].values[0]-13) ,arrowprops=dict(arrowstyle='->',facecolor='black'))
            '''
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
            axs_logp[0].set_ylabel("Rate [ct/s]", fontsize=15)
            axs_logp[1].set_ylabel(r"$\alpha$", fontsize=15)
            axs_logp[2].set_ylabel(r"$\beta$", fontsize=15)
            axs_logp[3].set_ylabel("$N_H$ [$10^{22}$ cm$^{-2}$]", fontsize=15)
            axs_logp[3].set_xlabel("Year", fontsize=15)
            axs_logp[0].margins(0)
            axs_logp[1].margins(0)
            axs_logp[2].margins(0)
            axs_logp[3].margins(0)
            plt.xticks(ticks=year_endpoints, labels=labels, rotation=50, fontsize=11)
            plt.tight_layout()
            plt.savefig(os.path.join(target_dir, "Products", "Plots_spectra", "panel_logpar.png"))
            

        if args.powerlaw:

            data_spec_zpowe = data_spec_zpowe[data_spec_zpowe['phoindex_up']!=0]
            df_plot_powerlaw = pd.DataFrame({'tbinstart': data_spec_zpowe['tbinstart'].values, 'tbinstop': data_spec_zpowe['tbinstop'].values , 'phoindex': data_spec_zpowe['phoindex'].values,
                                            "phoindex_top": data_spec_zpowe['phoindex_up'].values - data_spec_zpowe['phoindex'].values,
                                            "phoindex_bot": data_spec_zpowe['phoindex'].values - data_spec_zpowe['phoindex_low'].values, 
                                            "nH": data_spec_zpowe['nH'].values, "nH_up": data_spec_zpowe['nH_up'].values, "nH_low": data_spec_zpowe['nH_low'].values,
                                            "nH_top": data_spec_zpowe['nH_up'].values - data_spec_zpowe['nH'].values,
                                            "nH_bot": data_spec_zpowe['nH'].values - data_spec_zpowe['nH_low'].values,
                                            "phoindex_top68": data_spec_zpowe['phoindex_up68'].values - data_spec_zpowe['phoindex'].values,
                                            "phoindex_bot68": data_spec_zpowe['phoindex'].values - data_spec_zpowe['phoindex_low68'].values, 
                                            "nH_up68": data_spec_zpowe['nH_up68'].values, "nH_low68": data_spec_zpowe['nH_low68'].values,
                                            "nH_top68": data_spec_zpowe['nH_up68'].values - data_spec_zpowe['nH'].values,
                                            "nH_bot68": data_spec_zpowe['nH'].values - data_spec_zpowe['nH_low68'].values,
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
            fig_pw, axs_pw = plt.subplots(3, 1, figsize=(13,10), sharex=True, gridspec_kw={'hspace':0, 'wspace':0})
                #90% confidence intervals
            axs_pw[0].errorbar('index', 'RATE', 'ERROR', data=data_lc, ecolor='b', linestyle='', color='b')
            axs_pw[1].errorbar('index', 'phoindex', yerr=(df_plot_powerlaw['phoindex_bot'].values, df_plot_powerlaw['phoindex_top'].values),
                                data=df_plot_powerlaw, ecolor='lightgray', linestyle='', marker='.', markersize=3, color='b', capthick=1, elinewidth=1)
                #68% confidence intervals
            axs_pw[1].errorbar('index', 'phoindex', yerr=(df_plot_powerlaw['phoindex_bot68'].values, df_plot_powerlaw['phoindex_top68'].values),
                                data=df_plot_powerlaw, ecolor='cornflowerblue', linestyle='', color='b', capthick=1, elinewidth=1)
            
            axs_pw[2].errorbar('index', 'nH', yerr=(df_plot_powerlaw['nH_bot'].values, df_plot_powerlaw['nH_top'].values), data=df_plot_powerlaw,
                                 ecolor='lightgray', linestyle='', marker='.', markersize=3, color='b', capthick=1, elinewidth=1)
            axs_pw[2].errorbar('index', 'nH', yerr=(df_plot_powerlaw['nH_bot68'].values, df_plot_powerlaw['nH_top68'].values), data=df_plot_powerlaw,
                                ecolor='cornflowerblue', linestyle='', color='b', capthick=1, elinewidth=1)
                #upper limits on nH 90% confidence
            axs_pw[2].errorbar('index', 'nH_up', data=df_plot_nH_pegged, uplims=True, linestyle='', marker='v', markersize=3, capthick=0.5, elinewidth=1, color='blue', markeredgewidth=0.3, markeredgecolor='white')

            #Add vertical lines separating years
            axs_pw[0].vlines(year_endpoints, 0, 60, colors='r', linestyles='solid')
            axs_pw[1].vlines(year_endpoints, 1.5, 3.7, colors='r', linestyles='solid')
            axs_pw[2].vlines(year_endpoints, -0.015, 0.145, colors='r', linestyles='solid')
            
            print(f"Mean and standard dev of NH evolution: {np.mean(df_plot_powerlaw['nH'].values)}, {np.std(df_plot_powerlaw['nH'].values, ddof=1)}")

            #Bellurie
            axs_pw[0].grid()
            axs_pw[1].grid()
            axs_pw[2].grid()
            axs_pw[0].set_ylabel("Rate [ct/s]", fontsize=15)
            axs_pw[1].set_ylabel("$\Gamma$", fontsize=15)
            axs_pw[2].set_ylabel("$N_H$ [$10^{22}$ cm$^{-2}$]", fontsize=15)
            axs_pw[2].set_xlabel("Year", fontsize=15)
            axs_pw[0].margins(0)
            axs_pw[1].margins(0)
            axs_pw[2].margins(0)
            plt.xticks(ticks=year_endpoints, labels=labels, rotation=50, fontsize=11)
            plt.tight_layout()
            plt.savefig(os.path.join(target_dir, "Products", "Plots_spectra", "panel_powerlaw.png"))
    

    if args.distribution:

        if args.logpar:

            # Separate data into low and high state (rate >20 cts/s)
            data_spec_zlogp = data_spec_zlogp[data_spec_zlogp['phoindex_up']!=0]
            data_spec_zlogp_high = data_spec_zlogp[data_spec_zlogp['rate']>=20]
            data_spec_zlogp_low = data_spec_zlogp[data_spec_zlogp['rate']<20]

            # Plot distribution for alpha and beta 
            fig, axs = plt.subplots(2, 2, figsize=(7,7), sharex=True, gridspec_kw={'hspace':0.1})
            
            #Alpha, high rate
            weights = np.ones_like(data_spec_zlogp_high['phoindex'].values) / len(data_spec_zlogp_high['phoindex'].values)
            y_high_alpha, x_high_alpha, p = axs[0,0].hist(data_spec_zlogp_high['phoindex'].values, bins=30, color='red', weights=weights)
            binsize = (x_high_alpha[1]-x_high_alpha[0])/2
            axs[0,0].vlines(x=x_high_alpha[np.argmax(y_high_alpha)]+binsize, ymin=0, ymax=y_high_alpha.max(), linewidth=1, color='black')
            free_text = f"max alpha: {round(x_high_alpha[np.argmax(y_high_alpha)]+binsize, 2)} $\pm$ {round(binsize, 2)}"
            text_box = AnchoredText(free_text, frameon=False, loc='lower left', pad=0.5)
            plt.setp(text_box.patch, facecolor='white', alpha=0.5)
            axs[0,0].add_artist(text_box)

            #Alpha, low rate
            weights = np.ones_like(data_spec_zlogp_low['phoindex'].values) / len(data_spec_zlogp_low['phoindex'].values)
            y_low_alpha, x_low_alpha, _ = axs[1,0].hist(data_spec_zlogp_low['phoindex'].values, bins=30, hatch='/', color='blue', weights=weights)
            binsize = (x_low_alpha[1]-x_low_alpha[0])/2            
            axs[1,0].vlines(x=x_low_alpha[np.argmax(y_low_alpha)]+binsize, ymin=0, ymax=y_low_alpha.max(), linewidth=1, color='black')
            free_text = f"max alpha: {round(x_low_alpha[np.argmax(y_low_alpha)]+binsize, 2)} $\pm$ {round(binsize, 2)}"
            text_box = AnchoredText(free_text, frameon=False, loc='lower left', pad=0.5)
            plt.setp(text_box.patch, facecolor='white', alpha=0.5)
            axs[1,0].add_artist(text_box)   

            #K-S alpha
            sample1 = data_spec_zlogp_high['phoindex'].values
            sample2 = data_spec_zlogp_low['phoindex'].values
            print('K-S for alpha (logpar): ',ks_2samp(sample1, sample2))

            #Beta, high rate
            weights = np.ones_like(data_spec_zlogp_high['beta'].values) / len(data_spec_zlogp_high['beta'].values)
            y_high_beta, x_high_beta, _ = axs[0,1].hist(data_spec_zlogp_high['beta'].values, bins=30, color='red', weights=weights)
            binsize = (x_high_beta[1]-x_high_beta[0])/2
            axs[0,1].vlines(x=x_high_beta[np.argmax(y_high_beta)]+binsize, ymin=0, ymax=y_high_beta.max(), linewidth=1, color='black')
            free_text = f"max beta: {round(x_high_beta[np.argmax(y_high_beta)]+binsize, 2)} $\pm$ {round(binsize, 2)}"
            text_box = AnchoredText(free_text, frameon=False, loc='lower right', pad=0.5)
            plt.setp(text_box.patch, facecolor='white', alpha=0.5)
            axs[0,1].add_artist(text_box)

            #Beta, low rate
            weights = np.ones_like(data_spec_zlogp_low['beta'].values) / len(data_spec_zlogp_low['beta'].values)
            y_low_beta, x_low_beta, _ = axs[1,1].hist(data_spec_zlogp_low['beta'].values, bins=30, hatch='/', color='blue', weights=weights)
            binsize = (x_low_beta[1]-x_low_beta[0])/2
            axs[1,1].vlines(x=x_low_beta[np.argmax(y_low_beta)]+binsize, ymin=0, ymax=y_low_beta.max(), linewidth=1, color='black')
            free_text = f"max beta: {round(x_low_beta[np.argmax(y_low_beta)]+binsize, 2)} $\pm$ {round(binsize, 2)}"
            text_box = AnchoredText(free_text, frameon=False, loc='lower right', pad=0.5)
            plt.setp(text_box.patch, facecolor='white', alpha=0.5)
            axs[1,1].add_artist(text_box)
            
            #K-S for beta
            sample1 = data_spec_zlogp_high['beta'].values
            sample2 = data_spec_zlogp_low['beta'].values
            print('K-S for beta (logpar): ',ks_2samp(sample1, sample2))

            red_patch =  Patch(facecolor='red',  edgecolor='black', label='high state')
            blue_patch =  Patch(facecolor='b', hatch='/', edgecolor='black', label='low state')
            leg = axs[0,1].legend(handles=[red_patch, blue_patch], labelspacing=1.5, handlelength=4, loc='upper right',  fancybox=True, shadow=True)
            for patch in leg.get_patches():
                patch.set_height(16)
                patch.set_y(-6)
            axs[1,0].set_xlabel(r'$\alpha$ (logparabola)')
            axs[1,0].set_ylabel('Normalized counts')
            axs[0,0].set_ylabel('Normalized counts')
            axs[1,1].set_xlabel(r'$\beta$ (logparabola)')
            plt.savefig(os.path.join(target_dir, "Products", "Plots_spectra", "distribution_logpar.png"))

        if args.powerlaw:

            # Separate data into low and high state (rate >20)
            data_spec_zpowe = data_spec_zpowe[data_spec_zpowe['phoindex_up']!=0]
            data_spec_zpowe_high = data_spec_zpowe[data_spec_zpowe['rate']>=20]
            data_spec_zpowe_low = data_spec_zpowe[data_spec_zpowe['rate']<20]

            #Plot distribution only for photon index
            fig, axs = plt.subplots(2, 1, figsize=(6,7), sharex=True, gridspec_kw={'hspace':0.1})
            
            #High rate
            weights = np.ones_like(data_spec_zpowe_high['phoindex'].values) / len(data_spec_zpowe_high['phoindex'].values)
            y_high, x_high, _ = axs[0].hist(data_spec_zpowe_high['phoindex'].values, bins=30, color='red', weights=weights)
            binsize = (x_high[1]-x_high[0])/2
            axs[0].vlines(x=x_high[np.argmax(y_high)]+binsize, ymin=0, ymax=y_high.max(), color='black', linewidth=1)
            free_text = f"max $\Gamma$: {round(x_high[np.argmax(y_high)]+binsize, 2)} $\pm$ {round(binsize, 2)}"
            text_box = AnchoredText(free_text, frameon=False, loc='upper left', pad=0.5)
            plt.setp(text_box.patch, facecolor='white', alpha=0.5)
            axs[0].add_artist(text_box)

            #Low rate
            weights = np.ones_like(data_spec_zpowe_low['phoindex'].values) / len(data_spec_zpowe_low['phoindex'].values)
            y_low, x_low, _ = axs[1].hist(data_spec_zpowe_low['phoindex'].values, bins=30, color='blue',hatch='/', weights=weights)
            binsize = (x_low[1]-x_low[0])/2
            axs[1].vlines(x=x_low[np.argmax(y_low)]+binsize, ymin=0, ymax=y_low.max(), color='black', linewidth=1)
            free_text = f"max $\Gamma$: {round(x_low[np.argmax(y_low)]+binsize, 2)} $\pm$ {round(binsize, 2)}"
            text_box = AnchoredText(free_text, frameon=False, loc='upper left', pad=0.5)
            plt.setp(text_box.patch, facecolor='white', alpha=0.5)
            axs[1].add_artist(text_box)
            
            #K-S for gamma 
            sample1 = data_spec_zpowe_high['phoindex'].values
            sample2 = data_spec_zpowe_low['phoindex'].values
            print('K-S for gamma (powerlaw): ',ks_2samp(sample1, sample2))

            red_patch =  Patch(facecolor='red', edgecolor='black', label='high state')
            blue_patch =  Patch(facecolor='b', hatch='/', edgecolor='black', label='low state')
            leg = axs[0].legend(handles=[red_patch, blue_patch], labelspacing=1.5, handlelength=4, loc='upper right',  fancybox=True, shadow=True)
            for patch in leg.get_patches():
                patch.set_height(16)
                patch.set_y(-6)           
            axs[1].set_xlabel('$\Gamma$ (powerlaw)')
            axs[0].set_ylabel('Normalized counts')
            axs[1].set_ylabel('Normalized counts')
            plt.savefig(os.path.join(target_dir, "Products", "Plots_spectra", "distribution_powerlaw.png"))
