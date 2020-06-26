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
                    help="make phoindex  plot")
parser.add_argument('--beta', action='store_true',
                    help='make beta plot')
parser.add_argument('--fvar', action='store_true',
                    help='make fvar plot')
parser.add_argument('--rate', action='store_true',
                    help='make rate plot')
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
            fig, axs =plt.subplots(1, 1, figsize=(10,10), gridspec_kw={'hspace':0.4})
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
            df_plot_zlogp = pd.DataFrame({"fvar": data_lc['F_var'].values, "fvar_err": data_lc['F_var_sigma'].values, 'beta': data_spec_zlogp['beta'].values, 
                                          'beta_top': data_spec_zlogp['beta_up'].values - data_spec_zlogp['beta'].values,
                                          'beta_bot': data_spec_zlogp['beta'].values - data_spec_zlogp['beta_low'].values,
                                          "obsid": data_lc['ObsId'].values})
            df_plot_zlogp = df_plot_zlogp[df_plot_zlogp['fvar'] != -1.]

            figure = plt.figure(figsize=(10,10))
            plt.errorbar(df_plot_zlogp['beta'].values, df_plot_zlogp['fvar'].values, yerr=df_plot_zlogp['fvar_err'].values,
                        xerr = (df_plot_zlogp['beta_bot'].values, df_plot_zlogp['beta_top'].values), color='b', linestyle='', marker='.', markersize=0.5, ecolor='b', label='zlogpar')
            plt.grid()
            plt.xlabel('beta', fontsize=15)
            plt.ylabel('$F_{var}$', fontsize=15)
            plt.show()
    
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
    
    if not args.fvar and not args.rate:
        if args.phoindex and if args.beta:
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

        for directory in os.listdir(target_dir):
            if directory.startswith('0'):
                os.chdir(f"{target_dir}/{directory}/rgs")

                for filename in glob.glob('*_RGS_rates.ds'):
                    x, y, yerr, fracexp, y_bg, yerr_bg = mask_fracexp15(filename)
                    total_lightcurve_rates.extend(y)
                    total_lightcurve_errates.extend(yerr)
                    total_lightcurve_times.extend(x)
            
        total_lightcurve_rates = np.asarray(total_lightcurve_rates)
        total_lightcurve_errates = np.asarray(total_lightcurve_errates)
        total_lightcurve_times = np.asarray(total_lightcurve_times)
        
        data_lc = pd.DataFrame({"RATE":total_lightcurve_rates, "TIME":total_lightcurve_times, "ERROR":total_lightcurve_errates})
        data_lc = data_lc.sort_values(by=['TIME'])
        data_lc = data_lc.reset_index(drop=True)
        data_lc = data_lc.reset_index()

        if args.logpar:
            pass

        if args.powerlaw:
            df_plot_powerlaw = pd.DataFrame('tbinstart': data_spec_zpowe['tbinstart'].values, 'tbinstop:' data_spec_zpowe['tbinstop'].values , 'phoindex': data_spec_zpowe['phoindex'].values,
                                            "phoindex_top": data_spec_zpowe['phoindex_up'].values - data_spec_zpowe['phoindex'].values,
                                            "phoindex_bot": data_spec_zpowe['phoindex'].values - data_spec_zpowe['phoindex_low'].values, 
                                            "nH": data_spec_zpowe['nH'].values, "nH_up": data_spec_zpowe['nH_up'].values, "nH_low": data_spec_zpowe['nH_low'].values,
                                            "nH_top": data_spec_zpowe['nH_up'].values - data_spec_zpowe['nH'].values,
                                            "nH_bot": data_spec_zpowe['nH'].values - data_spec_zpowe['nH_low'].values)
            
            #Order dataframe and reset index
            df_plot_powerlaw = df_plot_powerlaw.sort_values(by=['tbinstart')
            df_plot_powerlaw = df_plot_powerlaw.reset_index(drop=True)
            df_plot_powerlaw = df_plot_powerlaw.reset_index()

            #Upper limits on nH
            df_plot_nH_pegged = df_plot_powerlaw[df_plot_powerlaw['nH_low']<=1e-4]
            df_plot_powerlaw = df_plot_powerlaw[df_plot_powerlaw['nH_low']>1e-4]

            #Plot panel
            fig_pw, axs_pw = plt.subplots(3, 1, figsize=(15,20), sharex=True, gridspec_kw={'hspace':0, 'wspace':0})
            
            axs_pw[0].errorbar('index', 'RATE', 'ERROR', data=data_lc, fmt='.', ecolor='gray', elinewidth=1, capsize=2, capthick=1)
            axs_pw[1].errorbar('index', 'phoindex', yerr=(df_plot_powerlaw['phoindex_bot'].values, df_plot_powerlaw['phoindex_top'].values),
                                data=df_plot_powerlaw, fmt='.', ecolor='gray', elinewidth=1, capsize=2, capthick=1)
            axs_pw[2].errorbar('index', 'nH', yerr=(df_plot_powerlaw['nH_bot'].values, df_plot_powerlaw['nH_top'].values), data=df_plot_powerlaw,
                                 fmt='.', ecolor='gray', elinewidth=1, capsize=2, capthick=1)
            axs_pw[2].errorbar('index', 'nH_up', data=df_plot_nH_pegged, uplims=True)

            axs_pw[0].grid()
            axs_pw[1].grid()
            axs_pw[2].grid()
            axs_pw[0].set_ylabel("Rate [ct/s]")
            axs_pw[1].set_ylabel("phoindex")
            axs_pw[2].set_ylabel("nH [$10^22$ g/cm$^2$]")
