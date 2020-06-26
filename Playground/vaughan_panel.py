import os
from config import CONFIG
from tools import *
import glob
from astropy.table import Table
import numpy as np
from itertools import compress
from scipy.optimize import curve_fit
import pandas as pd
import sys
logging.basicConfig(level=logging.ERROR)

target_dir = CONFIG['target_dir'] 
sas_dir = CONFIG['sas_dir']
ccf_dir = CONFIG['ccf_dir']
observations = [ '0658801801','0136541001', '0150498701', '0791782001', '0560980101'] 
observations = ['0136541001']
N = 10
M = 10
def linear_fit(x, q):
    return 0*x+q

def synchronous_times(tstart1, tstart2, tstop1, tstop2):
    """
    Function that compares tstart1 and tstart2, and tstop1 and tstop2. The output is a list of two values:
    the maximum between tstart1 and tstart2 and the minimum between tstop1 and stop2.
    """    
    final_start = max(tstart1, tstart2)
    final_stop = min(tstop1, tstop2)
    overlap = max(0., final_stop-final_start)
    try:
        if overlap>0:
            return final_start, final_stop
        else:
            raise Exception
    except Exception as e:
        logging.error("The given exposures do not overlap. Please check the if the input exposures are correct.")

def binning(N, bintime, dataframe, colx, coly):
    
    segment = N*bintime
    x_in = dataframe[colx].values[0]
    x_fin = dataframe[colx].values[-1]
    x = x_in
    mean_x = []
    mean_y = []
    mean_yerr = []
    mean_xerr = []

    while(x + segment < x_fin):
        segment_df = dataframe[(dataframe[colx]<x+segment) & (dataframe[colx]>x)]
        n_in_segment = len(segment_df)
        if n_in_segment <= N/3:
            
            x += segment
            continue
        else:
            mean_x.append(np.mean(segment_df[colx].values))
            mean_y.append(np.mean(segment_df[coly].values))
            mean_yerr.append(np.std(segment_df[coly].values)/np.sqrt(len(segment_df[coly].values)))
            mean_xerr.append((segment_df[colx].values[0]- segment_df[colx].values[-1])/2.)
            x += segment

    return mean_x, mean_y, mean_xerr, mean_yerr

            

for observation in observations:
    os.chdir(f'{target_dir}/{observation}')
    setupSAS(sas_dir, ccf_dir)
    os.chdir(f'{target_dir}/{observation}/rgs')
    rgsevlists = glob.glob('*EVENLI0000.FIT')
    rgssrclists = glob.glob('*SRCLI_0000.FIT')
    pairs_events = sort_rgs_list(rgsevlists, 'expo_number')
    pairs_srcli = sort_rgs_list(rgssrclists, "expo_number")

    with fits.open(pairs_events[0][0]) as hdul:
        expid0 = hdul[1].header['EXP_ID']
        tstart0 = hdul[1].header['TSTART']
        tstop0 = hdul[1].header['TSTOP']
    print(f'RGS1 exposure {expid0} start and stop times:', tstart0,'-', tstop0)

    # Open RGS2 evelist and save tstart and tstop
    with fits.open(pairs_events[0][1]) as hdul:
        expid1 = hdul[1].header['EXP_ID']
        tstart1 = hdul[1].header['TSTART']
        tstop1 = hdul[1].header['TSTOP']
    print(f'RGS2 exposure {expid1} start and stop times:', tstart1,'-', tstop1)

    start_time, stop_time = synchronous_times(tstart0, tstart1, tstop0, tstop1)
    print('Final start and stop time for rgslccorr:', start_time, stop_time)

    #Run rgslccorr
    timebinsize = 25 #s
    '''
    logging.info(f"Running rgslccorr SAS command for observation number {observation}.")
    rgslc_command = f"rgslccorr evlist='{pairs_events[0][0]} {pairs_events[0][1]}' srclist='{pairs_srcli[0][0]} {pairs_srcli[0][1]}' withbkgsubtraction=yes timemin={start_time} timemax={stop_time} timebinsize={timebinsize} orders='1' sourceid=3 outputsrcfilename={observation}_RGS_rates_{timebinsize}bin.ds outputbkgfilename={observation}_bkg_rates_{timebinsize}bin.ds"
    status_rgslc = run_command(rgslc_command)
    '''
    #Read LC data
    time, rate, erate, fracexp, backv, backe = mask_fracexp15(f"{target_dir}/{observation}/rgs/{observation}_RGS_rates_{timebinsize}bin.ds")
    data = pd.DataFrame({'RATE': rate, "TIME": time, "ERROR": erate, 'BACKV': backv, "BACKE": backe})
    
    ## PANEL 
    fig, axs = plt.subplots(6, 1, figsize=(15,20), sharex=True, gridspec_kw={'hspace':0})
    fig.suptitle(f'RGS Lightcurve ObsId {observation} binsize {timebinsize}s \n N ={N}, M = {M}', fontsize=15, y=0.92)
    
    #Subplot x (lightcurve)
    axs[0].errorbar(data['TIME'].values, data['RATE'].values, yerr=data['ERROR'].values, linestyle='', color='black', marker='.', ecolor='gray', 
                label=f'RGS Lightcurve ObsId {observation} binsize {timebinsize}s ')
    axs[0].grid(True)
    axs[0].set_ylabel('x', fontsize=10)

    #Subplot <x> (mean lightcurve)
    mean_time, mean_data, mean_time_err, mean_error_data = binning(M, 25, data, 'TIME', 'RATE' )

    axs[1].errorbar(mean_time, mean_data, yerr=mean_error_data, color='black', linestyle='', marker='.', ecolor='gray')
    axs[1].grid()
    axs[1].set_ylabel('<x>', fontsize=10)

    #Subplot excess variance
    i = 0
    xs_arr = []
    xs_err_arr = []
    while (i+M < len(data)):
        segment = data[i:i+M]
        xs,err_xs = excess_variance(segment['RATE'].values, segment['ERROR'].values, normalized=False)
        xs_arr.append(xs)
        xs_err_arr.append(err_xs)

        i=i+M

    mask_negative = []
    for el in xs_arr:
        if el<0:
            mask_negative.append(False)
        else:
            mask_negative.append(True)

    xs_arr = list(compress(xs_arr,mask_negative))
    xs_err_arr = list(compress(xs_err_arr, mask_negative))
    mean_time_nonneg = list(compress(mean_time, mask_negative))
    
    axs[2].errorbar(mean_time_nonneg, xs_arr, xs_err_arr, color='black', marker='.', linestyle='', ecolor='gray' )
    axs[2].grid()
    axs[2].set_ylabel('$\sigma_{XS}^2$', fontsize=10)


    #Subplot mean excess variance 
    df_mean_xs = pd.DataFrame({'time': mean_time_nonneg, 'xs': xs_arr, 'xs_err': xs_err_arr})
    meanx2_times, mean_xs, meanx2_times_err, mean_xs_err = binning(N, M*timebinsize, df_mean_xs, 'time', 'xs')
    
    axs[3].errorbar(meanx2_times, mean_xs, mean_xs_err, xerr=meanx2_times_err,  linestyle='', color='black', marker='.', ecolor='gray')
    axs[3].grid()
    axs[3].set_ylabel('$<\sigma_{XS}^2>$', fontsize=10)

    #Subplot F_var
    i = 0
    fvar_arr = []
    fvar_err_arr = []
    while(i+M<len(data)):
        segment = data[i:i+M]
        fvar, fvar_err = fractional_variability(segment['RATE'].values, segment['ERROR'].values, segment['BACKV'].values, segment['BACKE'].values, netlightcurve=True)
        fvar_arr.append(fvar)
        fvar_err_arr.append(fvar_err)
        i=i+M
    
    mask_fvar= []
    for el in fvar_arr:
        if el==-1.:
            mask_fvar.append(False)
        else:
            mask_fvar.append(True)
    fvar_arr = list(compress(fvar_arr, mask_fvar))
    fvar_err_arr = list(compress(fvar_err_arr, mask_fvar))
        
    axs[4].errorbar(mean_time_nonneg, fvar_arr, fvar_err_arr, linestyle='', color='black', marker='.', ecolor='gray')
    axs[4].grid()
    axs[4].set_ylabel('$F_{var}$', fontsize=10)

    # Subplot mean Fvar
    df_mean_fvar = pd.DataFrame({'time': mean_time_nonneg, 'fvar': fvar_arr, 'fvar_err': fvar_err_arr})
    meanx2_times, fvar_mean_arr, meanx2_times_err, fvar_err_mean_arr = binning(N, M*timebinsize, df_mean_fvar, 'time', 'fvar')

    axs[5].errorbar(meanx2_times, fvar_mean_arr, fvar_err_mean_arr, xerr=meanx2_times_err, linestyle='', color='black', marker='.', ecolor='gray')
    axs[5].grid()
    axs[5].set_xlabel('TIME [s]', fontsize=10)
    axs[5].set_ylabel('$<F_{var}>$', fontsize=10)

    plt.savefig(f'{target_dir}/Products/Plots_timeseries/{observation}_variability_panel.png')

    #Fit 
    fvar_mean_arr = np.array(fvar_mean_arr)
    meanx2_times = np.array(meanx2_times)
    meanx2_times_err = np.array(meanx2_times_err)
    fvar_err_mean_arr = np.array(fvar_err_mean_arr)
    
    #meanx2_times = meanx2_times - meanx2_times[0]
    # Fit constant
    #valori iniziali e fit
    initial_values =(0.03)
    pars, covm = curve_fit(linear_fit, meanx2_times, fvar_mean_arr, initial_values, fvar_err_mean_arr) 

    q0 = pars    #a0 e b0 sono i valori dei parametri della curva del fit
    dq = np.sqrt(covm.diagonal())   #e i loro rispettivi errori (dalla matrice di covaranza covm)

    #stampa risultati e plot 
    #print('m = %f +- %f' % (m0,dm))
    print('q = %f +- %f' % (q0, dq))
    
    #chi2
    chisq =(((fvar_mean_arr-linear_fit(meanx2_times, q0) )/fvar_err_mean_arr)**2).sum()
    ndof = len(meanx2_times) - 1
    print('Chisquare/ndof = %f/%d' % (chisq, ndof))


    '''
    #PLOT CONFIDENCE INTERVALS
    fig_confidence = plt.figure(figsize=(10,10))
    plt.errorbar(mean_time_nonneg, xs_arr, xs_err_arr, color='black', marker='.', linestyle='', ecolor='gray' )
    plt.grid()
    plt.ylabel('$<\sigma_{XS}^2>$', fontsize=10)
    plt.xlabel('Time (s)')
    
    upper_90 = 0.46
    lower_90 = 0.78
    print(np.mean(xs_arr))
    plt.hlines(np.mean(xs_arr), plt.xlim()[0], plt.xlim()[1],color='black', linestyle='-')

    plt.show()
    '''

    #PLOT RATE CORRELATION

    #Sort xs_arr by rate
    i = 0
    mean_rate = []
    mean_rate_err = []
    xs_arr_rate = []
    xs_arr_err_rate = []

    while(i+M<len(data)):
        mean_rate.append(np.mean(data[i:i+M]['RATE'].values))
        xs,err_xs = excess_variance(data[i:i+M]['RATE'].values, data[i:i+M]['ERROR'].values, normalized=False)
        xs_arr_rate.append(xs)
        xs_arr_err_rate.append(err_xs)
        i+=M

    xs_sorted = pd.DataFrame({'rate': mean_rate, 'xs': xs_arr_rate, 'xs_err': xs_arr_err_rate})
    xs_sorted = xs_sorted.dropna()
    xs_sorted = xs_sorted.sort_values(by=['rate'])
    xs_sorted = xs_sorted[xs_sorted['xs']>0.]    
    
    #Binning
    print(xs_sorted)
    meanx2_rate, meanx2_xs, meanx2_rate_err, meanx2_xs_err = binning(N, 0.01, xs_sorted, 'rate', 'xs')
    '''
    
    i=0
    meanx2_rate = []
    meanx2_xs = []
    meanx2_xs_err = []
    while(i+N<len(xs_sorted)):
        meanx2_rate.append(np.mean(xs_sorted[i:i+N]['rate'].values))
        meanx2_xs.append(np.mean(xs_sorted[i:i+N]['xs'].values))
        meanx2_xs_err.append(np.std(xs_sorted[i:i+N]['xs'].values)/np.sqrt(len(xs_sorted[i:i+N])))
        #meanx2_xs_err.append(np.sqrt(1/ (1/np.square(xs_sorted[i:i+M]['xs_err'].values)).sum()))
        i+=N
    '''
    fig_rate, ax = plt.subplots(1,1, figsize=(10,10))
    ax.errorbar(x=meanx2_rate, y=meanx2_xs, yerr=meanx2_xs_err, linestyle='', marker='.')
    ax.set_xlabel('x')
    ax.set_ylabel('<$\sigma_{XS}^2$>')
    ax.set_xlim(0, max(meanx2_rate)+1)
    plt.savefig(f'{target_dir}/Products/Plots_timeseries/{observation}_xs_rate.png')
    
    
