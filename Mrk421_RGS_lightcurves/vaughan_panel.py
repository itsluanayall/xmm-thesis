import os
from config import CONFIG
from tools import *
import glob
from astropy.table import Table
import numpy as np
from itertools import compress


target_dir = CONFIG['target_dir'] 
sas_dir = CONFIG['sas_dir']
ccf_dir = CONFIG['ccf_dir']
observations = ['0150498701']

for observation in observations:
    os.chdir(f'{target_dir}/{observation}')
    setupSAS(sas_dir, ccf_dir)
    os.chdir(f'{target_dir}/{observation}/rgs')
    rgsevlists = glob.glob('*EVENLI0000.FIT')
    rgssrclists = glob.glob('*SRCLI_0000.FIT')
    pairs_events = sort_rgs_list(rgsevlists, 'expo_number')
    pairs_srcli = sort_rgs_list(rgssrclists, "expo_number")

    #Run rgslccorr
    timebinsize = 25 #s
    '''
    logging.info(f"Running rgslccorr SAS command for observation number {observation}.")
    rgslc_command = f"rgslccorr evlist='{pairs_events[0][0]} {pairs_events[0][1]}' srclist='{pairs_srcli[0][0]} {pairs_srcli[0][1]}' withbkgsubtraction=yes timebinsize={timebinsize} orders='1' sourceid=3 outputsrcfilename={observation}_RGS_rates_{timebinsize}bin.ds outputbkgfilename={observation}_bkg_rates_{timebinsize}bin.ds"
    status_rgslc = run_command(rgslc_command)
    '''
    #Read LC data
    hdul = Table.read(f"{target_dir}/{observation}/rgs/{observation}_RGS_rates_{timebinsize}bin.ds", hdu=1)    
    data = hdul.to_pandas()
    data = data.dropna()
    
    ## PANEL 
    fig, axs = plt.subplots(6, 1, figsize=(15,20), sharex=True, gridspec_kw={'hspace':0})
    fig.suptitle(f'RGS Lightcurve ObsId {observation} binsize {timebinsize}s', fontsize=15, y=0.92)
    
    #Subplot x (lightcurve)
    axs[0].errorbar(data['TIME'].values, data['RATE'].values, yerr=data['ERROR'].values, linestyle='', color='black', marker='.', ecolor='gray', 
                label=f'RGS Lightcurve ObsId {observation} binsize {timebinsize}s ')
    axs[0].grid(True)
    axs[0].set_ylabel('x', fontsize=10)

    #Subplot <x> (mean lightcurve)
    mean_data = []
    mean_error_data = []
    mean_time = []
    mean_time_err = []
    i = 0
    M = 20
    while(i+M<len(data)):
        mean_data.append(np.mean(data[i:i+M]['RATE'].values))
        mean_error_data.append(np.sqrt(1/ (np.sum(1/np.square(data[i:i+M]['ERROR'].values)))))
        mean_time.append(np.mean([data['TIME'].loc[i], data['TIME'].loc[i+M]]))
        i+=M

    for i in range(len(mean_time)):
        mean_time_err.append(np.std(mean_time)/np.sqrt(len(mean_time)))
    
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
    print(len(xs_arr))

    #Subplot mean excess variance 
    i = 0 
    mean_xs = []
    mean_xs_err = []
    meanx2_times = []
    meanx2_times_err = []
    N=17
    while (i+N< len(xs_arr)):
        mean_xs.append(np.mean(xs_arr[i:i+N]))
        mean_xs_err.append(np.sqrt(1/ (np.sum(1/np.square(xs_err_arr[i:i+M])))))
        meanx2_times.append(np.mean([mean_time_nonneg[i], mean_time_nonneg[i+N]]))
        i+=N
    for i in range(len(meanx2_times)):
        meanx2_times_err.append(np.std(meanx2_times)/np.sqrt(len(meanx2_times)))
    print(len(meanx2_times_err))
    print(len(meanx2_times))

    axs[3].errorbar(meanx2_times, mean_xs, mean_xs_err, xerr=meanx2_times_err, linestyle='', color='black', marker='.', ecolor='gray')
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

    # Subplot maean Fvar
    i = 0
    fvar_mean_arr = []
    fvar_err_mean_arr = []
    while(i+N<len(fvar_arr)):
        fvar_mean_arr.append(np.mean(fvar_arr[i:i+N]))
        fvar_err_mean_arr.append(np.sqrt(1/ (np.sum(1/np.square(fvar_err_arr[i:i+M])))))
        i+=N
    
    axs[5].errorbar(meanx2_times, fvar_mean_arr, fvar_err_mean_arr, xerr=meanx2_times_err, linestyle='', color='black', marker='.', ecolor='gray')
    axs[5].grid()
    axs[5].set_xlabel('TIME [s]', fontsize=10)
    axs[5].set_ylabel('$<F_{var}>$', fontsize=10)

    plt.savefig(f'{target_dir}/Products/{observation}_variability_panel.png')
    plt.show()


