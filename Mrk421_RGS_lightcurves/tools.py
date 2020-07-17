import subprocess
import sys
import os
import numpy as np
import logging
from astropy.io import fits
import xspec
import matplotlib.pyplot as plt

def setupSAS(sas_dir, ccf_dir):
    """
    Set up the fundamental environment variables for a new session of SAS: the directory in which SAS is installed
    and the directory in which the CCF (Current Calibration Files) have been downloaded.
    To update the valid CCFs go to the following website: https://www.cosmos.esa.int/web/xmm-newton/current-calibration-files
    """
    logging.info('Setting up the environment variables for SAS...')
    os.environ["SAS_DIR"] = sas_dir
    os.environ["SAS_PATH"] = os.environ["SAS_DIR"]
    os.environ["SAS_CCFPATH"] = ccf_dir
    run_command('sasversion')


def run_command(command):
    """
    Execute a shell command. If an error occurs, it will be reported in the STDOUT and the exception will be handled,
    so that the program will keep on running.
    """
    try:
        result = subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        retcode = result.returncode
        logging.debug("Execution of {} returned {}, \n {}".format(command,retcode,result.stdout.decode()))
        return retcode

    except subprocess.CalledProcessError as e:
        print(f'\033[91m Execution of {command} failed: \033[0m')
        logging.error(e.output)
        return 1
    

def split_rgs_filename(filename):
    """
    Splits the RGS filname in substrings, each of which indicate a specific characteristic of the observation.
    Returns the dictionary of the substrings.
    """
    n = [1, 10, 2, 1, 3, 6, 1, 3]  #subsets of the filename string
    split_name = (([filename[sum(n[:i]):sum(n[:i+1])] for i in range(len(n))]))
    
    instr = split_name[2]
    exposure = split_name[3]
    expo_number = split_name[4]
    type_file = split_name[5]
    order = split_name[6]
    source_number = split_name[7]
    
    return {"instr":instr, "exposure": exposure, "expo_number":expo_number, "type_file":type_file,
            "order":order, "source_number":source_number, "filename": filename}


def sort_rgs_list(l, variable='expo_number'):
    """
    This method is useful when there are more than two RGS exposures in the observation.
    If this is the case, the exposures must be paired accordingly to their exposure number in order to
    correctly use the rgslccorr SAS command, which extracts the background-subtracted lightcurves. 
    Given a list of rgs products of the same type (e.g. 'EVENLI' or 'SRCLI_') this function sorts the list according
    to the 'variable' argument so that the list can be divided into pairs ready for the generation of the lightcurves.
    The variable argument is the parameter according to which the list will be sorted, e.g. "expo_number".
    """
    
    #Make a list of dictionaries associated to each event
    l_dict = []
    for ele in l:
        l_dict.append(split_rgs_filename(ele))
    
    #Sort the list according to the given variable
    sorted_l_dict = sorted(l_dict, key=lambda k: k[variable])
    
    #Reassemble the filenames and split the sorted list into pairs made up of a RGS1 and a RGS2
    #with consecutive exposure numbers
    sorted_l = []
    for ele in sorted_l_dict:
        sorted_l.append(ele['filename'])
    
    return [sorted_l[x:x+2] for x in range(0, len(sorted_l), 2)]


class RangeException(Exception):
    """
    Custom Exception to check if array lengths are the same.
    """
    def __init__(self, *args):
        if args:
            self.message = args[0]
        else:
            self.message = None

    def __str__(self):
        if self.message:
            return 'RangeException, {0} '.format(self.message)
        else:
            return 'RangeException has been raised.'


class NoDataException(Exception):
    """
    Custom Exception to check if there are enough datapoints.
    """
    def __init__(self, *args):
        if args:
            self.message = args[0]
        else:
            self.message = None

    def __str__(self):
        if self.message:
            return 'NoDataException, {0} '.format(self.message)
        else:
            return 'NoDataException has been raised.'


def excess_variance(rates, errrates, normalized=True):
    """
    Calculates the excess variance from the rates given as argument. 
    If normalized = True, the function returns the normalized excess variance and its error.
    
    :param rates: rate lightcurve [ct/s]
    :param errrates: error on rates
    :param normalized: if set to True, returns the normalized excess variance otherwise returns the excess variance.
    """

    try:
        mean = np.mean(rates)
        if mean<0:
            raise ValueError('Negative count rates in Input File.')
    except ValueError as e:
        loggng.error(e)

    #Calculation of excess variance beginning from the rates (definitions taken from Aggrawal et al., 2018)
    variance = np.var(rates, ddof=1)     #variance
    N = len(rates)
    mse = np.mean(np.square(errrates))   #mean square error
    xs = variance - mse                 #excess variance
    nxs = xs/(np.square(mean))          #normalized excess variance
    f_var = np.sqrt(nxs)                #fractional variability

    # Errors on excess variance (taken from Vaughan et al., 2003)
    err_xs = np.sqrt(2*np.square(mse)/N + 4*mse*np.square(xs)/N)
    if nxs>0:
        err_nxs = np.sqrt( np.square(np.sqrt(2/N)*mse/np.square(mean)) + np.square(np.sqrt(mse/N) *2*f_var/mean) )
    else:
        err_nxs = np.sqrt(2/N) * (mse/np.square(mean))

    if normalized:
        return nxs, err_nxs
    else:
        return xs, err_xs


def fractional_variability(rates, errrates, backv, backe, netlightcurve=True):
    """
    Returns the fractional variability and its error given the rates and error rates as arguments.
    The fractional variability is the root square of excess variance, so if the excess variance is negative,
    the functions sets the fractional variability and its error automatically to -1.

    :param rates: rate lightcurve [ct/s]
    :param errrates: error on rates
    :param backv: background rate lightcurve
    :param backe: errors on background  rate
    :param netlightcurve: if set to True, consider the net lightcurve, if set to False consider the total (source+background) lightcurve
    """
    if netlightcurve:
        nxs, err_nxs = excess_variance(rates, errrates, True)
        
        #A value of -1 indicates that the noise of the data is much greater than the scatter of the data.
        if nxs<0:
            nxs = nxs + err_nxs*1.64
            
            f_var = np.sqrt(nxs)
            err_fvar = 1/(2*f_var) * err_nxs
            #f_var = -1.
            #err_fvar = -1.
        else:
            f_var = np.sqrt(nxs)
            err_fvar = 1/(2*f_var) * err_nxs
    else:
        total_rate = rates + backv
        total_errates = np.sqrt( np.square(errrates) + np.square(backe) )
        nxs, err_nxs = excess_variance(total_rate, total_errates, True)
        f_var = np.sqrt(nxs)
        err_fvar = 1/(2*f_var) * err_nxs

    return f_var, err_fvar


def mask_fracexp15(fits_file):
    """
    """
    #Check FRACEXP column, the fractional exposure extension. If not at least 15%, discard the datapoint
    with fits.open(fits_file) as hdul:
        masknan = np.invert(np.isnan(hdul['RATE'].data['RATE']))
        fracexp = hdul['RATE'].data['FRACEXP'][masknan]
        time = hdul['RATE'].data['TIME'][masknan]
        rate = hdul['RATE'].data['RATE'][masknan]
        erate = hdul['RATE'].data['ERROR'][masknan]
        backv = hdul['RATE'].data['BACKV'][masknan]
        backe = hdul['RATE'].data['BACKE'][masknan]

        mask15 = fracexp >= 0.15
        
        time = time[mask15]
        rate = rate[mask15]
        erate = erate[mask15]
        fracexp = fracexp[mask15]
        backv = backv[mask15]
        backe = backe[mask15]
    return (time.byteswap().newbyteorder(), rate.byteswap().newbyteorder(), erate.byteswap().newbyteorder(),
             fracexp.byteswap().newbyteorder(), backv.byteswap().newbyteorder(), backe.byteswap().newbyteorder())


def spectrum_plot_xspec(observation, expid0, expid1, model, target_dir, i=0):
    """
    Plot spectrum in PyXspec.
    """
    xspec.Plot.device = '/null'
    xspec.Plot.xAxis = 'keV'
    xspec.Plot.setRebin(minSig=3, maxBins=4096) 
    xspec.Plot.background = True
    xspec.Plot('data')

    # Use matplotlib
    # Store x and y for RGS1 and RGS2
    chans1 = xspec.Plot.x(1)
    chans1_err = xspec.Plot.xErr(1)
    rates1 = xspec.Plot.y(1)
    rates1_err = xspec.Plot.yErr(1)
    chans2 = xspec.Plot.x(2)
    chans2_err = xspec.Plot.xErr(2)
    rates2 = xspec.Plot.y(2)
    rates2_err = xspec.Plot.yErr(2)

    fig, axs = plt.subplots(3, 1, sharex=True, gridspec_kw={'hspace':0.1}, figsize=(20,15))
    if i==0:
        fig.suptitle(f"Average spectra {observation}, expo {expid0}+{expid1}, {model} ", fontsize=25, y=0.92)
    else:
        fig.suptitle(f"{observation} expo {expid0}+{expid1} {model} Part {i}", fontsize=25, y=0.92)

    # First panel: Source Spectrum
    axs[0].errorbar(chans1, rates1, yerr=rates1_err, xerr=chans1_err, linestyle='', color='black', label='RGS1')
    axs[0].errorbar(chans2, rates2, yerr=rates2_err, xerr=chans2_err, linestyle='', color='red', label='RGS2')
    
    # First panel: Background Spectrum
    xspec.Plot('back')
    axs[0].errorbar(chans1, xspec.Plot.y(1), xerr=chans1_err, yerr=xspec.Plot.yErr(1), color='wheat', linestyle='', marker='.', markersize=0.2, label='Background RGS1')
    axs[0].errorbar(chans2, xspec.Plot.y(2), xerr=chans2_err, yerr=xspec.Plot.yErr(2), color='tan', linestyle='', marker='.', markersize=0.2, label='Background RGS2')
    
    #axs[0].set_yscale('log')
    axs[0].set_xscale('log')
    #axs[0].set_ylim(1e-7)
    axs[0].set_xlim(0.331, 2.001)
    axs[0].set_ylabel('Normalized counts [s-1 keV-1]', fontsize=17)
    axs[0].legend(loc='upper right', fontsize='x-large')

    # Second panel: Residuals
    xspec.Plot('resid')
    axs[1].errorbar(chans1, xspec.Plot.y(1), yerr=xspec.Plot.yErr(1), linestyle='', color='black', label='RGS1')
    axs[1].errorbar(chans2, xspec.Plot.y(2), yerr=xspec.Plot.yErr(2), linestyle='', color='red', label='RGS2')
    axs[1].hlines(0, plt.xlim()[0], plt.xlim()[1], color='lime')   
    axs[1].set_ylabel('Normalized counts [s-1 keV-1]', fontsize=17)

    # Third panel: delchi
    xspec.Plot('delchi')
    axs[2].errorbar(chans1, xspec.Plot.y(1), yerr=1., linestyle='', color='black', label='RGS1')
    axs[2].errorbar(chans2, xspec.Plot.y(2), yerr=1., linestyle='', color='red', label='RGS2')
    axs[2].hlines(0, plt.xlim()[0], plt.xlim()[1], color='lime')  
    axs[2].set_xlabel('Energy [keV]', fontsize=17)
    axs[2].set_ylabel('(data-model)/error', fontsize=17)

    if i==0:
        plt.savefig(f"{target_dir}/Products/RGS_Spectra/{observation}/average_{observation}_{expid0}+{expid1}_{model}.png")
    else:
        plt.savefig(f"{target_dir}/Products/RGS_Spectra/{observation}/{observation}_{expid0}+{expid1}_{model}_{i}.png")
    plt.close()


def epic_spectrum_plot_xspec(observation, expid, model, target_dir):
    """
    Plot spectrum in PyXspec for EPIC-pn.
    """

    xspec.Plot.device = '/null'
    xspec.Plot.xAxis = 'keV'
    xspec.Plot.setRebin(minSig=3, maxBins=4096) 
    xspec.Plot.background = True
    xspec.Plot('data')

    # Use matplotlib
    # Store x and y for EPIC-pn
    chans1 = xspec.Plot.x(1)
    chans1_err = xspec.Plot.xErr(1)
    rates1 = xspec.Plot.y(1)
    rates1_err = xspec.Plot.yErr(1)

    fig, axs = plt.subplots(3, 1, sharex=True, gridspec_kw={'hspace':0.1}, figsize=(11,9))
    fig.suptitle(f"Average spectra {observation}, expo {expid}, {model} ", fontsize=15, y=0.92)

    # First panel: Source Spectrum
    axs[0].errorbar(chans1, rates1, yerr=rates1_err, xerr=chans1_err, linestyle='', color='black', label='EPIC-pn')
    
    # First panel: Background Spectrum
    xspec.Plot('back')
    axs[0].errorbar(chans1, xspec.Plot.y(1), xerr=chans1_err, yerr=xspec.Plot.yErr(1), color='salmon', linestyle='', marker='.', markersize=0.2, label='Background')
    
    axs[0].set_yscale('log')
    axs[0].set_xscale('log')
    #axs[0].set_ylim(1e-7)
    axs[0].set_xlim(0.4, 12)
    axs[0].set_ylabel('Normalized counts [s-1 keV-1]', fontsize=10)
    axs[0].legend(loc='upper right', fontsize='x-large')

    # Second panel: Residuals
    xspec.Plot('resid')
    axs[1].errorbar(chans1, xspec.Plot.y(1), yerr=xspec.Plot.yErr(1), xerr=chans1_err, linestyle='', color='black', label='RGS1')
    axs[1].hlines(0, plt.xlim()[0], plt.xlim()[1], color='lime')   
    axs[1].set_ylabel('Normalized counts [s-1 keV-1]', fontsize=10)

    # Third panel: delchi
    xspec.Plot('delchi')
    axs[2].errorbar(chans1, xspec.Plot.y(1), yerr=1., xerr=chans1_err, linestyle='', color='black', label='RGS1')
    axs[2].hlines(0, plt.xlim()[0], plt.xlim()[1], color='lime')  
    axs[2].set_xlabel('Energy [keV]', fontsize=10)
    axs[2].set_ylabel('(data-model)/error', fontsize=10)

    plt.savefig(f"{target_dir}/Products/EPIC_Spectra/average_{observation}_{model}.png")
    plt.close()


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
        #print(n_in_segment)
        if n_in_segment <= 2:
            x += segment
            continue

        else:
            mean_x.append(np.mean(segment_df[colx].values))
            mean_y.append(np.mean(segment_df[coly].values))
            mean_yerr.append(np.std(segment_df[coly].values)/np.sqrt(len(segment_df[coly].values)))
            mean_xerr.append((segment_df[colx].values[-1]- segment_df[colx].values[0])/2.)
            x += segment

    return mean_x, mean_y, mean_xerr, mean_yerr


def constant(x, q):
    return 0*x+q