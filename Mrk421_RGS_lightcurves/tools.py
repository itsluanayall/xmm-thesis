import subprocess
import sys
import os
import numpy as np
import logging
from astropy.io import fits

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
    """
    try:
        mean = np.mean(rates)
        if mean<0:
            raise ValueError('Negative count rates in Input File.')
    except ValueError as e:
        loggng.error(e)

    variance = np.var(rates, ddof=1)
    N = len(rates)
    mse = np.mean(np.square(errrates))
    xs = variance - mse
    nxs = xs/(np.square(mean))
    f_var = np.sqrt(nxs)

    err_xs = np.sqrt(2*np.square(mse)/N + 4*mse*np.square(xs)/N)
    err_nxs = np.sqrt( np.square(np.sqrt(2/N)*mse/np.square(mean)) + np.square(np.sqrt(mse/N) *2*f_var/mean) )
    
    if normalized:
        return nxs, err_nxs
    else:
        return xs, err_xs

def fractional_variability(rates, errrates, backv, backe, netlightcurve=True):
    """
    Returns the fractional variability and its error given the rates and error rates as arguments.
    The fractional variability is the root square of excess variance, so if the excess variance is negative,
    the functions sets the fractional variability and its error automatically to -1.
    """
    if netlightcurve:
        nxs, err_nxs = excess_variance(rates, errrates, True)
        f_var = np.sqrt(nxs)
        err_fvar = 1/(2*f_var) * err_nxs
    else:
        total_rate = rates + backv
        total_errates = np.sqrt( np.square(errrates) + np.square(backe) )
        nxs, err_nxs = excess_variance(total_rate, total_errates, True)
        f_var = np.sqrt(nxs)
        err_fvar = 1/(2*f_var) * err_nxs

    #A value of -1 indicates that the noise of the data is much greater than the scatter of the data.
    if nxs<0:
        logging.warning("Excess variance is negative - the noise of the data is much greater than the scatter of the data. Fractional variability is its square root and will return -1.0")
        f_var = -1.
        err_fvar = -1.

    return f_var, err_fvar