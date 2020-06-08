import logging
import os
from datetime import datetime as dt
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
import glob
import pandas as pd
from config import CONFIG
from tools import *
import numpy as np
import matplotlib.pyplot as plt
import xspec

class Exposure:
    """
    Class for an Exposure of an Observation.
    """

    def __init__(self, evenli, srcli, specli=[], bkgli=[], respli=[]):
        """
        Constructor. evenli is the Eventlist of the exposure, which contains all the information
        that will be stored as attributes of the Exposure.
        """
        self.evenli = evenli
        self.srcli = srcli
        self.specli = specli
        self.bkgli = bkgli
        self.respli = respli

        with fits.open(evenli) as hdul:
            self.obsid = hdul[1].header['OBS_ID']
            self.expid = str(hdul[1].header['EXP_ID']).split(self.obsid)[1]
            self.fullid = hdul[1].header['EXP_ID']
            self.instrume = hdul[1].header['INSTRUME']
            self.tstart = hdul[1].header['TSTART']
            self.tstop = hdul[1].header['TSTOP']
            self.telapse = hdul[1].header['TELAPSE']
            self.start_date_str = hdul[1].header['DATE-OBS']
            self.end_date_str = hdul[1].header['DATE-END']

            '''
            start_exp = dt.strptime(start_date_str, "%Y-%m-%dT%H:%M:%S")
            end_exp = end_date = dt.strptime(end_date_str, "%Y-%m-%dT%H:%M:%S")
            duration_exp = end_exp - start_exp

            self.end_exp = end_exp.timestamp()
            self.start_exp = start_exp.timestamp()
            self.duration_exp = duration_exp.total_seconds()
            '''

    def synchronous_times(self, exp2):
        """
        """
        final_start = max(self.tstart, exp2.tstart)
        final_stop = min(self.tstop, exp2.tstop)
        overlap = max(0., final_stop-final_start)

        try:
            if overlap>0:
                return final_start, final_stop
            else:
                raise Exception
        except Exception as e:
            logging.error("The given exposures do not overlap. Please check the if the input exposures are correct.")


class Observation:
    """
    Observation class for a specific target.
    In order to avoid unnecessary complexity, the method names of the class recall
    the SAS commands used to analyse observations.
    """
    def __init__(self, obsid, target_dir):
        """
        Constructor. The target directory and the observation ID must be passed as arguments.
        The other attributes are then made based on these arguments.
        """
        self._target_dir = target_dir
        self._obsid = obsid
        self._obsdir = os.path.join(self.target_dir, self.obsid) 
        self._odfdir = os.path.join(self.obsdir, 'odf')
        self._rgsdir = os.path.join(self.obsdir, 'rgs')
        
        #The following attributes will be modified in the class methods. 
        # For now, we just inizialize them
        self.revolution = 0    
        self.starttime = 0.
        self.endtime = 0.
        self.duration = 0.
        self.rgsrate = []
        self.stdev = []
        self.fracvardict = []
        self.expoid = []
        self.longterm_lc_times = []

        self.duration_lc_ks = []

    @property
    def target_dir(self):
        return self._target_dir
    @property
    def obsid(self):
        return self._obsid
    @property
    def obsdir(self):
        return self._obsdir
    @property
    def odfdir(self):
        return self._odfdir
    @property
    def rgsdir(self):
        return self._rgsdir


    def cifbuild(self):
        """
        Generates the Calibration Index File (CIF) for the observation
        in the directory of the observation, and sets the 'SAS_CCF' variable pointing to it.
        """
        #Point SAS_ODF to the odf of the current observation 
        os.environ["SAS_ODF"] = self.odfdir
        os.chdir(self.obsdir)

        #Run the SAS command to make the ccf.cif file (output in logfile)
        if not glob.glob('ccf.cif'):    #check if the ccf file hasn't already been processed
            logging.info(f'Building CIF file for observation number {self.obsid}')
            ccf_command = "cifbuild > my_cifbuild_logfile"
            ccf_status = run_command(ccf_command)
            if (ccf_status != 0):
                raise Exception    
            else:
                logging.info(f"CIF file created for observation number {self.obsid}")
        else:
            logging.info(f'CIF file for observation number {self.obsid} already exists.')

        #Set SAS_CCF variable pointing to the product of 'cifbuild'
        os.environ["SAS_CCF"] = os.path.join(self.obsdir, 'ccf.cif')
        logging.info(f"SAS_CCF pointing to {os.path.join(self.obsdir, 'ccf.cif')}")


    def odfingest(self):
        """
        Generates the SUM.ASC file in the observation directory,
        that summarizes all the observational info
        and sets 'SAS_ODF' variable pointing to it.
        """
        #Point the SAS_ODF environment variable to the odf of the current observation
        os.environ["SAS_ODF"] = self.odfdir
        os.chdir(self.obsdir)

        #Run the SAS command odfingest (output in logfile)
        if not glob.glob('*SUM.SAS'):
            logging.info(f'Building *SUM.SAS file for observation number {self.obsid}')
            odf_command = "odfingest > my_odfingest_logfile"
            odf_status = run_command(odf_command)
            if (odf_status != 0):
                raise Exception
            else:
                logging.info(f"*SUM.SAS file created for observation number {self.obsid}")
        else:
            logging.info(f'*SUM.SAS file for observation number {self.obsid} already exists.')

        #Set SAS_ODF variable pointing to the product of 'odfingest'
        sum_odf_dir = glob.glob('*SUM.SAS')[0]
        os.environ["SAS_ODF"] = os.path.join(self.obsdir, sum_odf_dir)
        logging.info(f"SAS_ODF pointing to {os.path.join(self.obsdir, sum_odf_dir)}")

       # Mark info of the observation and assign Revolution Identifier, Observation Start time and End time
       # to attributes of the class. To move in a new method (?)
        with open(sum_odf_dir) as f:    
            sum_odf = f.read()

        sum_odf = sum_odf.split('\n') #divide text file into lines
        for line in sum_odf:
            if not line.startswith('//'): #skip commented lines
                if line.endswith('Revolution Identifier'):
                    self.revolution = line.split('/')[0]
                if line.endswith('Observation Start Time'):
                    self.starttime = Time(line.split('/')[0], format='isot', scale='utc')
                if line.endswith('Observation End Time'):
                    self.endtime = Time(line.split('/')[0], format='isot', scale='utc')
        self.duration = ((self.endtime - self.starttime)*86400).value    #duration observation in seconds


    def rgsproc(self):
        """
        Runs the rgsproc SAS command to process and reduce RGS data. 
        """
        os.chdir(self.rgsdir)

        #Check if the data has already been processed: if not, run the command.
        if not glob.glob('*EVENLI0000.FIT'):
            logging.info(f'Running rgsproc command for observation number {self.obsid}...')
            ra = CONFIG['target_RA']
            dec = CONFIG['target_DEC']
            rgs_command = f'rgsproc withsrc=yes srclabel=USER srcra={ra} srcdec=+{dec} > my_rgsproc_logfile'
            rgs_status = run_command(rgs_command)
            if (rgs_status != 0):
                print(f'\033[91m An error has occurred running rgsproc for observation {self.obsid}! \033[0m')
            else:
                logging.info(f'Done processing RGS data. The products are in {self.rgsdir}.')
        else:
            logging.info(f'RGS event lists for observation number {self.obsid} already exist.')
        
        #Define the products as new attributes of the Observation
        self.rgsevlists = glob.glob('*EVENLI0000.FIT')
        self.rgssrclists = glob.glob('*SRCLI_0000.FIT')

        
    def rgslccorr(self):
        """
        Runs the rgslccorr SAS command for each pair (RGS1+RGS2) of exposures present in the observation.
        The products are the lightcurves .ds files.
        """
        os.chdir(self.rgsdir)
        
        #Sort RGS eventlists according to exposure number
        self.pairs_events = sort_rgs_list(self.rgsevlists, 'expo_number')


        if self.obsid=='0510610101':
            self.pairs_events = [['P0510610101R1S004EVENLI0000.FIT', 'P0510610101R2S005EVENLI0000.FIT'], ['P0510610101R1S004EVENLI0000.FIT', 'P0510610101R2S013EVENLI0000.FIT']]
        if self.obsid=='0510610201':
            self.pairs_events = [['P0510610201R1S004EVENLI0000.FIT', 'P0510610201R2S005EVENLI0000.FIT'], ['P0510610201R1S015EVENLI0000.FIT', 'P0510610201R2S005EVENLI0000.FIT']]
        if self.obsid=='0136540701':
            self.pairs_events = [['P0136540701R1S001EVENLI0000.FIT','P0136540701R2S002EVENLI0000.FIT' ], ['P0136540701R1S001EVENLI0000.FIT', 'P0136540701R2S018EVENLI0000.FIT'],
                            ['P0136540701R1S011EVENLI0000.FIT','P0136540701R2S018EVENLI0000.FIT'], ['P0136540701R1S020EVENLI0000.FIT', 'P0136540701R2S018EVENLI0000.FIT'],
                            ['P0136540701R1S021EVENLI0000.FIT', 'P0136540701R2S019EVENLI0000.FIT'], ['P0136540701R1S022EVENLI0000.FIT','P0136540701R2S019EVENLI0000.FIT']]
        self.npairs = len(self.pairs_events)
        logging.info(f'There is(are) {self.npairs} set(s) of exposures for observation {self.obsid}.')
        print(self.pairs_events)

        #Sort RGS sourcelists according to exposure number
        self.pairs_srcli = sort_rgs_list(self.rgssrclists, 'expo_number')



        if self.obsid=='0510610101':
            self.pairs_srcli = [['P0510610101R1S004SRCLI_0000.FIT', 'P0510610101R2S005SRCLI_0000.FIT'], ['P0510610101R1S004SRCLI_0000.FIT', 'P0510610101R2S013SRCLI_0000.FIT']]
        if self.obsid=='0510610201':
            self.pairs_srcli = [['P0510610201R1S004SRCLI_0000.FIT', 'P0510610201R2S005SRCLI_0000.FIT'], ['P0510610201R1S015SRCLI_0000.FIT', 'P0510610201R2S005SRCLI_0000.FIT']]
        if self.obsid=='0136540701':
            self.pairs_srcli = [['P0136540701R1S001SRCLI_0000.FIT', 'P0136540701R2S002SRCLI_0000.FIT'], ['P0136540701R1S001SRCLI_0000.FIT', 'P0136540701R2S018SRCLI_0000.FIT'],
                           ['P0136540701R1S011SRCLI_0000.FIT', 'P0136540701R2S018SRCLI_0000.FIT'], ['P0136540701R1S020SRCLI_0000.FIT', 'P0136540701R2S018SRCLI_0000.FIT'],
                           ['P0136540701R1S021SRCLI_0000.FIT', 'P0136540701R2S019SRCLI_0000.FIT'], ['P0136540701R1S022SRCLI_0000.FIT', 'P0136540701R2S019SRCLI_0000.FIT']]
        
          
        for i in range(self.npairs):
            
            # Istance of the exposures
            expos0 = Exposure(self.pairs_events[i][0], self.pairs_srcli[i][0])
            expos1 = Exposure(self.pairs_events[i][1], self.pairs_srcli[i][1])
            self.expoid.append([expos0.expid, expos1.expid])
            
            # Make sure exposure times overlap
            start_time, stop_time = expos0.synchronous_times(expos1)

            if not glob.glob('*_RGS_rates.ds'): #If the lightcurves haven't already been generated, run rgslccorr
                
                logging.info(f"Running rgslccorr SAS command for observation number {self.obsid} and exposures {expos0.expid}, {expos1.expid} ...")
                rgslc_command = f"rgslccorr evlist='{expos0.evenli} {expos1.evenli}' srclist='{expos0.srcli} {expos1.srcli}' withbkgsubtraction=yes timebinsize=1000 timemin={start_time} timemax={stop_time} orders='1' sourceid=3 outputsrcfilename={self.obsid}_{expos0.expid}+{expos1.expid}_RGS_rates.ds outputbkgfilename={self.obsid}_{expos0.expid}+{expos1.expid}_bkg_rates.ds"
                status_rgslc = run_command(rgslc_command)
            
                #If an error occurred try running on separate exposures rgslccorr
                if status_rgslc!=0:
                    print(f'\033[91m An error has occurred running rgslccorr for observation {self.obsid}! \033[0m')
                
                #If no errors occurred, print to stdio success message
                else:
                    logging.info(f'RGS lightcurves successfully extracted.')
                
            else:
                logging.info(f'Lightcurves already extracted.')
            

    def lightcurve(self, mjdref, use_grace=False):
        """
        Makes the lightcurve plot and saves it in the rgs directory of the current observation
        The user here has two options: you can either save the lightcurve using the dsplot command
        with XmGrace as plotting package in '.ps' format, or you can plot it using python's
        module matplotlib and save it in the '.png' format. 
        You can choose the desired option setting the USE_GRACE boolean in the config.json file.
        """
        os.chdir(self.rgsdir)


        for filename in glob.glob('*_RGS_rates.ds'):
            output_name = filename
            try:
                if use_grace: 
                    logging.info(f'The RGS lightcurve {output_name} will be plotted with xmgrace.')
                    plot_lc_command = f'dsplot table={output_name} withx=yes x=TIME withy=yes y=RATE plotter="xmgrace -hardcopy -printfile {output_name}.ps"'
                    plot_status = run_command(plot_lc_command)
                    if (plot_status!=0):
                        raise Exception
                    else:
                        logging.info(f'Lightcurve {output_name} ready and saved.')

                else:  
                    logging.info(f'The lightcurve {output_name} will be plotted with matplotlib.')

                    #Extract data from the lightcurve fits file produced with rgslccorr
                    data = fits.open(output_name)
                    x = data[1].data['TIME']
                    y = data[1].data['RATE']               
                    yerr = data[1].data['ERROR']

                    #Drop NaN values by making a numpy mask
                    mask_nan = np.invert(np.isnan(y)) 
                    x = x[mask_nan]
                    y = y[mask_nan]
                    yerr = yerr[mask_nan]
                    
                    #Store average rate into Observation attribute
                    self.rgsrate.append(np.mean(y))
                    avg_rate = np.mean(y)
                    stdev_rate = np.sqrt(1/(np.sum(1/np.square(yerr))))  #weighted error of mean
                    self.stdev.append(stdev_rate)
                    avg_time = np.mean((x[0], x[-1]))
                    self.duration_lc_ks.append((x[-1] - x[0])/1000.)

                    #Conversion in MJD (note that 86400 are the seconds in one day)
                    avg_time_mjd = mjdref + (avg_time/86400.0)

                    self.longterm_lc_times.append(avg_time_mjd)


                    #Plot data and add labels and title
                    fig = plt.figure(figsize=(20,10))
                    ax = fig.add_subplot(1, 1, 1)
                    plt.errorbar(x, y, yerr=yerr, color='black', marker='.', ecolor='gray', 
                                label=f'RGS Lightcurve ObsId {self.obsid} ')
                    plt.grid(True)
                    plt.title(output_name, fontsize=30)
                    plt.xlabel('TIME [s]', fontsize=25)
                    plt.ylabel('RATE [count/s]', fontsize=25)
                    plt.xticks(fontsize=20)
                    plt.yticks(fontsize=20)

                    #Plot average rate and legend
                    plt.hlines(avg_rate, plt.xlim()[0], plt.xlim()[1], colors='red', label=f'Average rate: {avg_rate: .2f} +- {stdev_rate:.2f} [ct/s]')
                    ax.legend(loc='lower right', fontsize='x-large')

                    #Save figure in rgs directory of the current Observation
                    plt.savefig(f'{output_name}.png')
                    plt.savefig(f'{self.target_dir}/Products/RGS_Lightcurves/{output_name}.png')
                    plt.close()
                    logging.info(f'The lightcurve is saved as {output_name}.png')

            except Exception as e:
                logging.error(e)


    def bkg_lightcurve(self):
        
        logging.info('Generating Background lightcurve...')
        os.chdir(self.rgsdir)
        flat_evenli = [item for sublist in self.pairs_events for item in sublist]
        flat_srcli = [item for sublist in self.pairs_srcli for item in sublist]
        evenli_srcli = list(zip(flat_evenli, flat_srcli))
       
        for (evenli, srcli) in evenli_srcli:
        
            expos0 = Exposure(evenli, srcli)
            title_outputbkg0 = f"bkg{expos0.fullid}_check_rates.fit"
            if not glob.glob(f'{self.target_dir}/Products/Backgrnd_LC/{title_outputbkg0}.png'):
                select_bkg_cmmd0 = f"evselect table={expos0.evenli} timebinsize=1000 rateset={title_outputbkg0} makeratecolumn=yes maketimecolumn=yes expression='(CCDNR==9)&&(REGION({expos0.srcli}:{expos0.instrume}_BACKGROUND, M_LAMBDA, XDSP_CORR))'"
                status_cmmd0 = run_command(select_bkg_cmmd0)
            else:
                continue
            
            #back0_lc = f'dsplot table={title_outputbkg0} withx=yes x=TIME withy=yes y=RATE plotter="xmgrace -hardcopy -printfile {title_outputbkg0}.ps"'
            #plot0_status = run_command(back0_lc)
            
            data = fits.open(title_outputbkg0)
            x = data['RATE'].data['TIME']
            y = data['RATE'].data['RATE']               
            yerr = data['RATE'].data['ERROR']

            #Drop NaN values by making a numpy mask
            mask_nan = np.invert(np.isnan(y)) 
            x = x[mask_nan]
            y = y[mask_nan]
            yerr = yerr[mask_nan]
                
            #Plot data and add labels and title
            fig = plt.figure(figsize=(20,10))
            ax = fig.add_subplot(1, 1, 1)
            plt.errorbar(x, y, yerr=yerr, color='black', marker='.', ecolor='gray')
            plt.grid(True)
            plt.title(title_outputbkg0, fontsize=30)
            plt.xlabel('TIME [s]', fontsize=25)
            plt.ylabel('RATE [count/s]', fontsize=25)
            plt.xticks(fontsize=20)
            plt.yticks(fontsize=20)

            #Save figure in rgs directory of the current Observation
            plt.savefig(f'{self.target_dir}/Products/Backgrnd_LC/{title_outputbkg0}.png')
            plt.close()
            
        logging.info('Done generating background lightcurve!')


    def fracvartest(self, screen=True, netlightcurve=True, timescale=10):
        """
        Reads the FITS file containing the RGS source and background timeseries produced by rgslccorr. 
        It then calculates excess variance, normalized excess variance and fractional variability of the lightcurve,
        storing all these values into a dictionary that will then be reported in the final .csv file.
        """
        os.chdir(self.rgsdir)

        #Recover all lightcurves included in the rgs Observation directory
        i = 0
        for filename in glob.glob('*_RGS_rates.ds'):
            timeseries = filename
            try:
                #Get dataset(s) and table.
                with fits.open(timeseries) as hdul:
                    header = hdul['RATE'].header
                    astro_table = Table(hdul['RATE'].data)
                    dataset = astro_table.to_pandas()
                    dataset = dataset.sort_values(by=['TIME'])

                #Check important keyword consistency 
                if 'TIMEDEL' in header:
                    if header['TIMEDEL']<=0:
                        raise ValueError('\033[91m Null or negative Bin Width in Input File. This should be a positive value. \033[0m')
                else:
                    raise IOError('\033[91m The TIMEDEL keyword is missing in the timeseries FITS file, which is necessary to determine the binning factor of the data. \033[0m')

                if 'TSTART' not in header:
                    raise IOError('\033[91m Keyword TSTART missing in Input. There could be a problem with the input timeseries. \033[0m')
                
                if 'HDUCLASS' in header:
                    if not header['HDUCLASS']=='OGIP':
                        raise ValueError('\033[91m HDUCLASS is not equal to OGIP. \033[0m')
                else:
                    raise IOError('\033[91m Keyword HDUCLASS missing in Input File. There could be a problem with the input FITS timeseries. \033[0m')

                if 'HDUCLAS1' in header:
                    if not header['HDUCLAS1']=='LIGHTCURVE':
                        raise ValueError('\033[91m HDUCLAS1 is not equal to LIGHTCURVE. \033[0m')
                else:
                    raise IOError('\033[91m Keyword HDUCLAS1 missing in Input File. There could be a problem with the input FITS timeseries. \033[0m')
                
            except Exception as e:
                logging.error(e)

            try:
                #Delete gaps in data 
                dataset_cleaned = dataset.dropna()

                #Trim the lightcurve so to calculate the F_var for the chosen timescale
                duration_lc = (dataset_cleaned['TIME'].values[-1]-dataset_cleaned['TIME'].values[0])/1000.  #in kiloseconds
                
                if duration_lc < timescale:
                    logging.info(f'Lightcurve is not long enough to calculate fractional variability on a timescale of {timescale} kiloseconds. Moving on.') 
                    self.fracvardict.append({"Excess variance": np.nan, "Excess variance error": np.nan, "Normalized excess variance": np.nan,
                                "Normalized excess variance error": np.nan, "Fractional Variability": np.nan, 
                                "Fractional Variability Error": np.nan,
                                "Variability Amplitude": np.nan, "Variability amplitude error": np.nan,
                                "Number of non null data points": np.nan})
                    continue 

                elif duration_lc >= timescale:
                    logging.info(f'The lightcurve can be used to calculate fractional variability on a timescale of {timescale} kiloseconds.')
                    max_time = dataset_cleaned['TIME'].values[0] + (timescale*1000)
                    dataset_cleaned = dataset_cleaned[dataset_cleaned['TIME'] <= max_time]

                #Net source rates and errors and background rates and errors are recorded in arrays
                numnonnull = len(dataset_cleaned)
                rates = np.array(dataset_cleaned['RATE'])
                errrates = np.array(dataset_cleaned['ERROR'])
                backv = np.array(dataset_cleaned['BACKV'])
                backe = np.array(dataset_cleaned['BACKE'])
                time = np.array(dataset_cleaned['TIME'])

                #Sanity checks
                if numnonnull<2:
                    raise NoDataException('\033[91m Less than two good values in the timeseries FITS file.')
                if not len(rates)==len(errrates) and len(rates)==len(backv) and len(rates)==len(backe):
                    raise RangeException('\033[91m Different number of rows in columns between RATE, ERROR, BACKV and BACKE. \033[0m')
                if not (rates>0).all() and (backv>0).all():
                    raise ValueError('\033[91m Negative count rates in Input File. \033[0m')

            except NoDataException as e:
                logging.error(e)   
            except RangeException as e:
                logging.error(e)
            except ValueError as e:
                logging.error(e)

            #Calculate Variability parameters and respective errors
            xs, err_xs = excess_variance(rates, errrates, normalized=False)
            nxs, err_nxs = excess_variance(rates, errrates, normalized=True)
            f_var, err_fvar = fractional_variability(rates, errrates, backv, backe, netlightcurve=netlightcurve)
            va = (max(rates) - min(rates))/ (min(rates))
            err_va = va* ( (errrates[rates.argmax()] + errrates[rates.argmin()])/(max(rates) - min(rates)) +  (errrates[rates.argmin()] / min(rates) ))

            logging.info(f'Do you want to carry out the fractional varability amplitude test on the net lightcurve? {netlightcurve}.')
            self.fracvardict.append({"Excess variance": xs, "Excess variance error": err_xs, "Normalized excess variance": nxs,
                                "Normalized excess variance error": err_nxs, "Fractional Variability": f_var, 
                                "Fractional Variability Error": err_fvar,
                                "Variability Amplitude": va, "Variability amplitude error": err_va,
                                "Number of non null data points": numnonnull})
            
            #Write variability test results into header and/or to screen.
            if screen:
                for key, value in self.fracvardict[i].items():
                    print(key, ' : ', value)
            i+=1


    def divide_spectrum(self):
        """
        Splits the spectrum into pieces of 1000 seconds each. To do this, I use a while loop with 2 indices:
        - j represents the end_time of each piece, so after an iteration you just add 1000 seconds to it;
        - i represents the start_time of each piece, so after an iteration you just set it to j.
        At the end of the processing, we will have cut into pieces of 1000s each exposure. For instance, if there are 50 pieces for each exposure, the output will be 50 spectra for RGS1 and 50 spectra for RGS2. The nomenclature of the output spectra is the following:

        sourcespec{instrume}_gti{k}.fits for the source spectrum
        bkgspec{instrume}_gti{k}.fits for the background spectrum

        where instrume can be RGS1 or RGS2, and k=1,2,3...50 is the piece we are considering.
        """

        #Store the eventlists and sourcelists
        os.chdir(self.rgsdir)
        respli = glob.glob('*RSPMAT1*')
        total_spectra = glob.glob('*SRSPEC1*')
        total_bkgr = glob.glob('*BGSPEC1*')

        # Pair the events (RGS1+RGS2) and sort them according to the exposure number(first RGS1, second RGS2)
        self.pairs_respli = sort_rgs_list(respli, "expo_number")
        self.pairs_spectra = sort_rgs_list(total_spectra, "expo_number")
        self.pairs_bkg = sort_rgs_list(total_bkgr, "expo_number")


        if self.obsid=='0510610101':
            self.pairs_respli = [['P0510610101R1S004RSPMAT1003.FIT', 'P0510610101R2S005RSPMAT1003.FIT'], ['P0510610101R1S004RSPMAT1003.FIT', 'P0510610101R2S013RSPMAT1003.FIT']]
            self.pairs_spectra = [['P0510610101R1S004SRSPEC1003.FIT', 'P0510610101R2S005SRSPEC1003.FIT'], ['P0510610101R1S004SRSPEC1003.FIT', 'P0510610101R2S013SRSPEC1003.FIT']]
            self.pairs_bkg = [['P0510610101R1S004BGSPEC1003.FIT','P0510610101R2S005BGSPEC1003.FIT'], ['P0510610101R1S004BGSPEC1003.FIT', 'P0510610101R2S013BGSPEC1003.FIT']]

        if self.obsid=='0510610201':
            self.pairs_respli = [['P0510610201R1S004RSPMAT1003.FIT', 'P0510610201R2S005RSPMAT1003.FIT'], ['P0510610201R1S015RSPMAT1003.FIT', 'P0510610201R2S005RSPMAT1003.FIT']]
            self.pairs_spectra = [['P0510610201R1S004SRSPEC1003.FIT', 'P0510610201R2S005SRSPEC1003.FIT'], ['P0510610201R1S015SRSPEC1003.FIT', 'P0510610201R2S005SRSPEC1003.FIT']]
            self.pairs_bkg = [['P0510610201R1S004BGSPEC1003.FIT', 'P0510610201R2S005BGSPEC1003.FIT'], ['P0510610201R1S015BGSPEC1003.FIT', 'P0510610201R2S005BGSPEC1003.FIT']]

        if self.obsid=='0136540701':
            self.pairs_respli = [['P0136540701R1S001RSPMAT1003.FIT','P0136540701R2S002RSPMAT1003.FIT' ], ['P0136540701R1S001RSPMAT1003.FIT', 'P0136540701R2S018RSPMAT1003.FIT'],
                            ['P0136540701R1S011RSPMAT1003.FIT','P0136540701R2S018RSPMAT1003.FIT'], ['P0136540701R1S020RSPMAT1003.FIT', 'P0136540701R2S018RSPMAT1003.FIT'],
                            ['P0136540701R1S021RSPMAT1003.FIT', 'P0136540701R2S019RSPMAT1003.FIT'], ['P0136540701R1S022RSPMAT1003.FIT','P0136540701R2S019RSPMAT1003.FIT']]
            self.pairs_spectra = [['P0136540701R1S001SRSPEC1003.FIT','P0136540701R2S002SRSPEC1003.FIT' ], ['P0136540701R1S001SRSPEC1003.FIT', 'P0136540701R2S018SRSPEC1003.FIT'],
                            ['P0136540701R1S011SRSPEC1003.FIT','P0136540701R2S018SRSPEC1003.FIT'], ['P0136540701R1S020SRSPEC1003.FIT', 'P0136540701R2S018SRSPEC1003.FIT'],
                            ['P0136540701R1S021SRSPEC1003.FIT', 'P0136540701R2S019SRSPEC1003.FIT'], ['P0136540701R1S022SRSPEC1003.FIT','P0136540701R2S019SRSPEC1003.FIT']]
            self.pairs_bkg = [['P0136540701R1S001BGSPEC1003.FIT','P0136540701R2S002BGSPEC1003.FIT' ], ['P0136540701R1S001BGSPEC1003.FIT', 'P0136540701R2S018BGSPEC1003.FIT'],
                            ['P0136540701R1S011BGSPEC1003.FIT','P0136540701R2S018BGSPEC1003.FIT'], ['P0136540701R1S020BGSPEC1003.FIT', 'P0136540701R2S018BGSPEC1003.FIT'],
                            ['P0136540701R1S021BGSPEC1003.FIT', 'P0136540701R2S019BGSPEC1003.FIT'], ['P0136540701R1S022BGSPEC1003.FIT','P0136540701R2S019BGSPEC1003.FIT']]

        self.npairs = len(self.pairs_events)
        print("Event lists: ", self.pairs_events)
        print("Source lists: ", self.pairs_srcli)
        print("Response matrices: ", self.pairs_respli)
        print("Spectra: ", self.pairs_spectra)
        print("Bkgrounds: ", self.pairs_bkg)

        if not os.path.isdir(f'divided_spectra'):
            os.mkdir('divided_spectra')
        
        k = 1
        self.n_intervals_array =[]
        #Check if the exposures are synchronous
        for l in range(self.npairs):  

            expos0 = Exposure(self.pairs_events[l][0], self.pairs_srcli[l][0], self.pairs_spectra[l][0])
            expos1 = Exposure(self.pairs_events[l][1], self.pairs_srcli[l][1], self.pairs_spectra[l][1])
            start_time, stop_time = expos0.synchronous_times(expos1)
            print(f"Synchronous start and stop times for exposures {expos0.evenli}, {expos1.evenli}: {start_time} - {stop_time}")

        #Divide each eventlist into pieces of 1000 seconds

            if not glob.glob(f'divided_spectra/sourcespec{expos0.instrume}_{expos0.expid}_gti*') or not glob.glob(f'divided_spectra/sourcespec{expos1.instrume}_{expos1.expid}_gti*') :                

                print(f"Processing {expos0.evenli} and {expos1.evenli}...")

                # Initialize indices for loop
                step = 1000 #seconds
                i = start_time
                j = i + step
                if self.obsid not in ['0510610101', '0510610201', '0136540701']:
                    k = 1   #counter of pieces
                self.n_intervals_array.append(k)
                print(k)

                # Divide evenlist in parts of 1000s
                while (i<j) and (j<=stop_time):

                    # Cut and save each piece using SAS command tabgtigen and save output in rgs/divided_spectra/gti{instrume}_file{k}.fits
                    tabgtigen_cmd0 = f"tabgtigen table={expos0.evenli} expression='TIME in [{i}:{j}]' " \
                                    f"gtiset=divided_spectra/gti{expos0.instrume}_{expos0.expid}_file{k}.fits prefraction=0 postfraction=0" 
                    status_tabgtigen0 = run_command(tabgtigen_cmd0)
                    if status_tabgtigen0 !=0:    #Debug
                        print(f'Exception occured with piece n.{k}, exposure {expos0.evenli}, tabgtigen command')

                    tabgtigen_cmd1 = f"tabgtigen table={expos1.evenli} expression='TIME in [{i}:{j}]' " \
                                    f"gtiset=divided_spectra/gti{expos1.instrume}_{expos1.expid}_file{k}.fits prefraction=0 postfraction=0" 
                    status_tabgtigen1 = run_command(tabgtigen_cmd1)
                    if status_tabgtigen1 !=0:    #Debug
                        print(f'Exception occured with piece n.{k}, exposure {expos1.evenli}, tabgtigen command')


                    # Filter the piece that we just created
                    rgsfilter_cmd0 = f"rgsfilter mergedset={expos0.evenli} " \
                                    f"evlist=divided_spectra/event_gti{expos0.instrume}_{expos0.expid}_file{k}.fits " \
                                    f"auxgtitables=divided_spectra/gti{expos0.instrume}_{expos0.expid}_file{k}.fits" 
                    status_rgsfilter0 = run_command(rgsfilter_cmd0)
                    if status_rgsfilter0 !=0:    #Debug
                        print(f'Exception occured with piece n.{k}, exposure {expos0.evenli}, rgsfilter command')

                    rgsfilter_cmd1 = f"rgsfilter mergedset={expos1.evenli} " \
                                    f"evlist=divided_spectra/event_gti{expos1.instrume}_{expos1.expid}_file{k}.fits " \
                                    f"auxgtitables=divided_spectra/gti{expos1.instrume}_{expos1.expid}_file{k}.fits" 
                    status_rgsfilter1 = run_command(rgsfilter_cmd1)
                    if status_rgsfilter1 !=0:    #Debug
                        print(f'Exception occured with piece n.{k}, exposure {expos1.evenli}, rgsfilter command')

                    # Extract spectrum from the piece and save ouput as rgs/divided_spectra/sourcespec{instrume}_gti{k}.fits
                    rgsspectrum_cmd0 = f"rgsspectrum evlist=divided_spectra/event_gti{expos0.instrume}_{expos0.expid}_file{k}.fits " \
                                    f"srclist={expos0.srcli} withspectrum=yes bkgcorrect=no " \
                                    f"spectrumset=divided_spectra/sourcespec{expos0.instrume}_{expos0.expid}_gti{k}.fits withbkgset=yes " \
                                    f"bkgset=divided_spectra/bgspec{expos0.instrume}_{expos0.expid}_gti{k}.fits order=1 rebin=1 edgechannels=2 "\
                                    f"spectrumbinning=lambda withfracexp=no badquality=1"
                    status_rgsspectrum0 = run_command(rgsspectrum_cmd0)
                    if status_rgsspectrum0 !=0:    #Debug
                        print(f'Exception occured with piece n.{k}, exposure {expos0.evenli}, rgsspectrum command')
                    else:
                        print(f'Extracted piece number {k}, exposure {expos0.evenli}.')

                    rgsspectrum_cmd1 = f"rgsspectrum evlist=divided_spectra/event_gti{expos1.instrume}_{expos1.expid}_file{k}.fits " \
                                    f"srclist={expos1.srcli} withspectrum=yes bkgcorrect=no " \
                                    f"spectrumset=divided_spectra/sourcespec{expos1.instrume}_{expos1.expid}_gti{k}.fits withbkgset=yes " \
                                    f"bkgset=divided_spectra/bgspec{expos1.instrume}_{expos1.expid}_gti{k}.fits order=1 rebin=1 edgechannels=2 "\
                                    f"spectrumbinning=lambda withfracexp=no badquality=1"
                    status_rgsspectrum1 = run_command(rgsspectrum_cmd1)
                    if status_rgsspectrum1 !=0:    #Debug
                        print(f'Exception occured with piece n.{k}, exposure {expos1.evenli}, rgsspectrum command')
                    else:
                        print(f'Extracted piece number {k}, exposure {expos1.evenli}.')
                    
                    # Change index after each iteration
                    i = j
                    j = i + step
                    if j>stop_time:
                        j = stop_time
                    k += 1

                print(f"Done dividing {expos0.evenli} and {expos1.evenli}!") 
                print(k-1)
                self.n_intervals_array.append(k-1)         

            else:
                print(f"Divided spectra already extracted for {expos0.evenli} and {expos1.evenli}.")
                
                
    def xspec_divided_spectra_average(self, target_REDSHIFT):
        """
        XSPEC analysis of the total, average spectrum (RGS1+RGS2) of the observation. The steps are all logged into a file called XSPECLogFile_average_spectrum.txt.
        The fit is performed on two different models: logparabola and powerlaw. The plot of the spectra and residuals is done
        using matplotlib. The plotting function is written in tools.py
        The flux and luminosity are stored, given the target_REDSHIFT as argument.
        The fitted parameters and the spectrum counts are all stored into an astropy Table that is then saved as a FITS file.
        """             
        if not glob.glob(f'{self.target_dir}/Products/RGS_Spectra/{self.obsid}/*_1.png'):
                    
            logging.info(f"Starting spectral analysis with XSPEC for observation {self.obsid} (total, average spectra).")
            os.chdir(f"{self.target_dir}/{self.obsid}/rgs")
            


            if not os.path.isdir(f'{self.target_dir}/Products/RGS_Spectra/{self.obsid}'):
                os.mkdir(f'{self.target_dir}/Products/RGS_Spectra/{self.obsid}')

            #Set XSPEC verbosity and create xspec log file
            xspec.Xset.chatter = 10
            xspec.Xset.logChatter = 20
            logFile = xspec.Xset.openLog(f"{self.target_dir}/Products/RGS_Spectra/{self.obsid}/XSPECLogFile_average_spectrum.txt") 
            logFile = xspec.Xset.log

            # Create Table that we will fill with the output (parameters, flux, luminosity...)
            self.spectra_table = Table(names=('obsid', 'instrid', 'exposures_id',
                                        'tbinid', 'tbinstart', 'tbinstop', 'exposure', 'model',
                                        'nH', 'nH_low', 'nH_up', 'phoindex', 'phoindex_low', 'phoindex_up',
                                        'beta', 'beta_low', 'beta_up', 'norm', 'norm_low', 'norm_up',
                                        'cstat', 'chi2', 'dof', 
                                        'src_cts', 'esrc_cts', 'bkg_cts', 'ebkg_cts',
                                        'rate', 'erate', 
                                        'flux', 'flux_up', 'flux_low',
                                        'lumin', 'lumin_up', 'lumin_low'
                                        ),
                                        dtype=('i','U9', 'U9',
                                            'i', 'd', 'd', 'd','U20',
                                            'd','d', 'd','d', 'd', 'd',
                                            'd','d','d','d', 'd', 'd',
                                            'd','d','i', 
                                            'd', 'd', 'd', 'd',
                                            'd', 'd',
                                            'd', 'd', 'd', 
                                            'd', 'd', 'd'))

            
            # Xspec models we want to use for fitting
            model_list = ['const*tbabs*zlogpar', 'const*tbabs*zpowerlw']  

            for i in range(self.npairs):
                expos0 = Exposure(self.pairs_events[i][0], self.pairs_srcli[i][0], self.pairs_spectra[i][0], self.pairs_bkg[i][0], self.pairs_respli[i][0])
                expos1 = Exposure(self.pairs_events[i][1], self.pairs_srcli[i][1], self.pairs_spectra[i][1], self.pairs_bkg[i][1], self.pairs_respli[i][1])
                start_time, stop_time = expos0.synchronous_times(expos1)
                #Load RGS1 + RGS2 data                
                os.chdir(f"{self.target_dir}/{self.obsid}/rgs")
                xspec.AllData(f"1:1 {expos0.specli} 2:2 {expos1.specli}")

                spectrum1 = xspec.AllData(1)
                spectrum1.background = f"{expos0.bkgli}"
                spectrum1.response = f"{expos0.respli}"

                spectrum2 = xspec.AllData(2)
                spectrum2.background = f"{expos1.bkgli}"
                spectrum2.response = f"{expos1.respli}"

                # Ignore bad channels
                xspec.AllData.ignore("bad")
                xspec.AllData.ignore('**-0.331 2.001-**')

                for model in model_list:
                    m1 = xspec.Model(model)
                    m2 = xspec.AllModels(2) #Retrieve the model object assigned to data group 2

                    if m1.expression=='constant*TBabs*zlogpar':
                        m1.setPars("1.0 -1", {5:0.331, 6:0.0308})
                        m2.setPars("1.0", "/*")

                    if m1.expression=='constant*TBabs*zpowerlw':
                        m1.setPars("1.0 -1", {4:0.0308})
                        m2.setPars("1.0", "/*")

                    xspec.AllModels.show()

                    #Perform Fit using cstat as statistic for parameter estimation
                    xspec.Fit.statMethod = "cstat" 
                    xspec.Fit.statTest = "pchi"
                    xspec.Fit.renorm()    #renormalize model to minimize statistic with current parameters
                    xspec.Fit.query = 'yes'
                    xspec.Fit.perform() 

                    #Error calculation (confidence intervals)
                    if m1.expression=='constant*TBabs*zlogpar':
                        xspec.Fit.error("stopat 1000,, maximum 1000.0 2,3,4,7")    

                    if m1.expression=='constant*TBabs*zpowerlw':
                        xspec.Fit.error("stopat 1000,, maximum 1000.0 2,3,5")

                    #Plot
                    spectrum_plot_xspec(self.obsid, expos0.expid, expos1.expid, model, self.target_dir, 0)

                    #Calculate Flux and Luminosity and store their values 
                    xspec.AllModels.calcFlux('0.331 2.001 err 100 90')
                    xspec.AllModels.calcLumin(f'0.331 2.001 {target_REDSHIFT} err') 
                    flux = spectrum1.flux[0] #erg/cm2/s
                    lumin = spectrum1.lumin[0] #e+44 erg/s
                    flux_up = spectrum1.flux[2]
                    flux_low = spectrum1.flux[1]
                    lumin_up = spectrum1.lumin[2]
                    lumin_low = spectrum1.lumin[1]

                    #Store parameter results of fit
                    if m1.expression=='constant*TBabs*zpowerlw':
                        nH = m1(2).values[0]
                        phoindex = m1(3).values[0]
                        norm = m1(5).values[0]
                        #Confidence intervals
                        nH_low = m1(2).error[0]
                        nH_up = m1(2).error[1]
                        phoindex_low = m1(3).error[0]
                        phoindex_up = m1(3).error[1]
                        norm_low = m1(5).error[0]
                        norm_up = m1(5).error[1]
                        beta, beta_up, beta_low = (np.nan, np.nan, np.nan)

                    if m1.expression=='constant*TBabs*zlogpar':
                        nH = m1(2).values[0]
                        phoindex = m1(3).values[0]
                        beta = m1(4).values[0]
                        norm = m1(7).values[0]
                        #Confidence intervals
                        nH_low = m1(2).error[0]
                        nH_up = m1(2).error[1]
                        phoindex_low = m1(3).error[0]
                        phoindex_up = m1(3).error[1]
                        beta_low = m1(4).error[0]
                        beta_up = m1(4).error[1]
                        norm_low = m1(7).error[0]
                        norm_up = m1(7).error[1]

                    fit_statistic = xspec.Fit.statistic
                    test_statistic = xspec.Fit.testStatistic
                    dof = xspec.Fit.dof

                    # Retrieve rates and counts from xspec
                    exposure_time = max(spectrum1.exposure, spectrum2.exposure)
                    src_rate = spectrum1.rate[0] + spectrum2.rate[0]  #net rate
                    src_rate_std = spectrum1.rate[1] + spectrum2.rate[1]
                    frac = spectrum1.rate[3] + spectrum2.rate[3] #model rate

                    src_cts = spectrum1.exposure*spectrum1.rate[0] + spectrum2.exposure*spectrum2.rate[0]
                    src_ects = spectrum1.exposure*spectrum1.rate[1] + spectrum2.exposure*spectrum2.rate[1] 
                    bkg_cts = (1. - frac/100) *src_cts
                    bkg_ects = (1. - frac/100) *src_ects

                    #Save output table
                    self.spectra_table.add_row((self.obsid, "rgs12", f"{expos0.expid}+{expos1.expid}", 0, start_time, stop_time,  exposure_time, m1.expression, nH, nH_low, nH_up,
                                        phoindex, phoindex_low, phoindex_up, beta, beta_low, beta_up, norm, norm_low, norm_up, fit_statistic, test_statistic,
                                        dof, src_cts, src_ects, bkg_cts, bkg_ects, src_rate, src_rate_std, 
                                        flux, flux_up, flux_low, lumin, lumin_up, lumin_low))

                    # Set column units
                    self.spectra_table['tbinstart'].unit = 's'
                    self.spectra_table['tbinstop'].unit = 's'
                    self.spectra_table['exposure'].unit = 's'
                    self.spectra_table['nH'].unit = '10**(22) g / (cm2)'
                    self.spectra_table['nH_up'].unit = '10**(22) g/ (cm2)'
                    self.spectra_table['nH_low'].unit = '10**(22) g/ (cm2)'
                    self.spectra_table['norm'].unit = 'ph /(cm2 s)'
                    self.spectra_table['norm_low'].unit = 'ph/(cm2 s)'
                    self.spectra_table['norm_up'].unit = 'ph/(cm2 s)'
                    self.spectra_table['src_cts'].unit = 'ct'
                    self.spectra_table['bkg_cts'].unit = 'ct'
                    self.spectra_table['ebkg_cts'].unit = 'ct'
                    self.spectra_table['esrc_cts'].unit = 'ct'
                    self.spectra_table['rate'].unit = 'ct/s'
                    self.spectra_table['flux'].unit = 'erg/(cm2 s)'
                    self.spectra_table['flux_up'].unit = 'erg/(cm2 s)'
                    self.spectra_table['flux_low'].unit = 'erg/(cm2 s)'
                    self.spectra_table['lumin'].unit = '10**(44) erg /s'
                    self.spectra_table['lumin_up'].unit = '10**(44) erg/s'
                    self.spectra_table['lumin_low'].unit = '10**(44) erg/s'


            # Close XSPEC's currently opened log file.
            xspec.Xset.closeLog()

            # Write FITS output file 
            if glob.glob(f'{self.target_dir}/Products/RGS_Spectra/{self.obsid}/{self.obsid}_table.fits'):
                os.remove(glob.glob(f'{self.target_dir}/Products/RGS_Spectra/{self.obsid}/{self.obsid}_table.fits')[0])

            self.spectra_table.write(f'{self.target_dir}/Products/RGS_Spectra/{self.obsid}/{self.obsid}_table.fits', format='fits', overwrite=True)
            logging.info('Done average spectral analysis. Please check your Products directory.')
        else:
            logging.info('Spectral analysis on average spectrum already performed.')


    def xspec_divided_spectra(self, target_REDSHIFT):
        """
        If the divide_spectrum method has already been called, this method allows to perform analysis on all the pieces 
        into which we have divided the spectrum using a for loop. The results are stored in the same astropy Table
        of the average spectrum.
        """
        if not glob.glob(f'{self.target_dir}/Products/RGS_Spectra/{self.obsid}/{self.obsid}*_1.png'):
                
            logging.info(f"Starting spectral analysis with XSPEC for observation {self.obsid}, split spectrum.")
            os.chdir(f"{self.target_dir}/{self.obsid}/rgs")

            #Set XSPEC verbosity and create xspec log file
            xspec.Xset.chatter = 10
            xspec.Xset.logChatter = 20
            logFile = xspec.Xset.openLog(f"{self._target_dir}/Products/RGS_Spectra/{self.obsid}/XSPECLogFile_divided_spectra.txt") 
            logFile = xspec.Xset.log
            print(len(self.pairs_spectra))
            for s in range(len(self.pairs_spectra)):
                expos0 = Exposure(self.pairs_events[s][0], self.pairs_srcli[s][0], self.pairs_spectra[s][0], self.pairs_bkg[s][0], self.pairs_respli[s][0])
                expos1 = Exposure(self.pairs_events[s][1], self.pairs_srcli[s][1], self.pairs_spectra[s][1], self.pairs_bkg[s][1], self.pairs_respli[s][1])
                start_time, stop_time = expos0.synchronous_times(expos1)

                # Perform spectral analysis looping over the pieces
                print(self.n_intervals_array)
                self.n_intervals_pairs = [self.n_intervals_array[x:x+2] for x in range(0, len(self.n_intervals_array), 2)]
                print(self.n_intervals_pairs)
                for i in range(self.n_intervals_pairs[s][0], self.n_intervals_pairs[s][1]):
        
                    logging.info(f'Processing gtiRGS1_{expos0.expid}_file{i}.fits for observation {self.obsid} ')
                    if not glob.glob(f'divided_spectra/gtiRGS1_{expos0.expid}_file{i}.fits'):
                        continue

                    with fits.open(f'divided_spectra/gtiRGS1_{expos0.expid}_file{i}.fits') as hdul:
                        tstart = hdul['STDGTI'].header['TSTART']
                        tstop = hdul['STDGTI'].header['TSTOP']

                    #Load RGS1 + RGS2 data
                    xspec.AllData(f"1:1 divided_spectra/sourcespecRGS1_{expos0.expid}_gti{i}.fits 2:2 divided_spectra/sourcespecRGS2_{expos1.expid}_gti{i}.fits")

                    spectrum1 = xspec.AllData(1)
                    spectrum1.background = f"divided_spectra/bgspecRGS1_{expos0.expid}_gti{i}.fits"
                    spectrum1.response = f"{expos0.respli}"

                    spectrum2 = xspec.AllData(2)
                    spectrum2.background = f"divided_spectra/bgspecRGS2_{expos1.expid}_gti{i}.fits"
                    spectrum2.response = f"{expos1.respli}"

                    # Ignore bad channels
                    xspec.AllData.ignore("bad")
                    xspec.AllData.ignore('**-0.331 2.001-**')

                    #Calculate exposure time: if not at least 500 s, skip interval
                    exposure_time = max(spectrum1.exposure, spectrum2.exposure)
                    if exposure_time<500:
                        continue
                    #Loop over models
                    model_list = ['const*tbabs*zlogpar', 'const*tbabs*zpowerlw']
                    for model in model_list:
                        m1 = xspec.Model(model)
                        m2 = xspec.AllModels(2) #Retrieve the model object assigned to data group 2

                        if m1.expression=='constant*TBabs*zlogpar':
                            m1.setPars("1.0 -1", {5:0.331, 6:0.0308})
                            m2.setPars("1.0", "/*")

                        if m1.expression=='constant*TBabs*zpowerlw':
                            m1.setPars("1.0 -1", {4:0.0308})
                            m2.setPars("1.0", "/*")

                        #Perform Fit using cstat as statistic for parameter estimation
                        xspec.Fit.statMethod = "cstat" 
                        xspec.Fit.statTest = "pchi"
                        xspec.Fit.renorm()    #renormalize model to minimize statistic with current parameters
                        xspec.Fit.nIterations = 100
                        xspec.Fit.criticalDelta = 1e-1
                        xspec.Fit.query = 'no' 
                        try:
                            xspec.Fit.perform() 
                        except Exception as e:
                            logging.error(e)
                            continue

                        #Error calculation (confidence intervals)
                        try:
                            if m1.expression=='constant*TBabs*zlogpar':
                                xspec.Fit.error("stopat 1000,, maximum 1000.0 2,3,4,7")    

                            if m1.expression=='constant*TBabs*zpowerlw':
                                xspec.Fit.error("stopat 1000,, maximum 1000.0 2,3,5")
                        except Exception as e:
                            logging.error(e)
                            continue

                        #Plotting
                        spectrum_plot_xspec(self.obsid, expos0.expid, expos1.expid, model, self.target_dir, i)

                        #Calculate Flux and Luminosity and store their values 
                        try:
                            xspec.AllModels.calcFlux('0.331 2.001 err 100 90')
                            xspec.AllModels.calcLumin(f'0.331 2.001 {target_REDSHIFT} err') 
                            flux = spectrum1.flux[0] #erg/cm2/s
                            lumin = spectrum1.lumin[0] #e+44 erg/s
                            flux_up = spectrum1.flux[2]
                            flux_low = spectrum1.flux[1]
                            lumin_up = spectrum1.lumin[2]
                            lumin_low = spectrum1.lumin[1]
                        except Exception as e:
                            logging.error(e)
                            continue

                        #Store parameter results of fit
                        if m1.expression=='constant*TBabs*zpowerlw':
                            nH = m1(2).values[0]
                            phoindex = m1(3).values[0]
                            norm = m1(5).values[0]
                            #Confidence intervals
                            nH_low = m1(2).error[0]
                            nH_up = m1(2).error[1]
                            phoindex_low = m1(3).error[0]
                            phoindex_up = m1(3).error[1]
                            norm_low = m1(5).error[0]
                            norm_up = m1(5).error[1]
                            beta, beta_up, beta_low = (np.nan, np.nan, np.nan)

                        if m1.expression=='constant*TBabs*zlogpar':
                            nH = m1(2).values[0]
                            phoindex = m1(3).values[0]
                            beta = m1(4).values[0]
                            norm = m1(7).values[0]
                            #Confidence intervals
                            nH_low = m1(2).error[0]
                            nH_up = m1(2).error[1]
                            phoindex_low = m1(3).error[0]
                            phoindex_up = m1(3).error[1]
                            beta_low = m1(4).error[0]
                            beta_up = m1(4).error[1]
                            norm_low = m1(7).error[0]
                            norm_up = m1(7).error[1]

                        fit_statistic = xspec.Fit.statistic
                        test_statistic = xspec.Fit.testStatistic
                        dof = xspec.Fit.dof

                        # Retrieve rates and counts from xspec
                        src_rate = spectrum1.rate[0] + spectrum2.rate[0]  #net rate
                        src_rate_std = spectrum1.rate[1] + spectrum2.rate[1]
                        frac = spectrum1.rate[3] + spectrum2.rate[3] #model rate

                        src_cts = spectrum1.exposure*spectrum1.rate[0] + spectrum2.exposure*spectrum2.rate[0]
                        src_ects = spectrum1.exposure*spectrum1.rate[1] + spectrum2.exposure*spectrum2.rate[1] 
                        bkg_cts = (1. - frac/100) *src_cts
                        bkg_ects = (1. - frac/100) *src_ects

                        #Save output table
                        self.spectra_table.add_row((self.obsid,"rgs12", f"{expos0.expid}+{expos1.expid}", i, tstart, tstop,  exposure_time, m1.expression, nH, nH_low, nH_up,
                                        phoindex, phoindex_low, phoindex_up, beta, beta_low, beta_up, norm, norm_low, norm_up, fit_statistic, test_statistic,
                                        dof, src_cts, src_ects, bkg_cts, bkg_ects, src_rate, src_rate_std,
                                        flux, flux_up, flux_low, lumin, lumin_up, lumin_low))

            # Close XSPEC's currently opened log file.
            xspec.Xset.closeLog()
            self.spectra_table.write(f'{self.target_dir}/Products/RGS_Spectra/{self.obsid}/{self.obsid}_table.fits', format='fits', overwrite=True)
            logging.info('Done spectral analysis. Please check your Products directory.')
        else:
            logging.info('Spectral analysis on pieces of spectrum already performed.')