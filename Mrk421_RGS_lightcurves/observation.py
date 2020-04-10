import logging
import os
from statistics import mean
from datetime import datetime as dt
from astropy.io import fits
from astropy.time import Time
import glob
from config import CONFIG
from tools import run_command, split_rgs_filename, sort_rgs_list
import numpy as np
import matplotlib.pyplot as plt

class Exposure:
    """
    Class for an Exposure of an Observation.
    """

    def __init__(self, evenli, srcli):
        """
        Constructor. evenli is the Eventlist of the exposure, which contains all the information
        that will be stored as attributes of the Exposure.
        """
        self.evenli = evenli
        self.srcli = srcli
        with fits.open(evenli) as hdul:
            self.expid = hdul[1].header['EXP_ID']
            
            #Extrapolate start time and date. 
            start_date_str = hdul[1].header['DATE-OBS']
            end_date_str = hdul[1].header['DATE-END']

            start_exp = dt.strptime(start_date_str, "%Y-%m-%dT%H:%M:%S")
            end_exp = end_date = dt.strptime(end_date_str, "%Y-%m-%dT%H:%M:%S")
            duration_exp = end_exp - start_exp

            self.end_exp = end_exp.timestamp()
            self.start_exp = start_exp.timestamp()
            self.duration_exp = duration_exp.total_seconds()

    def overlaps_with(self, exp2):
        """
        Method that returns True if the Exposure's operating time overlaps
        at least 90% with another exposure exp2.
        """
        #Define overlapping interval of the two exposures
        overlap = max(0., min(self.end_exp, exp2.end_exp) - max(self.start_exp, exp2.start_exp))
        #Save the longest duration  
        max_duration = max(self.duration_exp, exp2.duration_exp) 
        #Calculate the non-overlapping interval 
        diff = abs(overlap - max_duration)
        diff_percent = (diff/max_duration)*100

        print(f'Non-overlapping time percentage: {diff_percent}%')
        if diff_percent>10: #Arbitrary value
            return False
        else:
            return True

    def longer_than(self, exp2):
        """
        Method that returns True if self.duration is longer than exposure2's duration.
        """
        max_duration = max(self.duration_exp, exp2.duration_exp)
        if max_duration == self.duration_exp:
            return True
        else: 
            return False


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
        
        #The following attributes will be modified in the odfingest method. 
        # For now, we just inizialize them
        self.revolution = 0    
        self.starttime = 0.
        self.endtime = 0.
        self.duration = 0.
        self.rgsrate = 0.

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
        self.duration = format( ((self.endtime - self.starttime)*86400).value/1000., '.3f')    #duration observation in kiloseconds


    def rgsproc(self):
        """
        Runs the rgsproc SAS command to process and reduce RGS data. 
        """
        os.chdir(self.rgsdir)

        #Check if the data has already been processed: if not, run the command.
        if not glob.glob('*EVENLI0000.FIT'):
            logging.info(f'Running rgsproc command for observation number {self.obsid}...')
            rgs_command = 'rgsproc > my_rgsproc_logfile'
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
        It can occur that the exposures do not come in pairs (RGS1+RGS2): in this case 
        the rgslccorr command takes as argument only one event list (either RGS1 or RGS2).
        """
        os.chdir(self.rgsdir)

        #Sort RGS eventlists according to exposure number
        pairs_events = sort_rgs_list(self.rgsevlists, 'expo_number')
        self.npairs = len(pairs_events)
        logging.info(f'There is(are) {self.npairs} set(s) of exposures for observation {self.obsid}.')
        print(pairs_events)

        #Sort RGS sourcelists according to exposure number
        pairs_srcli = sort_rgs_list(self.rgssrclists, 'expo_number')
        print(pairs_srcli)

        if not glob.glob('*_RGS_rates.ds'):    #If the lightcurves haven't already been generated
            
            for i in range(self.npairs):
                    
                if len(pairs_events[i])==2:
                    #Making the lightcurve for RGS1+RGS2 pairs
                    expos0 = Exposure(pairs_events[i][0], pairs_srcli[i][0])
                    expos1 = Exposure(pairs_events[i][1], pairs_srcli[i][1])
                    expo_num0 = split_rgs_filename(expos0.evenli)['expo_number']
                    expo_num1 = split_rgs_filename(expos1.evenli)['expo_number']

                    #Make sure exposure times overlap
                    if expos0.overlaps_with(expos1):   
                        logging.info(f"Running rgslccorr SAS command for observation number {self.obsid} and exposures {expo_num0}, {expo_num1} ...")
                        rgslc_command = f"rgslccorr evlist='{expos0.evenli} {expos1.evenli}' srclist='{expos0.srcli} {expos1.srcli}' timebinsize=1000 orders='1' sourceid=1 outputsrcfilename={self.obsid}_{expo_num0}+{expo_num1}_RGS_rates.ds"
                        status_rgslc = run_command(rgslc_command)
                    else:   #if the exposures do not overlap at least for 90% of their duration
                        logging.error(f"Exposures {expo_num0}, {expo_num1} do not overlap entirely.")
                        status_rgslc = 1
                        
                elif len(pairs_events[i])==1:
                    #Making the lightcurve for single RGS
                    logging.info(f"Running rgslccorr SAS command for observation number {self.obsid} and exposure {split_rgs_filename(pairs_events[i][0])['expo_number']} ...")
                    rgslc_command = f"rgslccorr evlist='{pairs_events[i][0]}' srclist='{pairs_srcli[i][0]}' timebinsize=1000 orders='1' sourceid=1 outputsrcfilename={self.obsid}_{split_rgs_filename(pairs_events[i][0])['expo_number']}_RGS_rates.ds"
                    status_rgslc = run_command(rgslc_command)

                #If an error occurred try running on separate exposures rgslccorr
                if status_rgslc!=0:
                    print(f'\033[91m An error has occurred running rgslccorr for observation {self.obsid}! \033[0m')

                    if len(pairs_events[i])==2:
                        print(f'\033[91m Will try to run rgslccorr with a single exposure (the longest). \033[0m')
                        if expos0.longer_than(expos1):
                            logging.info(f"Running rgslccorr SAS command for observation number {self.obsid} and exposure {expo_num0} ...")
                            rgslc_command2 = f"rgslccorr evlist='{expos0.evenli}' srclist='{expos0.srcli}' timebinsize=1000 orders='1' sourceid=1 outputsrcfilename={self.obsid}_{expo_num0}_RGS_rates.ds"
                            status_rgslc2 = run_command(rgslc_command2)
                        else:
                            logging.info(f"Running rgslccorr SAS command for observation number {self.obsid} and exposure {expo_num1} ...")
                            rgslc_command2 = f"rgslccorr evlist='{expos1.evenli}' srclist='{expos1.srcli}' timebinsize=1000 orders='1' sourceid=1 outputsrcfilename={self.obsid}_{expo_num1}_RGS_rates.ds"
                            status_rgslc2 = run_command(rgslc_command2)

                    if status_rgslc2==0:
                        print(f'\33[32m Error solved for {self.obsid}! \033[0m ')
                        logging.info(f"RGS lightcurves extracted.")
                #If no errors occurred, print to stdio success message
                else:
                    logging.info(f'RGS lightcurves successfully extracted.')
                
        else:   #if the Lightcurves have already beeen extracted
            logging.info(f'Lightcurves already extracted.')


    def lightcurve(self, use_grace=False):
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
                    avg_rate = mean(y)
                    if self.npairs!=0.:
                        self.rgsrate = format(avg_rate, '.3f')
                    else:
                        self.rgsrate = 'N/A'
                        
                    #Plot data and add labels and title
                    fig = plt.figure(figsize=(20,10))
                    ax = fig.add_subplot(1, 1, 1)
                    plt.errorbar(x, y, yerr=yerr, color='black', marker='.', ecolor='gray', 
                                label=f'RGS Lightcurve ObsId {self.obsid}')
                    plt.grid(True)
                    plt.title(output_name, fontsize=30)
                    plt.xlabel('TIME [s]', fontsize=25)
                    plt.ylabel('RATE [count/s]', fontsize=25)
                    plt.xticks(fontsize=20)
                    plt.yticks(fontsize=20)

                    #Plot average rate and legend
                    plt.hlines(self.rgsrate, plt.xlim()[0], plt.xlim()[1] ,colors='red', label=f'Average rate: {self.rgsrate: .2f} [count/s]')
                    ax.legend(loc='lower right', fontsize='x-large')

                    #Save figure in rgs directory of the current Observation
                    plt.savefig(f'{output_name}.png')
                    plt.savefig(f'{self.target_dir}/Products/RGS_Lightcurves/{output_name}.png')
                    plt.close()
                    logging.info(f'The lightcurve is saved as {output_name}.png')

            except Exception as e:
                logging.error(e)
