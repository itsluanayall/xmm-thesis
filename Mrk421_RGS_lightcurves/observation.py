import logging
import os
import glob
from config import CONFIG
from tools import run_command, split_rgs_filename, sort_rgs_list

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
                raise Exception
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

        for i in range(self.npairs):
            output_name = f'{self.obsid}_expo{i}_RGS_rates.ds'

            if not glob.glob(output_name):
                logging.info(f'Running rgslccorr SAS command for observation number {self.obsid}...')
                
                if len(pairs_events[i])==2:
                    #Making the lightcurve for RGS1+RGS2 pairs
                    rgslc_command = f"rgslccorr evlist='{pairs_events[i][0]} {pairs_events[i][1]}' srclist='{pairs_srcli[i][0]} {pairs_srcli[i][1]}' timebinsize=100 orders='1' sourceid=1 outputsrcfilename={output_name}"
                    status_rgslc = run_command(rgslc_command)

                elif len(pairs_events[i])==1:
                    #Making the lightcurve for single RGS
                    rgslc_command = f"rgslccorr evlist='{pairs_events[i][0]}' srclist='{pairs_srcli[i][0]}' timebinsize=100 orders='1' sourceid=1 outputsrcfilename={output_name}"
                    status_rgslc = run_command(rgslc_command)

                if (status_rgslc!=0):
                    raise Exception
                else:
                    logging.info(f'RGS lightcurves extracted in {output_name}.')
            else:
                logging.info(f'Lightcurves already extracted in {output_name}.')


    def lightcurve(self, use_grace=False):
        """
        Makes the lightcurve plot and saves it in the rgs directory of the current observation
        The user here has two options: you can either save the lightcurve using the dsplot command
        with XmGrace as plotting package in '.ps' format, or you can plot it using python's
        module matplotlib and save it in the '.png' format. 
        You can choose the desired option setting the USE_GRACE boolean in the config.json file.
        """
        os.chdir(self.rgsdir)

        for i in range(self.npairs):
            output_name = f'{self.obsid}_expo{i}_RGS_rates.ds'

            if use_grace: 
                logging.info(f'The RGS lightcurve {output_name} will be plotted with xmgrace.')
                plot_lc_command = f'dsplot table={output_name} withx=yes x=TIME withy=yes y=RATE plotter="xmgrace -hardcopy -printfile RGS_lightcurve{i}.ps"'
                plot_status = run_command(plot_lc_command)
                if (plot_status!=0):
                    raise Exception
                else:
                    logging.info(f'Lightcurve {output_name} ready and saved.')

            else:  
                logging.info(f'The lightcurve {output_name} will be plotted with matplotlib.')
                import matplotlib.pyplot as plt
                from astropy.io import fits

                #Extract data from the lightcurve fits file produced with rgslccorr
                data = fits.open(output_name) 
                x = data[1].data['TIME']
                y = data[1].data['RATE']
                yerr = data[1].data['ERROR']

                #Plot data and add labels and title
                fig = plt.figure(figsize=(20,10))
                plt.errorbar(x, y, yerr=yerr, color='black', marker='.', ecolor='gray')
                plt.grid(True)
                plt.title(output_name, fontsize=30)
                plt.xlabel('TIME [s]', fontsize=25)
                plt.ylabel('RATE [count/s]', fontsize=25)
                plt.xticks(fontsize=20)
                plt.yticks(fontsize=20)

                #Save figure in rgs directory of the current Observation
                plt.savefig(f'{output_name}.png')
                plt.close()
                logging.info(f'The lightcurve is saved as {output_name}.png')

