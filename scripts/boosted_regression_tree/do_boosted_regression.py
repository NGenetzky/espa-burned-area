#! /usr/bin/env python
import sys
import os
import re
import subprocess
import datetime
from argparse import ArgumentParser
import logging


#######################################################################
# Created on September 3, 2013 by Gail Schmidt, USGS/EROS
#     Created Python script to run the boosted regression tree algorithm.
# 
# History:
# 
# Usage: do_boosted_regression.py --help prints the help message
#######################################################################
class BoostedRegression():
    """Class for handling boosted regression tree processing.
    """

    def __init__(self):
        pass


    def runBoostedRegression (self, config_file=None, usebin=None):
        """Runs the boosted regression algorithm for the specified file.
        Description: runBoostedRegression will use the parameter passed for
        the input configuration file.  If input config file is None (i.e. not
        specified) then the command-line parameters will be parsed for this
        information.  The boosted regression tree application is then executed
        to run the regression on the specified input surface reflectance file
        (specified in the input configuration file).  If a log file was
        specified, then the output from this application will be logged to that
        file.
        
        History:
          Created in 2013 by Jodi Riegle and Todd Hawbaker, USGS Rocky Mountain
              Geographic Science Center
          Updated on Dec. 2, 2013 by Gail Schmidt, USGS/EROS LSRD Project
              Modified to use argparser vs. optionparser, since optionparser
              is deprecated.
        Args:
          config_file - name of the input configuration file to be processed
          usebin - this specifies if the boosted regression tree exe resides
              in the $BIN directory; if None then the boosted regression exe
              is expected to be in the PATH
        
        Returns:
            ERROR - error running the boosted regression tree application
            SUCCESS - successful processing
        
        Notes:
          1. The script changes directories to the directory of the
             configuration file.  If absolute paths are not provided in the
             configuration file, then the location of those input/output files
             will need to be the location of the configuration file.
        """
        logger = logging.getLogger(__name__)  # Obtain logger for this module.

        # if no parameters were passed then get the info from the command line
        if config_file is None:
            # get the command line argument for the reflectance file
            parser = ArgumentParser(  \
                description='Run boosted regression algorithm for the scene')
            parser.add_argument ('-c', '--config_file', type=str,
                dest='config_file',
                help='name of configuration file', metavar='FILE')
            parser.add_argument ('--usebin', dest='usebin', default=False,
                action='store_true',
                help='use BIN environment variable as the location of ' \
                     'boosted regression tree application')

            options = parser.parse_args()
    
            # validate the command-line options
            usebin = options.usebin          # should $BIN directory be used

            # surface reflectance file
            config_file = options.config_file
            if config_file is None:
                logger.error('missing configuration file command-line '
                             'argument')
                return ERROR

        # should we expect the boosted regression application to be in the PATH
        # or in the BIN directory?
        if usebin:
            # get the BIN dir environment variable
            bin_dir = os.environ.get('BIN')
            bin_dir = bin_dir + '/'
            logger.info('BIN environment variable: {0}'.format(bin_dir))
        else:
            # don't use a path to the boosted regression application
            bin_dir = ""
            logger.info('boosted regression executable expected to be in the PATH')

        # make sure the configuration file exists
        if not os.path.isfile(config_file):
            logger.error('configuration file does not exist or is not'
                        ' accessible: {0}'.format(config_file))
            return ERROR

        # get the path of the config file and change directory to that location
        # for running this script.  save the current working directory for
        # return to upon error or when processing is complete.  Note: use
        # abspath to handle the case when the filepath is just the filename
        # and doesn't really include a file path (i.e. the current working
        # directory).
        mydir = os.getcwd()
        configdir = os.path.dirname (os.path.abspath (config_file))
        if not os.access(configdir, os.W_OK):
            logger.error('Path of configuration file is not writable: {0}.'
                         '  Boosted regression may need write access to the'
                         ' configuration directory, depending on whether the'
                         ' output files in the configuration file have been'
                         ' specified.'.format(configdir))
            return ERROR
        logger.info('Changing directories for boosted regression processing:'
                    ' {0}'.format(configdir))
        os.chdir (configdir)

        # run boosted regression algorithm, checking the return status.  exit
        # if any errors occur.
        cmdstr = "%spredict_burned_area --config_file %s --verbose" %  \
            (bin_dir, config_file)
        cmdlist = cmdstr.split(' ')
        try:
            output = subprocess.check_output (cmdlist, stderr=None)
            logger.info(output)
        except subprocess.CalledProcessError, e:
            logger.error('Error running boosted regression. Processing will '
                         'terminate.\n ' + e.output)
            os.chdir (mydir)
            return ERROR

        # successful completion.  return to the original directory.
        logger.info('Completion of boosted regression.')

        os.chdir (mydir)
        return SUCCESS

######end of BoostedRegression class######

if __name__ == "__main__":
    sys.exit (BoostedRegression().runBoostedRegression())
