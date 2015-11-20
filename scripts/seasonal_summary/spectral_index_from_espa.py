from sys import path
import os, sys, time, datetime, getopt

from numpy import *
from osgeo import gdal
from osgeo import ogr
from osgeo import osr
from osgeo import gdal_array
from osgeo import gdalconst
from spectral_indices import *
import logging

### Error/Success codes ###
ERROR = 1
SUCCESS = 0

#############################################################################
# Created on May 1, 2013 by Gail Schmidt, USGS/EROS
# Created class to hold the methods and attributes for reading the ESPA
# file and generating spectral indices.
#
# History:
#   Updated on 5/21/2013 by Gail Schmidt, USGS/EROS LSRD Project
#       Modified to process all the indices one line  at a time (vs. the
#       entire band) since this is faster.
#   Updated on 3/17/2014 by Gail Schmidt, USGS/EROS LSRD Project
#       Modified to use the ESPA internal raw binary format
#   Updated on 7/8/2015 by Gail Schmidt, USGS/EROS LSRD Project
#       Thermal band is not used in burned area processing.
#
############################################################################
class spectralIndex:
    """Class for producing the spectral indices.
    """

    # datasets created by gdal.Open
    dataset1 = None
    dataset2 = None
    dataset3 = None
    dataset4 = None
    dataset5 = None
    dataset7 = None
    dataset_mask = None             # QA mask created from QA bands

    # bands
    band1 = None
    band2 = None
    band3 = None
    band4 = None
    band5 = None
    band7 = None
    band_mask = None                # QA mask created from QA bands

    def __init__ (self, band_dict):
        """Class constructor which opens the band files.
        Description: spectralIndex class constructor opens the input band
            files and obtains pointers to each of the desired bands.
        
        History:
          Created in 2013 by Jodi Riegle and Todd Hawbaker, USGS Rocky Mountain
              Geographic Science Center
          Updated on 5/2/2013 by Gail Schmidt, USGS/EROS LSRD Project
              Modified to utilize the class constructor to open and set up
              band pointers from the file, so they are available for use by
              other methods within the class.
        
        Args:
          band_dict - dictionary of bands to be opened and processed for this
              class
        
        Returns:
            None - error opening the file
            Object - successful processing
        """
        logger = logging.getLogger(__name__)  # Obtain logger for this module.
        # open connections to the individual bands
        self.dataset1 = gdal.Open(band_dict['band1'])
        if self.dataset1 is None:
            logger.error('GDAL could not open input file: {0}'
                         .format(band_dict['band1']))
            return None

        self.dataset2 = gdal.Open(band_dict['band2'])
        if self.dataset2 is None:
            logger.error('GDAL could not open input file: {0}'
                         .format(band_dict['band2']))
            return None

        self.dataset3 = gdal.Open(band_dict['band3'])
        if self.dataset3 is None:
            logger.error('GDAL could not open input file: {0}'
                         .format(band_dict['band3']))
            return None

        self.dataset4 = gdal.Open(band_dict['band4'])
        if self.dataset4 is None:
            logger.error('GDAL could not open input file: {0}'
                         .format(band_dict['band4']))
            return None

        self.dataset5 = gdal.Open(band_dict['band5'])
        if self.dataset5 is None:
            logger.error('GDAL could not open input file: {0}'
                         .format(band_dict['band5']))
            return None

        self.dataset7 = gdal.Open(band_dict['band7'])
        if self.dataset7 is None:
            logger.error('GDAL could not open input file: {0}'
                         .format(band_dict['band7']))
            return None

        self.dataset_mask = gdal.Open(band_dict['band_qa'])
        if self.dataset_mask is None:
            logger.error('GDAL could not open input mask file: {0}'
                         .format(band_dict['band_qa']))
            return None

        # create connections to the bands
        self.band1 = self.dataset1.GetRasterBand(1)
        self.band2 = self.dataset2.GetRasterBand(1)
        self.band3 = self.dataset3.GetRasterBand(1)
        self.band4 = self.dataset4.GetRasterBand(1)
        self.band5 = self.dataset5.GetRasterBand(1)
        self.band7 = self.dataset7.GetRasterBand(1)
        self.band_mask = self.dataset_mask.GetRasterBand(1)

        # verify the bands were actually accessed successfully
        if self.band1 is None:
            logger.error('Input band1 connection failed')
            return None
        if self.band2 is None:
            logger.error('Input band2 connection failed')
            return None
        if self.band3 is None:
            logger.error('Input band3 connection failed')
            return None
        if self.band4 is None:
            logger.error('Input band4 connection failed')
            return None
        if self.band5 is None:
            logger.error('Input band5 connection failed')
            return None
        if self.band7 is None:
            logger.error('Input band7 connection failed')
            return None
        if self.band_mask is None:
            logger.error('Input band_mask connection failed')
            return None


    def __del__ (self):
        """Class desctructor to clean up band pointers.
        Description: class destructor cleans up all the sub dataset and band
            pointers.
        
        History:
          Created in 2013 by Jodi Riegle and Todd Hawbaker, USGS Rocky Mountain
              Geographic Science Center
        
        Args: None
        
        Returns: Nothing
        """

        # cleanup
        self.band1 = None
        self.band2 = None
        self.band3 = None
        self.band4 = None
        self.band5 = None
        self.band7 = None
        self.band_mask = None

        self.dataset1 = None
        self.dataset2 = None
        self.dataset3 = None
        self.dataset4 = None
        self.dataset5 = None
        self.dataset7 = None
        self.dataset_mask = None


    def createSpectralIndices (self, index_dict):
        """Generates the specified spectral indices.
        Description: createSpectralIndices creates the desired spectral index
            products.  If mask is specified, then a combined mask file is
            generated using the various input masks.
        
        History:
          Created in 2013 by Jodi Riegle and Todd Hawbaker, USGS Rocky Mountain
              Geographic Science Center
          Updated on 5/2/2013 by Gail Schmidt, USGS/EROS LSRD Project
              Modified to utilize a class structure and only read the bands
              if they haven't been read already.  This saves from duplication
              of reading the same band over and over for different indices.
          Updated on 5/21/2013 by Gail Schmidt, USGS/EROS LSRD Project
              Modified to process all the indices one line  at a time (vs. the
              entire band) since this is faster.
        
        Args:
          index_dict - dictionary of index types (ndvi, nbr, nbr2, ndmi, mask)
              and the associated filename for the index file
        
        Returns:
            ERROR - error generating the spectral indices or mask
            SUCCESS - successful processing
        """
        logger = logging.getLogger(__name__)  # Obtain logger for this module.
        num_indices = len(index_dict)

        # ignore divide by zero and invalid (NaN) values when doing array
        # division.  these will be handled on our own.
        seterr(divide='ignore', invalid='ignore')

        # grab the noData value from the input band 1
        ncol = self.dataset1.RasterXSize
        nrow = self.dataset1.RasterYSize
        nodata = self.band1.GetNoDataValue()
        if nodata is None:
            nodata = -9999
            
        # loop through the indices specified to be processed and create the
        # output ENVI files
        output_ds = {}
        output_band = {}
        for index in index_dict.keys():
            # figure out which spectral index to generate
            if not (index in ['ndvi','ndmi','nbr','nbr2']):
                logger.error('Algorithm for {0} is not implemented'
                             .format(index))
                return ERROR

            # create the output folder if it does not exist
            output_dir = os.path.dirname(index_dict[index])
            if not os.path.exists(output_dir):
                logger.info('Creating output directory {0}'
                            .format(output_dir))
                os.makedirs(output_dir)

            # create the output file; spectral indices are multiplied by 1000.0
            # and the mask file is as-is.
            mydriver = gdal.GetDriverByName('ENVI')
            my_ds = mydriver.Create (index_dict[index], ncol, nrow, 1,  \
                gdal.GDT_Int16)
            my_ds.SetGeoTransform (self.dataset1.GetGeoTransform())
            my_ds.SetProjection (self.dataset1.GetProjection())
            output_ds[index] = my_ds
            my_band = my_ds.GetRasterBand(1)    
            my_band.SetNoDataValue(nodata)
            output_band[index] = my_band

        # loop through each line in the image and process
        for y in range (0, nrow):
            # read the QA data
            qa = self.band_mask.ReadAsArray(0, y, ncol, 1)

            # loop through the indices specified and process each index product
            # reusing the line from each band where possible
            b3 = b4 = b5 = b7 = None
            for index in index_dict.keys():
                # calculate the spectral index
                if index == 'nbr':
                    if b4 is None:
                        b4 = self.band4.ReadAsArray(0, y, ncol,1)
                    if b7 is None:
                        b7 = self.band7.ReadAsArray(0, y, ncol, 1)
                    newVals = 1000.0 * NBR(b4, b7, nodata)
                    newVals[qa < 0] = nodata

                    # write the output 
                    my_output_band = output_band[index]
                    my_output_band.WriteArray(newVals, 0, y)

                elif index == 'nbr2':
                    if b5 is None:
                        b5 = self.band5.ReadAsArray(0, y, ncol, 1)
                    if b7 is None:
                        b7 = self.band7.ReadAsArray(0, y, ncol, 1)
                    newVals = 1000.0 * NBR2(b5, b7, nodata)
                    newVals[qa < 0] = nodata

                    # write the output 
                    my_output_band = output_band[index]
                    my_output_band.WriteArray(newVals, 0, y)

                elif index == 'ndmi':
                    if b4 is None:
                        b4 = self.band4.ReadAsArray(0, y, ncol, 1)
                    if b5 is None:
                        b5 = self.band5.ReadAsArray(0, y, ncol, 1)
                    newVals = 1000.0 * NDMI(b4, b5, nodata)
                    newVals[qa < 0] = nodata

                    # write the output 
                    my_output_band = output_band[index]
                    my_output_band.WriteArray(newVals, 0, y)

                elif index == 'ndvi':
                    if b4 is None:
                        b4 = self.band4.ReadAsArray(0, y, ncol, 1)
                    if b3 is None:
                        b3 = self.band3.ReadAsArray(0, y, ncol, 1)
                    newVals = 1000.0 * NDVI(b3, b4, nodata)
                    newVals[qa < 0] = nodata

                    # write the output 
                    my_output_band = output_band[index]
                    my_output_band.WriteArray(newVals, 0, y)
                # end if
            # end for index
        # end for y

        # cleanup the bands
        b3 = b4 = b5 = b7 = None

        # cleanup
        del (output_band)
        del (output_ds)

        return SUCCESS
######end of spectralIndex class######
