#! /usr/bin/env python
from numpy import *
from osgeo import gdal
from osgeo import ogr
from osgeo import osr
from osgeo import gdal_array
from osgeo import gdalconst
import os
import time
import shutil
import logging


#############################################################################
# Created in 2014 by Gail Schmidt, USGS/EROS
# Created Python script to open and read the input ESPA XML file to obtain
# particular attributes from the file, including filenames as well as
# opening the image and QA surface reflectance data.
#
# History:
#   Updated on 3/26/2014 by Gail Schmidt, USGS/EROS
#       Removed metadata reads from the old HDF files that were not being
#       used
#   Updated on 7/8/2015 by Gail Schmidt, USGS/EROS LSRD Project
#       Thermal band is not used in burned area processing.
#
############################################################################
class XML_Scene:
    """Class for handling ESPA scene related functions.
    """

    band_dict = {}               # dictionary of bands in the current XML file
    xml_file = ""                # name of XML file
    dX = 0                       # pixel size in X direction
    dY = 0                       # pixel size in Y direction
    NRow = 0                     # number of lines in the scene
    NCol = 0                     # number of samples in the scene
    NorthBoundingCoordinate = 0  # northern bounding coord
    SouthBoundingCoordinate = 0  # southern bounding coord
    WestBoundingCoordinate = 0   # western bounding coord
    EastBoundingCoordinate = 0   # eastern bounding coord
    
    # datasets created by gdal.Open
    dataset1 = None
    dataset2 = None
    dataset3 = None
    dataset4 = None
    dataset5 = None
    dataset7 = None 
    dataset_fill_QA = None
    dataset_cloud_QA = None
    dataset_shadow_QA = None
    dataset_snow_QA = None
    dataset_land_water_QA = None
    dataset_adjacent_cloud_QA = None

    # bands
    band1 = None
    band2 = None
    band3 = None
    band4 = None
    band5 = None
    band7 = None
    band_fill_QA = None
    band_cloud_QA = None
    band_shadow_QA = None
    band_snow_QA = None
    band_land_water_QA = None
    band_adjacent_cloud_QA = None

    def __init__ (self, xml_file):
        """Class constructor.
        Description: class constructor verifies the input file exists, then
            opens it, reads the metadata, and establishes pointers to the
            various SDSs.
        
        History:
          Created on 3/17/2014 by Gail Schmidt, USGS EROS LSRD Project

        Args:
          xml_file - name of the input XML reflectance file to be processed
        
        Returns:
            None - error opening or reading the file via GDAL
            Object - successful processing
        """
        logger = logging.getLogger(__name__)  # Obtain logger for this module.

        # make sure the file exists then open it with GDAL
        self.xml_file = xml_file               # store the filename
        if not os.path.exists (xml_file):
            logger.error('Input XML file does not exist: ' + xml_file)
            return None

        # parse the XML file looking for the surface reflectance bands 1-7
        # and the QA bands.  then pass those files to GDAL for resampling.
        self.band_dict['band1'] = xml_file.replace ('.xml', '_sr_band1.img')
        self.band_dict['band2'] = xml_file.replace ('.xml', '_sr_band2.img')
        self.band_dict['band3'] = xml_file.replace ('.xml', '_sr_band3.img')
        self.band_dict['band4'] = xml_file.replace ('.xml', '_sr_band4.img')
        self.band_dict['band5'] = xml_file.replace ('.xml', '_sr_band5.img')
        self.band_dict['band7'] = xml_file.replace ('.xml', '_sr_band7.img')
        self.band_dict['band_fill'] = xml_file.replace ('.xml', \
            '_sr_fill_qa.img')
        self.band_dict['band_cloud'] = xml_file.replace ('.xml', \
            '_sr_cloud_qa.img')
        self.band_dict['band_cloud_shadow'] = xml_file.replace ('.xml', \
            '_sr_cloud_shadow_qa.img')
        self.band_dict['band_snow'] = xml_file.replace ('.xml', \
            '_sr_snow_qa.img')
        self.band_dict['band_land_water'] = xml_file.replace ('.xml', \
            '_sr_land_water_qa.img')
        self.band_dict['band_adjacent_cloud'] = xml_file.replace ('.xml', \
            '_sr_adjacent_cloud_qa.img')
        logger.info(self.band_dict)  # TODO: Should this be debug?
 
        # open connections to the individual bands
        self.dataset1 = gdal.Open(self.band_dict['band1'])
        if self.dataset1 is None:
            logger.error('GDAL could not open input file: ' +
                         self.band_dict['band1'])
            return None

        self.dataset2 = gdal.Open(self.band_dict['band2'])
        if self.dataset2 is None:
            logger.error('GDAL could not open input file: ' +
                         self.band_dict['band2'])
            return None

        self.dataset3 = gdal.Open(self.band_dict['band3'])
        if self.dataset3 is None:
            logger.error('GDAL could not open input file: ' +
                         self.band_dict['band3'])
            return None

        self.dataset4 = gdal.Open(self.band_dict['band4'])
        if self.dataset4 is None:
            logger.error('GDAL could not open input file: ' +
                         self.band_dict['band4'])
            return None

        self.dataset5 = gdal.Open(self.band_dict['band5'])
        if self.dataset5 is None:
            logger.error('GDAL could not open input file: ' +
                         self.band_dict['band5'])
            return None

        self.dataset7 = gdal.Open(self.band_dict['band7'])
        if self.dataset7 is None:
            logger.error('GDAL could not open input file: ' +
                         self.band_dict['band7'])
            return None

        self.dataset_fill_QA = gdal.Open(self.band_dict['band_fill'])
        if self.dataset_fill_QA is None:
            logger.error('GDAL could not open input file: ' +
                         self.band_dict['band_fill'])
            return None

        self.dataset_cloud_QA = gdal.Open(self.band_dict['band_cloud'])
        if self.dataset_cloud_QA is None:
            logger.error('GDAL could not open input file: ' +
                         self.band_dict['band_cloud'])
            return None

        self.dataset_shadow_QA = gdal.Open(self.band_dict['band_cloud_shadow'])
        if self.dataset_shadow_QA is None:
            logger.error('GDAL could not open input file: ' +
                         self.band_dict['band_cloud_shadow'])
            return None

        self.dataset_snow_QA = gdal.Open(self.band_dict['band_snow'])
        if self.dataset_snow_QA is None:
            logger.error('GDAL could not open input file: ' +
                         self.band_dict['band_snow'])
            return None

        self.dataset_land_water_QA =  \
            gdal.Open(self.band_dict['band_land_water'])
        if self.dataset_land_water_QA is None:
            logger.error('GDAL could not open input file: ' +
                         self.band_dict['band_land_water'])
            return None

        self.dataset_adjacent_cloud_QA =  \
            gdal.Open(self.band_dict['band_adjacent_cloud'])
        if self.dataset_adjacent_cloud_QA is None:
            logger.error('GDAL could not open input file: ' +
                         self.band_dict['band_adjacent_cloud'])

        # create connections to the bands
        self.band1 = self.dataset1.GetRasterBand(1)
        self.band2 = self.dataset2.GetRasterBand(1)
        self.band3 = self.dataset3.GetRasterBand(1)
        self.band4 = self.dataset4.GetRasterBand(1)
        self.band5 = self.dataset5.GetRasterBand(1)
        self.band7 = self.dataset7.GetRasterBand(1)
        self.band_fill_QA = self.dataset_fill_QA.GetRasterBand(1)
        self.band_cloud_QA = self.dataset_cloud_QA.GetRasterBand(1)
        self.band_shadow_QA = self.dataset_shadow_QA.GetRasterBand(1)
        self.band_snow_QA = self.dataset_snow_QA.GetRasterBand(1)
        self.band_land_water_QA = self.dataset_land_water_QA.GetRasterBand(1)
        self.band_adjacent_cloud_QA =  \
            self.dataset_adjacent_cloud_QA.GetRasterBand(1)

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
        if self.band_fill_QA is None:
            logger.error('Input band_fill_QA connection failed')
            return None
        if self.band_cloud_QA is None:
            logger.error('Input band_cloud_QA connection failed')
            return None
        if self.band_shadow_QA is None:
            logger.error('Input band_shadow_QA connection failed')
            return None
        if self.band_snow_QA is None:
            logger.error('Input band_snow_QA connection failed')
            return None
        if self.band_land_water_QA is None:
            logger.error('Input band_land_water_QA connection failed')
            return None
        if self.band_adjacent_cloud_QA is None:
            logger.error('Input band_adjacent_cloud_QA connection failed')
            return None

        # get coordinate information from band 1, get the number of rows
        # and columns, then determine the actual bounding coordinates and
        # the delta X and Y values (from west to east and north to south)
        gt = self.dataset1.GetGeoTransform()

        self.NCol = self.dataset1.RasterXSize
        self.dX = gt[1]
        self.WestBoundingCoordinate = gt[0]
        self.EastBoundingCoordinate = self.WestBoundingCoordinate +  \
            (self.NCol * self.dX)

        self.NRow = self.dataset1.RasterYSize        
        self.dY = gt[5]
        self.NorthBoundingCoordinate = gt[3]
        self.SouthBoundingCoordinate = self.NorthBoundingCoordinate +  \
            (self.NRow * self.dY)


    def __del__ (self):
        """Class destructor.
        Description: class destructor cleans up all the dataset and band
            pointers.
        
        History:
          Created in 2013 by Jodi Riegle and Todd Hawbaker, USGS Rocky Mountain
              Geographic Science Center
        
        Args: None
        
        Returns: Nothing
        """

        self.dataset1 = None
        self.dataset2 = None
        self.dataset3 = None
        self.dataset4 = None
        self.dataset5 = None
        self.dataset7 = None
        self.dataset_fill_QA = None
        self.dataset_cloud_QA = None
        self.dataset_shadow_QA = None
        self.dataset_snow_QA = None
        self.dataset_land_water_QA = None
        self.dataset_adjacent_cloud_QA = None

        self.band1 = None
        self.band2 = None
        self.band3 = None
        self.band4 = None
        self.band5 = None
        self.band7 = None
        self.band_fill_QA = None
        self.band_cloud_QA = None
        self.band_shadow_QA = None
        self.band_snow_QA = None
        self.band_land_water_QA = None
        self.band_adjacent_cloud_QA = None


    def xy2ij(self, x, y):
        """Convert from projection space to pixel space.
        Description: xy2ij converts projection coordinates (x,y) to pixel
            space (i,j)
        
        History:
          Created in 2013 by Jodi Riegle and Todd Hawbaker, USGS Rocky Mountain
              Geographic Science Center
        
        Args:
          x - x projection coordinate
          y - y projection coordinate
        
        Returns:
            (col, row) - i,j pixel referring to x,y coordinate
        """

        col = int((x - self.WestBoundingCoordinate) / self.dX)
        row = int((y - self.NorthBoundingCoordinate) / self.dY)
        return([col, row])


    def ij2xy(self, col, row):
        """Convert from pixel space to projection space.
        Description: ij2xy converts pixel space values (i,j) to projection
            coordinates (x,y)
        
        History:
          Created in 2013 by Jodi Riegle and Todd Hawbaker, USGS Rocky Mountain
              Geographic Science Center
        
        Args:
          row - row (line) pixel space coordinate
          col - column (sample) pixel space coordinate
        
        Returns:
            (x, y) - x,y projection coord referring to pixel i,j
        """

        x = self.WestBoundingCoordinate + (col * self.dX)
        y = self.NorthBoundingCoordinate + (row * self.dY)
        return([x,y])


    def getStackBandValues(self, x, y):
        """Reads a stack of pixels at the specified pixel location.
        Description: getStackBandValues reads a stack of pixel values (in the Z
            direction) at the specified x,y location
        
        History:
          Created in 2013 by Jodi Riegle and Todd Hawbaker, USGS Rocky Mountain
              Geographic Science Center
          Updated on 10/30/2013 by Gail Schmidt, USGS/EROS LSRD Project
              Modified to use fmask band for QA vs. LEDAPS QA bands.
          Updated on 12/8/2013 by Gail Schmidt, USGS/EROS LSRD Project
              Modified to back out the fmask band and return to the LEDAPS
              QA bands.
        
        Args:
          x - x projection coordinate
          y - y projection coordinate
        
        Returns:
            None - if pixel location does not fall within the coordinates of
                this scene
            pixel stack - stack of values for the x,y location as band1, band2,
                band3, band4, band5, band7, and QA band
        """

        if (x < self.EastBoundingCoordinate) and  \
           (x > self.WestBoundingCoordinate) and  \
           (y < self.NorthBoundingCoordinate) and  \
           (y > self.SouthBoundingCoordinate):

            # get the i,j pixel for the x,y projection coordinates
            ij = self.xy2ij(x,y)

            # read the pixel from each of the bands in the scene, including
            # the QA bands
            x1 = self.band1.ReadAsArray(ij[0], ij[1], 1, 1)[0,0]
            x2 = self.band2.ReadAsArray(ij[0], ij[1], 1, 1)[0,0]
            x3 = self.band3.ReadAsArray(ij[0], ij[1], 1, 1)[0,0]
            x4 = self.band4.ReadAsArray(ij[0], ij[1], 1, 1)[0,0]
            x5 = self.band5.ReadAsArray(ij[0], ij[1], 1, 1)[0,0]
            x7 = self.band7.ReadAsArray(ij[0], ij[1], 1, 1)[0,0]

            fill_QA = self.band_fill_QA.ReadAsArray(ij[0], ij[1], 1, 1)[0,0]
            cloud_QA = self.band_cloud_QA.ReadAsArray(ij[0], ij[1], 1, 1)[0,0] 
            shadow_QA = self.band_shadow_QA.ReadAsArray(ij[0], ij[1], 1, 1)[0,0]
            snow_QA = self.band_snow_QA.ReadAsArray(ij[0], ij[1], 1, 1)[0,0]
            land_water_QA =  \
                self.band_land_water_QA.ReadAsArray(ij[0], ij[1], 1, 1)[0,0]
            adjacent_cloud_QA =  \
                self.band_adjacent_cloud_QA.ReadAsArray(ij[0], ij[1], 1, 1)[0,0]

            # turn the QA values into one overall value with -9999 representing
            # the noData value
            QA = 0
            if land_water_QA > 0:
                QA = -3
            if snow_QA > 0:
                QA = -4
            if adjacent_cloud_QA > 0:
                QA = -5
            if shadow_QA > 0:
                QA = -6
            if cloud_QA > 0:
                QA = -7
            if fill_QA > 0:    # fill
                QA = -9999
                
            return( {'band1':x1, 'band2':x2, 'band3':x3, 'band4':x4, \
                'band5':x5, 'band7':x7, 'QA':QA} )
        else:
            return (None)


    def getRowOfBandValues(self, y):
        """Reads a line of band values at projection y location.
        Description: getRowOfBandValues reads a row of band values for each
            of the bands specified at the associated row for the projection
            y coordinate.  The projection y coordinate is converted to the
            actual row/line in the current scene, and then that row is read
            (reading all the samples in the line) and returned.
        
        History:
          Created in 2013 by Jodi Riegle and Todd Hawbaker, USGS Rocky Mountain
              Geographic Science Center
          Updated on 5/1/2013 by Gail Schmidt, USGS/EROS LSRD Project
              Modified to validate the y projection coordinate falls within
              the current scene.
          Updated on 10/30/2013 by Gail Schmidt, USGS/EROS LSRD Project
              Modified to use fmask band for QA vs. LEDAPS QA bands.
          Updated on 12/8/2013 by Gail Schmidt, USGS/EROS LSRD Project
              Modified to back out the fmask band and return to the LEDAPS
              QA bands.
        Args:
          y - y projection coordinate
        
        Returns:
            None - if pixel location does not fall within the coordinates of
                this scene
            Array - associated line of data y location as band1, band2,
                band3, band4, band5, band7, and QA band
        """

        if (y < self.NorthBoundingCoordinate) and  \
           (y > self.SouthBoundingCoordinate):

            # get the i,j pixel (line) for the y projection coordinates
            ij = self.xy2ij(0, y)

            # read the line from each of the bands in the scene, including
            # the QA bands
            x1 = self.band1.ReadAsArray(0, ij[1], self.NCol, 1)[0,]
            x2 = self.band2.ReadAsArray(0, ij[1], self.NCol, 1)[0,]
            x3 = self.band3.ReadAsArray(0, ij[1], self.NCol, 1)[0,]
            x4 = self.band4.ReadAsArray(0, ij[1], self.NCol, 1)[0,]
            x5 = self.band5.ReadAsArray(0, ij[1], self.NCol, 1)[0,]
            x7 = self.band7.ReadAsArray(0, ij[1], self.NCol, 1)[0,]
            
            fill_QA = self.band_fill_QA.ReadAsArray(0, ij[1], self.NCol, 1)[0,]
            snow_QA = self.band_snow_QA.ReadAsArray(0, ij[1], self.NCol, 1)[0,]
            land_water_QA = self.band_land_water_QA.ReadAsArray(  \
                0, ij[1], self.NCol, 1)[0,]
            adjacent_cloud_QA = self.band_adjacent_cloud_QA.ReadAsArray(  \
                0, ij[1], self.NCol, 1)[0,]
            shadow_QA = self.band_shadow_QA.ReadAsArray(  \
                0, ij[1], self.NCol, 1)[0,]
            cloud_QA = self.band_cloud_QA.ReadAsArray(  \
                0, ij[1], self.NCol, 1)[0,]

            # turn the QA values into one overall value with -9999 representing
            # the noData value
            QA = zeros( shape(fill_QA), dtype=int16)
            QA[land_water_QA > 0] = -3
            QA[snow_QA > 0] = -4
            QA[shadow_QA > 0] = -5
            QA[adjacent_cloud_QA > 0] = -6
            QA[cloud_QA > 0] = -7
            QA[fill_QA > 0] = -9999  # fill

            return( {'band1':x1, 'band2':x2, 'band3':x3, 'band4':x4,  \
                'band5':x5, 'band7':x7, 'QA':QA } )
        else:
            return (None)


    def getLineOfBandValues(self, j):
        """Reads a line of band values at line j.
        Description: getLineOfBandValues reads a row of band values for each
            of the bands at the specified row.
        
        History:
          Created on 5/1/2013 by Gail Schmidt, USGS/EROS LSRD Project
              Based on the getRowOfBandValues except it passes in the actual
                  line to be read vs. the projection coordinates and therefore
                  needing to convert proj x,y to line,sample (i,j).
          Updated on 10/30/2013 by Gail Schmidt, USGS/EROS LSRD Project
              Modified to use fmask band for QA vs. LEDAPS QA bands.
          Updated on 12/8/2013 by Gail Schmidt, USGS/EROS LSRD Project
              Modified to back out the fmask band and return to the LEDAPS
              QA bands.
        
        Args:
          j - row (line) pixel space coordinate
        
        Returns:
            None - if pixel location does not fall within the coordinates of
                this scene
            Array - associated line of data y location as band1, band2,
                band3, band4, band5, band7, and QA band
        """

        # read the line from each of the bands in the scene, including
        # the QA bands
        x1 = self.band1.ReadAsArray(0, j, self.NCol, 1)[0,]
        x2 = self.band2.ReadAsArray(0, j, self.NCol, 1)[0,]
        x3 = self.band3.ReadAsArray(0, j, self.NCol, 1)[0,]
        x4 = self.band4.ReadAsArray(0, j, self.NCol, 1)[0,]
        x5 = self.band5.ReadAsArray(0, j, self.NCol, 1)[0,]
        x7 = self.band7.ReadAsArray(0, j, self.NCol, 1)[0,]
        
        fill_QA = self.band_fill_QA.ReadAsArray(0, j, self.NCol, 1)[0,]
        snow_QA = self.band_snow_QA.ReadAsArray(0, j, self.NCol, 1)[0,]
        land_water_QA =   \
            self.band_land_water_QA.ReadAsArray(0, j, self.NCol, 1)[0,]
        adjacent_cloud_QA =   \
            self.band_adjacent_cloud_QA.ReadAsArray(0, j, self.NCol, 1)[0,]
        shadow_QA = self.band_shadow_QA.ReadAsArray(0, j, self.NCol, 1)[0,]
        cloud_QA = self.band_cloud_QA.ReadAsArray(0, j, self.NCol, 1)[0,]

        # combine all the QA bands to one output with negative values to
        # indicate the various types of QA values, with -9999 representing
        # the noData value
        QA = zeros (shape(fill_QA), dtype=int16)
        QA[land_water_QA > 0] = -3
        QA[snow_QA > 0] = -4
        QA[shadow_QA > 0] = -5
        QA[adjacent_cloud_QA > 0] = -6
        QA[cloud_QA > 0] = -7
        QA[fill_QA > 0] = -9999  # fill

        return( {'band1':x1, 'band2':x2, 'band3':x3, 'band4':x4,  \
            'band5':x5, 'band7':x7, 'QA':QA } )


    def getBandValues(self, band):
        """Reads the specified band.
        Description: getBandValues reads an entire band of data for the
            specified band.  If it's the QA band, then special processing will
            occur to combine the various QA band values into one representative
            band.
        
        History:
          Created on 5/1/2013 by Gail Schmidt, USGS/EROS LSRD Project
          Updated on 10/30/2013 by Gail Schmidt, USGS/EROS LSRD Project
              Modified to use fmask band for QA vs. LEDAPS QA bands.
          Updated on 12/8/2013 by Gail Schmidt, USGS/EROS LSRD Project
              Modified to back out the fmask band and return to the LEDAPS
              QA bands.
        
        Args:
          band - string representing which band to read (band1, band2, band3,
                 band4, band5, band7, band_qa)
        
        Returns:
            None - if error occurred trying to read the band or combine the
                QA bands
            Array - associated band of data
        """
        logger = logging.getLogger(__name__)  # Obtain logger for this module.
        if not (band in ['band1', 'band2', 'band3', 'band4', 'band5', \
            'band7', 'band_qa']):
            logger.error('Band {0} is not supported. Needs to be one of '
                         'band1, band2, band3, band4, band5, band7,'
                         ' or band_qa.'.format(band))
            return None

        # read and return the specified band
        if band == 'band1':
            return (self.band1.ReadAsArray())

        elif band == 'band2':
            return (self.band2.ReadAsArray())

        elif band == 'band3':
            return (self.band3.ReadAsArray())

        elif band == 'band4':
            return (self.band4.ReadAsArray())

        elif band == 'band5':
            return (self.band5.ReadAsArray())

        elif band == 'band7':
            return (self.band7.ReadAsArray())

        elif band == 'band_qa':
            # read all the QA-related bands
            fill_QA = self.band_fill_QA.ReadAsArray()
            snow_QA = self.band_snow_QA.ReadAsArray()
            land_water_QA = self.band_land_water_QA.ReadAsArray()
            adjacent_cloud_QA = self.band_adjacent_cloud_QA.ReadAsArray()
            shadow_QA = self.band_shadow_QA.ReadAsArray()
            cloud_QA = self.band_cloud_QA.ReadAsArray()
        
            # combine all the QA bands to one output with negative values to
            # indicate the various types of QA values, with -9999 representing
            # the noData value
            QA = zeros (shape(fill_QA), dtype=int16)
            QA[land_water_QA > 0] = -3
            QA[snow_QA > 0] = -4
            QA[shadow_QA > 0] = -5
            QA[adjacent_cloud_QA > 0] = -6
            QA[cloud_QA > 0] = -7
            QA[fill_QA > 0] = -9999  # fill
            return QA


    def createQaBand(self):
        """Creates a single QA band from the multiple QA bands in the surface
           reflectance file.
        Description: createQaBand will create a single QA band using the various
            QA bands from the surface reflectance product.  This will be an
            INT16 product and the noData value will be set to -9999.  The
            name of the band will be {xml_base_name}_mask.img.
        
        History:
          Created in 2013 by Jodi Riegle and Todd Hawbaker, USGS Rocky Mountain
              Geographic Science Center
          Updated on 3/18/2014 by Gail Schmidt, USGS/EROS LSRD Project
              Modified to work with the ESPA internal file format.
        Returns:  N/A
        """
        # create the name of the QA file from the XML filename
        qa_file = self.xml_file.replace ('.xml', '_mask.img')
        self.band_dict['band_qa'] = qa_file

        # create an output file with a single int16 band and get a pointer
        # to this band
        driver = gdal.GetDriverByName('ENVI')
        output_ds = driver.Create (qa_file, self.NCol, self.NRow, 1,  \
            gdal.GDT_Int16)
        output_ds.SetGeoTransform (self.dataset1.GetGeoTransform())
###        output_ds.SetProjection (self.dataset1.GetProjection())
        output_band_QA = output_ds.GetRasterBand(1)

        # initialize the noData value to -9999
        output_band_QA.SetNoDataValue(-9999)

        # read each surface reflectance QA band, then generate the overall QA
        # band, which is a combination of all the QA values (negative values
        # flag any non-clear pixels and -9999 represents the fill pixels).
        vals = self.getBandValues ('band_qa')
        output_band_QA.WriteArray(vals, 0, 0)

        # close the datasets so we can copy over the ENVI header file
        vals = None
        output_band_QA = None
        output_ds = None
        driver = None

        # the GDAL SetGeoTransform and SetProjection don't play completely
        # well with our ENVI header.  just copy the ENVI head for band1 to
        # the ENVI header for the mask band.
        qa_hdr = qa_file.replace ('.img', '.hdr')
        band_hdr = qa_file.replace ('_mask.img', '_sr_band1.hdr')
        shutil.copyfile (band_hdr, qa_hdr)

        return
#####end of XML_Scene class#####
