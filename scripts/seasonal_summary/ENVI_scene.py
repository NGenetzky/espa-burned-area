#! /usr/bin/env python
from numpy import *
from osgeo import gdal
from osgeo import ogr
from osgeo import osr
from osgeo import gdal_array
from osgeo import gdalconst
import os
import logging

#############################################################################
# Created on 3/24/2014 by Gail Schmidt, USGS/EROS LSRD Project
# Created Python script to open and read the input single band ENVI file
# to obtain particular attributes and data from the file.
#
# History:
#   Updated on 3/26/2014 by Gail Schmidt, USGS/EROS
#       Removed metadata reads from the old TIF files that were not being
#       used
############################################################################
class ENVI_Scene:
    """Class to work with the single band ENVI files.
    """

    filename = ""                 # name of input image file
    dX = 0                        # pixel size in X direction
    dY = 0                        # pixel size in Y direction
    NRow = 0                      # number of lines in the scene
    NCol = 0                      # number of samples in the scene
    NoData = 0                    # fill value
    NorthBoundingCoordinate = 0   # northern bounding coord
    SouthBoundingCoordinate = 0   # southern bounding coord
    WestBoundingCoordinate = 0    # western bounding coord
    EastBoundingCoordinate = 0    # eastern bounding coord
    
    # datasets created by gdal.Open
    dataset = None

    # bands
    band1 = None

    def __init__(self, fname, log_handler=None):
        """Class constructor to open and read metadata of the ENVI file.
        Description: Class constructor verifies the input file exists, then
            opens it, reads the metadata, and establishes pointers to the
            single band.

        History:
          Created on 3/24/2014 by Gail Schmidt, USGS/EROS LSRD Project

        Args:
          fname - name of the input ENVI file to be processed
          log_handler - open log file for logging or None for stdout

        Returns:
            None - error opening or reading the file via GDAL
            Object - successful processing
        """
        logger = logging.getLogger(__name__)
        # make sure the file exists
        if not os.path.exists(fname):
            logger.error('Input file does not exist: ' + fname)
            return None

        # store the filename for this class and open it
        self.filename = fname
        self.dataset = gdal.Open (fname)
        if self.dataset is None:
            logger.error('GDAL could not open input file: ' + fname)
            return None

        # create connections to the band
        self.band1 = self.dataset.GetRasterBand(1)

        # get scene and coordinate information
        gt = self.dataset.GetGeoTransform()
        self.NCol = self.dataset.RasterXSize
        self.dX = gt[1]
        self.WestBoundingCoordinate = gt[0]
        self.EastBoundingCoordinate = self.WestBoundingCoordinate +  \
            (self.NCol * self.dX)

        self.NRow = self.dataset.RasterYSize        
        self.dY = gt[5]
        self.NorthBoundingCoordinate = gt[3]
        self.SouthBoundingCoordinate = self.NorthBoundingCoordinate +  \
            (self.NRow * self.dY)

        self.NoData = self.band1.GetNoDataValue()


    def __del__(self):
        """Class constructor to clear the bands.
        Description: class destructor cleans up all the sub dataset and band
            pointers.
        
        History:
          Created in 2013 by Jodi Riegle and Todd Hawbaker, USGS Rocky Mountain
              Geographic Science Center
        
        Args: None
        
        Returns: Nothing
        """

        self.dataset = None
        self.band1 = None

    
    def xy2ij(self, x, y):
        """Converts projection space to pixel space.
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
        """Converts pixel space to projection space.
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


    def getRowOfBandValues(self, y):
        """Reads the specified line at specified projection y coordinate.
        Description: getRowOfBandValues reads a row of band values from band1
            at the associated row for the projection y coordinate.  The
            projection y coordinate is converted to the actual row/line in the
            current scene, and then that row is read (reading all the samples
            in the line) and returned.
        
        History:
          Created in 2013 by Jodi Riegle and Todd Hawbaker, USGS Rocky Mountain
              Geographic Science Center
          Updated on 5/1/2013 by Gail Schmidt, USGS/EROS LSRD Project
              Modified to validate the y projection coordinate falls within
              the current scene.
        
        Args:
          y - y projection coordinate
        
        Returns:
            None - if pixel location does not fall within the coordinates of
                this scene
            Array - associated line of data y location as band1
        """

        if (y < self.NorthBoundingCoordinate) and  \
           (y > self.SouthBoundingCoordinate):
            # get the i,j pixel (line) for the y projection coordinates
            ij = self.xy2ij(0, y)

            # read the line from band1 in the the scene
            x1 = self.band1.ReadAsArray(0, ij[1], self.NCol, 1)[0,]
        
            return (x1)
        else:
            return (None)


    def getBandValues(self):
        """Reads the specified band of data.
        Description: getBandValues reads an entire band of data for the
            band 1.
        
        History:
          Created on 5/1/2013 by Gail Schmidt, USGS/EROS LSRD Project
        
        Args: None
        
        Returns:
            Array - associated data for band 1
        """

        x = zeros ((self.NRow, self.NCol))
        return (self.band1.ReadAsArray())
        
#### end ENVI_Scene class ####

