#! /usr/bin/env python
#############################################################################
# Created in 2013 by Jodi Riegle and Todd Hawbaker, USGS Rocky Mountain
#       Geographic Science Center
# Created Python script to loop through all the individual burn probability
# (BP) images and then create a burned/unburned classified (BC) image.  The
# BP images are first filtered so only pixels with probabilities greater than
# or equal to the seed probability threshold are returned (as a binary result).
# Then the algorithm finds patches in the binary image and calculates the
# size of the burned areas.  The seed size threshold is used to eliminate
# small patches.  The remaining patches are grown by adding adjacent pixels
# with probabilities greater than or equal to the flood fill probability
# threshold.
#
# History:
#   Updated on 11/26/2013 by Gail Schmidt, USGS/EROS
#       Modified to incorporate into the ESPA environment
#   Updated on 4/14/2014 by Gail Schmidt, USGS/EROS
#       Modified to utilize the ESPA raw binary file format
#       Modified to run the scenes in parallel
#   Updated on 2/11/2015 by Gail Schmidt, USGS/EROS
#       Modified the recfromcsv calls to not specify the datatype and to
#       instead use the automatically-determined datatype from the read itself.
#############################################################################

import sys
import os
import time
import getopt
import multiprocessing, Queue
import logging

import numpy
import scipy.ndimage
import skimage.measure

from argparse import ArgumentParser
from osgeo import gdal
from osgeo import ogr
from osgeo import osr
from osgeo import gdal_array
from osgeo import gdalconst

ERROR = 1
SUCCESS = 0

#############################################################################
# Created on April 14, 2014 by Gail Schmidt, USGS/EROS
# Created Python class to handle the multiprocessing of a stack of scenes.
#
# History:
#
############################################################################
class parallelSceneThresholdWorker(multiprocessing.Process):
    """Runs the burn thresholding in parallel for a stack of scenes.
    """
 
    def __init__ (self, work_queue, result_queue, stackObject):
        # base class initialization
        multiprocessing.Process.__init__(self)
 
        # job management stuff
        self.work_queue = work_queue
        self.result_queue = result_queue
        self.stackObject = stackObject
        self.kill_received = False
 

    def run(self):
        logger = logging.getLogger(__name__)  # Obtain logger for this module.
        while not self.kill_received:
            # get a task
            try:
                xml_file = self.work_queue.get_nowait()
            except Queue.Empty:
                break
 
            # process the scene
            msg = 'Processing %s ...' % xml_file
            logger.info(msg)
            status = self.stackObject.sceneBurnThreshold (xml_file)
            if status != SUCCESS:
                msg = 'Error running burn thresholding on the XML file ' \
                    '(%s). Processing will terminate.' % xml_file
                logger.info(msg)
 
            # store the result
            self.result_queue.put(status)


#############################################################################
# Created on November 29, 2013 by Gail Schmidt, USGS/EROS
# Turned into a class to run the overall burn thresholds on the burn
#     probabilities.
#
# History:
#   Updated on 4/13/2014 by Gail Schmidt, USGS/EROS
#       Modified to utilize the ESPA internal file format.
#
# Usage: do_threshold_stack.py --help prints the help message
############################################################################
class BurnAreaThreshold():
    """Class for handling the burned area thresholding functions.
    """

    def __init__(self):
        pass


    def writeResults(self, outputData, outputFilename, geotrans, prj, nodata, \
                     outputRAT=None):
        """Writes an array of data to an output file.
        Description: simple function to write an array of data to an output
            image file
        
        History:
          Created in 2013 by Jodi Riegle and Todd Hawbaker, USGS Rocky Mountain
              Geographic Science Center
          Updated in April, 2013 by Gail Schmidt, USGE/EROS LSRD Project
              Modified to utilize the ESPA internal file format.
        
        Args:
          outputData - output data structure to be written
          outputFilename - name of output file to write the array of data
          geotrans - affine transform for mapping pixel coordinates into
              projection space (returned from GDAL GetGeoTransform())
          prj - projection coordinate system (returned from GDAL
              GetProjectionRef())
          nodata - fill or nodata data value for the output image
          outputRAT - name of output raster attribute table (RAT)
        
        Returns:
            Nothing
        """

        # create the ENVI driver for output data
        driver = gdal.GetDriverByName('ENVI')
        
        # create the output dataset
        bp_dataset = driver.Create(outputFilename, outputData.shape[1],  \
            outputData.shape[0], 1, gdal.GDT_Int16)
        bp_dataset.SetGeoTransform(geotrans)
        bp_dataset.SetProjection(prj)
        
        # get the output band
        bp_band = bp_dataset.GetRasterBand(1)
        bp_band.SetNoDataValue(nodata)
        bp_band.WriteArray(outputData)
        
        if outputRAT <> None:
            bp_band.SetDefaultRAT(outputRAT)


    def floodFill(self, input_image, row, col, output_image, output_label=1,
                  local_threshold=75, nodata=-9999):
        """Implements the flood fill algorithm for the identified burned areas.
        Description: routine to implement the flood fill of the burned areas
            using the lower threshold test.  Adjacent pixels are added to the
            burn patch if their thresholds are high enough.
        
        History:
          Created in 2013 by Jodi Riegle and Todd Hawbaker, USGS Rocky Mountain
              Geographic Science Center
   
        Args:
          input_image - input image of burn probabilities
          row, col - row and column (line and sample) to start at for the flood
              filling
          output_image - output image of burn/unburned pixels
          output_label - output value to be used to identify the burn pixels as
              part of the flood filling; default is 1
          local_threshold - threshold to be used to add burn pixels from the
              burn probability image to the burn classification; default is 75%
          nodata - pixel value used to identify nodata pixels in the input image
   
        Returns:
            nFill - number of pixels that were flood filled
   
        Notes:
          1. The default lower threshold for flood filling is 75% burn
             probability.
        """
    
        # returns number of pixels that were flood filled
        nFill = 0
        toFill = list()
        toFill.append((row,col))
    
        while len(toFill) > 0:
            (row,col) = toFill.pop()
            x = input_image[row,col]
            lab = output_image[row,col]
                    
            # flood fill pixels that are greater than the defined threshold
            if (x <> nodata) & (lab <> output_label) & (x > local_threshold):
                output_image[row,col] = output_label
                nFill += 1
    
                if row > 1:
                    if output_image[row-1,col] == 0:
                        toFill.append((row-1,col))
    
                if (row+1) < input_image.shape[0]:
                    if output_image[row+1,col] == 0:
                        toFill.append((row+1,col))
    
                if col > 1:
                    if output_image[row,col-1] == 0:
                        toFill.append((row,col-1))
    
                if (col+1) < input_image.shape[1]:
                    if output_image[row,col+1] == 0:
                        toFill.append((row,col+1))
    
        return nFill
        
        
    def findBurnScars(self, bp_image, seed_prob_thresh=97.5,
        seed_size_thresh=5, flood_fill_prob_thresh=75):
        """Identify the seeds for burn scars from the input burn probabilities.
        Description: routine to find burn scars using the flood-fill approach.
          Seed pixels are found by using the seed threshold.  Any pixels with
          a probability higher than the seed threshold are identified as seed
          pixels.  Any area with more than the seed size threshold is used as
          a seed area for growing the burn extent using the flood fill
          process.  Areas with fewer than the seed size are ignored and not
          flagged as burn areas.
        
        History:
          Created in 2013 by Jodi Riegle and Todd Hawbaker, USGS Rocky Mountain
              Geographic Science Center
          Updated on Nov. 26, 2013 by Gail Schmidt, USGS/EROS LSRD Project
              Modified to use int32 arrays vs. int64 arrays for the burn
              classification images
          Updated on Feb. 13, 2015 by Gail Schmidt, USGS/EROS LSRD Project
              Modified the skimage.measure.regionprops call to not use the
              deprecated 'properties' parameter.  Also changed the properties
              values to match the correct names of the dynamic list of props
              which is now created.
        
        Args:
          bp_image - input image of burn probabilities
          seed_prob_thresh - threshold to be used to identify burn pixels
              in the burn probability image as seed pixels for the burn area;
              default is 97.5%
          seed_size_thresh - threshold to be used to identify burn areas in the
              probability image. If the count of burn probability seed pixels is
              greater than this threshold for a certain area then the burn area
              is left and flood-filled.  If the count is less than the seed
              threshold, then this is a false positive and the area is not used
              in the burn classification.  Default is 5 pixels.
          flood_fill_prob_thresh - threshold to be used to add burn pixels
              from the burn probability image to the burn classification via
              flood filling; default is 75%
        
        Returns:
          nFill - number of pixels that were flood filled
        
        Notes:
          1. The default lower threshold for flood filling is 75% burn
             probability.
        """
        logger = logging.getLogger(__name__)  # Obtain logger for this module.
    
        # set up array to hold filled region labels, initialize to 0s which
        # means unburned
        bp_regions = numpy.zeros_like(bp_image, dtype=numpy.int32)
        
        # find regions to start the flood fill from; these regions are seed
        # pixels that are greater than the seed probability threshold
        bp_seeds = bp_image >= seed_prob_thresh    
        
        # group the seed pixels into regions of connected components; these
        # regions will be the start of the flood-fill algorithm
        bp_seed_regions = numpy.zeros_like(bp_seeds, dtype=numpy.int32)
        n_seed_labels = scipy.ndimage.label(bp_seeds, output=bp_seed_regions)
        msg = 'Found %d seeds to use for flood fill' % n_seed_labels
        logger.info(msg)
        
        # get list of region pixel coordinates, use the first pixel from each
        # as the seed for the region
        bp_region_coords = skimage.measure.regionprops(  \
            label_image=bp_seed_regions)

        # loop through regions and flood fill to expand them where they are of
        # an appropriate size
        for i in range(0, len(bp_region_coords)):
            temp_label = bp_region_coords[i]['label'] 
            temp_area = bp_region_coords[i]['area']
            temp_coords = bp_region_coords[i]['coords'][0]
            
            # if the number of pixels in this region exceeds the seed size
            # threshold then process the region by flood-filling to grow the
            # burn area
            if temp_area >= seed_size_thresh:
                col = temp_coords[1]
                row = temp_coords[0]
                
                nFilled = self.floodFill(input_image=bp_image, row=row,  \
                    col=col, output_image=bp_regions, output_label=temp_label, \
                    local_threshold=flood_fill_prob_thresh, nodata=-9999)
                
                if False:
                    # logger.info('#' * 80) # Creates header bar:############################################################
                    logger.debug('''{0}
Seed region label:{1}
Number of pixels in region:{2}
First coordinate:{3}
Filled pixels:{4}'''.format('80' * 60, temp_label, temp_area,
                            temp_coords, nFilled))
        
        
        # find region properties for the flood filled burn areas
        bc2 = bp_regions > 0
        bp_regions2 = numpy.zeros_like(bc2, dtype=numpy.int32)
        n_labels = scipy.ndimage.label(bc2, output=bp_regions2)
        prop_names = ['area','filled_area','max_intensity','mean_intensity',  \
            'min_intensity']
        bp_region2_props = skimage.measure.regionprops(  \
            label_image=bp_regions2, intensity_image=bp_image)
        
        # define the RAT (raster attribute table)
        #print 'Creating raster attribute table...'
        label_rat = gdal.RasterAttributeTable()
        label_rat.CreateColumn("Value", gdalconst.GFT_Integer,  \
            gdalconst.GFU_MinMax)
        
        for prop in prop_names:
            label_rat.CreateColumn(prop, gdalconst.GFT_Real,  \
                gdalconst.GFU_MinMax)
            
        # resize the RAT
        label_rat.SetRowCount(n_labels)
        
        # set values in the RAT
        #print 'Populating raster attribute table...'
        for i in range(0, n_labels):
            # label id
            label_rat.SetValueAsInt(i, 0, bp_region2_props[i]['label'])
            
            for j in range(0, len(prop_names)):
                temp_prop = prop_names[j]
                label_rat.SetValueAsDouble(i, j+1,  \
                    float(bp_region2_props[i][temp_prop]))
        
        return ([bp_regions2, label_rat])
    
    
    def sceneBurnThreshold(self, bp_file):
        """Runs the burn thresholding on the current scene.
        Description: sceneBurnThreshold will run the burn thresholding
            on the current burn probability file.

        History:
          Created in 2013 by Jodi Riegle and Todd Hawbaker, USGS Rocky Mountain
              Geographic Science Center
          Updated on 4/10/2014 by Gail Schmidt, USGS/EROS LSRD Project
              Modified to run as a multi-threaded process.
        
        Args:
          bp_file - name of burn probability file to process
        
        Returns:
            ERROR - error running the thresholding on this file
            SUCCESS - successful processing

        Notes: If the nodata value is not obtainable from the header data,
            then it will be set to -9999 which has been the common value
            used for the overall burned area applications.
        """
        logger = logging.getLogger(__name__)  # Obtain logger for this module.

        # determine the output classification name
        fname = os.path.basename(bp_file).replace('burn_probability.img', \
            'burn_class.img')
        bc_file_name = self.output_dir + '/' + fname
        
        # process the current burn probability file
        bp_dataset = gdal.Open(bp_file)
        if bp_dataset is None:
            msg = 'Failed to open bp file: ' + bp_file
            logger.error(msg)
            return ERROR
        
        # read the only band in the file, band 1
        bp_band = bp_dataset.GetRasterBand(1)
        if bp_band is None:
            msg = 'Failed to open bp band 1 from ' + bp_file
            logger.error(msg)
            return ERROR
        
        # get the projection and scene information
        geotrans = bp_dataset.GetGeoTransform()
        if geotrans is None:
            msg = 'Failed to obtain the GeoTransform info from ' + bp_file
            logger.error(msg)
            return ERROR

        prj = bp_dataset.GetProjectionRef()
        if prj is None:
            msg = 'Failed to obtain the ProjectionRef info from ' + bp_file
            logger.error(msg)
            return ERROR

        nrow = bp_dataset.RasterYSize
        ncol = bp_dataset.RasterXSize
        if (nrow is None) or (ncol is None):
            msg = 'Failed to obtain the RasterXSize and RasterYSize from ' +  \
                bp_file
            logger.error(msg)
            return ERROR

        nodata = bp_band.GetNoDataValue()
        if nodata is None:
            nodata = -9999
            msg = 'Failed to obtain the NoDataValue from %s.  Using %d.' % \
                (bp_file, nodata)
            logger.error(msg)
            
        # array to hold burn scars
        bp_scars = numpy.zeros((nrow, ncol))
        bp_rats = []
        
        # read the probabilities for the current scene
        bp_data = bp_band.ReadAsArray()
        
        # find the final burn scars from the burn probabilities
        bp_scar_results = self.findBurnScars(bp_data, self.seed_prob_thresh,
                                             self.seed_size_thresh,
                                             self.flood_fill_prob_thresh)
        bp_scars = bp_scar_results[0]
        bp_scars[ bp_data < 0 ] = bp_data[ bp_data < 0 ]
        bp_rats.append(bp_scar_results[1])
            
        # output the burn classifications for this scene
        msg = 'Writing output to %s ... ' % bc_file_name
        logger.info(msg)
        self.writeResults(outputData=bp_scar_results[0],
            outputFilename=bc_file_name, geotrans=geotrans, prj=prj,
            nodata=nodata, outputRAT=bp_scar_results[1])

        return SUCCESS

    def runBurnThreshold(self, stack_file=None, input_dir=None,
                         output_dir=None, start_year=None, end_year=None,
                         seed_prob_thresh=97.5, seed_size_thresh=5,
                         flood_fill_prob_thresh=75, num_processors=1):
        """Runs the burn thresholding algorithm to find the burn scars from the
           input burn probabilities.
        Description: routine to find the burn scars using the flood-fill
            approach.  Seed pixels are found by using the seed threshold.  Any
            pixels with a probability higher than the seed threshold are
            identified as seed pixels.  Any area with more than the seed size
            threshold is used as a seed area for growing the burn extent using
            the flood fill process.  Areas with fewer than the seed size are
            ignored and not flagged as burn areas.
        
        History:
          Created in 2013 by Jodi Riegle and Todd Hawbaker, USGS Rocky Mountain
              Geographic Science Center
          Updated on Nov. 26, 2013 by Gail Schmidt, USGS/EROS LSRD Project
              Modified to use int32 arrays vs. int64 arrays for the burn
              classification images
          Updated on Dec. 2, 2013 by Gail Schmidt, USGS/EROS LSRD Project
              Modified to use argparser vs. optionparser, since optionparser
              is deprecated.
          Updated on April 13, 2014 by Gail Schmidt, USGS/EROS LSRD Project
              Modified to utilize the ESPA internal file format.

        Args:
          stack_file - input CSV file with information about the files to be
              processed.  this is generated as part of the seasonal summaries
              application.
          input_dir - location of the burn probability files
          output_dir - location to write the output burn classifications
          start_year - starting year of the stack_file to process; default is
              to start with the lowest year + 1
          end_year - ending year of the stack_file to process; default is to
              end with the highest year
          seed_prob_thresh - threshold to be used to identify burn pixels
              in the burn probability image as seed pixels for the burn area;
              default is 97.5%
          seed_size_thresh - threshold to be used to identify burn areas in the
              probability image. If the count of burn probability seed pixels
              is greater than this threshold for a certain area then the burn
              area is left and flood-filled.  If the count is less than the
              seed threshold, then this is a false positive and the area is
              not used in the burn classification; default is 5 pixels
          flood_fill_prob_thresh - threshold to be used to add burn pixels
              from the burn probability image to the burn classification via
              flood filling; default is 75%
          num_processors - how many processors should be used for parallel
              processing sections of the application
        
        Returns:
            ERROR - error running the burn threshold application
            SUCCESS - successful processing
        """
        logger = logging.getLogger(__name__)  # Obtain logger for this module.

        # if no parameters were passed then get the info from the command line
        if stack_file is None:
            # get the command line argument for the input parameters
            parser = ArgumentParser(  \
                description='Process burn scars using thresholds')
            parser.add_argument ('-c', '--stack_file', type=str,
                dest='stack_file',
                help='input file, csv delimited, each row contains '  \
                     'information about a landsat image',
                metavar='FILE', required=True)
            parser.add_argument ('-i', '--input_dir', type=str,
                dest='input_dir',
                help='input directory, location to find input burn '  \
                     'probability files',
                metavar='DIR', required=True)
            parser.add_argument ('-o', '--output_dir', type=str,
                dest='output_dir',
                help='output directory, location to write burn '  \
                     'classification files',
                metavar='DIR', required=True)
            parser.add_argument ('-b', '--start_year', type=int,
                dest='start_year',
                help='starting year for processing the stack of burn '  \
                     'probabilities; default is to use the minimum year '  \
                     'of the files in the stack_file',
                metavar='YEAR')
            parser.add_argument ('-e', '--end_year', type=int, dest='end_year',
                help='ending year for processing the stack of burn '  \
                    'probabilities; default is to use the maximum year '  \
                    'of the files in the stack_file',
                metavar='YEAR')
            parser.add_argument ('-d', '--seed_prob_thresh', type=float,
                dest='seed_prob_thresh',
                help='probability threshold used to seed patch locations; '  \
                    'default is 97.5%%',
                metavar='THRESHOLD')
            parser.add_argument ('-s', '--seed_size_thresh', type=int,
                dest='seed_size_thresh',
                help='minimum size of seed patches, in pixels; default is 5',
                metavar='THRESHOLD')
            parser.add_argument ('-f', '--flood_fill_prob_thresh', type=float,
                dest='flood_fill_prob_thresh',
                help='probability threshold used to grow seeds, pixels '  \
                    'with probabilities above this value are added to '  \
                    'the patch; default is 75%%',
                metavar='THRESHOLD')
            parser.add_argument ('-p', '--num_processors', type=int,
                dest='num_processors',
                help='how many processors should be used for parallel '  \
                    'processing sections of the application '  \
                    '(default = 1, single threaded)')

            options = parser.parse_args()

            # validate command-line options and arguments
            stack_file = options.stack_file
            if stack_file is None:
                logger.error('missing CSV stack file cmd-line argument')
                return ERROR

            input_dir = options.input_dir
            if input_dir is None:
                logger.error('missing input directory cmd-line argument')
                return ERROR

            output_dir = options.output_dir
            if output_dir is None:
                logger.error('missing output directory cmd-line argument')
                return ERROR

            if options.start_year is not None:
                start_year = options.start_year

            if options.end_year is not None:
                end_year = options.end_year

            if options.seed_prob_thresh is not None:
                seed_prob_thresh = options.seed_prob_thresh
    
            if options.seed_size_thresh is not None:
                seed_size_thresh = options.seed_size_thresh
    
            if options.flood_fill_prob_thresh is not None:
                flood_fill_prob_thresh = options.flood_fill_prob_thresh

            # number of processors
            if options.num_processors is not None:
                num_processors = options.num_processors
        else:
            num_processors = num_processors

        self.seed_prob_thresh = seed_prob_thresh
        self.seed_size_thresh = seed_size_thresh
        self.flood_fill_prob_thresh = flood_fill_prob_thresh

        # validate options and arguments
        if start_year is not None:
            if (start_year < 1984):
                msg = 'start_year cannot begin before 1984: %d' % start_year
                logger.error(msg)
                return ERROR

        if end_year is not None:
            if (end_year < 1984):
                msg = 'end_year cannot begin before 1984: %d' % end_year
                logger.error(msg)
                return ERROR

        if (end_year is not None) & (start_year is not None):
            if end_year < start_year:
                msg = 'end_year (%d) is less than start_year (%d)' %  \
                    (end_year, start_year)
                logger.error(msg)
                return ERROR
            
        if not os.path.exists(stack_file):
            msg = 'CSV stack file does not exist: ' + stack_file
            logger.error(msg)
            return ERROR
    
        if not os.path.exists(input_dir):
            msg = 'Input directory does not exist: ' + input_dir
            logger.error(msg)
            return ERROR
    
        if not os.path.exists(output_dir):
            msg = 'Output directory does not exist: %s. Creating ...' % \
                output_dir
            logger.error(msg)
            os.makedirs(output_dir, 0755)
        self.output_dir = output_dir

        # save the current working directory for return to upon error or when
        # processing is complete
        mydir = os.getcwd()
        msg = 'Changing directories for burn threshold processing: ' +  \
            output_dir
        logger.info(msg)
        os.chdir (output_dir)

        # open the stack file
        stack = numpy.recfromcsv(stack_file, delimiter=",", names=True)
        
        # use the minimum and maximum years in the stack if the start year and
        # end year were not specified on the command line.  start year needs
        # to be one more than the actual starting year in the stack since the
        # burn products require a previous year to process.  so the start year
        # for the burn products is one year after the actual starting year in
        # the stack.
        if start_year is None:
            start_year = numpy.min(stack['year']) + 1
        
        if end_year is None:
            end_year = numpy.max(stack['year'])
        
        stack_mask = (stack['year'] >= start_year) & (stack['year'] <= end_year)
        stack2 = stack[stack_mask]
        
        # read the input data from the stack, for the years specified
        msg = 'Processing burn probabilities for %d-%d' % (start_year, end_year)
        logger.info(msg)
        msg = '  seed probability threshold: %f' % seed_prob_thresh
        logger.info(msg)
        msg = '  seed size threshold: %d' % seed_size_thresh
        logger.info(msg)
        msg = '  flood fill probability threshold: %f' %  \
            flood_fill_prob_thresh
        logger.info(msg)

        # load up the work queue for processing scenes in parallel for burn
        # thresholding
        work_queue = multiprocessing.Queue()
        num_scenes = stack2.shape[0]
        for i in range(num_scenes):
            # use the XML filename in the CSV file to obtain the burn
            # probability filename to be thresholded
            xml_file = stack2['file_'][i]
            bp_file_name = xml_file.replace('.xml','_burn_probability.img')
            if not os.path.exists(bp_file_name):
                msg = 'burn probability file does not exist: ' +  bp_file_name
                logger.info(msg)
                os.chdir (mydir)
                return ERROR

            # add this file to the queue to be processed
            print 'Pushing on the queue ... ' + bp_file_name
            work_queue.put(bp_file_name)

        # create a queue to pass to workers to store the processing status
        result_queue = multiprocessing.Queue()
 
        # spawn workers to process each scene in the stack - run the burn
        # thresholding on each scene in the stack
        msg = 'Spawning %d scenes for burn thresholding via %d '  \
            'processors ....' % (num_scenes, num_processors)
        logger.info(msg)
        for i in range(num_processors):
            worker = parallelSceneThresholdWorker(work_queue, result_queue,
                self)
            worker.start()
 
        # collect the burn threshold results off the queue
        for i in range(num_scenes):
            status = result_queue.get()
            if status != SUCCESS:
                msg = 'Error in burn threshold for %d file in the list '  \
                    '(associated XML file is %s).' % (i, stack2['file_'][i])
                logger.info(msg)
                return ERROR

        # successful completion.  return to the original directory.
        msg = 'Completion of burn threshold.'
        logger.info(msg)
        os.chdir (mydir)
        return SUCCESS

######end of BurnAreaThreshold class######

if __name__ == "__main__":
    sys.exit (BurnAreaThreshold().runBurnThreshold())
