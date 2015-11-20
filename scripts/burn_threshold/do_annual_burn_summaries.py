#! /usr/bin/env python
#############################################################################
# Created in 2013 by Jodi Riegle and Todd Hawbaker, USGS Rocky Mountain
#       Geographic Science Center
# Created Python script to generate a single image of burned/unburned,
# maximum burned probability, first julian day a pixel was labeled as burned,
# and the number of good looks for each pixel (non cloud, non water, etc.) on
# an annual basis.
#
# History:
#   Updated on 12/2/2013 by Gail Schmidt, USGS/EROS LSRD Project
#       Modified to incorporate into the ESPA environment
#   Updated on 5/15/2014 by Gail Schmidt, USGS/EROS LSRD Project
#       Modified to clean up the .img.aux.xml files created by GDAL for the
#       burned area products
#   Updated on 5/19/2014 by Gail Schmimdt, USGS/EROS LSRD Project
#       Changed the use of burn scar to burned area
#############################################################################

import sys
import glob
import os
import time
import datetime as datetime_
import getopt
import csv
import logging

import numpy

from argparse import ArgumentParser
from osgeo import gdal
from osgeo import ogr
from osgeo import osr
from osgeo import gdal_array
from osgeo import gdalconst

import metadata_api

### Error/Success codes ###
ERROR = 1
SUCCESS = 0

def convert_imageXY_to_mapXY (image_x, image_y, transform):
    '''
    Description:
      Translate image coordinates into mapp coordinates
    '''
    """Translate the image coordinates into map coordinates.
    Description: routine to convert the image coordinates into map coordinates
        using the geotransform information obtained via GDAL.
    
    History:
      Created in 2014 by Ron Dilley, USGS/EROS LSRD Project
        Pulled from warp.py in the ESPA project to use for this script

    Args:
      image_x - x-coordinate from the image to be converted to map coordinates,
          floating point
      image_y - y-coordinate from the image to be converted to map coordinates,
          floating point
      transform - geotransform array from GDAL GetGeoTransform()

    Returns:
        map coordinate pair (map_x, map_y)
    """

    map_x = transform[0] + image_x * transform[1] + image_y * transform[2]
    map_y = transform[3] + image_x * transform[4] + image_y * transform[5]

    return (map_x, map_y)



#############################################################################
# Created on December 2, 2013 by Gail Schmidt, USGS/EROS LSRD Project
# Turned into a class to run the overall annual burn summaries.
#
# History:
#   Updated on 7/9/2015 by Gail Schmimdt, USGS/EROS LSRD Project
#       Updated the radiance and reflectance info for the output bands in the
#       XML product
#
# Usage: do_annual_burn_summaries.py --help prints the help message
############################################################################
class AnnualBurnSummary():
    """Class for processing the annual burn summaries.
    Usage: do_annual_burn_summaries.py --help prints the help message
    """
    burned_area_version = "version 1.0.0"

    def __init__(self):
        pass


    def createXML(self, scene_xml_file=None, output_xml_file=None,
        start_year=None, end_year=None, fill_value=None, imgfile=None):
        """Creates an XML file for the products produced by
           runAnnualBurnSummaries.
        Description: routine to create the XML file for the burned area summary
            bands.  The sample scene-based XML file will be used as the basis
            for the projection information for the output XML file.  The image
            size, extents, etc. will need to be updated, as will the band
            information.
        
        History:
          Created on May 12, 2014 by Gail Schmidt, USGS/EROS LSRD Project

        Args:
          scene_xml_file - scene-based XML file to be used as the base XML
              information for the projection metadata.
          output_xml_file - name of the XML file to be written
          start_year - starting year of the scenes to process
          end_year - ending year of the scenes to process
          fill_value - fill or nodata value for this dataset
          imgfile - name of burned area image file with associated ENVI header
              which can be used to obtain the extents and geographic
              information for these products
   
        Returns:
            ERROR - error creating the XML file
            SUCCESS - successful creation of the XML file
        """
        logger = logging.getLogger(__name__)  # Obtain logger for this module.

        # parse the scene-based XML file, just as a basis for the output XML
        # file.  the global attributes will be similar, but the extents and
        # size of the image will be different.  the bands will be based on the
        # bands that are output from this routine.
        xml = metadata_api.parse (scene_xml_file, silence=True)
        meta_bands = xml.get_bands()
        meta_global = xml.get_global_metadata()

        # update the global information
        meta_global.set_data_provider("USGS/EROS")
        meta_global.set_satellite("LANDSAT")
        meta_global.set_instrument("combination")
        del (meta_global.acquisition_date)
        meta_global.set_acquisition_date(None)

        # open the image file to obtain the geospatial and spatial reference
        # information
        ds = gdal.Open (imgfile)
        if ds is None:
            logger.error('GDAL failed to open {0}'.format(imgfile))
            return ERROR

        ds_band = ds.GetRasterBand (1)
        if ds_band is None:
            logger.error('GDAL failed to get the first band in {0}'
                         .format(imgfile))
            return ERROR
        nlines = float(ds_band.YSize)
        nsamps = float(ds_band.XSize)
        nlines_int = ds_band.YSize
        nsamps_int = ds_band.XSize
        del (ds_band)

        ds_transform = ds.GetGeoTransform()
        if ds_transform is None:
            logger.error('GDAL failed to get the geographic transform'
                         ' information from {0}'.format(imgfile))
            return ERROR

        ds_srs = osr.SpatialReference()
        if ds_srs is None:
            logger.error('GDAL failed to get the spatial reference information'
                         ' from {0}'.format(imgfile))
            return ERROR
        ds_srs.ImportFromWkt (ds.GetProjection())
        del (ds)

        # get the UL and LR center of pixel map coordinates
        (map_ul_x, map_ul_y) = convert_imageXY_to_mapXY (0.5, 0.5,
            ds_transform)
        (map_lr_x, map_lr_y) = convert_imageXY_to_mapXY (
            nsamps - 0.5, nlines - 0.5, ds_transform)

        # update the UL and LR projection corners along with the origin of the
        # corners, for the center of the pixel (global projection information)
        for mycorner in meta_global.projection_information.corner_point:
            if mycorner.location == 'UL':
                mycorner.set_x (map_ul_x)
                mycorner.set_y (map_ul_y)
            if mycorner.location == 'LR':
                mycorner.set_x (map_lr_x)
                mycorner.set_y (map_lr_y)
        meta_global.projection_information.set_grid_origin("CENTER")

        # update the UL and LR latitude and longitude coordinates, using the
        # center of the pixel
        srs_lat_lon = ds_srs.CloneGeogCS()
        coord_tf = osr.CoordinateTransformation (ds_srs, srs_lat_lon)
        for mycorner in meta_global.corner:
            if mycorner.location == 'UL':
                (lon, lat, height) = \
                    coord_tf.TransformPoint (map_ul_x, map_ul_y)
                mycorner.set_longitude (lon)
                mycorner.set_latitude (lat)
            if mycorner.location == 'LR':
                (lon, lat, height) = \
                    coord_tf.TransformPoint (map_lr_x, map_lr_y)
                mycorner.set_longitude (lon)
                mycorner.set_latitude (lat)

        # determine the bounding coordinates; initialize using the UL and LR
        # then work around the scene edges
        # UL
        (map_x, map_y) = convert_imageXY_to_mapXY (0.0, 0.0, ds_transform)
        (ul_lon, ul_lat, height) = coord_tf.TransformPoint (map_x, map_y)
        # LR
        (map_x, map_y) = convert_imageXY_to_mapXY (nsamps, nlines, ds_transform)
        (lr_lon, lr_lat, height) = coord_tf.TransformPoint (map_x, map_y)

        # find the min and max values accordingly, for initialization
        west_lon = min (ul_lon, lr_lon)
        east_lon = max (ul_lon, lr_lon)
        north_lat = max (ul_lat, lr_lat)
        south_lat = min (ul_lat, lr_lat)

        # traverse the boundaries of the image to determine the bounding
        # coords; traverse one extra line and sample to get the the outer
        # extents of the image vs. just the UL of the outer edge.
        # top and bottom edges
        for samp in range (0, nsamps_int+1):
            # top edge
            (map_x, map_y) = convert_imageXY_to_mapXY (samp, 0.0, ds_transform)
            (top_lon, top_lat, height) = coord_tf.TransformPoint (map_x, map_y)

            # lower edge
            (map_x, map_y) = convert_imageXY_to_mapXY (samp, nlines,
                ds_transform)
            (low_lon, low_lat, height) = coord_tf.TransformPoint (map_x, map_y)

            # update the min and max values
            west_lon = min (top_lon, low_lon, west_lon)
            east_lon = max (top_lon, low_lon, east_lon)
            north_lat = max (top_lat, low_lat, north_lat)
            south_lat = min (top_lat, low_lat, south_lat)

        # left and right edges
        for line in range (0, nlines_int+1):
            # left edge
            (map_x, map_y) = convert_imageXY_to_mapXY (0.0, line, ds_transform)
            (left_lon, left_lat, height) = coord_tf.TransformPoint (map_x,
                map_y)

            # right edge
            (map_x, map_y) = convert_imageXY_to_mapXY (nsamps, line,
                ds_transform)
            (right_lon, right_lat, height) = coord_tf.TransformPoint (map_x,
                map_y)

            # update the min and max values
            west_lon = min (left_lon, right_lon, west_lon)
            east_lon = max (left_lon, right_lon, east_lon)
            north_lat = max (left_lat, right_lat, north_lat)
            south_lat = min (left_lat, right_lat, south_lat)

        # update the XML
        bounding_coords = meta_global.get_bounding_coordinates()
        bounding_coords.set_west (west_lon)
        bounding_coords.set_east (east_lon)
        bounding_coords.set_north (north_lat)
        bounding_coords.set_south (south_lat)

        del (ds_transform)
        del (ds_srs)

        # clear some of the global information that doesn't apply for these
        # products
        del (meta_global.scene_center_time)
        meta_global.set_scene_center_time(None)
        del (meta_global.lpgs_metadata_file)
        meta_global.set_lpgs_metadata_file(None)
        del (meta_global.orientation_angle)
        meta_global.set_orientation_angle(None)
        del (meta_global.level1_production_date)
        meta_global.set_level1_production_date(None)

        # clear the solar angles
        del (meta_global.solar_angles)
        meta_global.set_solar_angles(None)

        # save the first band and then wipe the bands out so that new bands
        # can be added for the burned area bands
        myband_save = meta_bands.band[0]
        del (meta_bands.band)
        meta_bands.band = []

        # create the band information; there are 4 output products per year
        # for the burned area dataset; add enough bands to cover the products
        # and years
        #    1. first date a burned area was observed (burned_area)
        #    2. number of times burn was observed (burn_count)
        #    3. number of good looks (good_looks_count)
        #    4. maximum probability for burned area (max_burn_prob)
        nproducts = 4
        nyears = end_year - start_year + 1
        nbands = nproducts * nyears
        for i in range (0, nbands):
            # add the new band
            myband = metadata_api.band()
            meta_bands.band.append(myband)

        # how many bands are there in the new XML file
        num_scene_bands =  len (meta_bands.band)
        logger.info('New XML file has {0} bands'.format(num_scene_bands))

        # loop through the products and years to create the band metadata
        band_count = 0
        for product in range (1, nproducts+1):
            for year in range (start_year, end_year+1):
                myband = meta_bands.band[band_count]
                myband.set_product("burned_area")
                myband.set_short_name("LNDBA")
                myband.set_data_type("INT16")
                myband.set_pixel_size(myband_save.get_pixel_size())
                myband.set_fill_value(fill_value)
                myband.set_nlines(nlines)
                myband.set_nsamps(nsamps)
                myband.set_app_version(self.burned_area_version)
                production_date = time.strftime("%Y-%m-%dT%H:%M:%S",
                    time.gmtime())
                myband.set_production_date (  \
                    datetime_.datetime.strptime(production_date,
                    '%Y-%m-%dT%H:%M:%S'))

                # clear some of the band-specific fields that don't apply for
                # this product (see metadata_api.py for the full list of
                # band-related values under class band)
                del (myband.source)
                myband.set_source(None)
                del (myband.saturate_value)
                myband.set_saturate_value(None)
                del (myband.scale_factor)
                myband.set_scale_factor(None)
                del (myband.add_offset)
                myband.set_add_offset(None)
                del (myband.radiance)
                myband.set_radiance(None)
                del (myband.reflectance)
                myband.set_reflectance(None)
                del (myband.thermal_const)
                myband.set_thermal_const(None)
                del (myband.bitmap_description)
                myband.set_bitmap_description(None)
                del (myband.class_values)
                myband.set_class_values(None)
                del (myband.qa_description)
                myband.set_qa_description(None)
                del (myband.calibrated_nt)
                myband.set_calibrated_nt(None)

                # handle the band-specific differences
                valid_range = metadata_api.valid_range()
                if product == 1:
                    name = "burned_area_%d" % year
                    long_name = "first DOY a burn was observed"
                    file_name = "burned_area_%d.img" % year
                    category = "image"
                    data_units = "day of year"
                    valid_range.min = 0
                    valid_range.max = 366
                    qa_description = "0: no burn observed" 

                elif product == 2:
                    name = "burn_count_%d" % year
                    long_name = "number of times a burn was observed"
                    file_name = "burn_count_%d.img" % year
                    category = "image"
                    data_units = "count"
                    valid_range.min = 0
                    valid_range.max = 366
                    qa_description = "0: no burn observed" 

                elif product == 3:
                    name = "good_looks_count_%d" % year
                    long_name = "number of good looks (pixels with good QA)"
                    file_name = "good_looks_count_%d.img" % year
                    category = "qa"
                    data_units = "count"
                    valid_range.min = 0
                    valid_range.max = 366
                    qa_description = "0: no valid pixels (water, cloud, " \
                        "snow, etc.)"

                elif product == 4:
                    name = "max_burn_prob_%d" % year
                    long_name = "maximum probability for burned area"
                    file_name = "max_burn_prob_%d.img" % year
                    category = "image"
                    data_units = "probability"
                    valid_range.min = 0
                    valid_range.max = 100
                    qa_description = "-9998: bad QA (water, cloud, snow, etc.)"

                myband.set_name(name)
                myband.set_long_name(long_name)
                myband.set_file_name(file_name)
                myband.set_category(category)
                myband.set_data_units(data_units)
                myband.set_valid_range(valid_range)
                myband.set_qa_description(qa_description)

                # increment the band counter
                band_count += 1

            # end for year
        # end for nproducts

        # write out a the XML file after validation
        # call the export with validation
        fd = open (output_xml_file, 'w')
        if fd == None:
            logger.error('Unable to open the output XML file ({0}) for writing.'
                         .format(output_xml_file))
            return ERROR

        metadata_api.export (fd, xml)
        fd.flush()
        fd.close()

        return SUCCESS


    def runAnnualBurnSummaries(self, stack_file=None, bp_dir=None, bc_dir=None,
        output_dir=None, start_year=None, end_year=None):
        """Processes the annual burn summaries for each year in the stack.
        Description: routine to process the annual burn summaries for each
            pixel.
        
        History:
          Created in 2013 by Jodi Riegle and Todd Hawbaker, USGS Rocky Mountain
              Geographic Science Center
          Updated on Dec. 2, 2013 by Gail Schmidt, USGS/EROS LSRD Project
              Modified to use argparser vs. optionparser, since optionparser
              is deprecated.
          Updated on Dec. 4, 2013 by Gail Schmidt, USGS/EROS LSRD Project
              Modified to resize the scene-based burn probabilities and
              classifications to the maximum geographic extent for stacking
              annual summaries.
          Updated on April 13, 2014 by Gail Schmidt, USGS/EROS LSRD Project
              Modified to utilize the ESPA raw binary file format.  And, the
              resizing of inputs will no longer occur, since the scene-based
              probabilities are already processed at the maximum geographic
              extent of the stack and therefore are one common extent.
              Removed the writing of data in blocks and instead process data
              as a line at a time.
          Updated on Feb. 11, 2015 by Gail Schmidt, USGS/EROS LSRD Project
              Modified the recfromcsv calls to not specify the datatype and to
              instead use the automatically-determined datatype from the read
              itself.

        Args:
          stack_file - input CSV file with information about the files to be
              processed.  this is generated as part of the seasonal summaries
              application.
          bp_dir - location of the burn probability files
          bc_dir - location of the burn classification files
          output_dir - location to write the output burn classifications
          start_year - starting year of the stack_file to process; default is
              to start with the lowest year + 1
          end_year - ending year of the stack_file to process; default is to end
              with the highest year
   
        Returns:
            ERROR - error running the annual burn summary application
            SUCCESS - successful processing
        """
        logger = logging.getLogger(__name__)  # Obtain logger for this module.

        # if no parameters were passed then get the info from the command line
        if stack_file is None:
            # get the command line argument for the input parameters
            parser = ArgumentParser(description='Run the annual burn summaries')
            parser.add_argument ('-f', '--stack_file', type=str,
                dest='stack_file',
                help='input file, csv delimited, each row contains '  \
                     'information about a landsat image',
                metavar='FILE', required=True)
            parser.add_argument ('-p', '--bp_dir', type=str, dest='bp_dir',
                help='input directory, location to find input burn '  \
                     'probability files',
                metavar='DIR', required=True)
            parser.add_argument ('-c', '--bc_dir', type=str, dest='bc_dir',
                help='input directory, location to find input burn '  \
                     'classification files',
                metavar='DIR', required=True)
            parser.add_argument ('-o', '--output_dir', type=str,
                dest='output_dir',
                help='output directory, location to write burn '  \
                     'classification files',
                metavar='DIR', required=True)
            parser.add_argument ('-s', '--start_year', type=int,
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

            options = parser.parse_args()

            # validate command-line options and arguments
            stack_file = options.stack_file
            if stack_file is None:
                logger.error('missing CSV stack file cmd-line argument')
                return ERROR

            bp_dir = options.bp_dir
            if input_dir is None:
                logger.error('missing input directory for the burn '
                             'probabilities cmd-line argument')
                return ERROR

            bc_dir = options.bc_dir
            if input_dir is None:
                logger.error('missing input directory for the burn '
                             'classifications cmd-line argument')
                return ERROR

            output_dir = options.output_dir
            if output_dir is None:
                logger.error('missing output directory cmd-line argument')
                return ERROR

            if options.start_year is not None:
                start_year = options.start_year

            if options.end_year is not None:
                end_year = options.end_year

        # validate options and arguments
        if start_year is not None:
            if (start_year < 1984):
                logger.error('start_year cannot begin before 1984: {0}'
                             .format(start_year))
                return ERROR

        if end_year is not None:
            if (end_year < 1984):
                logger.error('end_year cannot begin before 1984: {0}'
                             .format(end_year))
                return ERROR

        if (end_year is not None) & (start_year is not None):
            if end_year < start_year:
                logger.error('end_year ({0}) is less than start_year ({0})'
                             .format(end_year, start_year))
                return ERROR

        if not os.path.exists(stack_file):
            logger.error('CSV stack file does not exist: {0}'
                         .format(stack_file))
            return ERROR

        if not os.path.exists(bp_dir):
            logger.error('Burn probability directory does not exist: {0}'
                         .format(bp_dir))
            return ERROR

        if not os.path.exists(bc_dir):
            logger.error('Burn classificaton directory does not exist: {0}'
                         .format(bc_dir))
            return ERROR

        if not os.path.exists(output_dir):
            logger.info('Output directory does not exist: {0}. Creating ...'
                        .format(output_dir))
            os.makedirs(output_dir, 0755)

        # save the current working directory for return to upon error or when
        # processing is complete
        mydir = os.getcwd()
        logger.info('Changing directories for burn threshold processing: {0}'
                    .format(output_dir))
        os.chdir (output_dir)

        # start of threshold processing
        start_time0 = time.time()

        # open the stack file
        stack = numpy.recfromcsv(stack_file, delimiter=',', names=True)

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

        # given that all burn products in this temporal stack have the same
        # scene extents and projection information, just obtain that
        # information from the first file and use it for all of the files.
        # use the XML filename in the CSV file to obtain the burn probability
        # filename
        xml_file = stack2['file_'][0]
        fname = os.path.basename(xml_file).replace  \
            ('.xml','_burn_probability.img')
        bp_file = bp_dir + '/' + fname
        if not os.path.exists(bp_file):
            logger.error('burn probability file does not exist: {0}'
                         .format(bp_file))
            os.chdir (mydir)
            return ERROR

        bp_dataset = gdal.Open(bp_file)
        if bp_dataset is None:
            logger.error('Failed to open bp file: {0}'.format(bp_file))
            os.chdir (mydir)
            return ERROR
        
        bp_band = bp_dataset.GetRasterBand(1)
        if bp_band is None:
            logger.error('Failed to open bp band 1 from {0}'.format(bp_file))
            os.chdir (mydir)
            return ERROR
        
        geotrans = bp_dataset.GetGeoTransform()
        if geotrans is None:
            logger.error('Failed to obtain the GeoTransform info from {0}'
                         .format(bp_file))
            return ERROR

        prj = bp_dataset.GetProjectionRef()
        if prj is None:
            logger.error('Failed to obtain the ProjectionRef info from {0}'
                         .format(bp_file))
            return ERROR

        nrow = bp_dataset.RasterYSize
        ncol = bp_dataset.RasterXSize
        if (nrow is None) or (ncol is None):
            logger.error('Failed to obtain the RasterXSize and RasterYSize'
                         ' from {0}'.format(bp_file))
            return ERROR

        nodata = bp_band.GetNoDataValue()
        if nodata is None:
            nodata = -9999
            logger.warn('Failed to obtain the NoDataValue from {0}.'
                        '  Using {1}'.format(bp_file, nodata))

        # close the file
        bp_band = None
        bp_dataset = None

        # create the ENVI driver for output data
        driver = gdal.GetDriverByName('ENVI')
    
        # process the data for the years specified
        # create images for:
        #    1. first date a burned area was observed (burned_area)
        #    2. number of times burn was observed (burn_count)
        #    3. number of good looks (good_looks_count)
        #    4. maximum probability for burned area (max_burn_prob)
    
        # loop through the years in the stack
        logger.info('Processing burn files for {0}-{1}'
                    .format(start_year, end_year))
        for year in range(start_year, end_year+1):
            logger.info('########################################################')
            logger.info('Processing {0} ...'.format(year))

            stack_mask = stack2['year'] == year
            stack3 = stack2[stack_mask]

            # initialize the input and output datasets
            input_datasets = numpy.empty( (stack2.shape[0],2), dtype=object )
            input_bands = numpy.empty( (stack2.shape[0],2), dtype=object )
            
            output_datasets = numpy.empty((4), dtype=object)
            output_bands = numpy.empty((4), dtype=object)
        
            # open the input datasets - 1st band is burn probability,
            # 2nd band is burn classification
            for i in range(0, stack3.shape[0]):
                xml_file = stack3['file_'][i]
                
                # construct the burn probability and classification filenames
                # from the XML filenames in the CSV
                fname = os.path.basename(xml_file).replace  \
                    ('.xml','_burn_probability.img')
                bp_file = bp_dir + '/' + fname
                if not os.path.exists(bp_file):
                    logger.error('burn probability file does not exist: {0}'
                                 .format(bp_file))
                    os.chdir (mydir)
                    return ERROR

                logger.info('    Reading {0} ...'.format(bp_file))
                input_datasets[i,0] = gdal.Open(bp_file)
                input_bands[i,0] = input_datasets[i,0].GetRasterBand(1)

                fname = os.path.basename(xml_file).replace  \
                    ('.xml','_burn_class.img')
                bc_name = bc_dir + '/' + fname
                if not os.path.exists(bc_name):
                    logger.error('burn classification file does not exist: {0}'
                                 .format(bc_name))
                    os.chdir (mydir)
                    return ERROR

                logger.info('    Reading {0} ...'.format(bc_name))
                input_datasets[i,1] = gdal.Open(bc_name)
                input_bands[i,1] = input_datasets[i,1].GetRasterBand(1)

            # open the output datasets
            # first date of burned area (burned_area)
            fname = output_dir + '/burned_area_' + str(year) + '.img'
            output_datasets[0] = driver.Create(fname, ncol, nrow, 1, \
                gdal.GDT_Int16)
            output_datasets[0].SetGeoTransform(geotrans)
            output_datasets[0].SetProjection(prj)
            output_bands[0] = output_datasets[0].GetRasterBand(1)
            output_bands[0].SetNoDataValue(nodata)
            
            # count of times a pixel was burned (burn_count)
            fname = output_dir + '/burn_count_' + str(year) + '.img'
            output_datasets[1] = driver.Create(fname, ncol, nrow, 1,  \
                gdal.GDT_Int16)
            output_datasets[1].SetGeoTransform(geotrans)
            output_datasets[1].SetProjection(prj)
            output_bands[1] = output_datasets[1].GetRasterBand(1)
            output_bands[1].SetNoDataValue(nodata)
            
            # count of good looks (good_looks_count)
            fname = output_dir + '/good_looks_count_' + str(year) + '.img'
            output_datasets[2] = driver.Create(fname, ncol, nrow, 1,  \
                gdal.GDT_Int16)
            output_datasets[2].SetGeoTransform(geotrans)
            output_datasets[2].SetProjection(prj)
            output_bands[2] = output_datasets[2].GetRasterBand(1)
            output_bands[2].SetNoDataValue(nodata)
            
            # maximum burn probability (max_burn_prob)
            fname = output_dir + '/max_burn_prob_' + str(year) + '.img'
            output_datasets[3] = driver.Create(fname, ncol, nrow, 1,  \
                gdal.GDT_Int16)
            output_datasets[3].SetGeoTransform(geotrans)
            output_datasets[3].SetProjection(prj)
            output_bands[3] = output_datasets[3].GetRasterBand(1)
            output_bands[3].SetNoDataValue(nodata)

            # loop through the lines in the images
            for y in range (0, nrow):
                # create the arrays to hold input and output data (one line)
                input_data = numpy.empty((stack3.shape[0], 2, 1, ncol),  \
                    dtype=numpy.int16)
                input_data.fill(nodata)

                output_data = numpy.empty((1, 4, 1, ncol), dtype=numpy.int16)
                output_data.fill(nodata)

                # read input data for burn probs and burn classes
                for i in range(0, stack3.shape[0]):
                    input_data[i,0,:,:] = input_bands[i,0].ReadAsArray(  \
                        0, y, ncol, 1)
                    input_data[i,1,:,:] = input_bands[i,1].ReadAsArray(  \
                        0, y, ncol, 1)

                # find the maximum burn probability (using burn prob)
                bp_max = numpy.apply_over_axes(numpy.max, input_data[:,0,:,:], \
                    axes=[0])[0,:,:]

                # find the count of burns - how many times a pixel burned
                # (using burn class)
                bc = numpy.apply_over_axes(numpy.sum,  \
                    input_data[:,1,:,:] >= 1, axes=[0])[0,:,:]
                bc[bp_max == nodata] = nodata
    
                # find the first date of burn (using burn class)
                bdi = numpy.apply_over_axes(numpy.argmax,  \
                    input_data[:,1,:,:] >= 1, axes=[0])[0,:,:]
                
                # convert bdi to julian date
                bd = stack3['julian'][bdi]
                bd[bc == 0] = 0
                bd[bp_max == nodata] = nodata
                
                # find the number of good looks (using burn class)
                gc = numpy.apply_over_axes(numpy.sum,  \
                    input_data[:,1,:,:] >= 0, axes=[0])[0,:,:]
                gc[bp_max == nodata] = nodata
            
                # write output data for the burned area DOY, burn count, good
                # looks count, and the maximum burn probability
                output_bands[0].WriteArray(bd, xoff=0, yoff=y)
                output_bands[1].WriteArray(bc, xoff=0, yoff=y)
                output_bands[2].WriteArray(gc, xoff=0, yoff=y)
                output_bands[3].WriteArray(bp_max, xoff=0, yoff=y)

            # close the input datasets 
            for i in range(0, stack3.shape[0]):
                input_datasets[i,0] = None
                input_bands[i,0] = None
                input_datasets[i,1] = None
                input_bands[i,1] = None

            # close the output datasets 
            output_datasets[0] = None
            output_datasets[1] = None
            output_datasets[2] = None
            output_datasets[3] = None
            output_bands[0] = None
            output_bands[1] = None
            output_bands[2] = None
            output_bands[3] = None
        # end for year

        # remove the .img.aux.xml files that are generated by GDAL as these
        # won't be delivered to the user
        rm_files = glob.glob (output_dir + '/burned_area_*.img.aux.xml')
        for file in rm_files:
            logger.info('Remove: {0}'.format(file))
            os.remove (os.path.join (file))

        rm_files = glob.glob (output_dir + '/burn_count_*.img.aux.xml')
        for file in rm_files:
            logger.info('Remove: {0}'.format(file))
            os.remove (os.path.join (file))

        rm_files = glob.glob (output_dir + '/good_looks_count_*.img.aux.xml')
        for file in rm_files:
            logger.info('Remove: {0}'.format(file))
            os.remove (os.path.join (file))

        rm_files = glob.glob (output_dir + '/max_burn_prob_*.img.aux.xml')
        for file in rm_files:
            logger.info('Remove: {0}'.format(file))
            os.remove (os.path.join (file))

        # create the output XML file which contains information for each of
        # the bands: burned area date, burn count, good looks count, and the
        # maximum burn probability
        logger.info('Creating output XML file for burned area ...')
        xml_file = stack2['file_'][0]
        fname = os.path.basename(xml_file).replace  \
            ('.xml','_burn_probability.img')
        output_xml_file = "burned_area_%d_%d.xml" % (start_year, end_year)
        status = self.createXML (xml_file, output_xml_file, start_year,
            end_year, nodata, fname)
        if status != SUCCESS:
            logger.error('Failed to write the output XML file: {0}'
                         .format(output_xml_file))
            return ERROR

        # successful completion.  return to the original directory.
        logger.info('Completion of annual burn summaries.')

        os.chdir (mydir)
        return SUCCESS

######end of AnnualBurnSummary class######

if __name__ == "__main__":
    # setup the default logger format and level. log to STDOUT.
    logging.basicConfig(format=('%(asctime)s.%(msecs)03d %(process)d'
                                ' %(levelname)-8s'
                                ' %(filename)s:%(lineno)d:'
                                '%(funcName)s -- %(message)s'),
                        datefmt='%Y-%m-%d %H:%M:%S',
                        level=logging.INFO)

    sys.exit (AnnualBurnSummary().runAnnualBurnSummaries())
