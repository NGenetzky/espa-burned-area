#! /usr/bin/env python
import multiprocessing, Queue
import time
import logging

#if temporalBAStack is already imported from a higher level script, then
#this import is not needed
#from process_temporal_ba_stack import temporalBAStack
 
class parallelSceneWorker(multiprocessing.Process):
    """Runs the scene resampling in parallel for a stack of scenes.
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
        logger = logging.getLogger(__name__)
        while not self.kill_received:
            # get a task
            try:
                xml_file = self.work_queue.get_nowait()
            except Queue.Empty:
                break
 
            # process the scene
            logger.info('Processing {0} ...'.format(xml_file))
            status = SUCCESS
            status = self.stackObject.sceneResample (xml_file)
            if status != SUCCESS:
                logger.error('Error resampling the surface reflectance bands '
                            'in the XML file ({0}). Processing will terminate'
                            ..format(xml_file))
 
            # store the result
            self.result_queue.put(status)


class parallelSummaryWorker(multiprocessing.Process):
    """Runs the seasonal summaries in parallel for a temporal stack.
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
        logger = logging.getLogger(__name__)
        while not self.kill_received:
            # get a task
            try:
                year_season = self.work_queue.get_nowait()
            except Queue.Empty:
                break
 
            # process the scene
            year = int (year_season[0])
            season = year_season[1]
            logger.info('Processing year {0}, season {1} ...'
                        .format(year, season))
            status = SUCCESS
            status = self.stackObject.generateYearSeasonalSummaries (year,
                season)
            if status != SUCCESS:
                logger.error('Error processing seasonal summaries for year'
                             ' {0}, season {1}. Processing will terminate.'
                             .format(year, season))
 
            # store the result
            self.result_queue.put(status)


class parallelMaxWorker(multiprocessing.Process):
    """Runs the annual maximums in parallel for a temporal stack.
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
        logger = logging.getLogger(__name__)
        while not self.kill_received:
            # get a task
            try:
                year = self.work_queue.get_nowait()
            except Queue.Empty:
                break
 
            # process the scene
            logger.info('Processing year {0} ...'.format(year))
            status = SUCCESS
            status = self.stackObject.generateYearMaximums (year)
            if status != SUCCESS:
                logger.info('Error processing maximums for year {0}.'
                            ' Processing will terminate.'.format(year))
 
            # store the result
            self.result_queue.put(status)

