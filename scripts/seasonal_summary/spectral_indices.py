# use true division so we don't have to worry about scalar divided by scalar
# not being a floating point
from __future__ import division
from numpy import *

########################################################################
# Description: The various modules calculate spectral indices utilized as
# part of the burned area mapping processing.
#
# Returns:
#    Spectral index array
#
# Notes:
#   1. Many of the calculations kick off a divide by zero error, if the
#      pixel values are both zero.  The divide by zero warning is shut
#      off for this processing.
#   2. To deal with infinity and NaN values, the nan_to_num method is used
#      after the spectral index calculation to handle any associated NaNs
#      due to the denomerator being zero.
#######################################################################

# Normalized Burn Ratios
# From Key and Benson 1999, Measuring and remote sensing of burn severity. In Proceedings of the Joint Fire Science Conference and Workshop, vol. II, Boise, ID, 15-17 June 1999. University of Idaho and International Association of Wildland Fire
def NBR(b4, b7, nodata=-9999):
    """Computes the normalized burn index.
    """
    x = (b4-b7) / (b4+b7)
    nan_to_num(x)
    x[(b4==nodata) | (b7==nodata)] = nodata
    return(x)

def NBR2(b5, b7, nodata=-9999):
    """Computes the normalized burn index 2.
    """
    seterr(divide= 'ignore')  # ignore divide by zero
    x = (b5-b7) / (b5+b7)
    nan_to_num(x)
    x[(b5==nodata) | (b7==nodata)] = nodata
    return(x)

# Normalized Difference Moisture Index
def NDMI(b4, b5, nodata=-9999):
    """Computes the normalized difference moisture index.
    """
    seterr(divide= 'ignore')  # ignore divide by zero
    x = (b4-b5) / (b4+b5)
    nan_to_num(x)
    x[(b4==nodata) | (b5==nodata)] = nodata
    return(x)

# Normalized Difference Vegetation Index
# From Rouse et al. 1973, Monitoring vegetation systems in the Great Plains with ERTS. In: Proc. Third ERTS Symposium, NASA, SP-351, vol. 1, pp. 309-317
def NDVI(b3, b4, nodata=-9999):
    """Computes the normalized difference vegetation index.
    """
    seterr(divide= 'ignore')  # ignore divide by zero
    x = (b4-b3) / (b4+b3)
    nan_to_num(x)
    x[(b3==nodata) | (b4==nodata)] = nodata
    return(x)

# Char Soil Index
# Smith et al. 2005, Testing the potential of multi-spectral remote sensing for retrospectively estimating fire severity in African savanna environments. Remote Sensing of Environment 97 (1):92-115
def CSI(b4, b5, nodata=-9999):
    """Computes the char soil index.
    """
    seterr(divide= 'ignore')  # ignore divide by zero
    x = b4/b5
    nan_to_num(x)
    x[(b4==nodata) | (b5==nodata)] = nodata
    return(x)

# Mid-Infrared Burn Index 
# Trigg and Flasse 2001, An evaluation of different bi-spectral spaces for discriminating burned shrub savanna.  International Journal of Remote Sensing 22(13):2641-2647
def MIRBI(b5, b7, nodata=-9999):
    """Computes the mid-infrared burn index.
    """
    x = ( (10*b7) - (9.5*b5) + 2 )
    x[(b5==nodata) | (b7==nodata)] = nodata
    return(x)

# Burned Area Index
# Martin et al. 2005, Performanec of a burned-area index (BAIM) for mapping Mediterranean burned scars from MODIS data. In J. Ria, F. Perez-Cabello, and E. Chuvieco (Editors), Proceedings of the 5th International Workshop on Remote Sensing and GIS Applications to Forest Fire Management: Fire Effects Assessment (pp. 193-198). Paris: Universidad de Zaragoza, GOFC-GOLD, EARSel.
def BAI(b3, b4, nodata=-9999):
    """Computes the burned area index.
    """
    x = 1 / ( pow(b4-0.06,2) + pow(b3-0.1,2) )
    x[(b4==nodata) | (b3==nodata)] = nodata
    return(x)

# Martin et al. 2005, 
def BAIM(b4, b5, nodata=-9999):
    """Computes the burned area index for mapping Mediterranean burn scars.
    """
    x = 1 / ( pow(b4-0.05,2) + pow(b5-0.2,2) )
    x[(b4==nodata) | (b5==nodata)] = nodata
    return(x)

def BAIM2(b4, b7, nodata=-9999):
    """Computes the burned area index for mapping Mediterranean burn scars 2.
    """
    x = ( 1 / ( pow(b4-0.05,2) + pow(b7-0.2,2) ) )
    x[(b4==nodata) | (b7==nodata)] = nodata
    return(x)

# Soil-Adjusted Vegetation Index
# Huete 1998, A soil adjusted vegetation index (SAVI). Remote Sensing of Environment, 25(3):295-309
def SAVI(b3, b4, nodata=-9999):
    """Computes the soil adjusted vegetation index.
    """
    seterr(divide= 'ignore')  # ignore divide by zero
    x = ( 1.5 * (b4-b3) / (b4 + b3 + 0.5) )
    nan_to_num(x)
    x[(b3==nodata) | (b4==nodata)] = nodata
    return(x)

# Enhanced Vegetaion Index
# Huete et al. 2002, Overview of the radiometric and biophysical performance of the MODIS vegetation indices. Remote Sensing of Environment, 83(1-2):195-213
# Enhanced Vegetation Index 2
def EVI(b1, b3, b4, nodata=-9999):
    """Computes the enhanced vegetation index.
    """
    seterr(divide= 'ignore')  # ignore divide by zero
    x = ( 2.5 * (b4-b3) / (b4 + 6*b3 - 7.5*b1 + 1) )
    nan_to_num(x)
    x[(b1==nodata) | (b3==nodata)| (b4==nodata)] = nodata
    return(x)

# Jiang et al. 2008, Development of a two-band enhanced vegetation index without a blue band. Remote Sensing of Environment 112(10):3833-3845
def EVI2(b3, b4, nodata=-9999):
    """Computes the enhanced vegetation index without a blue band.
    """
    seterr(divide= 'ignore')  # ignore divide by zero
    x = ( 2.5 * (b4-b3) / (b4 + 2.4*b3 + 1) ) 
    nan_to_num(x)
    x[(b3==nodata) | (b4==nodata)] = nodata
    return(x)

