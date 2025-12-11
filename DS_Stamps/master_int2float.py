#!/usr/bin/env python
import os,sys,time
import math
import numpy as np
from numpy import *

import scipy as Sci
import scipy.linalg
import scipy.io as sio
#import pandas as pd
#import matplotlib.pyplot as plt

import gdal,gdalconst
from gdalconst import *

def usage():
    print '\nUsage: python  master_int2float.py  resFilename '
    print '  where                                           '
    print '        resFilename    is the .res file of burst  ' 
    print '                                                  '
    print '                                                  '
    print ' This python applies cpxint16 transform  to  cpxcomplex64'
    print '                                                          '
    print ' for example                                              '
    print ' python   master_int2float.py  master.res                 '
    print '                                                          '
    print ' Python: Wu Wenhao   Wuhan university   QQ:460249274'
try:
    resFilename          = sys.argv[1]
     
except:
    print 'Unrecognized input'
    usage()
    sys.exit(1)
if len(sys.argv) == 2:
    print 'Expected input'
else:
    print 'Unrecognized input'
    usage()
    sys.exit(1) 


################################################################################
#function
def get_parameter(First_param,file_name,format_flag=1,Second_param=None,Third_param=None):
    Read_contine_flag=0
    orbit_info=""
    value=None
    for line in open(file_name):
        if format_flag==1:
            if not (line.find(First_param)):
                index=line.find(':')
                value=(line[(index+1):].strip(' \n\t'))
                return value 

        if format_flag==2:
            if not(line.find(Second_param)):
                Read_contine_flag=1             
            if (Read_contine_flag==1) and (not (line.find(First_param))):  #Be careful
                index=line.find(':')
                value=(line[(index+1):].strip(' \n\t'))
                continue               
            if Read_contine_flag==1  and (not(line.find(Third_param))):  #Be careful               
                Read_contine_flag=0         
                return value
               

        if format_flag==3:            
            if not (line.find(First_param)):
                index=line.find(':')
                pixel_time=(line[(index+1):].strip(' \n\t')).split(' ')[1].split(':')
                return pixel_time  


        if format_flag==4:
            
            if not (line.find(First_param)):
                index=line.find(':')
                
                value=int(line[(index+1):].strip(' \n\t'))
                Read_contine_flag=1;
                #print line
                continue                 
            if (Read_contine_flag>=1) :
               
                orbit_info=orbit_info+line
                Read_contine_flag=Read_contine_flag+1
                if (Read_contine_flag==(value+1)):
                    return orbit_info  
################################################################################


###############################################################################
def freadbk(path_file,line_start=1, Pixels_start=1,nofLines1=None,nofPixels1=None):
    #Driver
    driver=gdal.GetDriverByName('MFF')
    driver.Register()
    gdal.AllRegister()
    thisBurstData_file=gdal.Open(path_file,GA_ReadOnly)
    if thisBurstData_file is None:
        print 'Could not open'+Path_MFF_HDR
        sys.exit(1)
    #print 'Driver: ', thisBurstData_file.GetDriver().ShortName,'/', \
    #      thisBurstData_file.GetDriver().LongName
    #print 'Size is ',thisBurstData_file.RasterXSize,'x',thisBurstData_file.RasterYSize, \
    #      'x',thisBurstData_file.RasterCount
    #print 'Projection is ',thisBurstData_file.GetProjection()
    #geotransform = thisBurstData_file.GetGeoTransform()
    #if not geotransform is None:
    #    print 'Origin = (',geotransform[0], ',',geotransform[3],')'
    #    print 'Pixel Size = (',geotransform[1], ',',geotransform[5],')'

    cint_srd=thisBurstData_file.GetRasterBand(1)
    print 'Band Type=',gdal.GetDataTypeName(cint_srd.DataType)

    if cint_srd.GetOverviewCount() > 0:
            print 'Band has ', cint_srd.GetOverviewCount(), ' overviews.'
    thisBurstData= cint_srd.ReadAsArray(int(Pixels_start-1),int(line_start-1),nofPixels1,nofLines1)
    return thisBurstData
##################################################################################
dataFilename=get_parameter('Data_output_file',resFilename,1)
print "dataFilename =",dataFilename
masterDataFormat = get_parameter('Data_output_format',resFilename,2,'*_Start_crop','* End_crop:_NORMAL')
print "masterDataFormat =",masterDataFormat

if masterDataFormat=='complex_real4':
    dataFormat = 'cpxfloat32'
    print 'dataFormat =',dataFormat
else:
    dataFormat = 'cpxint16'


l0 = int(get_parameter('First_line (w.r.t. original_image)',resFilename,2,'*_Start_crop','* End_crop:_NORMAL'))

lN = int(get_parameter('Last_line (w.r.t. original_image)',resFilename,2,'*_Start_crop','* End_crop:_NORMAL'))

p0 = int(get_parameter('First_pixel (w.r.t. original_image)',resFilename,2,'*_Start_crop','* End_crop:_NORMAL'))

pN = int(get_parameter('Last_pixel (w.r.t. original_image)',resFilename,2,'*_Start_crop','* End_crop:_NORMAL'))

#% Get resampled Slv size
Naz_res = lN-l0+1
Nrg_res = pN-p0+1
############################################################################################
#%% Read the  image 

Path_MFF_HDR   =dataFilename.split('.')[0]+'.hdr'

if dataFormat == 'cpxfloat32':
    Link_DATA  =dataFilename.split('.')[0]+'.x00'
elif dataFormat == 'cpxint16':
    Link_DATA  =dataFilename.split('.')[0]+'.j00'
else:
    print 'Unrecognized input'
    usage()
    sys.exit(1)
outStream      = open(Path_MFF_HDR,'w')
outStream.write('IMAGE_FILE_FORMAT = MFF\n')
outStream.write('FILE_TYPE = IMAGE\n')
outStream.write('IMAGE_LINES = %d\n' % int(Naz_res))
outStream.write('LINE_SAMPLES = %d\n'% int(Nrg_res))
outStream.write('BYTE_ORDER = LSB\n')
outStream.write('END\n')
outStream.close()

if (os.path.isfile(Link_DATA)):
    os.remove(Link_DATA)

os.symlink(dataFilename,Link_DATA)

slc=freadbk(Path_MFF_HDR,1, 1,int(Naz_res),int(Nrg_res))

slc_complex64=slc.astype(np.complex64)

###################################################################
dataFilename=dataFilename+'.cr4'
masterDataFormat='complex_real4';
cols = slc_complex64.shape[1]
rows = slc_complex64.shape[0]
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = open(dataFilename,'wb')
if masterDataFormat=='complex_real4':
      #It's not useful ,may be a bug
    if slc_complex64.dtype=='complex64':
        print 'data format is right\n!!!!'
    else:
        print "Waring !The data format may be wrong !"
    for i_temp in arange(0,rows):
        fid.write(slc_complex64[i_temp,:])
else:
     save_row_sl_deramped=zeros((2*cols,1),dtype='int16')
     for i_temp in arange(0,rows):        
         save_row_sl_deramped[0::2]= (slc_complex64[i_temp,:].real.astype(np.int16)).reshape(-1,1,order='F').copy()
         save_row_sl_deramped[1::2]= (slc_complex64[i_temp,:].imag.astype(np.int16)).reshape(-1,1,order='F').copy() 
         fid.write(save_row_sl_deramped)
     if save_row_sl_deramped.dtype=='int16':
         print '\n Good ! data format is right!!!!\n'
     else:
         print 'Waring !The data format may be wrong !'
     del save_row_sl_deramped

fid.close()
print "\nThe data_format transform operation completed!\n"
