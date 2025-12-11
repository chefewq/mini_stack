#!/usr/bin/env python
import os,sys,time
import math
from scipy.stats import ks_2samp
import numpy as np
from numpy import *
import scipy as Sci
import scipy.linalg
import scipy.io as sio
from  scipy  import ndimage
import matplotlib.pyplot as plt
import struct,shutil
import gdal,gdalconst
from gdalconst import *
from scipy import stats
from scipy import signal
import multiprocessing
import re
import scipy.linalg
import scipy.io as sio
import linecache
#from osgeo import gdal,gdalconst
# NOT FETOOLS

#Version:20181101
def usage():
    print '\nUsage: python concatenateIfg.py <burstList> [Degree] [do_ESD_flag] [Plotflag] '
    print '  where burstList  is the list of bursts you want to merge        '
    print '        Degree     is the  order of the polynomial interpolation  '
    print '        do_ESD_flag is a boolean to save the result of azimuth shift by ESD '
    print '        Plotflag    is a boolean to plot, if true, plot the ESD  interferogram  '
    print ' For example (default instructon)                                '
    print '    python KS_filter.py 1:9  1(default)                                '
    print ' or\n    python KS_filter.py 1:9  2 True True                           '
    print ' \n                                                              '
    print " Freek's Suggestion                                               "
    print ' Python: Wu Wenhao    QQ:460249274              '
try:
    burstList  = sys.argv[0] 
    
except:
    print 'Unrecognized input'
    usage()
    sys.exit(1)
    value=None

#******************************************************************************************#
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
            if  not (line.find(Second_param)):
                Read_contine_flag=1             
            if (Read_contine_flag==1) and (not (line.find(First_param))):  #Be careful
                index=line.find(':')
                value=(line[(index+1):].strip(' \n\t'))                             
            if  (not (line.find(Third_param))):  #Be careful                      
                Read_contine_flag=0              
                return value
               

        if format_flag==3:            
            if not (line.find(First_param)):
                index=line.find(':')
                pixel_time=(line[(index+1):].strip(' \n\t')).split(' ')[1].split(':')
                return pixel_time  


        if format_flag==4:   # for orbit

            
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
        if format_flag==5:   # for corgistration

            
            if not (line.find(First_param)):
                index=line.find(':')
                
                value=int(line[(index+1):].strip(' \n\t'))
                Read_contine_flag=1;
                print line
                continue                 
            if (Read_contine_flag>=1) :
               
                orbit_info=orbit_info+line
                Read_contine_flag=Read_contine_flag+1
                if (Read_contine_flag==(value+2)):
                    return orbit_info 
        if format_flag==6:   # for cpm

            
            if not (line.find(First_param)):
             
                Read_contine_flag=1;
                #print line
                continue                 
            if (Read_contine_flag>=1) :
               
                orbit_info=orbit_info+line
                Read_contine_flag=Read_contine_flag+1
                nofCoef= 0.5*((int(Second_param)+1)**2+int(Second_param)+1)
                if (Read_contine_flag==int(nofCoef+1)):
                    return orbit_info
 
#---------------------------------------------------------------#
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
    geotransform = thisBurstData_file.GetGeoTransform()
    #if not geotransform is None:
    #    print 'Origin = (',geotransform[0], ',',geotransform[3],')'
    #    print 'Pixel Size = (',geotransform[1], ',',geotransform[5],')'
    cint_srd=thisBurstData_file.GetRasterBand(1)
    print 'Band Type=',gdal.GetDataTypeName(cint_srd.DataType)

    if cint_srd.GetOverviewCount() > 0:
            print 'Band has ', cint_srd.GetOverviewCount(), ' overviews.'
    thisBurstData= cint_srd.ReadAsArray(int(Pixels_start-1),int(line_start-1),nofPixels1,nofLines1)
    return thisBurstData
###############################################################################
#******************************************************************************
def Linux_cmd(python_instruction):
    #python_instruction='source '+work_dir+'/'+Configfile
    print 'python_instruction=',python_instruction
    with os.popen(python_instruction,"r") as p:
        cmd_profile_copy=p.read()
    #print 'cmd_profile_copy=',cmd_profile_copy       
    return cmd_profile_copy

#******************************************************************************
################################################################################################
def main():
     
    cores=multiprocessing.cpu_count()
    print 'Total cores is ', cores
    #pool=multiprocessing.Pool(processes=(cores-1))
  
    
    resFilename='master.res'
    l0 = int(get_parameter('First_line (w.r.t. original_image)',resFilename,1))
    lN = int(get_parameter('Last_line (w.r.t. original_image)',resFilename,1))
    p0 = int(get_parameter('First_pixel (w.r.t. original_image)',resFilename,1))
    pN = int(get_parameter('Last_pixel (w.r.t. original_image)',resFilename,1))
    nofLines=lN-l0+1
    nofPixels=pN-p0+1
    dataFormat = 'cpxfloat32'    
    ##############################################################
    ##Parameter of filter Windows
    #############################################################
   ##############################################################
    ##Parameter of filter Windows
    #############################################################
    work_dir=os.getcwd();
    print 'work_dir=',work_dir.split('/')[-1].split('_')[-1]
    Master_date=work_dir.split('/')[-1].split('_')[-1]
    
    f=open('day.1.in')
    cint_day_file=f.readlines()
   
    cint_day_file.append(Master_date)
    for file_cint  in cint_day_file:

        Output_file_path    =os.getcwd()+'/'+'mini_stack_compress'+'/'+file_cint.strip()
        if not os.path.exists(Output_file_path):
            continue
        SLC_path_file  =os.getcwd()+'/'+file_cint.strip()+'/'+'slave_res.slc' 
        """
        if os.path.isfile(SLC_path_file):
    
            Path_MFF_HDR   =os.getcwd()+'/'+file_cint.strip()+'/'+'slave_rsmp'+'.hdr'
            Link_DATA      =os.getcwd()+'/'+file_cint.strip()+'/'+'slave_rsmp'+'.x00'
      
            print 'Path_MFF_HDR=',Path_MFF_HDR 
      
            if os.path.exists(Path_MFF_HDR)==True:
            #os.remove(Path_MFF_HDR)
                print 'Path_MFF_HDR=',Path_MFF_HDR
            else:       
                outStream      = open(Path_MFF_HDR,'w')
                outStream.write('IMAGE_FILE_FORMAT = MFF\n')
                outStream.write('FILE_TYPE = IMAGE\n')
                outStream.write('IMAGE_LINES = %d\n' % int(nofLines))
                outStream.write('LINE_SAMPLES = %d\n'% int(nofPixels))
                outStream.write('BYTE_ORDER = LSB\n')
                outStream.write('END\n')
                outStream.close()

            if (os.path.isfile(Link_DATA)):
                 #os.remove(Link_DATA)
                print  'Link_DATA=' ,Link_DATA
            else:       
                os.symlink(SLC_path_file ,Link_DATA)     
            Single_SLC_Data=freadbk(Path_MFF_HDR,1,1,int(nofLines),int(nofPixels)).copy() 
        else:
        """
        Single_SLC_Data= np.zeros((int(nofLines),int(nofPixels)),dtype=np.complex64)
             #Single_SLC_Data=
        print  'Single_SLC_Data=',Single_SLC_Data.shape  

        #patch_Name='PATCH_1'
        #if os.path.exists(patch_Name)==True:
        for patch_Name in open('patch.list'):
            patch_Data_Dimension=patch_Name.strip()+'/'+'patch.in'
            file_patch=open(patch_Data_Dimension);
            Patch_First_Pixel= int(file_patch.readline().strip())-1
            Patch_Last_Pixel = int(file_patch.readline().strip())-1
            Patch_First_Line = int(file_patch.readline().strip())-1
            Patch_Last_Line  = int(file_patch.readline().strip())-1
            Patch_nofPixels  = Patch_Last_Pixel- Patch_First_Pixel+1
            Patch_nofLines   = Patch_Last_Line - Patch_First_Line +1

            print Patch_First_Pixel
            print Patch_Last_Pixel
            print Patch_First_Line
            print Patch_Last_Line
            print Patch_nofPixels
            print Patch_nofLines
            patch_Data_Dimension_noover=patch_Name.strip()+'/'+'patch_noover.in'
            file_patch=open(patch_Data_Dimension_noover);
            Patch_noover_First_Pixel= int(file_patch.readline().strip())-1
            Patch_noover_Last_Pixel = int(file_patch.readline().strip())-1
            Patch_noover_First_Line = int(file_patch.readline().strip())-1
            Patch_noover_Last_Line  = int(file_patch.readline().strip())-1
            Patch_noover_nofPixels  = Patch_noover_Last_Pixel- Patch_noover_First_Pixel+1
            Patch_noover_nofLines   = Patch_noover_Last_Line - Patch_noover_First_Line +1
        
           
            Read_file      =os.getcwd()+'/'+'mini_stack_compress'+'/'+file_cint.strip()+'/'+'resample_ministack'+'_'+patch_Name.strip()+'.raw'            
            print 'Read_file=',Read_file

            if os.path.exists(Read_file)==True:
 
                Path_MFF_HDR   =os.getcwd()+'/'+'mini_stack_compress'+'/'+file_cint.strip()+'/'+'resample_ministack'+'_'+patch_Name.strip()+'.hdr'
                Link_DATA      =os.getcwd()+'/'+'mini_stack_compress'+'/'+file_cint.strip()+'/'+'resample_ministack'+'_'+patch_Name.strip()+'.x00'

                #if os.path.exists(Path_MFF_HDR)==True:
                #    #os.remove(Path_MFF_HDR)
                    
                if os.path.exists(Path_MFF_HDR)==False:       
                    outStream      = open(Path_MFF_HDR,'w')
                    outStream.write('IMAGE_FILE_FORMAT = MFF\n')
                    outStream.write('FILE_TYPE = IMAGE\n')
                    outStream.write('IMAGE_LINES = %d\n' % int(Patch_nofLines))
                    outStream.write('LINE_SAMPLES = %d\n'% int(Patch_nofPixels))
                    outStream.write('BYTE_ORDER = LSB\n')
                    outStream.write('END\n')
                    outStream.close()
                          
                if (os.path.isfile(Link_DATA)):
                    #os.remove(Link_DATA)
                    print  'Link_DATA=' ,Link_DATA
                else:       
                    os.symlink(Read_file ,Link_DATA)     
                KS_SLC_Data=freadbk(Path_MFF_HDR,1,1,int(Patch_nofLines),int(Patch_nofPixels)).copy() 
                print 'KS_SLC_Data=',KS_SLC_Data.shape
                print Patch_noover_First_Line-Patch_First_Line
                print Patch_noover_First_Line-Patch_First_Line+Patch_noover_nofLines
                slc_cut=KS_SLC_Data[(Patch_noover_First_Line-Patch_First_Line):(Patch_noover_Last_Line-Patch_First_Line+1),\
                         (Patch_noover_First_Pixel-Patch_First_Pixel):(Patch_noover_Last_Pixel-Patch_First_Pixel+1)]  
                #print 'slc_cut=',slc_cut.shape
                #print 'Patch_noover_nofPixels=',Patch_noover_nofPixels
                #print 'Patch_noover_nofLines=',Patch_noover_nofLines
                #print 'slc_cut=',slc_cut
                #print 'slc_cut=',np.angle(slc_cut[600:610,1400:1410])
                #print 'Single_SLC_Data=',Single_SLC_Data[600:610,1400:1410]
                #test_slc=Single_SLC_Data[int(Patch_noover_First_Line-1):int(Patch_noover_Last_Line),int(Patch_noover_First_Pixel-1):int(Patch_noover_Last_Pixel)]
                #print  'test_slc=',test_slc.shape
                #print '\n\n'
                Single_SLC_Data[int(Patch_noover_First_Line):int(Patch_noover_Last_Line+1),int(Patch_noover_First_Pixel):int(Patch_noover_Last_Pixel+1)]=slc_cut.copy()
            else:
                sys.exit(1)
        Single_SLC_Data           =Single_SLC_Data.astype(np.complex64)
        Single_SLC_Data_minus_sign=np.abs(Single_SLC_Data)*np.exp(-1.0J*np.angle(Single_SLC_Data))
        Output_file_path    =os.getcwd()+'/'+'mini_stack_compress'+'/'+file_cint.strip()
        if os.path.exists(Output_file_path):
            Output_file         =Output_file_path+'/'+'resample_ministack.raw'

            print 'Output_file=',Output_file
            
            fid = open(Output_file,'wb');
            if fid < 0:
                print "ERROR :cannot open file for writing"
                sys.exit(1)
            else: 
                for i_temp in np.arange(0,int(nofLines)):
                    fid.write(Single_SLC_Data[i_temp,:])
                fid.close() 
                print "File Output complete !!!!" 
            
            Output_minus_file         =Output_file_path+'/'+'resample_ministack_minus_sign.raw'

            print 'Output_minus_file=',Output_minus_file

            fid = open(Output_minus_file,'wb')
            if fid < 0:
                print "ERROR :cannot open file for writing"
                sys.exit(1)
            else:
                for i_temp in np.arange(0,int(nofLines)):
                    fid.write(Single_SLC_Data_minus_sign[i_temp,:])
                fid.close()
                print "File Output complete !!!!"


            if os.path.isfile(Output_file):
                Crop_Image_width=nofPixels
                print  'Crop_Image_width=',Crop_Image_width
                python_instruction='cpxfiddle -w '+str(Crop_Image_width)+' -e 0.2 -s 1.8 -q mixed -o sunraster  -c hsv -M 40/8 -f cr4 -l1 -p1 '+\
        ' -P'+str(Crop_Image_width)+' '+Output_file+' '+'  > '+Output_file_path+'/interferogram_mix.ras'
                cmd_profile_copy=Linux_cmd(python_instruction)
                print 'cmd_profile_copy=',cmd_profile_copy

                jpg_name=os.getcwd()+'/'+'mini_stack_compress'+'/'+file_cint.strip()+'_crop.jpg'
                python_instruction='convert '+Output_file_path+'/interferogram_mix.ras '+' '+jpg_name
                cmd_profile_copy=Linux_cmd(python_instruction)
                print 'cmd_profile_copy=',cmd_profile_copy
                
                python_instruction='cpxfiddle -w '+str(Crop_Image_width)+' -e 0.2 -s 1.8 -q mixed -o sunraster  -c hsv -M 40/8 -f cr4 -l1 -p1 '+\
        ' -P'+str(Crop_Image_width)+' '+Output_minus_file+' '+'  > '+Output_file_path+'/interferogram_mix.ras'
                cmd_profile_copy=Linux_cmd(python_instruction)
                print 'cmd_profile_copy=',cmd_profile_copy

                jpg_name=os.getcwd()+'/'+'mini_stack_compress'+'/'+file_cint.strip()+'_minus_crop.jpg'
                python_instruction='convert '+Output_file_path+'/interferogram_mix.ras '+' '+jpg_name
                cmd_profile_copy=Linux_cmd(python_instruction)
                print 'cmd_profile_copy=',cmd_profile_copy


     
        print '\n\n' 
        
    
if __name__== '__main__':
    main()
