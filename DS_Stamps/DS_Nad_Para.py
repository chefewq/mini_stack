#!/usr/bin/env python
import multiprocessing
import numpy as np
from scipy.optimize import minimize
import os,sys,time
import numpy as np
from scipy.stats import ks_2samp
import math
import numpy as np
from numpy import *
import numpy
#import matplotlib.pyplot as plt
import scipy
import gdal,gdalconst
from gdalconst import *
from struct import *
import datetime
#version 2018 04 09
#Version 2018 10 25
#Version 2018 11 01 Wu Wenhao QQ:460249274
#Version 2019 07 30 Wu Wenhao QQ:460249274
#Version 2019 08 03 repair the bug: CPU CORE and len(win_mask)>0
#Version 2019 08 23 output the Parameter of filter Windows and amplitude dispersion
#Version 2019 08 25 Fixed the bug
#Version 20210805 Fixed the bug
###############################################################################
#
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
    #print 'Band Type=',gdal.GetDataTypeName(cint_srd.DataType)

    if cint_srd.GetOverviewCount() > 0:
            print 'Band has ', cint_srd.GetOverviewCount(), ' overviews.'
    thisBurstData= cint_srd.ReadAsArray(int(Pixels_start-1),int(line_start-1),nofPixels1,nofLines1)
    return thisBurstData
###############################################################################
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
 
###############################################################################################
def homofilter_KS(patch_Name):
    SLC_file='calamp.in'
    resFilename='master.res'
    l0 = int(get_parameter('First_line (w.r.t. original_image)',resFilename,1))
    lN = int(get_parameter('Last_line (w.r.t. original_image)',resFilename,1))
    p0 = int(get_parameter('First_pixel (w.r.t. original_image)',resFilename,1))
    pN = int(get_parameter('Last_pixel (w.r.t. original_image)',resFilename,1))
    nofLines=lN-l0+1
    nofPixels=pN-p0+1
    dataFormat = 'cpxfloat32'
    SLC_file_count=len(open(SLC_file).readlines())

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
    KS_folder='KS_result'

    if not os.path.exists(KS_folder):
        print "Can not find  KS_result!!!!"
        print "Can not find  KS_result!!!!"
        print "Can not find  KS_result!!!!"
        print "Can not find  KS_result!!!!"
        print "Can not find  KS_result!!!!"
        print "Can not find  KS_result!!!!"
        sys.exit(1)
    #    os.makedirs(KS_folder)
    #     
    #Homogenous_param = KS_folder+'/'+'homo_wins.par'
    
    ##############################################################
    ##Parameter of filter Windows
    #############################################################

    win_nofLines   =13
    win_nofPixels  =19
    amplitude_dispersion_param_value=0.25
   ##############################################################
    ##Parameter of filter Windows
    
    #file_para= open(Homogenous_param,'w')
    #file_para.write('win_nofLines:	%d\n'%win_nofLines)
    #file_para.write('win_nofPixels:	%d\n'%win_nofPixels)
    #file_para.write('amplitude_dispersion_param_value:		%f\n'%amplitude_dispersion_param_value)
    #file_para.close()
    Homogenous_param = KS_folder+'/'+'homo_wins.par'
    if os.path.exists(Homogenous_param):

        win_nofLines = int(get_parameter('win_nofLines',Homogenous_param,1))
        win_nofPixels = int(get_parameter('win_nofPixels',Homogenous_param,1))
        amplitude_dispersion_param_value=float(get_parameter('amplitude_dispersion_param_value',Homogenous_param,1))
        print 'win_nofLines =',win_nofLines 
        print 'win_nofPixels=',win_nofPixels
        print 'amplitude_dispersion_param_value=',amplitude_dispersion_param_value
    else:
        print "I can not find the homo_wins.par"
        sys.exit(1)
    #############################################################
    start_time = datetime.datetime.now()    

    Half_win_nofLines   =(win_nofLines-1)/2
    Half_win_nofPixels  =(win_nofPixels-1)/2
    
    Cen_Sta_Lines=(win_nofLines-1)/2
    Cen_Sta_Pixels=(win_nofPixels-1)/2

    Cen_End_Lines =nofLines -(win_nofLines+1)/2
    Cen_End_Pixels=nofPixels-(win_nofPixels+1)/2
    Total_win_pix=win_nofLines*win_nofPixels

    #win_mask=np.zeros((Total_win_pix),dtype=np.int8)
    #win_mask='0'*int(Total_win_pix)
    
    #for patch_Name in open('patch.list'):
    print 'patch_Name.strip()=',patch_Name.strip()
    patch_Data_Dimension=patch_Name.strip()+'/'+'patch.in'
    file_patch=open(patch_Data_Dimension)
    Patch_First_Pixel= int(file_patch.readline().strip())-1
    Patch_Last_Pixel = int(file_patch.readline().strip())-1
    Patch_First_Line = int(file_patch.readline().strip())-1
    Patch_Last_Line  = int(file_patch.readline().strip())-1
    Patch_nofPixels  = Patch_Last_Pixel- Patch_First_Pixel+1
    Patch_nofLines   = Patch_Last_Line - Patch_First_Line +1

    Patch_Sta_Lines  =Patch_First_Line +(win_nofLines+1)/2
    Patch_Sta_Pixels =Patch_First_Pixel+(win_nofPixels+1)/2
    Patch_End_Lines  =Patch_Last_Line  -(win_nofLines+1)/2
    Patch_End_Pixels =Patch_Last_Pixel -(win_nofPixels+1)/2

    patch_Data_Dimension_noover=patch_Name.strip()+'/'+'patch_noover.in'
    file_patch=open(patch_Data_Dimension_noover);
    Patch_noover_First_Pixel= int(file_patch.readline().strip())-1
    Patch_noover_Last_Pixel = int(file_patch.readline().strip())-1
    Patch_noover_First_Line = int(file_patch.readline().strip())-1
    Patch_noover_Last_Line  = int(file_patch.readline().strip())-1
    
    Patch_noover_Sta_Lines  =Patch_noover_First_Line 
    Patch_noover_Sta_Pixels =Patch_noover_First_Pixel
    Patch_noover_End_Lines  =Patch_noover_Last_Line 
    Patch_noover_End_Pixels =Patch_noover_Last_Pixel

    Limited_square_First_Line  =max(Patch_Sta_Lines,Patch_noover_First_Line)
    Limited_square_Last_Line   =min(Patch_End_Lines,Patch_noover_Last_Line)
    Limited_square_First_Pixel =max(Patch_Sta_Pixels,Patch_noover_First_Pixel)
    Limited_square_Last_Pixel  =min(Patch_End_Pixels,Patch_noover_Last_Pixel)

    file_count=0

    SLC_Data_Stack=scipy.zeros((int(SLC_file_count),int(Patch_nofLines),int(Patch_nofPixels)),dtype=np.float32)

    for dataFilename in open(SLC_file):
        SLCFilename=dataFilename.split('/')[-1].strip()
        SLC_Path,SLC_Name   =os.path.split(dataFilename)
        SLC_path_file       =SLC_Path+'/'+SLC_Name.strip()
        #print SLC_path_file
        Path_MFF_HDR   =SLC_Path+'/'+SLCFilename.split('.')[0]+'.hdr'
        Link_DATA      =SLC_Path+'/'+SLCFilename.split('.')[0]+'.x00'

        print 'Path_MFF_HDR=',Path_MFF_HDR 
        SLC_Data_Stack[file_count,:,:]=np.abs(freadbk(Path_MFF_HDR,Patch_First_Line+1,Patch_First_Pixel+1,Patch_nofLines,Patch_nofPixels))
        file_count=file_count+1
        print 'file_count=',file_count

    Total_line=0L
    psc_da=patch_Name.strip()+'/'+'pscands.1.da'
    psc_ij=patch_Name.strip()+'/'+'pscands.1.ij'
    

    KS_Path_file=KS_folder+'/'+patch_Name.strip()+'_data.dat'
    
    
    
    win_mask='' 
    print 'len(win_mask)=',len(win_mask)   
    Output_ks=open(KS_Path_file,'wb')
    with open(psc_da,'r') as fp1, open(psc_ij,'r') as fp2:
        for line_da in fp1:
            line_ij=fp2.readline() 
            string_line_ij=line_ij.strip('\n').split()  
            #print (float(line_da.strip('\n')),int(string_line_ij[0]),int(string_line_ij[1]),int(string_line_ij[2]))
            temp_da        =float(line_da.strip('\n'))
            temp_line_ind  =int(string_line_ij[1])
            temp_pixel_ind =int(string_line_ij[2])
            if (temp_da > amplitude_dispersion_param_value) and (temp_line_ind > Limited_square_First_Line) and (temp_line_ind < Limited_square_Last_Line) and (temp_pixel_ind > Limited_square_First_Pixel) and (temp_pixel_ind < Limited_square_Last_Pixel):                   
                    
                
                Relative_line_ind =int(temp_line_ind-Patch_First_Line)
                Relative_pixel_ind=int(temp_pixel_ind-Patch_First_Pixel)
                #for temp_Line in xrange(-Half_win_nofLines,Half_win_nofLines+1) :
                    #for temp_Pixel in xrange(-Half_win_nofPixels,Half_win_nofPixels+1):                  
                        #D,P = ks_2samp(SLC_Data_Stack[Relative_line_ind,Relative_pixel_ind,:],SLC_Data_Stack[Relative_line_ind+temp_Line,Relative_pixel_ind+temp_Pixel,:])
                Array_position=[(ks_2samp(SLC_Data_Stack[:,Relative_line_ind,Relative_pixel_ind], \
                                          SLC_Data_Stack[:,Relative_line_ind+temp_Line,Relative_pixel_ind+temp_Pixel])) \
                                          for temp_Line in xrange(-Half_win_nofLines,Half_win_nofLines+1) \
                                          for temp_Pixel in xrange(-Half_win_nofPixels,Half_win_nofPixels+1)]
                #for i_temp in 
                win_mask_ind=[ 1 if P_value > 0.05 else 0 for Statistic_KS,P_value in Array_position  ]
                #print 'win_mask_ind=',win_mask_ind
                win_mask_string = ''.join('%s' % id for id in win_mask_ind) 
                #win_mask_string=[]
                #print 'win_mask_string=',win_mask_string
                #print '\n'
                win_mask=win_mask+win_mask_string
     
                #print 'win_mask=',win_mask
                #print type(win_mask_ind)
                #print win_mask_ind
                       # if P > 0.05 :
                        #    win_mask=win_mask+'1'
                                #Output_ks.write('1')
                        #else:
                         #   win_mask=win_mask+'0'
                                #win_mask[temp_mask]=0
                                
                                #temp_mask=temp_mask+1                                         
                if int(string_line_ij[0])%25000==0:
                    if len(win_mask) > 0: 
                        Output_ks.write(win_mask)                   
                        win_mask=''
                    print 'processing the patch row:',patch_Name.strip()
                    print 'processing the data row:',string_line_ij[0]
                    end_time = datetime.datetime.now()
                    print "time-consuming=",end_time-start_time
    Output_ks.write(win_mask)    
    Output_ks.close()
    end_time = datetime.datetime.now()
    print "time-consuming=",end_time-start_time
    return                   




################################################################################################
def main():
     
    cores=multiprocessing.cpu_count()
    print 'Total cores is ', cores
    Patch_count=len(open('patch.list','r').readlines())
    multithreading_cpu=min(cores,Patch_count)
    pool=multiprocessing.Pool(processes=(multithreading_cpu))
  
    SLC_file='calamp.in'
    resFilename='master.res'
    l0 = int(get_parameter('First_line (w.r.t. original_image)',resFilename,1))
    lN = int(get_parameter('Last_line (w.r.t. original_image)',resFilename,1))
    p0 = int(get_parameter('First_pixel (w.r.t. original_image)',resFilename,1))
    pN = int(get_parameter('Last_pixel (w.r.t. original_image)',resFilename,1))
    nofLines=lN-l0+1
    nofPixels=pN-p0+1
    dataFormat = 'cpxfloat32'
    SLC_file_count=len(open(SLC_file).readlines())
    ##############################################################
    ##Parameter of filter Windows
    #############################################################
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
    KS_folder='KS_result'

    if not os.path.exists(KS_folder):
        os.makedirs(KS_folder)
        print "mkdir KS_result!!!!" 
    Homogenous_param = KS_folder+'/'+'homo_wins.par'
    
    ##############################################################
    ##Parameter of filter Windows
    #############################################################

    win_nofLines   =19
    win_nofPixels  =45
    amplitude_dispersion_param_value=0.1250001
    sum_filter_wins       =19
    ##############################################################
    ##Parameter of filter Windows
    
    file_para= open(Homogenous_param,'w')
    file_para.write('win_nofLines:	%d\n'%win_nofLines)
    file_para.write('win_nofPixels:	%d\n'%win_nofPixels)
    file_para.write('amplitude_dispersion_param_value:		%f\n'%amplitude_dispersion_param_value)
    file_para.write('sum_filter_wins:	%d\n'%sum_filter_wins)
    file_para.close()
    #############################################################
   ##############################################################
    ##Parameter of filter Windows
    #############################################################

    for dataFilename in open(SLC_file):
        SLCFilename=dataFilename.split('/')[-1].strip()
        SLC_Path,SLC_Name   =os.path.split(dataFilename)
        SLC_path_file       =SLC_Path+'/'+SLC_Name.strip()
        #print SLC_path_file
        Path_MFF_HDR   =SLC_Path+'/'+SLCFilename.split('.')[0]+'.hdr'
        Link_DATA      =SLC_Path+'/'+SLCFilename.split('.')[0]+'.x00'
        print Link_DATA
        if os.path.exists(Path_MFF_HDR)==False:
       
            outStream      = open(Path_MFF_HDR,'w')
            outStream.write('IMAGE_FILE_FORMAT = MFF\n')
            outStream.write('FILE_TYPE = IMAGE\n')
            outStream.write('IMAGE_LINES = %d\n' % int(nofLines))
            outStream.write('LINE_SAMPLES = %d\n'% int(nofPixels))
            outStream.write('BYTE_ORDER = LSB\n')
            outStream.write('END\n')
            outStream.close()
        
        if (os.path.isfile(Link_DATA)):
            os.remove(Link_DATA)        
        os.symlink(SLC_path_file ,Link_DATA)

    patch_Name=[]
    #f=open('patch.list')
    #lines=f.readlines()
    #print lines
    for lines in open('patch.list'):
    #for line_num in range(Patch_count-1,-1,-1):
        patch_string=lines.strip()
        if len(patch_string)>0:
            #print patch_string
            patch_Name.append(patch_string)
        #print lines[line_num]
    #    patch_Name=lines[line_num].strip()
    print 'patch_Name=',patch_Name
    pool.map(homofilter_KS,patch_Name)
    #patch_Name='PATCH_4'
    #homofilter_KS(patch_Name)
        
    
if __name__== '__main__':
    main()
