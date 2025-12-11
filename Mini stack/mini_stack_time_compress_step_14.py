#!/usr/bin/env python
import numpy as np
from scipy.optimize import minimize
#from scipy.optimize import Bounds
import scipy
print scipy.version.version
import multiprocessing
from scipy.optimize import minimize
import os,sys,time
from scipy.stats import ks_2samp
import math
import numpy as np
from numpy import *
from numpy import dot
from numpy import cos
import numpy
#import matplotlib.pyplot as plt
import gdal,gdalconst
from gdalconst import *
from struct import *
import re
import scipy.linalg
import scipy.io as sio
import linecache
import struct,shutil
import datetime
import copy
import gc
import fcntl
import threading
#Version 20181031
#Version 20190831
#Version 20190902 # correct the Bug:100 101 for three digital
#Version 20210207 # add the operation of chunksize
#####################################################################################
R=threading.Lock()
#global Observation_data_array
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
    #R.acquire()
    driver=gdal.GetDriverByName('MFF')
    driver.Register()
    gdal.AllRegister()
    thisBurstData_file=gdal.Open(path_file,GA_ReadOnly)
    if thisBurstData_file is None:
        #print 'Could not open'+Path_MFF_HDR
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
    #R.release()
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
###############################################################################
def SAR_write(output_file,Outdata,outformat):
    cols = Outdata.shape[1]
    rows = Outdata.shape[0]
    #print "cols = ",cols
    #print "rows = ", rows
    fid = open(output_file,'wb')
    if fid < 0:
        print "ERROR :cannot open file for writing\n"
        sys.exit(1)
    else: 
        if outformat=='complex_real4':
            Outdata=Outdata.astype(complex64)
            for i_temp in range(0,rows):
                fid.write(Outdata[i_temp,:])
            fid.close()
        elif outformat=='complex_short':
            data_short=zeros((2*cols,1),dtype='int16')
            for i_temp in arange(0,rows):        
                data_short[0::2]= (Outdata[i_temp,:].real.astype(np.int16)).reshape(-1,1,order='F').copy()
                data_short[1::2]= (Outdata[i_temp,:].imag.astype(np.int16)).reshape(-1,1,order='F').copy() 
                fid.write(data_short)
            fid.close() 
        elif outformat=='real4':
            Outdata=Outdata.real.astype(float32)
            for i_temp in range(0,rows):
                fid.write(Outdata[i_temp,:])
            fid.close()
        else:
             print "Sorry! not yet completed !"
             usage()
             sys.exit(1)
#****************************************************************************************
def SAR_data_Deci(SLC_path_file,data_format,nofLines,nofPixels,Patch_First_Line,Patch_First_Pixel,Patch_nofLines,Patch_nofPixels,sort_win_nofLines,sort_win_nofPixels=10):
    R.acquire()
    print SLC_path_file
    [fname,fename]=os.path.splitext(SLC_path_file)
    Path_MFF_HDR   =fname+'.hdr'
    if data_format=='complex_real4':
        Link_DATA      =fname+'.x00'
    elif data_format=='complex_short':
        Link_DATA      =fname+'.j00'
    elif data_format=='real4':
        Link_DATA      =fname+'.r00'
    else:
        print "Sorry! not yet completed !"
        usage()
        sys.exit(1)

    print 'Path_MFF_HDR=',Path_MFF_HDR 
    #if os.path.exists(Path_MFF_HDR)==True:
    #   os.remove(Path_MFF_HDR)
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
            #os.remove(Link_DATA)
        print  'Link_DATA=' ,Link_DATA
    else:        
        os.symlink(SLC_path_file,Link_DATA)
    Cint_Data=freadbk(Path_MFF_HDR,1,1,int(nofLines),int(nofPixels)) 
    Cint_Data=Cint_Data[::sort_win_nofLines,::sort_win_nofPixels].copy()
           
    #Cint_Data=freadbk(Path_MFF_HDR,Patch_First_Line+1,Patch_First_Pixel+1,Patch_nofLines,Patch_nofPixels)
    Deci_Patch_First_Line=int(Patch_First_Line/sort_win_nofLines)
    Deci_Patch_Last_Line=int(Patch_First_Line/sort_win_nofLines)+int((Patch_nofLines+sort_win_nofLines-1)/sort_win_nofLines)

    Deci_Patch_First_Pixel=int(Patch_First_Pixel/sort_win_nofPixels)
    Deci_Patch_Last_Pixel =int(Patch_First_Pixel/sort_win_nofPixels)+int((Patch_nofPixels+sort_win_nofPixels-1)/sort_win_nofPixels)

    Cint_Data =Cint_Data[Deci_Patch_First_Line:Deci_Patch_Last_Line,Deci_Patch_First_Pixel:Deci_Patch_Last_Pixel].copy()
                        
    R.release()
    return Cint_Data

#############################################################################################
#****************************************************************************************
def SAR_data_multilook(SLC_path_file,data_format,nofLines,nofPixels,Patch_First_Line,Patch_First_Pixel,Patch_nofLines,Patch_nofPixels,sort_win_nofLines,sort_win_nofPixels=10):
    R.acquire()
    print SLC_path_file
    [fname,fename]=os.path.splitext(SLC_path_file)
    Path_MFF_HDR   =fname+'.hdr'
    if data_format=='complex_real4':
        Link_DATA      =fname+'.x00'
    elif data_format=='complex_short':
        Link_DATA      =fname+'.j00'
    elif data_format=='real4':
        Link_DATA      =fname+'.r00'
    else:
        print "Sorry! not yet completed !"
        usage()
        sys.exit(1)

    print 'Path_MFF_HDR=',Path_MFF_HDR 
    #if os.path.exists(Path_MFF_HDR)==True:
    #   os.remove(Path_MFF_HDR)
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
            #os.remove(Link_DATA)
        print  'Link_DATA=' ,Link_DATA
    else:        
        os.symlink(SLC_path_file,Link_DATA)
    Cint_Data=freadbk(Path_MFF_HDR,1,1,int(nofLines),int(nofPixels)) 
    Cint_Data_Deci=Cint_Data[::sort_win_nofLines,::sort_win_nofPixels].copy()
    Cint_Data_temp=np.zeros((Cint_Data_Deci.shape[0],Cint_Data_Deci.shape[1]),dtype=np.complex64)
    Cint_Data_temp_edge=np.zeros((Cint_Data_Deci.shape[0],Cint_Data_Deci.shape[1]),dtype=np.complex64) 
    for i_temp in np.arange(0,sort_win_nofLines):
        for j_temp in np.arange(0,sort_win_nofPixels):
            a_temp=Cint_Data[i_temp::sort_win_nofLines,j_temp::sort_win_nofPixels].copy()
            Cint_Data_temp_edge[0:a_temp.shape[0],0:a_temp.shape[1]]=a_temp 
            Cint_Data_temp=Cint_Data_temp+Cint_Data_temp_edge
    Cint_Data_temp=Cint_Data_temp/sort_win_nofLines/sort_win_nofPixels
                 
    #Cint_Data=freadbk(Path_MFF_HDR,Patch_First_Line+1,Patch_First_Pixel+1,Patch_nofLines,Patch_nofPixels)
    Deci_Patch_First_Line=int(Patch_First_Line/sort_win_nofLines)
    Deci_Patch_Last_Line=int(Patch_First_Line/sort_win_nofLines)+int((Patch_nofLines+sort_win_nofLines-1)/sort_win_nofLines)

    Deci_Patch_First_Pixel=int(Patch_First_Pixel/sort_win_nofPixels)
    Deci_Patch_Last_Pixel =int(Patch_First_Pixel/sort_win_nofPixels)+int((Patch_nofPixels+sort_win_nofPixels-1)/sort_win_nofPixels)

    Cint_Data =Cint_Data_temp[Deci_Patch_First_Line:Deci_Patch_Last_Line,Deci_Patch_First_Pixel:Deci_Patch_Last_Pixel].copy()
                        
    R.release()
    return Cint_Data

#############################################################################################
#****************************************************************************************
def SAR_realdata_multilook(SLC_path_file,data_format,nofLines,nofPixels,Patch_First_Line,Patch_First_Pixel,Patch_nofLines,Patch_nofPixels,sort_win_nofLines,sort_win_nofPixels=10):
    R.acquire()
    print SLC_path_file
    [fname,fename]=os.path.splitext(SLC_path_file)
    Path_MFF_HDR   =fname+'.hdr'
    if data_format=='complex_real4':
        Link_DATA      =fname+'.x00'
    elif data_format=='complex_short':
        Link_DATA      =fname+'.j00'
    elif data_format=='real4':
        Link_DATA      =fname+'.r00'
    else:
        print "Sorry! not yet completed !"
        usage()
        sys.exit(1)

    print 'Path_MFF_HDR=',Path_MFF_HDR 
    #if os.path.exists(Path_MFF_HDR)==True:
    #   os.remove(Path_MFF_HDR)
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
            #os.remove(Link_DATA)
        print  'Link_DATA=' ,Link_DATA
    else:        
        os.symlink(SLC_path_file,Link_DATA)
    Cint_Data=np.abs(freadbk(Path_MFF_HDR,1,1,int(nofLines),int(nofPixels))) 
    Cint_Data_Deci=Cint_Data[::sort_win_nofLines,::sort_win_nofPixels].copy()
    Cint_Data_temp=np.zeros((Cint_Data_Deci.shape[0],Cint_Data_Deci.shape[1]),dtype=np.float32)
    Cint_Data_temp_edge=np.zeros((Cint_Data_Deci.shape[0],Cint_Data_Deci.shape[1]),dtype=np.float32) 
    for i_temp in np.arange(0,sort_win_nofLines):
        for j_temp in np.arange(0,sort_win_nofPixels):
            a_temp=Cint_Data[i_temp::sort_win_nofLines,j_temp::sort_win_nofPixels].copy()
            Cint_Data_temp_edge[0:a_temp.shape[0],0:a_temp.shape[1]]=a_temp
            Cint_Data_temp=Cint_Data_temp+Cint_Data_temp_edge
    Cint_Data_temp=Cint_Data_temp/sort_win_nofLines/sort_win_nofPixels
                 
    #Cint_Data=freadbk(Path_MFF_HDR,Patch_First_Line+1,Patch_First_Pixel+1,Patch_nofLines,Patch_nofPixels)
    Deci_Patch_First_Line=int(Patch_First_Line/sort_win_nofLines)
    Deci_Patch_Last_Line=int(Patch_First_Line/sort_win_nofLines)+int((Patch_nofLines+sort_win_nofLines-1)/sort_win_nofLines)

    Deci_Patch_First_Pixel=int(Patch_First_Pixel/sort_win_nofPixels)
    Deci_Patch_Last_Pixel =int(Patch_First_Pixel/sort_win_nofPixels)+int((Patch_nofPixels+sort_win_nofPixels-1)/sort_win_nofPixels)

    Cint_Data =Cint_Data_temp[Deci_Patch_First_Line:Deci_Patch_Last_Line,Deci_Patch_First_Pixel:Deci_Patch_Last_Pixel].copy()
                        
    R.release()
    return Cint_Data

#############################################################################################
#****************************************************************************************
def SAR_data_read(SLC_path_file,data_format,nofLines,nofPixels,Patch_First_Line,Patch_First_Pixel,Patch_nofLines,Patch_nofPixels):
    R.acquire()
    print SLC_path_file
    [fname,fename]=os.path.splitext(SLC_path_file)
    Path_MFF_HDR   =fname+'.hdr'
    if data_format=='complex_real4':
        Link_DATA      =fname+'.x00'
    elif data_format=='complex_short':
        Link_DATA      =fname+'.j00'
    elif data_format=='real4':
        Link_DATA      =fname+'.r00'
    else:
        print "Sorry! not yet completed !"
        usage()
        sys.exit(1)

    print 'Path_MFF_HDR=',Path_MFF_HDR 
    #if os.path.exists(Path_MFF_HDR)==True:
    #   os.remove(Path_MFF_HDR)
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
            #os.remove(Link_DATA)
        print  'Link_DATA=' ,Link_DATA
    else:        
        os.symlink(SLC_path_file,Link_DATA)
            
    Cint_Data=freadbk(Path_MFF_HDR,Patch_First_Line+1,Patch_First_Pixel+1,Patch_nofLines,Patch_nofPixels)
    R.release()
    return Cint_Data

#############################################################################################
def Optimization_wrap_phase(para_insar,Observation_data_array,Observation_coherence_array):
    #global A_matrix
    #print  'para_insar = ',para_insar.shape 
    #print  'Observation_data_array = ',Observation_data_array.shape
    #print  'para_insar = ',para_insar.shape
    #temp_master_phase = np.hstack((temp_master_phase,para_insar))
    #temp_master_phase = temp_master_phase.reshape(-1,1)
    #print  'temp_master_phase = ',para_insar.shape
    #print 'temp_master_phase = ',temp_master_phase.shape
    #print 'A_matrix = ',A_matrix.shape
    #para_InSAR_phase  = np.dot(A_matrix,para_insar)-Observation_data_array
    #.reshape(-1,1)
    #Ind=np.where(np.cos(para_InSAR_phase)<0)
    #print 'para_InSAR_phase = ',np.cos(para_InSAR_phase)
    #print 'para_InSAR_phase = ',np.where(np.cos(para_InSAR_phase)<0)
    #print 'para_InSAR_phase = ',np.sum(np.cos(para_InSAR_phase))
    #print 'para_InSAR_phase=',para_InSAR_phase[Ind]
    #print 'para_InSAR_phase=',np.cos(para_InSAR_phase[Ind])
    #print 'Observation_data_array=',Observation_data_array.shape

    #para_InSAR_phase=Observation_data_array+para_InSAR_phase
    #print np.sum(Observation_coherence_array)
    #coherence_array=-1*Observation_coherence_array.copy().reshape(1,-1)
    #coherence_test=-1*np.dot(Observation_coherence_array.ravel(),np.cos(para_InSAR_phase))
    #print 'coherence_test=', coherence_test
    #print 'coherence_test=', coherence_test.shape
    return -1*np.dot(Observation_coherence_array.ravel(),np.cos(np.dot(A_matrix,para_insar)-Observation_data_array))

#############################################################################################
def Optimization_wrap_phase_der(para_insar,Observation_data_array,Observation_coherence_array):
    #print  'Observation_data_array=',Observation_data_array
    #print 'A_matrix=',A_matrix.shape
    #print 'para_insar=',para_insar.shape
    #print '\n\n\n'
    global A_matrix
    der=np.zeros_like(para_insar)
    #temp_master_phase= np.array([0])
    #para_insar=np.hstack((temp_master_phase,para_insar))
    
    for temp_i in xrange(0,len(der)):
        Ind=np.where(A_matrix[:,temp_i] !=0)
        A_matrix_der=A_matrix[Ind,:] 
        A_matrix_array=A_matrix_der[0,:,temp_i] 
        para_InSAR_phase=np.dot(A_matrix_der[0,:,:],para_insar)-Observation_data_array[Ind]
        A_matrix_array=A_matrix_array*Observation_coherence_array[Ind]
        der[temp_i]=np.dot(A_matrix_array.reshape(1,-1) ,np.sin(para_InSAR_phase))
        
    #sum_A_Matrix=np.dot(A_matrix,axis=0)

    #Internel_Varible=np.dot(sum_A_Matrix,para_insar)
    #Internel_Result=np.sin(para_InSAR_phase)
    
    #para_insar=np.hstack((temp_master_phase,para_insar))
    #para_InSAR_phase=np.dot(A_matrix,para_insar)
    #print 'para_InSAR_phase=',para_InSAR_phase.shape
    #para_InSAR_phase=Observation_data_array+para_InSAR_phase
    return der
   
#################################################################################################
#def phase_linking_paral() 
#################################################################################################
def PCA_insar(input_fun_cmd):
    work_dir=input_fun_cmd[0]
    target_run_path=input_fun_cmd[0]
    Master_time=work_dir.split('/')[-1]
    Master_time=Master_time[-8:]
    print 'Master_time=', Master_time
    print 'Master_time=', type(Master_time)

    SBAS_work_dir=work_dir+'/SMALL_BASELINES'
    
    start_time = datetime.datetime.now()
    if os.path.exists(SBAS_work_dir):    
        print 'SBAS_work_dir=',SBAS_work_dir    
    else:
        print "This path is wrong !!!"
        sys.exit(1)
   
    KS_folder=work_dir+'/KS_result'

    ############################################################
    cint_day_file=input_fun_cmd[1] 
    SLC_list_Index=input_fun_cmd[1] 
    Resample_day_file=copy.deepcopy(input_fun_cmd[1])
   
    cint_minrefdem_file=SLC_list_Index=input_fun_cmd[2]
    Relative_MASTER_DATE_INDEX=int(input_fun_cmd[3])
    print 'Relative_MASTER_DATE_INDEX=',Relative_MASTER_DATE_INDEX
    #print 'Total_Numer_Mini_Stack=',Total_Numer_Mini_Stack
    #f=open(work_dir+'/'+'patch.list')
    #for patch_file in f.readlines():
    patch_file=input_fun_cmd[4]

    if os.path.exists(KS_folder):
        patch_Name=patch_file.strip()
    #for Mini_Stack_Index in range(0,Total_Numer_Mini_Stack):
     
        Num_Int =len(cint_minrefdem_file) 
        print 'Num_Int=',Num_Int 
        Num_SLC= len(cint_day_file)
        print 'Num_SLC=',Num_SLC
        if Num_SLC < 3:
            print "The number of mini stack image is not enough!!!!!!!!!"
            #break
            sys.exit(1)
        
   
        
        ##############################################################
        ##Parameter of filter Windows
        #############################################################

        win_nofLines   =13
        win_nofPixels  =19
        amplitude_dispersion_param_value=0.25
        sum_filter_wins=19
           ###
        #############
        #window
        #Sort_Window_para=
        sort_win_nofLines=2
        sort_win_nofPixels=10
        ################################################
        #MASTER_DATE_INDEX=SLC_list_Index.index(Master_time)
        #print 'MASTER_DATE_INDEX=',MASTER_DATE_INDEX
        ##############################################################
        ##Parameter of filter Windows
        #############################################################
        KS_folder=target_run_path+'/KS_result'
        Homogenous_param = KS_folder+'/'+'homo_wins.par'
        print 'os.path.exists(Homogenous_param)=',os.path.exists(Homogenous_param)
        if os.path.exists(Homogenous_param):
            win_nofLines = int(get_parameter('win_nofLines',Homogenous_param,1))
            win_nofPixels = int(get_parameter('win_nofPixels',Homogenous_param,1))
            amplitude_dispersion_param_value=float(get_parameter('amplitude_dispersion_param_value',Homogenous_param,1))
            sum_filter_wins = int(get_parameter('sum_filter_wins',Homogenous_param,1))
        else:
            print "I can not find the homo_wins.par"
            #sys.exit(1)
        ######################################################
        Half_win_nofLines   =(win_nofLines-1)/2
        Half_win_nofPixels  =(win_nofPixels-1)/2
    
        Cen_Sta_Lines=(win_nofLines-1)/2
        Cen_Sta_Pixels=(win_nofPixels-1)/2
        Total_win_pix=win_nofLines*win_nofPixels
        win_mask=np.zeros((Total_win_pix),dtype=np.int8)
        Number_win_mask=200000
        win_mask_total=np.zeros((Total_win_pix*Number_win_mask),dtype=np.int8)
        
        
        #i_win=0

        KS_Estimation = target_run_path+'/KS_result'
        KS_folder    =target_run_path+'/KS_result'

        print 'patch_Name.strip()=',patch_Name.strip()
        #print "patch_Name.split('_')=",patch_Name.split('_')
        patch_Data_Dimension=target_run_path+'/'+patch_Name.strip()+'/'+'patch.in'
        print 'patch_Data_Dimension=',patch_Data_Dimension
        
        if os.path.exists(patch_Data_Dimension)== True and\
             len(cint_minrefdem_file) and \
             os.path.isdir(KS_folder) >0 and \
             os.path.exists(Homogenous_param):

            file_patch=open(patch_Data_Dimension)
            fcntl.flock(file_patch.fileno(), fcntl.LOCK_EX)
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
            fcntl.flock(file_patch.fileno(),fcntl.LOCK_UN)
            file_patch.close()
            patch_Data_Dimension_noover=target_run_path+'/'+patch_Name.strip()+'/'+'patch_noover.in'
            

            file_patch=open(patch_Data_Dimension_noover)
            fcntl.flock(file_patch.fileno(), fcntl.LOCK_EX)
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
            fcntl.flock(file_patch.fileno(),fcntl.LOCK_UN)
            file_patch.close()

            file_count=0   
            Int_Data_Stack=np.zeros((int(Patch_nofPixels*Patch_nofLines),int(Num_SLC),int(Num_SLC)),dtype=np.complex64)
            #Decompostion_Int_Data_Stack=np.zeros((int(Num_SLC),int(Num_SLC),sort_win_nofLines*sort_win_nofPixels),dtype=np.complex64)
            #coh_mx = np.zeros((int(Patch_nofPixels*Patch_nofLines),int(Num_SLC),int(Num_SLC)),dtype=float32)
            for i_j_temp in np.arange(0,int(Num_SLC)):
                Int_Data_Stack[:,i_j_temp,i_j_temp]=1
                #coh_mx[:,i_j_temp,i_j_temp] = 1
            for dataFilename in cint_minrefdem_file:
                SLCFilename=SBAS_work_dir+'/'+dataFilename.strip()
                #####################################################################################
                Imag_Master=dataFilename.strip()[0:8]
                Imag_Slave =dataFilename.strip()[9:17]
                print 'Imag_Master=' ,Imag_Master
                print 'Imag_Slave=' ,Imag_Slave
                #print 'SLC_list_Index=',SLC_list_Index               
                Matrix_subscript_Master=Resample_day_file.index(Imag_Master)
                Matrix_subscript_Slave =Resample_day_file.index(Imag_Slave) 
                ####################################################################################
                resFilename=SLCFilename+'/'+'master.res'
                l0 = int(get_parameter('First_line (w.r.t. original_image)',resFilename,1))
                lN = int(get_parameter('Last_line (w.r.t. original_image)',resFilename,1))
                p0 = int(get_parameter('First_pixel (w.r.t. original_image)',resFilename,1))
                pN = int(get_parameter('Last_pixel (w.r.t. original_image)',resFilename,1))
                Total_nofLines=lN-l0+1
                Total_nofPixels=pN-p0+1

                SLC_path_file       =SLCFilename+'/'+'cint_filt_relative_master'+'_'+patch_Name.strip()+'.raw'
                if os.path.isfile(SLC_path_file):
                    print 'SLC_path_file =',SLC_path_file
                else:
                 
                    SLC_path_file       =SLCFilename+'/'+'cint.minrefpha'+'_'+patch_Name.strip()+'.raw'
                data_format='complex_real4'

                Cint_Dem=np.angle(SAR_data_read(SLC_path_file,data_format,\
                         Patch_nofLines,Patch_nofPixels,0,0,Patch_nofLines,Patch_nofPixels))
                       
                if (Matrix_subscript_Master > Matrix_subscript_Slave): 
                    Int_Data_Stack[:,Matrix_subscript_Master,Matrix_subscript_Slave]=np.exp(-1j*Cint_Dem).reshape(1,-1)
                    Int_Data_Stack[:,Matrix_subscript_Slave,Matrix_subscript_Master]=np.exp( 1j*Cint_Dem).reshape(1,-1)                       
                else:
                    Int_Data_Stack[:,Matrix_subscript_Master,Matrix_subscript_Slave]=np.exp( 1j*Cint_Dem).reshape(1,-1)                        
                    Int_Data_Stack[:,Matrix_subscript_Slave,Matrix_subscript_Master]=np.exp(-1j*Cint_Dem).reshape(1,-1)

                file_count=file_count+1
                print 'file_count=',file_count
                print '\n\n\n'
            del Cint_Dem
            gc.collect()
        ####################################################################################
            file_count=0
        
            for dataFilename in cint_minrefdem_file:
                SLCFilename=SBAS_work_dir+'/'+dataFilename.strip()
                #####################################################################################
                Imag_Master=dataFilename.strip()[0:8]
                Imag_Slave =dataFilename.strip()[9:17]
                print 'Imag_Master=' ,Imag_Master
                print 'Imag_Slave=' ,Imag_Slave
                                
                Matrix_subscript_Master=Resample_day_file.index(Imag_Master)
                Matrix_subscript_Slave =Resample_day_file.index(Imag_Slave) 
                ####################################################################################
                SLC_path_file       =SLCFilename+'/'+'coherence.raw'
                print SLC_path_file
              
                data_format='real4'
                Coherence_Data=SAR_data_read(SLC_path_file,data_format,\
                               Total_nofLines,Total_nofPixels,Patch_First_Line,\
                                Patch_First_Pixel,Patch_nofLines,Patch_nofPixels)

                #Coherence_Data=(freadbk(Path_MFF_HDR,Patch_First_Line,Patch_First_Pixel,Patch_nofLines,Patch_nofPixels)).copy()
                #print Coherence_Data_Stack[file_count,91:93,108:110]
                #if (Matrix_subscript_Master > Matrix_subscript_Slave): 
                Int_Data_Stack[:,Matrix_subscript_Master,Matrix_subscript_Slave]=\
                Int_Data_Stack[:,Matrix_subscript_Master,Matrix_subscript_Slave]*Coherence_Data.reshape(1,-1)
                Int_Data_Stack[:,Matrix_subscript_Slave,Matrix_subscript_Master]=\
                Int_Data_Stack[:,Matrix_subscript_Slave,Matrix_subscript_Master]*Coherence_Data.reshape(1,-1)                       
               

                file_count=file_count+1
            
                print 'file_count=',file_count 
                print '\n\n\n'
            del Coherence_Data
            gc.collect()   
            ########################################################################################
            # Read_Master_Slave
            ########################################################################################
            #
            resFilename=work_dir+'/'+'master.res'
            l0 = int(get_parameter('First_line (w.r.t. original_image)',resFilename,1))
            lN = int(get_parameter('Last_line (w.r.t. original_image)',resFilename,1))
            p0 = int(get_parameter('First_pixel (w.r.t. original_image)',resFilename,1))
            pN = int(get_parameter('Last_pixel (w.r.t. original_image)',resFilename,1))
            Total_nofLines=lN-l0+1
            Total_nofPixels=pN-p0+1
            file_count=0
            print 'cint_day_file=',cint_day_file
            
            #Relative_MASTER_DATE_INDEX=int((len(cint_day_file)-1)/2-1)
            print 'Relative_MASTER_DATE_INDEX=',Relative_MASTER_DATE_INDEX

            Relative_master_mini_stack=cint_day_file[Relative_MASTER_DATE_INDEX]
            print 'Relative_master_mini_stack=',Relative_master_mini_stack
            #del cint_day_file[Relative_MASTER_DATE_INDEX]
            Single_SLC_Data_Stack=np.zeros((len(cint_day_file)-1,int(Patch_nofLines),\
                                  int(Patch_nofPixels)),dtype=np.complex64)
            print 'cint_day_file=',cint_day_file
            for file_cint  in cint_day_file:
                print 'file_cint=',file_cint
                if int(file_cint) > int(Relative_master_mini_stack.strip()): 
                    temp_cint_file_name= Relative_master_mini_stack.strip()+'_'+file_cint.strip()
                if int(file_cint) < int(Relative_master_mini_stack.strip()):
                    temp_cint_file_name= file_cint.strip()+'_'+Relative_master_mini_stack.strip()
                if int(file_cint) == int(Relative_master_mini_stack.strip()):
                    continue
                print 'Relative_master_mini_stack.strip()=',Relative_master_mini_stack.strip()
                print 'temp_cint_file_name=',temp_cint_file_name 
                SLC_path_file  =SBAS_work_dir+'/'+temp_cint_file_name+'/'+'cint.mini_track_time_compress.raw'
              
                #cint_data_init=freadbk(Path_MFF_HDR,Patch_First_Line,Patch_First_Pixel,Patch_nofLines,Patch_nofPixels)
                data_format='complex_real4'
                cint_data_init=SAR_data_read(SLC_path_file,data_format,\
                               Total_nofLines,Total_nofPixels,\
                               Patch_First_Line,Patch_First_Pixel,Patch_nofLines,Patch_nofPixels)

                if int(file_cint) > int(Relative_master_mini_stack.strip()):
                    Single_SLC_Data_Stack[file_count,:,:]=cint_data_init.copy()########################
                if int(file_cint) < int(Relative_master_mini_stack.strip()):
                    Single_SLC_Data_Stack[file_count,:,:]=np.exp(-1j*np.angle(cint_data_init))*np.abs(cint_data_init)   
                file_count=file_count+1
                print 'file_count=',file_count 
                print '\n\n'  
            Single_SLC_KS =Single_SLC_Data_Stack.copy()
             
    #######################################################################################
           
    ######################################################################################
            #Total_coherence=0    
            #print 'Num_Int=',Num_Int
            #Expected_coherence  = Num_Int/4

            #Vari_Num_PS_Cint       =len(cint_day_file)-1
            #para_insar             =np.zeros((Vari_Num_PS_Cint),dtype=float32)
            #Observation_data_array =np.zeros((Num_Int),dtype=float32)
            #Observation_coherence_array=np.zeros((Num_Int),dtype=float32)
            Coherence_Ind=np.append(xrange(Relative_MASTER_DATE_INDEX),xrange(Relative_MASTER_DATE_INDEX+1,Num_SLC))
            Coherence_Ind=Coherence_Ind.astype(np.int16)
            print 'Coherence_Ind=',Coherence_Ind

            PCA_Energy_thres=15/Num_SLC
            temp_Relative_line_ind =-1
            temp_Relative_pixel_ind=-1
            temp_Decompostion_Int_Data_Stack_index=0
            total_count_tensor=0 
            #PCA_Energy_thres=0.35
            #PCA_Energy_thres=2.0/7.0
            #Total_line=0L
            psc_da=work_dir+'/'+patch_Name.strip()+'/'+'pscands.1.da'
            psc_ij=work_dir+'/'+patch_Name.strip()+'/'+'pscands.1.ij'
            if os.path.exists(KS_folder)==False:
                print "KS_result is not exit"
            else:
                file_KS_Estimation=KS_folder+'/'+patch_Name.strip()+'_data'+'.dat'
                print file_KS_Estimation
 
                Line_index_ks=0
                f=open(file_KS_Estimation,'rb')
                fcntl.flock(f.fileno(), fcntl.LOCK_EX)
                with open(psc_da,'r') as fp1, open(psc_ij,'r') as fp2:
                    fcntl.flock(fp1.fileno(), fcntl.LOCK_EX)
                    fcntl.flock(fp2.fileno(), fcntl.LOCK_EX)
                    for line_da in fp1:
                        line_ij=fp2.readline() 
                        string_line_ij=line_ij.strip('\n').split()  
                        #print (float(line_da.strip('\n')),int(string_line_ij[0]),int(string_line_ij[1]),int(string_line_ij[2]))
                        temp_da        =float(line_da.strip('\n'))
                        temp_line_ind  =int(string_line_ij[1])
                        temp_pixel_ind =int(string_line_ij[2])
                        if (temp_da > amplitude_dispersion_param_value) and (temp_line_ind > Limited_square_First_Line) \
                                       and (temp_line_ind < Limited_square_Last_Line)  \
                                       and (temp_pixel_ind > Limited_square_First_Pixel) \
                                       and (temp_pixel_ind < Limited_square_Last_Pixel):
                
                            #print 'temp_line_ind=',temp_line_ind
                            #print 'temp_pixel_ind=',temp_pixel_ind
                            Relative_line_ind =int(temp_line_ind-Patch_First_Line)
                            Relative_pixel_ind=int(temp_pixel_ind-Patch_First_Pixel)
                            
                            if int(Line_index_ks)%Number_win_mask==0:
                                total_length_win=Total_win_pix*Number_win_mask 
                                #print 'ks_win_total_length=',total_length_win
                                temp=f.read(total_length_win)
                                #position_tell_file=f.tell()
                                #print 'position_tell_file=',position_tell_file
                                temp_length=len(temp)
                                #print 'temp_length=',len(temp) 
                                if temp_length == total_length_win:
                                    print "Processing now !!!"
                                else:
                                    print 'Please Notice!!! The prcessing will end soon !!!!'
                                    #sys.exit(1)   
                                unpack_format=str(temp_length)+'c'
                                #print 'unpack_format=',unpack_format                           
                                win_mask_total=int8(struct.unpack(unpack_format,temp))
                                #print 'Line_index_ks=',Line_index_ks
                                Line_index_ks=0
                                 
                                print 'processing the data row:',string_line_ij[0]
                                print '\n\n'
                            #temp=f.read(Total_win_pix) 
                            #print 'Line_index_ks=',Line_index_ks
                            Line_index_ks+=1                         
                            win_mask=win_mask_total[((Line_index_ks-1)*Total_win_pix):(Line_index_ks*Total_win_pix)]
                            #int8(struct.unpack(unpack_format,temp))
                            #print  'win_mask= ',win_mask
                            sum_filter=win_mask.sum()   
                            #print 'sum_filter=',sum_filter
                            #sum_filter_wins=20000
                            if sum_filter > sum_filter_wins:                            
                                #Total_coherence=np.sum(Observation_coherence_array)
                                
                                # Initival Value                           
                                coh_mx_inv = np.linalg.pinv(np.abs(Int_Data_Stack[Relative_line_ind*Patch_nofPixels+Relative_pixel_ind,:,:]))
                                #print 'coh_mx_inv =',coh_mx_inv 
                                cov_mx = Int_Data_Stack[Relative_line_ind*Patch_nofPixels+Relative_pixel_ind,:,:] # covariance matrix
                                #print  'cov_mx=',cov_mx 
                                cov_EMI = coh_mx_inv*cov_mx # Eq.13 of EMI
                                #print 'cov_EMI =',cov_EMI                             
                                eig_val, eig_vec = np.linalg.eig(cov_EMI) # eigen value and vector                          
                                #print 'eig_val=',eig_val                                                                        
                                                                                                                                            
                                if np.min(np.real(eig_val)) > 1.0: # EMI. To avoid small negative eigenvalue.
                                                               # The minimum eigenvalue should be larger than 1 (By Eq.12 & 14 of EMI).
                                    abs_eig_val      = np.abs(eig_val)
             
                                    red_eig_vec=np.angle(eig_vec[:,np.argmin(abs_eig_val)])
                                    #print 'red_eig_vec=',red_eig_vec
                                    Optimization_pca=red_eig_vec[Relative_MASTER_DATE_INDEX]-red_eig_vec
                                    #Optimization_pca=np.delete(Optimization_pca,Relative_MASTER_DATE_INDEX)
                              
                                    Single_SLC_KS[:,Relative_line_ind,Relative_pixel_ind] = np.exp(1j*Optimization_pca[Coherence_Ind])*np.abs(Single_SLC_Data_Stack[:,Relative_line_ind,Relative_pixel_ind])

                                
                                else:

                                    eig_val, eig_vec = np.linalg.eig(Int_Data_Stack[Relative_line_ind*Patch_nofPixels+Relative_pixel_ind,:,:])
                                    abs_eig_val      = np.abs(eig_val)
                                    max_eig_val      = np.max(abs_eig_val)                            
                                
                                    if (max_eig_val/np.sum(abs_eig_val)) > PCA_Energy_thres:
                                    
                                        max_eig_vec=eig_vec[:,np.argmax(abs_eig_val)]
                                        max_eig_vec=max_eig_vec*conj(max_eig_vec[Relative_MASTER_DATE_INDEX])
                                
                                        max_eig_vec=max_eig_vec/np.linalg.norm(max_eig_vec)                                    
                                    
                                   
                                    ###########mini_stack#################################
                                        red_eig_vec=np.angle(eig_vec[:,np.argmax(abs_eig_val)])
                                        #print 'red_eig_vec=',red_eig_vec
                                        Optimization_pca=red_eig_vec[Relative_MASTER_DATE_INDEX]-red_eig_vec
                                    
                                        Optimization_pca=np.delete(Optimization_pca,Relative_MASTER_DATE_INDEX)
                                        Single_SLC_KS[:,Relative_line_ind,Relative_pixel_ind] = np.exp(1j*Optimization_pca)*np.abs(Single_SLC_Data_Stack[:,Relative_line_ind,Relative_pixel_ind])
                                    
                                 
                                
                                if int(string_line_ij[0])%25000==0:
                                       
                                    print 'processing the patch row:',patch_Name.strip()
                                    print 'processing the data row:',string_line_ij[0]
                                    print 'Relative_master_mini_stack=',Relative_master_mini_stack
                                    end_time = datetime.datetime.now()
                                    print "time-consuming=",end_time-start_time
                                    print '\n'
                    fcntl.flock(fp1.fileno(),fcntl.LOCK_UN)
                    fcntl.flock(fp2.fileno(),fcntl.LOCK_UN)
                fcntl.flock(f.fileno(),fcntl.LOCK_UN)       
                f.close() #close the file
            
            del Int_Data_Stack
            gc.collect()
            del Single_SLC_Data_Stack
            gc.collect()
            file_count=0 
            ###################################################
            
            ###################################################
            print 'cint_day_file=',cint_day_file
           
            for file_cint  in cint_day_file:
                print 'file_cint=',file_cint
                if int(file_cint) > int(Relative_master_mini_stack.strip()): 
                    temp_cint_file_name= Relative_master_mini_stack.strip()+'_'+file_cint.strip()
                if int(file_cint) < int(Relative_master_mini_stack.strip()):
                    temp_cint_file_name= file_cint.strip()+'_'+Relative_master_mini_stack.strip()
                if int(file_cint) == int(Relative_master_mini_stack.strip()):                    
                    continue 
                
                
                Output_fold    =work_dir+'/'+'Relative_master'+'/'+file_cint
                if not os.path.exists(Output_fold):
                    os.makedirs(Output_fold)
                Output_file    =Output_fold+'/'+'cint.relative_master_'+patch_Name.strip()+'.raw'     
                
                print 'Output_file=',Output_file
                outformat= 'complex_real4'
                print 'outformat=',outformat
                SAR_write(Output_file,Single_SLC_KS[file_count,:,:],outformat)
                file_count=file_count+1
                
            del Single_SLC_KS
            gc.collect()
            ###################################################


               
        else:
            print "Can not find the file or file fold !!!!!!!!!!!!!!!!!!!!!!!!!"
            print "Can not find the file or file fold !!!!!!!!!!!!!!!!!!!!!!!!!"
            print "Can not find the file or file fold !!!!!!!!!!!!!!!!!!!!!!!!!"
            print "Can not find the file or file fold !!!!!!!!!!!!!!!!!!!!!!!!!"
            print "Can not find the file or file fold !!!!!!!!!!!!!!!!!!!!!!!!!" 
            #sys.exit(1)
    for x in list(locals().keys())[:]:
        del locals()[x]
    gc.collect()
################################################################################################
def main():
     
    cores=multiprocessing.cpu_count()
    #print 'Total cores is ', cores
    pool=multiprocessing.Pool(processes=(int(cores/2)))

    work_dir=os.getcwd();
    print 'work_dir=' ,work_dir
    ##################################################################
    #N1datadir="/home/wuwenhao/data_disk/kunming/cangyuan_sentinel/data_processing"
    #####provide convenience for users
    ##################################################################

    os.chdir(work_dir)

    filename_list=[]
    large_mini_stack=[]
    for filename in open('mini_stack_combine.txt'):
        print 'filename=',filename.split('\t')
        master_date=filename.strip().split('\t')[2].strip()
        filename_list.append(filename.strip()) 
        large_mini_stack.append('stamps_'+filename.split('\t')[5].strip())
    print 'filename_list=',filename_list
    #print 'large_mini_stack=',large_mini_stack

    large_mini_stack =list(set(large_mini_stack))
    large_mini_stack=sorted(large_mini_stack)
    print 'large_mini_stack=',large_mini_stack

    #******************************************************************************
    geo_initialization_index=0

    for ks_stack_num in large_mini_stack:
        large_mini_stack_fold_name = ks_stack_num
        if  os.path.exists(large_mini_stack_fold_name):
            #python_instruction='rm -fr '+os.path.abspath(large_mini_stack_fold_name)
            #os.system(python_instruction)
            print "This is the target fold"
        else:
            #os.makedirs(large_mini_stack_fold_name) 
            print "can not find the fold"
 
        Master_Date_fold=large_mini_stack_fold_name+'/'+'INSAR_'+master_date
        if not os.path.exists(Master_Date_fold):
            #os.makedirs(Master_Date_fold)
            print "can not find the fold"
        print 'ks_stack_num=',ks_stack_num
        geo_initialization_index=0
        Insar_master_check=0
        ks_stack_index=ks_stack_num.split('_')[-1].strip()
##################################################################
    
        mini_stack_index_list=[]
        for insar_file_name in filename_list:
        
            #print '\n\n\n'
            #print 'insar_file_name=',insar_file_name
            file_time    =insar_file_name.split('\t')[3].strip()
        
            master_date=insar_file_name.split('\t')[2].strip()        
            #print 'ks_stack_index=',ks_stack_index
            if (insar_file_name.split('\t')[1].strip()!='-1') and\
            int(insar_file_name.split('\t')[5].strip())==int(ks_stack_index):
                #print "int(insar_file_name.split('\t')[-1].strip())=",int(insar_file_name.split('\t')[-1].strip())
                #if (file_time.strip()!=master_date.strip()) and\
                #int(insar_file_name.split('\t')[-1].strip())==int(ks_stack_index):
                mini_stack_index_list.append(insar_file_name.split('\t')[1].strip())
        mini_stack_index_list =list(set(mini_stack_index_list))
        mini_stack_index_list=sorted(mini_stack_index_list)
        print 'mini_stack_index_list=',mini_stack_index_list 
        
        #phase_linking_paral()

        mini_stack_combine_list=[]
        input_cmd_list=[]

        master_mini_stack_stack=[]
        for mini_stack_index_num in mini_stack_index_list:
            input_cmd=[]
            mini_stack_file_name_list=[]
            Relative_master_mini_stack=[]
            for insar_file_name in filename_list:
              
                file_time    =insar_file_name.split('\t')[3].strip()
        
                master_date=insar_file_name.split('\t')[2].strip()
                

                #print ' mini_stack_index_num=', mini_stack_index_num
                if (insar_file_name.split('\t')[1].strip()!='-1') and\
                    int(insar_file_name.split('\t')[1].strip())==int(mini_stack_index_num)\
                    and int(insar_file_name.split('\t')[5].strip())==int(ks_stack_index):

                    mini_stack_file_name_list.append(file_time) 
                    Relative_master_mini_stack.append(insar_file_name.split('\t')[4].strip())           
            #print 'mini_stack_file_name_list=',mini_stack_file_name_list
            mini_stack_combine_list=[]
            for i_index in mini_stack_file_name_list:     

                for j_index in mini_stack_file_name_list:
                    if int(j_index) > int(i_index):
                        master_time=i_index
                        slave_time=j_index 
                        temp_combine_file_name=master_time.strip()+'_'+slave_time.strip()
                        mini_stack_combine_list.append(temp_combine_file_name)
            Relative_master_mini_stack=list(set(Relative_master_mini_stack))
            print '\n\n\n'
            #print ' mini_stack_index_num=', list(mini_stack_index_num)
            #print 'Relative_master_mini_stack=',Relative_master_mini_stack
            
            #print 'mini_stack_combine_list=', mini_stack_combine_list
            #print 'mini_stack_file_name_list=',mini_stack_file_name_list

            if int(Relative_master_mini_stack[0])==int(master_date):
                Relative_MASTER_DATE_INDEX=mini_stack_file_name_list.index(Relative_master_mini_stack[0])
                
            else:
                Relative_MASTER_DATE_INDEX=mini_stack_file_name_list.index(Relative_master_mini_stack[0])
            #print 'Relative_MASTER_DATE_INDEX=',Relative_MASTER_DATE_INDEX
            print 'Relative_master_mini_stack=',mini_stack_file_name_list[Relative_MASTER_DATE_INDEX]
            #print mini_stack_file_name_list
            master_mini_stack_stack.append(int(mini_stack_file_name_list[Relative_MASTER_DATE_INDEX]))
        master_mini_stack_stack =list(set(master_mini_stack_stack))
        master_mini_stack_stack=sorted(master_mini_stack_stack)
        
        #print  'master_mini_stack_stack=',master_mini_stack_stack
        master_mini_stack_list=[]
        for i_temp_master in master_mini_stack_stack:
            master_mini_stack_list.append(str(i_temp_master))
        print 'master_mini_stack_list=',master_mini_stack_list

        target_run_path= work_dir+ '/'+ks_stack_num +'/'+'INSAR_'+insar_file_name.split('\t')[2].strip()
        #print 'target_run_path=',target_run_path
        mini_stack_para=target_run_path+'/'+'master_mini_stack.txt'
        file_master_mini_stack=open(mini_stack_para,'w')
        mini_stack_combine_list=[]
        for master_i_temp in master_mini_stack_stack:
            for master_j_temp in master_mini_stack_stack:
                if master_i_temp < master_j_temp:
                    print master_i_temp,master_j_temp
                    line_output=str(master_i_temp)+' '+str(master_j_temp)+'\n'
                    file_master_mini_stack.write(line_output)
                    line_output=str(master_i_temp)+'_'+str(master_j_temp)
                    mini_stack_combine_list.append(line_output.strip())                     
        file_master_mini_stack.close()
                    #file_mini_stack.write('\n')

        mini_stack_para=target_run_path+'/'+'small_baselines.list'
        file_master_mini_stack=open(mini_stack_para,'w')
        for master_i_temp in master_mini_stack_stack:
            for master_j_temp in master_mini_stack_stack:
                if master_i_temp < master_j_temp:
                    print master_i_temp,master_j_temp
                    line_output=str(master_i_temp)+' '+str(master_j_temp)+'\n'
                    #line_output=str(master_i_temp)+' '+str(master_j_temp)+'\n'
                    file_master_mini_stack.write(line_output)
        file_master_mini_stack.close()

        input_cmd=[]
        input_cmd.append(target_run_path)
        input_cmd.append(master_mini_stack_list)
        input_cmd.append(mini_stack_combine_list)
        #if int(Relative_master_mini_stack[0])==int(master_date):
        Relative_MASTER_DATE_INDEX=master_mini_stack_list.index(master_date)    
        input_cmd.append(Relative_MASTER_DATE_INDEX)
        
        if int(master_mini_stack_list[Relative_MASTER_DATE_INDEX])==int(master_date):
            #print input_cmd
            input_cmd_list.append(input_cmd)
            #PCA_insar(input_cmd)

        patch_list_name=target_run_path+'/'+'patch.list'
        
        if os.path.exists(patch_list_name):
            f=open(patch_list_name)
            fcntl.flock(f.fileno(), fcntl.LOCK_EX)
            patch_list_content=f.readlines()
            fcntl.flock(f.fileno(),fcntl.LOCK_UN)
            f.close()
            for patch_file in (patch_list_content):
                patch_Name=patch_file.strip()
                
                Every_patch_input_cmd=[]
                for input_fun_cmd in input_cmd_list:
                    function_input_cmd=copy.deepcopy(input_fun_cmd)
                    function_input_cmd.append(patch_Name)
                    print function_input_cmd
                    Every_patch_input_cmd.append(function_input_cmd)
                    #print 'Every_patch_input_cmd=',Every_patch_input_cmd
                    #print 'function_input_cmd=',function_input_cmd
                    #print 'function_input_cmd=',type(function_input_cmd[1][0])
                    #if function_input_cmd[1][0]=='20190329':
                    PCA_insar(function_input_cmd) 
                #pool.map(PCA_insar,Every_patch_input_cmd)
         

        ###################################################
if __name__== '__main__':
    main()
