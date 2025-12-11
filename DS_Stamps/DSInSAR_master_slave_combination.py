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

Rapid_Stamps_SCR=os.environ['DORIS_SCR']
#print os.environ['DORIS_SCR']
#print os.environ#['STAMPS']
InSAR_cmd='doris '
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
        os.symlink(SLC_path_file ,Link_DATA)
            
    Cint_Data=freadbk(Path_MFF_HDR,Patch_First_Line,Patch_First_Pixel,Patch_nofLines,Patch_nofPixels)
    R.release()
    return Cint_Data

#############################################################################################
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
################################################################################
def process_control(ifg_para_file,process_control_name,rmstep_flag):
    ifg_process_control= int(get_parameter(process_control_name, ifg_para_file,2,'Start_process_control','End_process_control'))
    if ifg_process_control==1:
        if rmstep_flag==1:
            python_instruction='doris.rmstep.sh '+' '+process_control_name+' '+ifg_para_file
            with os.popen(python_instruction,"r") as p:
                cmd_profile_copy=p.read()
            print 'cmd_profile_copy=',cmd_profile_copy
        return 1 
    else:
        return 0
             
#***********************************************################################
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


#############################################################################################


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
#################################################################################################
#def phase_linking_paral() 
#******************************************************************************
def Doris_coherence_operation(cmd_path):
    global Rapid_Stamps_SCR
    if os.path.exists(cmd_path):
        os.chdir(cmd_path)
        print "cmd_path=",cmd_path
        ifg_para_file=os.getcwd()+'/'+'interferogram.out'
        RapidSAR_stamps_coherence  =Rapid_Stamps_SCR+'/'+'input.coherence'
        process_control_flag=process_control(ifg_para_file,'coherence',1)



        python_instruction=InSAR_cmd+' '+RapidSAR_stamps_coherence
        cmd_profile_copy=Linux_cmd(python_instruction)
        print 'cmd_profile_copy=',cmd_profile_copy
        print "cmd_path=",cmd_path
    
#####################################################################################
#******************************************************************************
def Doris_filter_operation(cmd_path):
    global Rapid_Stamps_SCR
    if os.path.exists(cmd_path):
        os.chdir(cmd_path)
        print "cmd_path=",cmd_path
        ifg_para_file         =os.getcwd()+'/'+'interferogram.out'
        RapidSAR_stamps_dem   =Rapid_Stamps_SCR+'/'+'input.filtphase'
        process_control_flag  =process_control(ifg_para_file,'filtphase',1)

        python_instruction    =InSAR_cmd+' '+RapidSAR_stamps_dem  # ../dem.dorisin'
        cmd_profile_copy      =Linux_cmd(python_instruction)
        print 'cmd_profile_copy=',cmd_profile_copy
        print "cmd_path=",cmd_path
        
#******************************************************************************
#******************************************************************************
def Doris_dem_operation(cmd_path):
    global Rapid_Stamps_SCR
    if os.path.exists(cmd_path):
        os.chdir(cmd_path)
        print "cmd_path=",cmd_path
        ifg_para_file=os.getcwd()+'/'+'interferogram.out'
        RapidSAR_stamps_dem   ='input.dem'
        
        process_control_flag=process_control(ifg_para_file,'comp_refdem',1)

        python_instruction=InSAR_cmd+' '+RapidSAR_stamps_dem  # ../dem.dorisin'
        cmd_profile_copy=Linux_cmd(python_instruction)
        print 'cmd_profile_copy=',cmd_profile_copy
        print "cmd_path=",cmd_path
        
#******************************************************************************
def Doris_interferogram_operation(cmd_path):

    global Rapid_Stamps_SCR
    if os.path.exists(cmd_path):
        os.chdir(cmd_path)
        print "cmd_path=",cmd_path
        ifg_para_file=os.getcwd()+'/'+'interferogram.out'

        RapidSAR_stamps_interferogram   =Rapid_Stamps_SCR+'/'+'interferogram.dorisin'
        process_control_flag=process_control(ifg_para_file,'interfero',1) 
        #print "\nDo  comp_refphase Operation!!!!\n"
        #if process_control_flag=        

        python_instruction=InSAR_cmd+' '+RapidSAR_stamps_interferogram  # ../dem.dorisin'
        cmd_profile_copy=Linux_cmd(python_instruction)
        print 'cmd_profile_copy=',cmd_profile_copy
        print "cmd_path=",cmd_path
        
#####################################################################################
def Interferogram_combination(input_fun_cmd):
# by wu wen hao 
# 2021 12 29
    work_dir=input_fun_cmd[0]
    original_work_dir=input_fun_cmd[-1]
    target_run_path=input_fun_cmd[0]
    #print 'input_fun_cmd=',input_fun_cmd
    #print 'input_fun_cmd[0]=',input_fun_cmd[0]
    print 'work_dir=',work_dir
    Master_date_time=work_dir.split('/')[-1]
    Master_date_time=Master_date_time[-8:]
    print 'Master_date_time=', Master_date_time
    #print 'Master_time=', type(Master_time)
   
    resFilename=work_dir+'/'+'master.res'
    l0 = int(get_parameter('First_line (w.r.t. original_image)',resFilename,1))
    lN = int(get_parameter('Last_line (w.r.t. original_image)',resFilename,1))
    p0 = int(get_parameter('First_pixel (w.r.t. original_image)',resFilename,1))
    pN = int(get_parameter('Last_pixel (w.r.t. original_image)',resFilename,1))
    Total_nofLines=lN-l0+1
    Total_nofPixels=pN-p0+1

    
    SBAS_work_dir=work_dir+'/SMALL_BASELINES'
    if os.path.exists(SBAS_work_dir):    
        print 'SBAS_work_dir=',SBAS_work_dir    
    else:
        print "This path is wrong !!!"
        sys.exit(1)
    master_slave_combination=input_fun_cmd[2]
    print 'master_slave_combination=',master_slave_combination
    print '\n\n\n'
    SBAS_master_slave_list=[]
    for temp_master_slave in master_slave_combination:
        print 'temp_master_slave=',temp_master_slave
        master_time=temp_master_slave.split('_')[0].strip()
        slave_time=temp_master_slave.split('_')[1].strip()
        print 'master_time=',master_time
        print 'slave_time=',slave_time
        res_file_name=SBAS_work_dir+'/'+temp_master_slave+'/'+'interferogram.out'
        if os.path.isfile(res_file_name):
            if int(master_time)==int(Master_date_time):
                SLC_path_file=work_dir+'/'+slave_time+'/'+'cint.minrefdem.raw'
                data_format='complex_real4'
                cint_srd_data=SAR_data_read(SLC_path_file,data_format,\
                        Total_nofLines,Total_nofPixels,1,1,Total_nofLines,Total_nofPixels)
                #Virtual_master_data=SAR_data_read(SLC_path_file,data_format,\
                #         Patch_nofLines,Patch_nofPixels,1,1,Patch_nofLines,Patch_nofPixels)  
            elif int(slave_time)==int(Master_date_time):
                SLC_path_file=work_dir+'/'+master_time+'/'+'cint.minrefdem.raw'
                data_format='complex_real4'
                cint_srd_data=np.conj(SAR_data_read(SLC_path_file,data_format,\
                        Total_nofLines,Total_nofPixels,1,1,Total_nofLines,Total_nofPixels))
                
            else: 
                SLC_path_file=work_dir+'/'+master_time+'/'+'cint.minrefdem.raw'
                data_format='complex_real4'
                master_data=SAR_data_read(SLC_path_file,data_format,\
                        Total_nofLines,Total_nofPixels,1,1,Total_nofLines,Total_nofPixels)
                
                SLC_path_file=work_dir+'/'+slave_time+'/'+'cint.minrefdem.raw'
                data_format='complex_real4'
                slave_data=SAR_data_read(SLC_path_file,data_format,\
                        Total_nofLines,Total_nofPixels,1,1,Total_nofLines,Total_nofPixels)
                cint_srd_data=slave_data*np.exp(-1.0j*np.angle(master_data))*np.abs(master_data)
 
            Output_file    =SBAS_work_dir+'/'+temp_master_slave+'/'+'cint.minrefdem.raw' 
            outformat= 'complex_real4'
            #print 'outformat=',outformat 
            print 'Output_file=',Output_file 
            print '\n\n\n'    
            SAR_write(Output_file,cint_srd_data,outformat)
#******************************************************************************        
#####################################################################################
def coherence_operation(input_fun_cmd):
    work_dir=input_fun_cmd[0]
    original_work_dir=input_fun_cmd[-1]
    target_run_path=input_fun_cmd[0]
    #print 'input_fun_cmd=',input_fun_cmd
    #print 'input_fun_cmd[0]=',input_fun_cmd[0]
    print 'work_dir=',work_dir
    Master_time=work_dir.split('/')[-1]
    Master_time=Master_time[-8:]
    #print 'Master_time=', Master_time
    #print 'Master_time=', type(Master_time)
    Master_date_time=work_dir.split('/')[-1]
    Master_date_time=Master_date_time[-8:]
    print 'Master_date_time=', Master_date_time
    SBAS_work_dir         = work_dir+'/SMALL_BASELINES'
    
    """
    original_source_input = original_work_dir+'/'+'input.coherence'
    target_source_res=SBAS_work_dir+'/'+'input.coherence'
    shutil.copy(original_source_input,target_source_res)   
     
    original_source_input=original_work_dir+'/'+'dem.dorisin'
    target_source_res=SBAS_work_dir+'/'+'dem.dorisin'
    shutil.copy(original_source_input,target_source_res) 
    """
    python_instruction='dem_dorisin_replace.py  '+work_dir+'/'+'dem.dorisin'+' '+work_dir+'/'+'input.dem'
    cmd_profile_copy=Linux_cmd(python_instruction)
    print 'cmd_profile_copy=',cmd_profile_copy


    start_time = datetime.datetime.now()
    if os.path.exists(SBAS_work_dir):    
        print 'SBAS_work_dir=',SBAS_work_dir    
    else:
        os.makedirs(SBAS_work_dir)
        #print "This path is wrong !!!"
        #sys.exit(1)
    master_slave_combination=input_fun_cmd[2]
    print 'master_slave_combination=',master_slave_combination
    SBAS_master_slave_list=[]
    for temp_master_slave in master_slave_combination:
        #print 'temp_master_slave=',temp_master_slave
        master_time=temp_master_slave.split('_')[0]
        slave_time=temp_master_slave.split('_')[1]
        #print 'master_time=',master_time
        #print 'slave_time=',slave_time
        #original_source_input=work_dir+'/'+master_time.strip()+'/'+'interferogram.out'
        #target_source_res    =SBAS_work_dir+'/'+'interferogram.res'
        
        #if os.path.isfile(original_source_input):
        #    shutil.copy(original_source_input,target_source_res)
        SBAS_file_fold=SBAS_work_dir+'/'+temp_master_slave.strip()
        print 'SBAS_file_fold=',SBAS_file_fold
        print os.path.exists(SBAS_file_fold)
        if os.path.exists(SBAS_file_fold):    
            print 'SBAS_work_dir=',SBAS_work_dir    
        else:
            os.makedirs(SBAS_file_fold)


        if os.path.exists(SBAS_file_fold):
            #os.chdir(SBAS_file_fold)
            print 'SBAS_file_fold=',SBAS_file_fold
            print '\n'

                       
            if int(master_time)==int(Master_date_time):
                SLC_res_file=work_dir+'/'+'master.res'
            else:
                SLC_res_file=work_dir+'/'+master_time.strip()+'/'+'slave.res'
            target_source_res    =SBAS_work_dir+'/'+temp_master_slave+'/'+'master.res'
            if os.path.isfile(SLC_res_file):
                #shutil.copy(SLC_res_file,target_source_res)
                print "wuwenhao"
            if int(slave_time)==int(Master_date_time):
                SLC_res_file=work_dir+'/'+'master.res'
            else:
                SLC_res_file=work_dir+'/'+master_time.strip()+'/'+'slave.res'
            target_source_res    =SBAS_work_dir+'/'+temp_master_slave+'/'+'slave.res'
            if os.path.isfile(SLC_res_file):
                #shutil.copy(SLC_res_file,target_source_res)
                print "wuwenhao"

            original_dem_input  =work_dir+'/'+'input.dem'
            target_dem_input    =SBAS_work_dir+'/'+temp_master_slave+'/'+'input.dem'

            original_source_input=work_dir+'/'+master_time.strip()+'/'+'interferogram.out'
            target_source_res    =SBAS_work_dir+'/'+temp_master_slave+'/'+'interferogram.out'
            if os.path.isfile(original_source_input):
                shutil.copy(original_dem_input,target_dem_input)
                SBAS_master_slave_list.append(SBAS_file_fold)
            else:
                original_source_input=work_dir+'/'+slave_time.strip()+'/'+'interferogram.out'
                target_source_res    =SBAS_work_dir+'/'+temp_master_slave+'/'+'interferogram.out'
                if os.path.isfile(original_source_input):
                    shutil.copy(original_dem_input,target_dem_input)
                    SBAS_master_slave_list.append(SBAS_file_fold)
                #Doris_operation(SBAS_file_fold)
            #os.chdir(SBAS_work_dir)
            
    print 'SBAS_master_slave_list=',SBAS_master_slave_list
    cores=multiprocessing.cpu_count()
    #print 'Total cores is ', cores
    pool=multiprocessing.Pool(processes=(int(cores-1)))
    pool.map(Doris_dem_operation,SBAS_master_slave_list)
    pool.map(Doris_interferogram_operation,SBAS_master_slave_list)
    pool.map(Doris_coherence_operation,SBAS_master_slave_list)
   
    #pool.map(Doris_filter_operation,SBAS_master_slave_list)
####################################################################################
#################################################################################################
#################################################################################################
   
     
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
    original_large_mini_stack=[]
    Interfergrom_combination_info=[]

    for filename in open('Int_Data_merge.list'):
    #print 'filename=',filename
        master_date=os.path.basename(filename.strip()).split('_')[0].strip()
        slave_date=os.path.basename(filename.strip()).split('_')[1].strip()
        filename_list.append(slave_date.strip()) 
        Interfergrom_combination_info.append(filename.strip())
        filename_list.append(master_date)
    print 'filename_list=',filename_list

    filename_list =list(set(filename_list))
    filename_list=sorted(filename_list)
    #time.sleep(1000000000000)
    #******************************************************************************
    

    if len(filename_list) >0:
       
        Master_Date_fold=work_dir+'/'+'INSAR_'+master_date
        if not os.path.exists(Master_Date_fold):
            #os.makedirs(Master_Date_fold)
            print "can not find the fold"
        
##################################################################
        else:
        
            mini_stack_combine_list=[]
            
            print '\n\n\n'

            Relative_MASTER_DATE_INDEX=filename_list.index(master_date.strip())
                           
            #print 'Relative_MASTER_DATE_INDEX=',Relative_MASTER_DATE_INDEX           
            #print mini_stack_file_name_list
          
            target_run_path= work_dir+ '/'+'INSAR_'+master_date
            print 'target_run_path=',target_run_path
        

            mini_stack_para=target_run_path+'/'+'small_baselines.list'
            file_master_mini_stack=open(mini_stack_para,'w')
            for master_i_temp in filename_list:
                for master_j_temp in filename_list:
                    if int(master_i_temp) < int(master_j_temp):
                        print master_i_temp,master_j_temp
                        line_output=str(master_i_temp)+' '+str(master_j_temp)+'\n'
                        file_master_mini_stack.write(line_output)
                        line_output=str(master_i_temp)+'_'+str(master_j_temp)
                        mini_stack_combine_list.append(line_output.strip())
            file_master_mini_stack.close()

            if os.path.exists(target_run_path):

                os.chdir(target_run_path)

                python_instruction='make_small_baselines 1 1'
                cmd_profile_copy=Linux_cmd(python_instruction)
                print 'cmd_profile_copy=',cmd_profile_copy
 
                os.chdir(work_dir)


                input_cmd=[]
                input_cmd.append(target_run_path)
                input_cmd.append(filename_list)
                input_cmd.append(mini_stack_combine_list)
        
                Relative_MASTER_DATE_INDEX=filename_list.index(master_date)    
                input_cmd.append(Relative_MASTER_DATE_INDEX)
        
                if int(filename_list[Relative_MASTER_DATE_INDEX])==int(master_date):
           
                    input_cmd.append(work_dir)
                    print 'input_cmd=',input_cmd
           
                    coherence_operation(input_cmd)
                    Interferogram_combination(input_cmd)
    else:            
        print "can not find the fold" 

                    
        ###################################################
if __name__== '__main__':
    main()
