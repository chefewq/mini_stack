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
#def phase_linking_paral() 
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
                #Relative_MASTER_DATE_INDEX=mini_stack_file_name_list.index(Relative_master_mini_stack[0])
                if Relative_master_mini_stack[0] in mini_stack_file_name_list:
                     Relative_MASTER_DATE_INDEX=mini_stack_file_name_list.index(Relative_master_mini_stack[0])
                else:
                     Relative_MASTER_DATE_INDEX=0 #mini_stack_file_name_list.index(Relative_master_mini_stack[0])
            #print 'Relative_MASTER_DATE_INDEX=',Relative_MASTER_DATE_INDEX
            print 'Relative_master_mini_stack=',mini_stack_file_name_list[Relative_MASTER_DATE_INDEX]
            #print mini_stack_file_name_list
            master_mini_stack_stack.append(int(mini_stack_file_name_list[Relative_MASTER_DATE_INDEX]))
        master_mini_stack_stack =list(set(master_mini_stack_stack))
        master_mini_stack_stack=sorted(master_mini_stack_stack)
        
        print  'master_mini_stack_stack=',master_mini_stack_stack
        target_run_path= work_dir+ '/'+ks_stack_num +'/'+'INSAR_'+insar_file_name.split('\t')[2].strip()
        print 'target_run_path=',target_run_path
        mini_stack_para=target_run_path+'/'+'master_mini_stack.txt'
        file_master_mini_stack=open(mini_stack_para,'w')
        for master_i_temp in master_mini_stack_stack:
            for master_j_temp in master_mini_stack_stack:
                if master_i_temp < master_j_temp:
                    print master_i_temp,master_j_temp
                    line_output=str(master_i_temp)+' '+str(master_j_temp)+'\n'
                    file_master_mini_stack.write(line_output)
        file_master_mini_stack.close()
                    #file_mini_stack.write('\n')

        mini_stack_para=target_run_path+'/'+'small_baselines.list'
        file_master_mini_stack=open(mini_stack_para,'w')
        for master_i_temp in master_mini_stack_stack:
            for master_j_temp in master_mini_stack_stack:
                if master_i_temp < master_j_temp:
                    print master_i_temp,master_j_temp
                    line_output=str(master_i_temp)+' '+str(master_j_temp)+'\n'
                    file_master_mini_stack.write(line_output)
        file_master_mini_stack.close()

        if os.path.exists(target_run_path):

            os.chdir(target_run_path)

            python_instruction='make_small_baselines 1 1'
            cmd_profile_copy=Linux_cmd(python_instruction)
            print 'cmd_profile_copy=',cmd_profile_copy
 
            os.chdir(work_dir)

        else:
            #os.makedirs(Master_Date_fold)
            print "can not find the fold" 

                    
        ###################################################
if __name__== '__main__':
    main()
