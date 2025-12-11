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

#Version 20181031
#Version 20190831
#Version 20190902 # correct the Bug:100 101 for three digital
#####################################################################################

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
    return thisBurstData
###############################################################################
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
def PCA_insar(patch_Name):
    work_dir=os.getcwd()
    Master_time=work_dir.split('/')[-1]
    Master_time=Master_time[-8:]
    print 'Master_time=', Master_time
    print 'Master_time=', type(Master_time)

    SBAS_work_dir=work_dir+'/SMALL_BASELINES'
    print 'SBAS_work_dir=',SBAS_work_dir
    start_time = datetime.datetime.now()
    if os.path.exists(SBAS_work_dir):    
        file_dir=os.listdir(SBAS_work_dir)
    else:
        print "This path is wrong !!!"
        sys.exit(1)

    filename_list=[]
    for filename in file_dir:    
        m_index=re.search('\d\d\d\d_\d\d\d\d',filename)      
        if m_index is not None:
            print 'filename of Interferometry folder =',filename 
            filename_list.append(filename)

    filename_list.sort()
     
    cint_minrefdem_file=filename_list
    Num_Int =len(cint_minrefdem_file)
    print 'Num_Int=',Num_Int
    ###Next part
    #
    SLC_list_filename='slcs.list'   #not include master

    SLC_list_Index=[]
    Num_SLC=0
    for line in open(SLC_list_filename):
        SLC_list_Index.append(line.strip())
        #print line.strip()
        Num_SLC=Num_SLC+1
    SLC_list_Index.append(Master_time.strip())
    print SLC_list_Index
    print type(SLC_list_Index)
    SLC_list_Index.sort()
    print SLC_list_Index
    Num_SLC=Num_SLC+1
    print SLC_list_Index
    print 'Num_SLC=',Num_SLC

    MASTER_DATE_INDEX=SLC_list_Index.index(Master_time)
    print 'MASTER_DATE_INDEX=',MASTER_DATE_INDEX
    ##############################################################
    ##Parameter of filter Windows
    #############################################################

    #win_nofLines   =13
    #win_nofPixels  =19
    ##############################################################
    ##Parameter of filter Windows
    #############################################################
    KS_folder=os.getcwd()+'/KS_result'
    Homogenous_param = KS_folder+'/'+'homo_wins.par'
    if os.path.exists(Homogenous_param):
        win_nofLines = int(get_parameter('win_nofLines',Homogenous_param,1))
        win_nofPixels = int(get_parameter('win_nofPixels',Homogenous_param,1))
        amplitude_dispersion_param_value=float(get_parameter('amplitude_dispersion_param_value',Homogenous_param,1))
        sum_filter_wins = int(get_parameter('sum_filter_wins',Homogenous_param,1))
    else:
        print "I can not find the homo_wins.par"
        sys.exit(1)
    ######################################################
    Half_win_nofLines   =(win_nofLines-1)/2
    Half_win_nofPixels  =(win_nofPixels-1)/2
    
    Cen_Sta_Lines=(win_nofLines-1)/2
    Cen_Sta_Pixels=(win_nofPixels-1)/2
    Total_win_pix=win_nofLines*win_nofPixels
    win_mask=np.zeros((Total_win_pix),dtype=np.int8)
    Number_win_mask=50000
    win_mask_total=np.zeros((Total_win_pix*Number_win_mask),dtype=np.int8)
    work_dir=os.getcwd();
    unpack_format=str(Total_win_pix)+'c'
    #i_win=0

    KS_Estimation = os.getcwd()+'/KS_result'
    KS_folder=os.getcwd()+'/KS_result'

    #f=open('patch.list')
    #patch_file=f.readlines()
    #for patch_Name in patch_file:
    
    print 'patch_Name.strip()=',patch_Name.strip()
    #print "patch_Name.split('_')=",patch_Name.split('_')
    patch_Data_Dimension=patch_Name.strip()+'/'+'patch.in'
    print 'patch_Data_Dimension=',patch_Data_Dimension
    if os.path.exists(patch_Data_Dimension)== True:

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
        Int_Data_Stack=np.zeros((int(Patch_nofPixels*Patch_nofLines),int(Num_SLC),int(Num_SLC)),dtype=np.complex64)
        coh_mx = np.zeros((int(Patch_nofPixels*Patch_nofLines),int(Num_SLC),int(Num_SLC)),dtype=float32)
        #np.fill_diagonal(coh_mx, 1.0) # setting diag(coh_fisher) = 1.0
        #--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#
        
        for i_j_temp in np.arange(0,int(Num_SLC)):
            Int_Data_Stack[:,i_j_temp,i_j_temp] = 1
            
            coh_mx[:,i_j_temp,i_j_temp] = 1 # setting diag(coh_fisher) = 1.0
        for dataFilename in cint_minrefdem_file:
            SLCFilename=SBAS_work_dir+'/'+dataFilename.strip()[0:8]+'_'+dataFilename.strip()[9:17]
            #####################################################################################
            Imag_Master=dataFilename.strip()[0:8]
            Imag_Slave =dataFilename.strip()[9:17]
            print 'Imag_Master=' ,Imag_Master
            print 'Imag_Slave=' ,Imag_Slave
                                
            Matrix_subscript_Master=SLC_list_Index.index(Imag_Master)
            Matrix_subscript_Slave =SLC_list_Index.index(Imag_Slave) 
            ####################################################################################
            resFilename=SLCFilename+'/'+'master.res'
            l0 = int(get_parameter('First_line (w.r.t. original_image)',resFilename,1))
            lN = int(get_parameter('Last_line (w.r.t. original_image)',resFilename,1))
            p0 = int(get_parameter('First_pixel (w.r.t. original_image)',resFilename,1))
            pN = int(get_parameter('Last_pixel (w.r.t. original_image)',resFilename,1))
            Total_nofLines=lN-l0+1
            Total_nofPixels=pN-p0+1

            SLC_path_file       =SLCFilename+'/'+'cint.minrefpha'+'_'+patch_Name.strip()+'.raw'
            print SLC_path_file
            Path_MFF_HDR   =SLCFilename+'/'+'cint.minrefpha'+'_'+patch_Name.strip()+'.hdr'
            Link_DATA      =SLCFilename+'/'+'cint.minrefpha'+'_'+patch_Name.strip()+'.x00'

            print 'Path_MFF_HDR=',Path_MFF_HDR 
            if os.path.exists(Path_MFF_HDR)==True:
                os.remove(Path_MFF_HDR)
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
                os.symlink(SLC_path_file ,Link_DATA)
            
            Cint_Dem=(np.angle(freadbk(Path_MFF_HDR,1,1,Patch_nofLines,Patch_nofPixels))).copy()
                       
            if (Matrix_subscript_Master > Matrix_subscript_Slave): 
                Int_Data_Stack[:,Matrix_subscript_Master,Matrix_subscript_Slave]=np.exp(-1j*Cint_Dem).reshape(1,-1)
                Int_Data_Stack[:,Matrix_subscript_Slave,Matrix_subscript_Master]=np.exp( 1j*Cint_Dem).reshape(1,-1)                       
            else:
                Int_Data_Stack[:,Matrix_subscript_Master,Matrix_subscript_Slave]=np.exp( 1j*Cint_Dem).reshape(1,-1)                        
                Int_Data_Stack[:,Matrix_subscript_Slave,Matrix_subscript_Master]=np.exp(-1j*Cint_Dem).reshape(1,-1)

            file_count=file_count+1
            print 'file_count=',file_count
            print '\n\n\n'
        ####################################################################################
        file_count=0
        
        for dataFilename in cint_minrefdem_file:
            SLCFilename=SBAS_work_dir+'/'+dataFilename.strip()[0:8]+'_'+dataFilename.strip()[9:17]
            #####################################################################################
            Imag_Master=dataFilename.strip()[0:8]
            Imag_Slave =dataFilename.strip()[9:17]
            print 'Imag_Master=' ,Imag_Master
            print 'Imag_Slave=' ,Imag_Slave
                                
            Matrix_subscript_Master=SLC_list_Index.index(Imag_Master)
            Matrix_subscript_Slave =SLC_list_Index.index(Imag_Slave) 
            ####################################################################################
            SLC_path_file       =SLCFilename+'/'+'coherence.raw'
            print SLC_path_file
            Path_MFF_HDR   =SLCFilename+'/'+'coherence'+'.hdr'
            Link_DATA      =SLCFilename+'/'+'coherence'+'.r00'

            print 'Path_MFF_HDR=',Path_MFF_HDR 
            #if os.path.exists(Path_MFF_HDR)==True:
            #    os.remove(Path_MFF_HDR)
            if os.path.exists(Path_MFF_HDR)==False:
       
                outStream      = open(Path_MFF_HDR,'w')
                outStream.write('IMAGE_FILE_FORMAT = MFF\n')
                outStream.write('FILE_TYPE = IMAGE\n')
                outStream.write('IMAGE_LINES = %d\n' % int(Total_nofLines))
                outStream.write('LINE_SAMPLES = %d\n'% int(Total_nofPixels))
                outStream.write('BYTE_ORDER = LSB\n')
                outStream.write('END\n')
                outStream.close()

            if (os.path.isfile(Link_DATA)):
                #os.remove(Link_DATA) 
                print  'Link_DATA=' ,Link_DATA 
            else:      
                os.symlink(SLC_path_file ,Link_DATA)
        
            Coherence_Data=(freadbk(Path_MFF_HDR,Patch_First_Line+1,Patch_First_Pixel+1,Patch_nofLines,Patch_nofPixels)).copy()
            #print Coherence_Data_Stack[file_count,91:93,108:110]
            
            Int_Data_Stack[:,Matrix_subscript_Master,Matrix_subscript_Slave]=Int_Data_Stack[:,Matrix_subscript_Master,Matrix_subscript_Slave]*Coherence_Data.reshape(1,-1)
            Int_Data_Stack[:,Matrix_subscript_Slave,Matrix_subscript_Master]=Int_Data_Stack[:,Matrix_subscript_Slave,Matrix_subscript_Master]*Coherence_Data.reshape(1,-1)                       
            
            #--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#
            #--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#
            # coherence matrix, used for EMI, Homa Ansari
            coh_mx[:,Matrix_subscript_Master,Matrix_subscript_Slave] = Coherence_Data.reshape(1,-1)
            coh_mx[:,Matrix_subscript_Slave,Matrix_subscript_Master] = Coherence_Data.reshape(1,-1)
            #--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#
            #--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#
            file_count=file_count+1
            print 'file_count=',file_count 
            print '\n\n\n'   
        ########################################################################################
        # Read_Master_Slave
        ########################################################################################
        #
        resFilename='master.res'
        l0 = int(get_parameter('First_line (w.r.t. original_image)',resFilename,1))
        lN = int(get_parameter('Last_line (w.r.t. original_image)',resFilename,1))
        p0 = int(get_parameter('First_pixel (w.r.t. original_image)',resFilename,1))
        pN = int(get_parameter('Last_pixel (w.r.t. original_image)',resFilename,1))
        Total_nofLines=lN-l0+1
        Total_nofPixels=pN-p0+1
        file_count=0
        f=open('day.1.in')
        cint_day_file=f.readlines()
        Single_SLC_Data_Stack=np.zeros((len(cint_day_file),int(Patch_nofLines),int(Patch_nofPixels)),dtype=np.complex64)
        for file_cint  in cint_day_file:    
            SLC_path_file  =os.getcwd()+'/'+file_cint.strip()+'/'+'cint.minrefdem.raw'
            Output_file    =os.getcwd()+'/'+file_cint.strip()+'/'+'cint.minrefpha.raw'
     
            Path_MFF_HDR   =os.getcwd()+'/'+file_cint.strip()+'/'+'cint.minrefdem'+'.hdr'
            Link_DATA      =os.getcwd()+'/'+file_cint.strip()+'/'+'cint.minrefdem'+'.x00'
      
            print 'Path_MFF_HDR=',Path_MFF_HDR 
      
            #if os.path.exists(Path_MFF_HDR)==True:
            #    os.remove(Path_MFF_HDR)
            if os.path.exists(Path_MFF_HDR)==False:
       
                outStream      = open(Path_MFF_HDR,'w')
                outStream.write('IMAGE_FILE_FORMAT = MFF\n')
                outStream.write('FILE_TYPE = IMAGE\n')
                outStream.write('IMAGE_LINES = %d\n' % int(Total_nofLines))
                outStream.write('LINE_SAMPLES = %d\n'% int(Total_nofPixels))
                outStream.write('BYTE_ORDER = LSB\n')
                outStream.write('END\n')
                outStream.close()

            if (os.path.isfile(Link_DATA)):
                 #os.remove(Link_DATA)
                print  'Link_DATA=' ,Link_DATA
            else:       
                os.symlink(SLC_path_file ,Link_DATA)

     
            Single_SLC_Data_Stack[file_count,:,:]=freadbk(Path_MFF_HDR,Patch_First_Line+1,Patch_First_Pixel+1,Patch_nofLines,Patch_nofPixels).copy()
       
            file_count=file_count+1
            print 'file_count=',file_count 
            print '\n\n'  
        Single_SLC_KS =Single_SLC_Data_Stack.copy()  
    
    ######################################################################################
        Total_coherence=0    
        #print 'Num_Int=',Num_Int
        #Expected_coherence  = Num_Int/4

        #Vari_Num_PS_Cint       =len(cint_day_file)
        #para_insar             =np.zeros((Vari_Num_PS_Cint),dtype=float32)
        #Observation_data_array =np.zeros((Num_Int),dtype=float32)
        #Observation_coherence_array=np.zeros((Num_Int),dtype=float32)
        #PCA_Energy_thres=3.0/Num_SLC
        Coherence_Ind=np.append(xrange(MASTER_DATE_INDEX),xrange(MASTER_DATE_INDEX+1,Num_SLC))
        print 'Coherence_Ind=',Coherence_Ind
        PCA_Energy_thres=0.03
        Total_line=0L
        psc_da=patch_Name.strip()+'/'+'pscands.1.da'
        psc_ij=patch_Name.strip()+'/'+'pscands.1.ij'
        if os.path.exists(KS_folder)==False:
            print "KS_result is not exit"
        else:
            file_KS_Estimation=KS_folder+'/'+patch_Name.strip()+'_data.dat'
            print file_KS_Estimation
            Line_index_ks=0
            f=open(file_KS_Estimation,'rb')
            with open(psc_da,'r') as fp1, open(psc_ij,'r') as fp2:
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
                        #temp=f.read(Total_win_pix)
                         #print 'temp=',temp
                        #win_mask=int8(struct.unpack(unpack_format,temp))
                        sum_filter=win_mask.sum()   
                    
                        if sum_filter > sum_filter_wins:                            
                            #Total_coherence=np.sum(Observation_coherence_array)
                            #para_insar=np.angle(Single_SLC_Data_Stack[:,Relative_line_ind,Relative_pixel_ind])
                            #print 'para_insar=',para_insar
                            # Initival Value                                                                                                             
                            #if ( (Optimization_res.fun/Total_coherence) < -0.79) and (Total_coherence  > Expected_coherence):
                                                 
                            #print  'max_eig_val/np.sum(abs_eig_val)=',max_eig_val/np.sum(abs_eig_val)
                            #if (max_eig_val/np.sum(abs_eig_val)) > 0.169:
                            ############################################################
                            #coh_mx_inv = np.linalg.pinv(np.abs(Int_Data_Stack[Relative_line_ind*Patch_nofPixels+Relative_pixel_ind,:,:])) # seudo inverse of coherence matrix
                            
                            
                            eig_val, eig_vec = np.linalg.eig(Int_Data_Stack[Relative_line_ind*Patch_nofPixels+Relative_pixel_ind,:,:])
                            abs_eig_val      = np.abs(eig_val)
                            max_eig_val      = np.max(abs_eig_val) 
                            if (max_eig_val/np.sum(abs_eig_val)) > PCA_Energy_thres:
                                    #print  'max(eig_val)/np.sum(eig_val)=',np.max(eig_val)/np.sum(eig_val)
                                    #eig_val_Ind=np.argmax(eig_val) 
                                    #print  'eig_val_Ind=',eig_val_Ind 
                                    #print  'eig_val_Ind=',eig_val_Ind.shape  

                                                                
                                red_eig_vec=np.angle(eig_vec[:,np.argmax(abs_eig_val)])
                                red_eig_vec=red_eig_vec[MASTER_DATE_INDEX]-red_eig_vec
                                #Optimization_pca=np.delete(Optimization_pca,MASTER_DATE_INDEX)
                                #Optimization_pca=np.delete(Optimization_pca,MASTER_DATE_INDEX)
                                #print 'Optimization_pca=',Optimization_pca.shape   
                                Single_SLC_KS[:,Relative_line_ind,Relative_pixel_ind] = np.exp(1j*red_eig_vec[Coherence_Ind])*np.abs(Single_SLC_Data_Stack[:,Relative_line_ind,Relative_pixel_ind])                                   
                            else:
                                 coh_mx_inv = np.linalg.pinv(coh_mx[Relative_line_ind*Patch_nofPixels+Relative_pixel_ind,:,:])
                                 cov_mx = Int_Data_Stack[Relative_line_ind*Patch_nofPixels+Relative_pixel_ind,:,:] # covariance matrix
                            
                                 cov_EMI = coh_mx_inv*cov_mx # Eq.13 of EMI
                            
                                 eig_val, eig_vec = np.linalg.eig(cov_EMI) # eigen value and vector                                                                                                    
                                                                                                                                            
                                 if np.min(np.real(eig_val)) > 1.0: # EMI. To avoid small negative eigenvalue.
                                                               # The minimum eigenvalue should be larger than 1 (By Eq.12 & 14 of EMI).
                                     abs_eig_val      = np.abs(eig_val)
                            
                                     red_eig_vec=np.angle(eig_vec[:,np.argmin(abs_eig_val)])
                                     Optimization_pca=red_eig_vec[MASTER_DATE_INDEX]-red_eig_vec
                                     Optimization_pca=np.delete(Optimization_pca,MASTER_DATE_INDEX)
                                #Single_SLC_KS[:,Relative_line_ind,Relative_pixel_ind] = np.exp(1j*red_eig_vec[Coherence_Ind])*np.abs(Single_SLC_Data_Stack[:,Relative_line_ind,Relative_pixel_ind])
                                     Single_SLC_KS[:,Relative_line_ind,Relative_pixel_ind] = np.exp(1j*Optimization_pca)*np.abs(Single_SLC_Data_Stack[:,Relative_line_ind,Relative_pixel_ind])
                            ###########################################################
                            #else:
                            #    eig_val, eig_vec = np.linalg.eig(Int_Data_Stack[Relative_line_ind*Patch_nofPixels+Relative_pixel_ind,:,:])
                            #    abs_eig_val      = np.abs(eig_val)
                            #    max_eig_val      = np.max(abs_eig_val) 
                            #    if (max_eig_val/np.sum(abs_eig_val)) > PCA_Energy_thres:
                                    #print  'max(eig_val)/np.sum(eig_val)=',np.max(eig_val)/np.sum(eig_val)
                                    #eig_val_Ind=np.argmax(eig_val) 
                                    #print  'eig_val_Ind=',eig_val_Ind 
                                    #print  'eig_val_Ind=',eig_val_Ind.shape  

                                                                
                            #        red_eig_vec=np.angle(eig_vec[:,np.argmax(abs_eig_val)])
                            #        red_eig_vec=red_eig_vec[MASTER_DATE_INDEX]-red_eig_vec
                                #Optimization_pca=np.delete(Optimization_pca,MASTER_DATE_INDEX)
                                #Optimization_pca=np.delete(Optimization_pca,MASTER_DATE_INDEX)
                                #print 'Optimization_pca=',Optimization_pca.shape   
                            #        Single_SLC_KS[:,Relative_line_ind,Relative_pixel_ind] = np.exp(1j*red_eig_vec[Coherence_Ind])*np.abs(Single_SLC_Data_Stack[:,Relative_line_ind,Relative_pixel_ind])
                            if int(string_line_ij[0])%5000==0:
                                
                                #print 'Optimization_pca=',Optimization_pca
                                #print 'para_insar=',para_insar
                                print 'processing the patch row:',patch_Name.strip()
                                print 'processing the data row:',string_line_ij[0]
                                end_time = datetime.datetime.now()
                                print "time-consuming=",end_time-start_time
                                print '\n'
                       
            f.close() #close the file
        file_count=0 
        for file_cint  in cint_day_file:    
            Input_file     =os.getcwd()+'/'+file_cint.strip()+'/'+'cint.minrefdem.raw'
            Output_file    =os.getcwd()+'/'+file_cint.strip()+'/'+'cint.minrefpha'+'_'+patch_Name.strip()+'.raw'
            print 'Output_file=',Output_file 
            fid = open(Output_file,'wb');
            if fid < 0:
                print "ERROR :cannot open file for writing"
                sys.exit(1)
            else:
                print 'Patch_nofLines=',Patch_nofLines
                print 'Patch_nofPixels=',Patch_nofPixels
                print 'Single_SLC_KS.shape', Single_SLC_KS.shape
                for i_temp in np.arange(0,Patch_nofLines):
                    fid.write(Single_SLC_KS[file_count,i_temp,:])
            file_count=file_count+1
            fid.close()
    else:
        print "Can not find the file or file fold !!!!!!!!!!!!!!!!!!!!!!!!!"
        print "Can not find the file or file fold !!!!!!!!!!!!!!!!!!!!!!!!!"
        print "Can not find the file or file fold !!!!!!!!!!!!!!!!!!!!!!!!!"
        print "Can not find the file or file fold !!!!!!!!!!!!!!!!!!!!!!!!!"
        print "Can not find the file or file fold !!!!!!!!!!!!!!!!!!!!!!!!!" 
        sys.exit(1)
################################################################################################
def main():
     
    #cores=multiprocessing.cpu_count()
    #print 'Total cores is ', cores
    #pool=multiprocessing.Pool(processes=(2))
  
    
    resFilename='master.res'
    l0 = int(get_parameter('First_line (w.r.t. original_image)',resFilename,1))
    lN = int(get_parameter('Last_line (w.r.t. original_image)',resFilename,1))
    p0 = int(get_parameter('First_pixel (w.r.t. original_image)',resFilename,1))
    pN = int(get_parameter('Last_pixel (w.r.t. original_image)',resFilename,1))
    nofLines=lN-l0+1
    nofPixels=pN-p0+1
    dataFormat = 'cpxfloat32'
    #SLC_file_count=len(open(SLC_file).readlines())
    ##############################################################
    ##Parameter of filter Windows
    #############################################################
    ##############################################################
    ##Parameter of filter Windows
    #############################################################
    
    ###################################################################

    f=open('patch.list')
    for patch_file in f.readlines():
        PCA_insar(patch_file.strip()) 
    #patch_file=f.readlines()
    #cint_minrefdem_file=filename_list
    #print cint_minrefdem_file[0]
    #Num_SLC =len(cint_minrefdem_file)
    #print 'Num_SLC =',Num_SLC
    #pool.map(PCA_insar ,patch_file)
    #homofilter_KS(patch_Name)
    #PCA_insar(patch_file[0])    
    
if __name__== '__main__':
    main()
