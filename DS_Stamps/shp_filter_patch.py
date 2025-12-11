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
#Version:20180706
#Version:20180923
#Version:20181019
#Version"20181020
#Version"20181024
#Version:20181024 afternoon
#Version:20181024 evening
#version 20190903 repair the bug 
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
####
#-------------------------------------------------------------------------------------#
def homofilter_KS():


    #meminfo
    mem = {}
    f = open("/proc/meminfo")
    lines = f.readlines()
    f.close()
    for line in lines:
        if len(line) < 2: continue
        name = line.split(':')[0]
        var = line.split(':')[1].split()[0]
        mem[name] = long(var) * 1024.0
    #
    print round(mem['MemTotal']/1024/1024/1024)
    masterRes = os.getcwd()+'/'+'master.res'
    l0 = int(get_parameter('First_line (w.r.t. original_image)',masterRes,1))
    lN = int(get_parameter('Last_line (w.r.t. original_image)',masterRes,1))
    p0 = int(get_parameter('First_pixel (w.r.t. original_image)',masterRes,1))
    pN = int(get_parameter('Last_pixel (w.r.t. original_image)',masterRes,1))
    nofLines=lN-l0+1
    nofPixels=pN-p0+1
    patch_lines=len(open('patch.list').readlines())
    file_size_byte=nofPixels*nofLines*8/patch_lines
    print 'file_size_byte=',file_size_byte/1024/1024
    #
    print '\n'
    memory_file_count=int(mem['MemTotal']/3*1.15/file_size_byte)
    print 'memory_file_count=',memory_file_count
    
    SBAS_Int_file='small_baselines.list'
    myfile  =open(SBAS_Int_file)
    print myfile
    cint_minrefdem_file=myfile.readlines()
    
    filename_list=[]
    for cint_file in cint_minrefdem_file:
        print  cint_file.strip()
    file_strip= [x.strip() for x in cint_minrefdem_file if x.strip()!='' ] 
    cint_minrefdem_file=file_strip

    #print len(cint_minrefdem_file)
    Num_minrefdem_file =len(cint_minrefdem_file)
    print "There are %d files" % (Num_minrefdem_file)
    
    file_segment= math.ceil(float(Num_minrefdem_file) /memory_file_count)
    print "There are %d file segment" % (file_segment)

    
    ######################################################################################
    for file_segment_Index in np.arange(file_segment):
        
        
        #print file_segment_Index
        Start_Index=int(file_segment_Index*memory_file_count)
        End_Index  =int((file_segment_Index+1)*memory_file_count)
        print ' Start_Index,End_Index=',Start_Index,End_Index
        if End_Index > Num_minrefdem_file:
            End_Index = Num_minrefdem_file
        myfile=file_strip[Start_Index:End_Index]
        print myfile
        print ("\n\nProcessing Data\n");
        start_time=time.clock()
        print "Start time:",start_time;
        print "\n"
        ###############################
        
        #############################################################
        sum_filter_wins       =13
        #

        KS_folder=os.getcwd()+'/KS_result'
        Homogenous_param = KS_folder+'/'+'homo_wins.par'
        if os.path.exists(Homogenous_param):

            win_nofLines = int(get_parameter('win_nofLines',Homogenous_param,1))
            win_nofPixels = int(get_parameter('win_nofPixels',Homogenous_param,1))
            amplitude_dispersion_param_value=float(get_parameter('amplitude_dispersion_param_value',Homogenous_param,1))
            sum_filter_wins = int(get_parameter('sum_filter_wins',Homogenous_param,1))
            print 'win_nofLines =',win_nofLines 
            print 'win_nofPixels=',win_nofPixels
            print 'amplitude_dispersion_param_value=',amplitude_dispersion_param_value
            print 'sum_filter_wins=',sum_filter_wins
        else:
            print "I can not find the homo_wins.par"
            sys.exit(1)
        
        #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        #file_para= open(Homogenous_param,'a')
        #file_para.write('sum_filter_wins:	%d\n'%sum_filter_wins)
        #file_para.close()
        #win_nofLines   =13
        #win_nofPixels  =19
        ##############################################################
        ##Parameter of filter Windows
        #############################################################
        Half_win_nofLines   =(win_nofLines-1)/2
        Half_win_nofPixels  =(win_nofPixels-1)/2
    
        Cen_Sta_Lines=(win_nofLines-1)/2
        Cen_Sta_Pixels=(win_nofPixels-1)/2

        Cen_End_Lines =nofLines -(win_nofLines+1)/2
        Cen_End_Pixels=nofPixels-(win_nofPixels+1)/2
        Total_win_pix=win_nofLines*win_nofPixels
        win_mask=np.zeros((Total_win_pix),dtype=np.int8)
        Number_win_mask=50000
        win_mask_total=np.zeros((Total_win_pix*Number_win_mask),dtype=np.int8)        
 
        work_dir=os.getcwd();
        unpack_format=str(Total_win_pix)+'c'
        i_win=0
        
        


        KS_Estimation = os.getcwd()+'/KS_result'
        

        for patch_Name in open('patch.list'):
            patch_Data_Dimension=patch_Name.strip()+'/'+'patch.in'
            file_patch=open(patch_Data_Dimension);
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
            Patch_noover_nofPixels  = Patch_Last_Pixel- Patch_First_Pixel+1
            Patch_noover_nofLines   = Patch_Last_Line - Patch_First_Line +1

   
            Patch_noover_Sta_Lines  =Patch_noover_First_Line 
            Patch_noover_Sta_Pixels =Patch_noover_First_Pixel
            Patch_noover_End_Lines  =Patch_noover_Last_Line  
            Patch_noover_End_Pixels =Patch_noover_Last_Pixel 

            Limited_square_First_Line  =max(Patch_Sta_Lines,Patch_noover_First_Line)
            Limited_square_Last_Line   =min(Patch_End_Lines,Patch_noover_Last_Line)
            Limited_square_First_Pixel =max(Patch_Sta_Pixels,Patch_noover_First_Pixel)
            Limited_square_Last_Pixel  =min(Patch_End_Pixels,Patch_noover_Last_Pixel)
            

            ###############################
            file_count=0
            SBAS_Int_Data_Stack=scipy.zeros((int(len(myfile)),int(Patch_nofLines),int(Patch_nofPixels)),dtype=complex64)
            for cint_minrefdem_file in myfile:
                print 'cint_minrefdem_file=',cint_minrefdem_file
                cint_minrefdem_name=cint_minrefdem_file.strip()
                cint_minrefdem_name=cint_minrefdem_name[0:8]+'_'+cint_minrefdem_name[9:17]
                print 'cint_minrefdem_name=',cint_minrefdem_name
    
                #master res and ifg res files for current burst
                masterRes = os.getcwd()+'/SMALL_BASELINES/'+cint_minrefdem_name+'/'+'master.res'
                ifgRes    = os.getcwd()+'/SMALL_BASELINES/'+cint_minrefdem_name+'/'+'interferogram.out'
                print 'masterRes1=',masterRes 
                #print 'ifgRes=',ifgRes     
                #reads time and transform it to seconds
                #TODO for Taz_start near midnight this code will not work, because days are not used

                #First_pixel_azimuth_time=get_parameter('First_pixel_azimuth_time (UTC)',masterRes,3)
                #Taz_start=float(First_pixel_azimuth_time[0])*3600+float(First_pixel_azimuth_time[1])*60+float(First_pixel_azimuth_time[2])
                #Taz_start= Taz_start+ (l1_1-1)/azimFreq   
    	
                l0 = int(get_parameter('First_line (w.r.t. original_image)',masterRes,1))
                lN = int(get_parameter('Last_line (w.r.t. original_image)',masterRes,1))
                p0 = int(get_parameter('First_pixel (w.r.t. original_image)',masterRes,1))
                pN = int(get_parameter('Last_pixel (w.r.t. original_image)',masterRes,1))
                nofLines=lN-l0+1
                nofPixels=pN-p0+1
       
                #**************************#
                #    read/operate data     #
                #**************************#

                RAW_CINT_SRD =os.getcwd()+'/SMALL_BASELINES/'+cint_minrefdem_name+'/'+'cint.minrefdem.raw'
                Path_MFF_HDR =os.getcwd()+'/SMALL_BASELINES/'+cint_minrefdem_name+'/'+'cint.minrefdem.hdr'
                Link_CINT_SRD=os.getcwd()+'/SMALL_BASELINES/'+cint_minrefdem_name+'/'+'cint.minrefdem.x00'

                outStream      = open(Path_MFF_HDR,'w')
                outStream.write('IMAGE_FILE_FORMAT = MFF\n')
                outStream.write('FILE_TYPE = IMAGE\n')
                outStream.write('IMAGE_LINES = %d\n' % int(nofLines))
                outStream.write('LINE_SAMPLES = %d\n'% int(nofPixels))
                outStream.write('BYTE_ORDER = LSB\n')
                outStream.write('END\n')
                outStream.close()
    
                if (os.path.exists(Link_CINT_SRD)):
                    os.remove(Link_CINT_SRD)
                RAW_CINT_SRD_ABSOLUTE_PATH=os.path.abspath(RAW_CINT_SRD)
                print "RAW_CINT_SRD_ABSOLUTE_PATH=", RAW_CINT_SRD_ABSOLUTE_PATH
                os.symlink(RAW_CINT_SRD_ABSOLUTE_PATH,Link_CINT_SRD)
                #**************************************************************************************************# 
                #SLC_Data=freadbk(Path_MFF_HDR,int(nofLines1-line2), 1,int(line2),int(nofPixels1))
                SBAS_Int_Data_Stack[file_count,:,:]=freadbk(Path_MFF_HDR,Patch_First_Line+1,\
                                                           Patch_First_Pixel+1,Patch_nofLines,Patch_nofPixels)
                #**************************************************************************************************#
                if (os.path.exists(Path_MFF_HDR)):
                    os.remove(Path_MFF_HDR)
                if (os.path.exists(Link_CINT_SRD)):
                    os.remove(Link_CINT_SRD)
                #keeps track of total nuber of lines
                #*******************************************************#
                file_count=file_count+1
                print 'file_count=',file_count
                print '\n\n\n'
                dataFormat = 'cpxfloat32'

            KS_SBAS_Data=SBAS_Int_Data_Stack.copy()

            psc_da=patch_Name.strip()+'/'+'pscands.1.da'
            psc_ij=patch_Name.strip()+'/'+'pscands.1.ij'
            if os.path.exists(KS_folder)==False:
                print "KS_result is not exit"
            else:
                file_KS_Estimation=KS_folder+'/'+patch_Name.strip()+'_data.dat'
            print file_KS_Estimation
            Line_index_ks=0
            f=open(file_KS_Estimation,'rb')
            if  f < 0:
                print "ERROR :cannot open file for writing"
                sys.exit(1)
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
                
                    #temp_mask=0
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
                        Line_index_ks+=1                         
                        win_mask=win_mask_total[((Line_index_ks-1)*Total_win_pix):(Line_index_ks*Total_win_pix)]          
                        #temp=f.read(Total_win_pix)
                        #win_mask=int8(struct.unpack(unpack_format,temp))
                        sum_filter=win_mask.sum()                      
                        if sum_filter > sum_filter_wins: 
                            SLC_cut=SBAS_Int_Data_Stack[:,Relative_line_ind-Half_win_nofLines:Relative_line_ind+Half_win_nofLines+1,Relative_pixel_ind-Half_win_nofPixels:Relative_pixel_ind+Half_win_nofPixels+1].reshape(file_count,1,-1)
                            #print SLC_cut.shape       
                            #print win_mask
                            #KS_real_pixels=np.sum(win_mask*real(SLC_cut))
                            #KS_imag_pixels=np.sum(win_mask*imag(SLC_cut))
                            #angle_SLC_Data_temp =angle(np.sum(win_mask*SLC_cut))
                            #print '\nangle_SLC_Data_temp=',angle_SLC_Data_temp
                            #SLC_Data_KS[temp_line_ind-1,temp_pixel_ind-1]  =abs(SLC_Data[temp_line_ind-1,temp_pixel_ind-1])*(cos(angle_SLC_Data_temp)+1j*sin(angle_SLC_Data_temp))
                            #print 'SLC_Data_KS=',angle(SLC_Data_KS_origin[temp_line_ind-1,temp_pixel_ind-1])
                            SLC_Data_KS = np.dot(SLC_cut,win_mask)/sum_filter
                            #print 'SLC_Data_KS=',SLC_Data_KS.shape
                            #print SLC_Data_KS
                            KS_SBAS_Data[:,Relative_line_ind,Relative_pixel_ind]=SLC_Data_KS[:,0]
                            if int(string_line_ij[0])%250000==0:
                                print 'processing the patch row:',patch_Name.strip()
                                print 'processing the file name:',cint_minrefdem_name
                                print 'processing the data index:',string_line_ij[0]    #print win_mask
                                print win_mask[1:10]
                           
            f.close()
            file_count=0
            for cint_minrefdem_file in myfile:
                print 'cint_minrefdem_file=',cint_minrefdem_file
                cint_minrefdem_name=cint_minrefdem_file.strip()
                cint_minrefdem_name=cint_minrefdem_name[0:8]+'_'+cint_minrefdem_name[9:17]
                print 'cint_minrefdem_name=',cint_minrefdem_name           
                outfile=os.getcwd()+'/SMALL_BASELINES/'+cint_minrefdem_name+'/'+'cint.minrefpha'+'_'+patch_Name.strip()+'.raw'
                print 'outfile=',outfile
                fid = open(outfile,'wb');
                if fid < 0:
                    print "ERROR :cannot open file for writing"
                    sys.exit(1)
                else: 
                    for i_temp in np.arange(0,Patch_nofLines):
                        fid.write(KS_SBAS_Data[file_count,i_temp,:])
                fid.close()
                file_count=file_count+1

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
    #SLC_file_count=len(open(SLC_file).readlines())
    ##############################################################
    ##Parameter of filter Windows
    #############################################################
   ##############################################################
    ##Parameter of filter Windows
    #############################################################
    #win_mask=np.zeros((Total_win_pix),dtype=np.int8)
    #win_mask='0'*int(Total_win_pix)
    #Patch_count=len(open(Slave_SLC_File,'r').readlines())
    work_dir=os.getcwd();
    SBAS_work_dir=work_dir+'/SMALL_BASELINES'
    print 'SBAS_work_dir=',SBAS_work_dir

    if os.path.exists(SBAS_work_dir):    
        file_dir=os.listdir(SBAS_work_dir)
    else:
        print "This path is wrong !!!"
        sys.exit(1)

    filename_list=[]
    for filename in file_dir:    
        m_index=re.search('\d\d\d\d_\d\d\d\d',filename)
        #M_index=re.search('\d_\d',filename)
        if m_index is not None:
        #print 'filename of Interferometry folder =',filename
            filename_list.append(filename)

    filename_list.sort()

    filename_list=[]
    f=open('patch.list')
    patch_file_list=f.readlines()
    for patch_file in  patch_file_list:   
        print patch_file.strip()
        filename_list.append(patch_file.strip())
    #filename_list.sort()
    Num_patch =len(patch_file_list)
    #print 'Num_SLC =',Num_patch
    #pool.map(homofilter_KS,cint_minrefdem_file)
        #homofilter_KS(patch_Name)
    homofilter_KS()
   
        
    
if __name__== '__main__':
    main()
