#!/usr/bin/env python
import os,sys,time
import math
import numpy as np
from numpy import *
import numpy
#import matplotlib.pyplot as plt

import gdal,gdalconst
from gdalconst import *
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


#from osgeo import gdal,gdalconst
# NOT FETOOLS
#Version 20190824 updaqte the Midpoint 


#******************************************************************************#
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
#*******************************************************************************#
#******************************************************************************
def Linux_cmd(python_instruction):
    #python_instruction='source '+work_dir+'/'+Configfile
    print 'python_instruction=',python_instruction
    with os.popen(python_instruction,"r") as p:
        cmd_profile_copy=p.read()
    #print 'cmd_profile_copy=',cmd_profile_copy       
    return cmd_profile_copy

#******************************************************************************
#_______________________________________________________________________________________
def hms2sec(hmsString,convertFlag='float'):
    # input hmsString syntax: XX:XX:XX.xxxxxx
    print 'hmsString=',hmsString
    secString = int(hmsString[0:2])*3600 + int(hmsString[3:5])*60 + float(hmsString[6:])
    if convertFlag == 'int' :
        return int(secString)
    elif convertFlag == 'float' :
        return float(secString)
    else:
        return int(secString)
#_______________________________________________________________________________________ 
def lph2xyz(line,pixel,norm_orbit_line,resfile,height):
    #$$ initialization
    MAXITER=19900
    CRITERPOS=1.0e-16
    #WGS84 Elliposid:
    ellipsoid=[6378137.0 , 6356752.3141]
    num_points = np.array([line]).shape[0]
    xyz       = np.zeros((num_points,3))
    #$$$ parameter of the image
    tr1=float(get_parameter('Range_time_to_first_pixel (2way) (ms)',resfile,1))/2000
    #tr1  = float(container['rangeTimePix'][0])/2 # one way range time [sec]
    RSR  = float(get_parameter('Range_sampling_rate (computed, MHz)',resfile,1))*2*1000000 # one way in [HZ]
    #print RSR
    centerphi   =float(get_parameter('Scene_centre_latitude',resfile,1))
    centerlambda=float(get_parameter('Scene_centre_longitude',resfile,1))

    #centerphi   =4.102289556600350e+01
    #centerlambda=1.177526113418940e+02
    #print 'centerphi   =',centerphi
    #print 'centerlambda=',centerlambda
    #centerphi = float(centroid_lat)
    #centerlambda= float(centroid_lon)
    SOL=299792458
    # $$$ reference surace: WGS84
    
    ell_a  = ellipsoid[0]                  # semimajor of the ellipsoid
    ell_b  = ellipsoid[1]                  # semiminor of the ellipsoid
    ell_e2 = (ell_a**2-ell_b**2)/ell_a**2    # squared first eccentricity(derived)
    # $$$ ell_e2b=(ell_a^2-ell_b^2)/ell_b^2;  % squared second eccentricity(derived)

    # $$$ [lat,long,h] of scene center to [x,y,z]
    h            = height        # this is only for initial values
    centerphi    = centerphi*numpy.pi/180 
    centerlambda = centerlambda*numpy.pi/180
 
    Ncenter      = ell_a/np.sqrt(1-ell_e2*(numpy.sin(centerphi)**2)); # eccentricity
    scenecenterx = (Ncenter+h)*numpy.cos(centerphi)*numpy.cos(centerlambda);
    scenecentery = (Ncenter+h)*numpy.cos(centerphi)*numpy.sin(centerlambda);
    scenecenterz = (Ncenter+h-ell_e2*Ncenter)*numpy.sin(centerphi);
    #print 'scenecenterx =',scenecenterx
    #print 'scenecentery =',scenecentery
    #print 'scenecenterz =',scenecenterz
    for n in range(0,num_points):   # loop through points

        posonellx = scenecenterx
        posonelly = scenecentery
        posonellz = scenecenterz
        ratime = tr1 + (pixel-1.0)/RSR
      
    #get position and velocity of the satellite
        possatx = norm_orbit_line[n,1]
        possaty = norm_orbit_line[n,2]
        possatz = norm_orbit_line[n,3]
        velsatx = norm_orbit_line[n,4]
        velsaty = norm_orbit_line[n,5]
        velsatz = norm_orbit_line[n,6]
        
        equationset=np.zeros((3,1))
        partialsxyz=np.zeros((3,3))
        for iter in range(1, MAXITER+1):
        
            #update equations and slove system
            dsat_Px = posonellx - possatx    #vector of 'satellite to P on ellipsoid'
            dsat_Py = posonelly - possaty
            dsat_Pz = posonellz - possatz
        
            equationset[0,0] = -(velsatx*dsat_Px+velsaty*dsat_Py+velsatz* dsat_Pz)
        
            equationset[1,0] = -(dsat_Px*dsat_Px+dsat_Py*dsat_Py+dsat_Pz*dsat_Pz-(SOL*ratime)**2)
        
            equationset[2,0] = -((posonellx*posonellx+posonelly*posonelly)/((ell_a+height)**2)+(posonellz/(ell_b+height))**2-1.0)

            partialsxyz[0,0] = velsatx
            partialsxyz[0,1] = velsaty
            partialsxyz[0,2] = velsatz
            partialsxyz[1,0] = 2*dsat_Px
            partialsxyz[1,1] = 2*dsat_Py
            partialsxyz[1,2] = 2*dsat_Pz
            partialsxyz[2,0] = (2*posonellx)/((ell_a+height)**2)
            partialsxyz[2,1] = (2*posonelly)/((ell_a+height)**2)
            partialsxyz[2,2] = (2*posonellz)/((ell_b+height)**2)
        
        # solve system [NOTE] orbit_normalized, otherwise close to
        # singular
            solpos = numpy.linalg.solve(partialsxyz,equationset)
       
            solx = solpos[0,0]
            soly = solpos[1,0]
            solz = solpos[2,0]
           
        # update solution
            posonellx = posonellx + solx;
            posonelly = posonelly + soly;
            posonellz = posonellz + solz;
        
        # check convergence
            if (abs(solx)<CRITERPOS and abs(soly)<CRITERPOS and abs(solz)<CRITERPOS):
                break
            elif(iter>=MAXITER):
                MAXITER=MAXITER+1    
    # final solution: array of XYZ coordinates
        #print 'solx=',solx
        # print 'soly=',soly
        xyz[n,:]=np.array([posonellx, posonelly, posonellz]).copy()
         
        return xyz
#_____________________________________________________________________________________
def intrp_orbit(resfile,line):

    #intrpOrder = 'spline'
    orbit_number = int(get_parameter('NUMBER_OF_DATAPOINTS',resfile,1))
    orbit_time=np.zeros(orbit_number,dtype=float64)
    orbit_x   =np.zeros(orbit_number,dtype=float64)
    orbit_y   =np.zeros(orbit_number,dtype=float64)
    orbit_z   =np.zeros(orbit_number,dtype=float64) 
    orbit_info = get_parameter('NUMBER_OF_DATAPOINTS',resfile,4)
    orbit_info=orbit_info.split('\n')

    for row in range(orbit_number):

        orbit_time_position=orbit_info[row]
        orbit_time[row]=float64(orbit_time_position.strip().split()[0])
        orbit_x[row]=float64(orbit_time_position.strip().split()[1])
        orbit_y[row]=float64(orbit_time_position.strip().split()[2])
        orbit_z[row]=float64(orbit_time_position.strip().split()[3])
        #print orbit_x[row]
    # compute normalization factors
    px = orbit_time # time
    
    f  = min(px)
    g  = (max(px)-min(px))
    px = (px-f)/g
    # polyDegree
    polyDegree = 2

    coef_x1 = (np.polyfit(px,orbit_x,polyDegree));
    a = coef_x1[2]
    b = coef_x1[1]
    c = coef_x1[0]
    coef_x = [c/(g**2), b/g-(2*c*f)/(g**2), a-b*f/g+c*(f/g)**2]

    coef_y1 = (np.polyfit(px,orbit_y,polyDegree))
    a = coef_y1[2]
    b = coef_y1[1]
    c = coef_y1[0]
    coef_y = [c/(g**2), b/g-(2*c*f)/(g**2), a-b*f/g+c*(f/g)**2]

    coef_z1 = (np.polyfit(px,orbit_z,polyDegree));
    a = coef_z1[2]
    b = coef_z1[1]
    c = coef_z1[0]
    coef_z = [c/(g**2), b/g-(2*c*f)/(g**2), a-b*f/g+c*(f/g)**2]

    vel_x = numpy.polyval(numpy.polyder(coef_x),orbit_time)
    vel_y = numpy.polyval(numpy.polyder(coef_y),orbit_time)
    vel_z = numpy.polyval(numpy.polyder(coef_z),orbit_time)
    
    acc_x = numpy.kron(numpy.ones(orbit_number),numpy.polyder(numpy.polyder(coef_x)))
    acc_y = numpy.kron(numpy.ones(orbit_number),numpy.polyder(numpy.polyder(coef_y)))
    acc_z = numpy.kron(numpy.ones(orbit_number),numpy.polyder(numpy.polyder(coef_z)))
    #print 'acc_x.shape=',acc_x.shape

    # interpolated orbit
    norm_orbit = np.array([orbit_time, orbit_x,orbit_y,orbit_z,vel_x,  vel_y,  vel_z,acc_x,  acc_y,  acc_z])
    #print get_parameter('First_pixel_azimuth_time (UTC)',resfile,1).split(' ')[1]
    Taz_start    = hms2sec(get_parameter('First_pixel_azimuth_time (UTC)',resfile,1).split(' ')[1])
    #print 'Taz_start=',Taz_start
    ta1          = Taz_start  # start time in UTC [sec]
    PRF       = float(get_parameter('Pulse_Repetition_Frequency (computed, Hz)',resfile,1))   #[Hz]     # Pulse repeition frequency
    l_aztime     = (line-1)/PRF + ta1 
    #print 'PRF=',PRF
    pos_orbit_x = numpy.interp(l_aztime,orbit_time,orbit_x)
    pos_orbit_y = numpy.interp(l_aztime,orbit_time,orbit_y)
    pos_orbit_z = numpy.interp(l_aztime,orbit_time,orbit_z)


    vel_orbit_x = numpy.interp(l_aztime,orbit_time,vel_x)
    vel_orbit_y = numpy.interp(l_aztime,orbit_time,vel_y)
    vel_orbit_z = numpy.interp(l_aztime,orbit_time,vel_z)

    #acc = numpy.interp(orbit_time,[acc_x, acc_y, acc_z],    l_aztime,intrpOrder)
    acc_orbit_x = numpy.interp(l_aztime,orbit_time,acc_x)
    acc_orbit_y = numpy.interp(l_aztime,orbit_time,acc_y)
    acc_orbit_z = numpy.interp(l_aztime,orbit_time,acc_z)
    
    
    norm_orbit_line = numpy.array([l_aztime, pos_orbit_x,pos_orbit_y,pos_orbit_z,vel_orbit_x,vel_orbit_y,vel_orbit_z,acc_orbit_x,acc_orbit_y,acc_orbit_z])
    
    return norm_orbit.transpose(),norm_orbit_line.reshape(1,-1,order='F') 
  
####################################################################################################

def xyz2ell(position):

    
    ellipsoid=[6378137.0 , 6356752.3141]
    ell_a = ellipsoid[0]
    ell_b = ellipsoid[1]

    ell_e2  = (ell_a**2-ell_b**2)/ell_a**2;    #squared first eccentricity(derived)
    ell_e2b = (ell_a**2-ell_b**2)/ell_b**2;    # squared second eccentricity(derived)


    posx = position[:,0]
    posy = position[:,1]
    posz = position[:,2]
    print 'posx=',posx
    r    = math.sqrt(posx**2 + posy**2)
    mu   = math.atan2(posz*ell_a, r*ell_b)

    sin3 = (math.sin(mu))**3;
    cos3 = (math.cos(mu))**3;
    phi  = math.atan2((posz + ell_e2b * ell_b * sin3),(r - ell_e2 * ell_a* cos3));

    Radar_lambda = math.atan2(posy,posx)
    N = ell_a / math.sqrt(1 - ell_e2 * (math.sin(phi))**2) #for every point no
                                                # approx with scene.center
    height = (r/math.cos(phi)) - N
    
    phi_lam_height = np.zeros(3)
    phi_lam_height[0] = phi*180/math.pi
    phi_lam_height[1] = Radar_lambda*180/math.pi 
    phi_lam_height[2] = height
    return phi_lam_height

#-------------------------------------------------------------------------------------#
################################################################################   
###############################################################################
#thisBurstData = freadbk(['burst' num2str(nBurst)   '/cint_srd.raw'],nofLines1,formatData1, line1, nofLines1,1,nofPixels1);
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
def SAR_data_read(SLC_path_file,data_format,nofLines,nofPixels,Patch_First_Line,Patch_First_Pixel,Patch_nofLines,Patch_nofPixels):
    
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
    
    return Cint_Data

#############################################################################################
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
def InSAR_combination(input_fun_cmd):
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

    resFilename=work_dir+'/'+'master.res'
    l0 = int(get_parameter('First_line (w.r.t. original_image)',resFilename,2,'*_Start_crop','* End_crop:_NORMAL'))
    lN = int(get_parameter('Last_line (w.r.t. original_image)',resFilename,2,'*_Start_crop','* End_crop:_NORMAL'))
    p0 = int(get_parameter('First_pixel (w.r.t. original_image)',resFilename,2,'*_Start_crop','* End_crop:_NORMAL'))
    pN = int(get_parameter('Last_pixel (w.r.t. original_image)',resFilename,2,'*_Start_crop','* End_crop:_NORMAL'))

    #% Get resampled Slv size
    Naz_res = lN-l0+1
    Nrg_res = pN-p0+1



    SBAS_work_dir=work_dir+'/SMALL_BASELINES'
    
    start_time = datetime.datetime.now()
    if os.path.exists(SBAS_work_dir):    
        print 'SBAS_work_dir=',SBAS_work_dir    
    else:
        print "This path is wrong !!!"
        sys.exit(1)
    master_slave_combination=input_fun_cmd[2]
    print 'master_slave_combination=',master_slave_combination
    #SBAS_master_slave_list=[]
    data_format='complex_real4'

    for temp_master_slave in master_slave_combination:
        print 'temp_master_slave=',temp_master_slave
        master_time=temp_master_slave.split('_')[0]
        slave_time=temp_master_slave.split('_')[1]

        SBAS_file_fold=SBAS_work_dir+'/'+temp_master_slave.strip()
        print 'SBAS_file_fold=',SBAS_file_fold
        print os.path.exists(SBAS_file_fold)
        if os.path.exists(SBAS_file_fold):
            #os.chdir(SBAS_file_fold)
            print 'SBAS_file_fold=',SBAS_file_fold
            print '\n'
            original_source_input=work_dir+'/'+'mini_stack_compress'+'/'+master_time.strip()+'/'+'resample_ministack.raw'
            
            if os.path.isfile(original_source_input):
                master_SLC_data=SAR_data_read(original_source_input,data_format,Naz_res,Nrg_res,1,\
                              1,Naz_res,Nrg_res)               
            else:
                print "can not find the fold !!!"
                sys.exit(1)

            original_source_input=work_dir+'/'+'mini_stack_compress'+'/'+slave_time.strip()+'/'+'resample_ministack.raw'
            if os.path.isfile(original_source_input):
                slave_SLC_data=SAR_data_read(original_source_input,data_format,Naz_res,Nrg_res,1,\
                              1,Naz_res,Nrg_res)               
            else:
                print "can not find the fold !!!"
                sys.exit(1)

            Cint_data=master_SLC_data*np.conj(slave_SLC_data)
            Output_file=SBAS_work_dir+'/'+temp_master_slave.strip()+'/'+'cint.mini_track_time_compress.raw'
            #Output_file=SBAS_work_dir+'/'+temp_master_slave.strip()+'/'+'cint.minrefdem.raw' 
            print 'Output_file=',Output_file
            outformat= 'complex_real4'
            print 'outformat=',outformat
            SAR_write(Output_file,Cint_data,outformat)
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
    original_master_time=0
    master_data=np.zeros((Total_nofLines,Total_nofPixels),dtype=np.complex64)
    for temp_master_slave in master_slave_combination:
        print 'temp_master_slave=',temp_master_slave
        master_time=temp_master_slave.split('_')[0].strip()
        slave_time=temp_master_slave.split('_')[1].strip()
        print 'master_time=',master_time
        print 'slave_time=',slave_time
        #original_slave_time=
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
                if int(original_master_time)==int(master_time.strip()):
                    print 'original_master_time=',original_master_time 

                else:
                    SLC_path_file=work_dir+'/'+master_time+'/'+'cint.minrefdem.raw'
                    original_master_time=master_time.strip()
        
                    data_format='complex_real4'
                    master_data=SAR_data_read(SLC_path_file,data_format,\
                        Total_nofLines,Total_nofPixels,1,1,Total_nofLines,Total_nofPixels)
                
                SLC_path_file=work_dir+'/'+slave_time+'/'+'cint.minrefdem.raw'
                data_format='complex_real4'
                slave_data=SAR_data_read(SLC_path_file,data_format,\
                        Total_nofLines,Total_nofPixels,1,1,Total_nofLines,Total_nofPixels)
                cint_srd_data=slave_data*np.conj(master_data)
 
            Output_file    =SBAS_work_dir+'/'+temp_master_slave+'/'+'cint.minrefdem.raw' 
            outformat= 'complex_real4'
            #print 'outformat=',outformat 
            print 'Output_file=',Output_file 
            print '\n\n\n'    
            SAR_write(Output_file,cint_srd_data,outformat)
#******************************************************************************           

################################################################################################
def main():
     
    cores=multiprocessing.cpu_count()
    #print 'Total cores is ', cores
    pool=multiprocessing.Pool(processes=(int(cores/2)))

    work_dir=os.getcwd()
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
    #print 'filename_list=',filename_list
    #print 'large_mini_stack=',large_mini_stack

    large_mini_stack =list(set(large_mini_stack))
    large_mini_stack=sorted(large_mini_stack)
    #print 'large_mini_stack=',large_mini_stack

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
        #print 'ks_stack_num=',ks_stack_num
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
        #print 'mini_stack_index_list=',mini_stack_index_list 
        
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
                 if Relative_master_mini_stack[0] in mini_stack_file_name_list:
                     Relative_MASTER_DATE_INDEX=mini_stack_file_name_list.index(Relative_master_mini_stack[0])
                 else:
                     Relative_MASTER_DATE_INDEX=0 #mini_stack_file_name_list.index(Relative_master_mini_stack[0])
            print 'Relative_MASTER_DATE_INDEX=',Relative_MASTER_DATE_INDEX
            #print 'Relative_MASTER_DATE_INDEX=',Relative_MASTER_DATE_INDEX
            #print 'Relative_master_mini_stack=',mini_stack_file_name_list[Relative_MASTER_DATE_INDEX]
            #print mini_stack_file_name_list
            master_mini_stack_stack.append(int(mini_stack_file_name_list[Relative_MASTER_DATE_INDEX]))
        master_mini_stack_stack =list(set(master_mini_stack_stack))
        master_mini_stack_stack=sorted(master_mini_stack_stack)
        
        #print  'master_mini_stack_stack=',master_mini_stack_stack
        master_mini_stack_list=[]
        for i_temp_master in master_mini_stack_stack:
            master_mini_stack_list.append(str(i_temp_master))
        #print 'master_mini_stack_list=',master_mini_stack_list

        target_run_path= work_dir+ '/'+ks_stack_num +'/'+'INSAR_'+insar_file_name.split('\t')[2].strip()
        #print 'target_run_path=',target_run_path

       
        mini_stack_combine_list=[]
        for master_i_temp in master_mini_stack_stack:
            for master_j_temp in master_mini_stack_stack:
                if master_i_temp < master_j_temp:
                    print master_i_temp,master_j_temp
                    line_output=str(master_i_temp)+'_'+str(master_j_temp)     
                    mini_stack_combine_list.append(line_output.strip())                            
                    #file_mini_stack.write('\n')
        print 'mini_stack_combine_list=',mini_stack_combine_list

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
            
             os.chdir(target_run_path)

             original_source_input = work_dir+'/'+'concatenate_mini_stack_compress.py'
             target_source_res     = target_run_path+'/'+'concatenate_mini_stack_compress.py'
             if os.path.isfile(original_source_input):
                 shutil.copy(original_source_input,target_source_res) 

                 python_instruction='python concatenate_mini_stack_compress.py'
                 #cmd_profile_copy=Linux_cmd(python_instruction)
                 #print 'cmd_profile_copy=',cmd_profile_copy
             os.chdir(work_dir)


             input_cmd.append(work_dir)
             print 'input_cmd=',input_cmd
            
             #print 'input_cmd_list=',input_cmd_list
             InSAR_combination(input_cmd)

             os.chdir(target_run_path)

             python_instruction='shp_filter_relative_master_patch.py'
             cmd_profile_copy=Linux_cmd(python_instruction)
             print 'cmd_profile_copy=',cmd_profile_copy

             python_instruction='shp_filter_patch.py'
             #cmd_profile_copy=Linux_cmd(python_instruction)
             #print 'cmd_profile_copy=',cmd_profile_copy

             os.chdir(work_dir)
        
if (__name__ == "__main__"):
    main()          

              
