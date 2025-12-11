%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% created by wuwenhao
%  20120911
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
diary PS_INSAR_log.txt
getparm 
%%%%%%%%%creat the log file%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define the parmeter
setparm('max_topo_err',35)    
setparm('select_method','percent')
setparm('gamma_stdev_reject',0.02)
setparm('percent_rand',0.00000000000000001)
%setparm('drop_ifg_index',[4])
%setparm('drop_ifg_index',[])
setparm('gamma_change_convergence',0.002)
setparm('scla_deramp','y')
setparm('scn_deramp_ifg','y')
setparm('scn_time_win',66)
%setparm('n_cores',10)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%STEP 7: Estimate spatilly-correalted llok angle error%%%%%%%%%%
%setparm('ref_lat',[30.82 30.84])
%setparm('ref_lon',[110.98  111.03])
%setparm('ref_lat',[22.8 22.805])
%setparm('ref_lon',[108.34 108.35])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time calculate
Time_Step=now;
disp(['*******************Create by Wu Wen hao  20120912********************'])
disp([' Step 1: Load data  The start time is:' datestr(Time_Step,0)])
disp(['*******************Create by Wu Wen hao  20120912********************'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
stamps(1,1)
%quit
%stamps_mc_header(1,1)
disp(['*******************************************************************'])
toc 

disp([' Step 1: Load data is finished '])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time calculate
Time_Step=now;
disp(['*******************Create by Wu Wen hao  20120912********************'])
disp([' Step 2: Estimate Phase  The start time is:' datestr(Time_Step,0)])
disp(['*******************Create by Wu Wen hao  20120912********************'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic 
setparm('gamma_max_iteration',6)
%stamps_mc_header(2,2)
stamps(2,2)
%quit
disp(['*******************************************************************'])
toc 
disp([' Step 2: Estimate Phase  is finished!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!']) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time calculate
Time_Step=now;
disp(['*******************Create by Wu Wen hao  20120912********************'])
disp([' Step 3: PS selection  The start time is:' datestr(Time_Step,0)])
disp(['Step 3 need so much long time ,please have a cup of tea and have a sleep'])
disp(['*******************Create by Wu Wen hao  20120912********************'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic 
%stamps_mc_header(3,3)
stamps(3,3)
%quit
disp(['*******************************************************************'])
toc 
disp([' Step 3: PS selection  is finished!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['*******************Create by Wu Wen hao  20120912********************'])
disp([' Step 4: PS weeding  The start time is:' datestr(Time_Step,0)])
disp(['*******************Create by Wu Wen hao  20120912********************'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
setparm('weed_standard_dev',0.55)
setparm('weed_max_noise',0.55)
setparm('weed_time_win',380)
setparm('weed_neighbours','n')
setparm('weed_zero_elevation','y')
%stamps_mc_header(4,4)
stamps(4,4,'y')
%quit
toc
disp([' Step 4: PS selection  is finished!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['*******************Create by Wu Wen hao  20120912********************'])
disp([' Step 5: Phase weeding  The start time is:' datestr(Time_Step,0)])
disp(['*******************Create by Wu Wen hao  20120912********************'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic 
%stamps_mc_header(5,5)
stamps(5,5)
disp(['*******************************************************************'])
toc 
disp([' Step 5: PS weeding  is finished !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['*******************Create by Wu Wen hao  20120912********************'])
disp([' Step 6: Phase unwraping  The start time is:' datestr(Time_Step,0)])
disp(['*******************Create by Wu Wen hao  20120912********************'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic 
setparm('unwrap_prefilter_flag','y')
setparm('unwrap_grid_size',100)
setparm('unwrap_patch_phase','n')
setparm('unwrap_gold_alpha',0.8)
setparm('unwrap_gold_n_win',64)
setparm('unwrap_time_win',128)
setparm('unwrap_method','3D')
setparm('subtr_tropo','n')
setparm('tropo_method','a_l')
stamps(6,6)
toc 
disp([' Step 6: Phase unwrapping  is finished  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
tic
disp(['**************************!!!!!!!!!!!!!!!!!!!***************************'])
disp(['Then Step 7: Estimate spatially-correlated llok angle error  AND Step 8: Filter spatially correlated nose' ])
stamps(7,7,'y')

stamps(8,8,'y')
disp(['*******************************************************************'])
toc
disp(['Conculate ! All step finished !!!!!!!!!!!!!!!!!!!!!!!!!!!!'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

%#################plot figure################%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%ps_output
%ps_plot('a','a_linear',-1)

ps_plot('vs',-1)
%ps_plot('vs')
%print(gcf,'-djpeg','plot_vs.jpg')
%close 
ps_plot('m',-1)
ps_plot('d',-1)
ps_plot('w',-1)
ps_plot('o',-1)
ps_plot('s',-1)
ps_plot('w-d',-1)
ps_plot('w-o',-1)
ps_plot('w-dm',-1)
ps_plot('w-do',-1)
ps_plot('w-dmo',-1)
ps_plot('u',-1)
ps_plot('u-d',-1)
ps_plot('u-m',-1)
ps_plot('u-dm',-1)
ps_plot('u-dms',-1)
ps_plot('u-dmo',-1)
ps_plot('u-dmos',-1)
ps_plot('v',-1)
ps_plot('v-d',-1)
ps_plot('v-o',-1)
ps_plot('v-do',-1)

quit
%{
ps_plot('v-d',-1)
ps_plot('v-d')
print(gcf,'-djpeg','plot_v-d.jpg')
close
ps_plot('v-o',-1)
ps_plot('v-o')
print(gcf,'-djpeg','plot_v-o.jpg')
close
ps_plot('v-do',-1)
ps_plot('v-do')
print(gcf,'-djpeg','plot_v-do.jpg')
close
ps_plot('vs',-1)
ps_plot('vs')
print(gcf,'-djpeg','plot_vs.jpg')
close 
ps_plot('vs-d',-1)
ps_plot('vs-d')
print(gcf,'-djpeg','plot_vs-d.jpg')
close
%ps_plot('vs-d',4)
%print(gcf,'-djpeg','plot_vs-d-4.jpg')
%ps_plot('vs-d',5)
%print(gcf,'-djpeg','plot_vs-d-5.jpg')
ps_plot('vs-o',-1)
ps_plot('vs-o')
print(gcf,'-djpeg','plot_vs-o.jpg')
close
ps_plot('vs-do',-1)
ps_plot('vs-do')
print(gcf,'-djpeg','plot_vs-do.jpg')
close
%ps_plot('vs-do',4)
%print(gcf,'-djpeg','plot_vs-do-4.jpg')
%ps_plot('vs-do',5)
%print(gcf,'-djpeg','plot_vs-do-5.jpg')
%ps_plot('vdrop')
%print(gcf,'-djpeg','plot_vdrop.jpg')
%ps_plot('vdrop-d')
%print(gcf,'-djpeg','plot_vdrop-d.jpg')
%ps_plot('vdrop-o')
%print(gcf,'-djpeg','plot_vdrop-o.jpg')
%ps_plot('vdrop-do')
%print(gcf,'-djpeg','plot_vdrop-do.jpg')
%ps_plot('usb')
%print(gcf,'-djpeg','plot_usb.jpg')
%ps_plot('usb-d')
%print(gcf,'-djpeg','plot_usb-d.jpg')
%ps_plot('usb-o')
%print(gcf,'-djpeg','plot_usb-o.jpg')
%ps_plot('usb-do')
%print(gcf,'-djpeg','plot_usb-do.jpg')

%ps_plot('rsb')
%print(gcf,'-djpeg','plot_rsb.jpg')
%}
%print(gcf,'-djpeg','plot_v.jpg')
%print -djpeg70  -r600  plot_v.jpg
%close
%ps_plot('v',4)
%print(gcf,'-djpeg','plot_v-4.jpg')
%close


%ps_plot('v-a',-1)
%
%print(gcf,'-djpeg','plot_v-a.jpg')
%close
%{
ps_plot('v-da',-1)
ps_plot('v-da')
print(gcf,'-djpeg','plot_v-da.jpg')
close
ps_plot('v-dao',-1)
ps_plot('v-dao')
print(gcf,'-djpeg','plot_v-dao.jpg')
close
ps_plot('v-ds',-1)
ps_plot('v-ds')
print(gcf,'-djpeg','plot_v-ds.jpg')
close
}%
%ps_plot('v-das',-1)
%ps_plot('v-das')
%print(gcf,'-djpeg','plot_v-das.jpg')
%close
}%
%{
%ps_plot('v-d',1)
%print(gcf,'-djpeg','plot_v-d-1.jpg')
%close
%ps_plot('v-d',2)
%print(gcf,'-djpeg','plot_v-d-2.jpg')
%close
%ps_plot('v-d',3)
%print(gcf,'-djpeg','plot_v-d-3.jpg')
%close
%ps_plot('v-d',4)
%print(gcf,'-djpeg','plot_v-d-4.jpg')
%close
%ps_plot('v-d',5)
%print(gcf,'-djpeg','plot_v-d-5.jpg')
%ps_plot('v-d',5)
%print(gcf,'-djpeg','plot_v-d-5.jpg')
%ps_plot('v-d',6)
%print(gcf,'-djpeg','plot_v-d-6.jpg')
%}
%}
%}
ps_plot('v-d')
print(gcf,'-djpeg','plot_v-d.jpg')
close
ps_plot('u')
print(gcf,'-djpeg','plot_u.jpg')
close
ps_plot('u-dmo')
print(gcf,'-djpeg','plot_u-dmo.jpg')
close
ps_plot('u-dmos')
print(gcf,'-djpeg','plot_u-dmos.jpg')
close
ps_plot('w-d')
print(gcf,'-djpeg','plot_w-d.jpg')
close
ps_plot('v-do')
print(gcf,'-djpeg','plot_v-do.jpg')
close
ps_plot('v-o')
print(gcf,'-djpeg','plot_v-o.jpg')
close
ps_plot('v')
print(gcf,'-djpeg','plot_v.jpg')
close 
ps_plot('o')
print(gcf,'-djpeg','plot_o.jpg')
close 
ps_plot('w-dm')
print(gcf,'-djpeg','plot_w-dm.jpg')
close
%ps_plot('a','a_linear')
%print(gcf,'-djpeg','plot_a_linear.jpg')
%close
ps_plot('w')
print(gcf,'-djpeg','plot_w.jpg')
close 
ps_plot('s')
print(gcf,'-djpeg','plot_s.jpg')
close 
ps_plot('w-o')
print(gcf,'-djpeg','plot_w-o.jpg')
close
ps_plot('u-dms')
print(gcf,'-djpeg','plot_u-dms.jpg')
close
ps_plot('v-d')
print(gcf,'-djpeg','plot_v-d.jpg')
close
ps_plot('w-dmo')
print(gcf,'-djpeg','plot_w-dmo.jpg')
close
ps_plot('u-m')
print(gcf,'-djpeg','plot_u-m.jpg')
close
ps_plot('u-dm')
print(gcf,'-djpeg','plot_u-dm.jpg')
close
ps_plot('m')
print(gcf,'-djpeg','plot_m.jpg')
close
ps_plot('w-do')
print(gcf,'-djpeg','plot_w-do.jpg')
close

ps_plot('d')
print(gcf,'-djpeg','plot_d.jpg')
close 
ps_plot('u-d')
print(gcf,'-djpeg','plot_u-d.jpg')
close


%load ps_plot_v-ds ph_disp 
%ps_gescatter('project_velo.kml',ph_disp,10,0.4)
ps_output
%plot_v.gmt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
diary off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all 
quit
