
clear;



%% Working file
%{
Path='C:\Users\bridg\Desktop\Work New\TMS ToolBox Development\Working_Folder\Test Data\Example Data\Signal\';
cfs2mat.Convert('E021_T1_RC_IAD_FDI000.cfs', Path); %convert .cfs file to .mat file, saves to folder
%}


%% not working due to different number of samples
%
Path='C:\Users\bridg\Desktop\Work New\TMS ToolBox Development\Working_Folder\Test Data\Example Data\Signal\';
cfs2mat.Convert('MIRANDA_SICI_1mV.cfs', Path); %convert .cfs file to .mat file, saves to folder

%}