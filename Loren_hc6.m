%%
close all
clear
clc

%%
% load('C:\Tamir\work\Courses\MBL\project\Datasets\Loren_Frank\hc6\Ten\tencellinfo.mat');
% load('C:\Tamir\work\Courses\MBL\project\Datasets\Loren_Frank\hc6\Ten\tentetinfo.mat')
% load('C:\Tamir\work\Courses\MBL\project\Datasets\Loren_Frank\hc6\Ten\tenpos01.mat')
% load('C:\Tamir\work\Courses\MBL\project\Datasets\Loren_Frank\hc6\Ten\tenspikes01.mat')

%%
main_dir = 'C:\Tamir\work\Courses\MBL\project\Datasets\Loren_Frank\hc6';
animal_name = 'ten';
day = 1;
session = 2;
TT=2;
unit=1;
load(fullfile(main_dir,animal_name,sprintf('%scellinfo',animal_name)));
load(fullfile(main_dir,animal_name,sprintf('%stetinfo',animal_name)));
load(fullfile(main_dir,animal_name,sprintf('%spos%.2d',animal_name,day)));
load(fullfile(main_dir,animal_name,sprintf('%sspikes%.2d',animal_name,day)));
pos_data = pos{day}{session}.data;
spikes_data = spikes{day}{session}{TT}{unit}.data;

%%
figure
hold on
plot(pos_data(:,1),pos_data(:,5),'.')
plot(pos_data(:,1),pos_data(:,6),'.')
% plot(sdf(:,1),sdf(:,7),'.')

%%
figure
hold on
IX = find(pos_data(:,7)==0);
plot(pos_data(IX,2),pos_data(IX,3),'.-')
IX = find(pos_data(:,7)~=0);
plot(pos_data(IX,2),pos_data(IX,3),'.-')

%%
figure
hold on
plot(pos_data(:,2),pos_data(:,3),'.k')
plot(spikes_data(:,2),spikes_data(:,3),'.r')


%%










%%