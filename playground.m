clc,
clear 
close all
%% determine which foci shall be used
carbonBaseDataHIT=matRad_getDDDfromTxt('p','E:\TRiP98DATA_HIT-20131120\DDD\p\RF0MM');

%% get dEdx* alpha from CNAO files or get them from the 37 spc files
pathCNAO = 'C:\Users\wieserh\Documents\matRad\CNAO_baseData';
pathHIT = '\\psf\Home\Documents\Heidelberg\matRad\baseDataHIT2'; 

[ sDataHIT ] = matRad_ParseBioDataHIT(pathHIT,0);
%[ sDataCNAO ] = matRad_ParseBioDataCNAO(pathCNAO,'C',1);

%% merge ddd and dose averaged depth alpha curves

 load('C:\Users\wieserh\Documents\matRad\baseDataHIT2\sDataHIT.mat');
 load('C:\Users\wieserh\Documents\matRad\baseDataHIT2\carbonBaseDataHIT.mat');

[ baseData ] = matRad_interpDoseAvgBioData(carbonBaseDataHIT, sDataHIT ,1);

%% interpolate to deeper depths

