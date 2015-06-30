clc,
clear 
close all

pathCNAO = 'C:\Users\wieserh\Documents\matRad\CNAO_baseData';
pathHIT = 'C:\Users\wieserh\Documents\matRad\baseDataHIT2'; 

%[ sDataHIT ] = matRad_ParseBioDataHIT(pathHIT,0);
%[ sDataCNAO ] = matRad_ParseBioDataCNAO(pathCNAO,'C',1);

%% 
load('\\psf\Home\Documents\Heidelberg\matRad\baseDataHIT2\sDataHIT.mat');
load('\\psf\Home\Documents\Heidelberg\matRad\baseDataHIT2\carbonBaseDataHIT.mat');
[ carbonBaseDataHITBIO ] = matRad_interpDoseAvgBioData(carbonBaseDataHIT, sDataHIT ,0);