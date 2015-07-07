% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad_BioDataGeneration script
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015, Hans-Peter Wieser
%
% h.wieser@dkfz.de
%
% This file is NOT part of the official matRad release. 
% This file has to be used only for internal purposes! 
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script can be used to generate baseData sets from files stemming from monte carlo
% simulation . Each section in this file contains one specific step in the 
% base data generation procedure. First, depth dose distributions for
% protons or carbon ions can be extracted (including laterals sigmas and weigths).
% Once spc files are converted to *.xls files the function XlsSpectra2Mat
% can used to parse the spectra data in mat files.

%% extract carbon and proton depth dose profiles. Please note that the first foci/width will be extracted from the sis file
clc,
clear 
close all
%e.g. rename protonBaseDataHIT into baseData and save it as protonBaseData
%in the matRad root directory in order to use within matRad - analog for
%carbons
addpath([pwd filesep 'BioDataGeneration']);
protonBaseDataHIT=matRad_getDDDfromTxt('p','E:\TRiP98DATA_HIT-20131120',0);
carbonBaseDataHIT=matRad_getDDDfromTxt('C','E:\TRiP98DATA_HIT-20131120',0);


%% each *xls file in the provided folder will be converted to *.mat files to use
% them further on for the bio base data generation.

clc,
clear 
close all
addpath([pwd filesep 'BioDataGeneration']);
% provide the folder of the *.xls spectra files and the destination path
matRad_XlsSpectra2Mat('E:\TRiP98DATA_HIT-20131120\SPC\12C\RF3MM',[pwd filesep 'baseDataHIT'])


%% get dEdx* alpha from CNAO files or get them from the 37 spc files
pathCNAO = 'C:\Users\wieserh\Documents\matRad\CNAO_baseData';
pathHIT = 'C:\Users\wieserh\Documents\matRad\baseDataHIT'; 

% get dose averaged alpha and beta depth curves from the 37spc files, RBE
% inital file and dEdx file using the rapidScholz algorithm
[ sDataHIT ] = matRad_ParseBioDataHIT(pathHIT,0);

% get dEdx* alpha and dEdx*sqBeta curves from CNAO files. Please keep in
% mind that these curves needs to be divided by the ddd in order to use
% them in matRad. Takes about 10s per alphaBetaRatio. The fields
% sDataCNAO.alphaX and sDataCNAO.betaX have to be added manually for each
% cellline (requirement by matRad).
[ sDataCNAO ] = matRad_ParseBioDataCNAO(pathCNAO,'C',0);

%% merge ddd and dose averaged depth alpha curves

 load('C:\Users\wieserh\Documents\matRad\sDataHIT.mat');
 load('C:\Users\wieserh\Documents\matRad\carbonBaseDataHIT.mat');
 CNAOisUsed = true;
[ baseData ] = matRad_interpDoseAvgBioData(baseData, sDataCNAO ,CNAOisUsed,1);

%% interpolate to deeper depths
clc
clear
close all

load('\\psf\Home\Documents\Heidelberg\matRad\carbonBaseDataHITBio.mat');
[ baseData ] = extrapDeeper( baseData,0 );



