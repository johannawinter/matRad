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
% simulations stored as *.ddd. Each section in this file contains one specific step in the 
% base data generation procedure. First, depth dose distributions for
% protons or carbon ions can be extracted (including single and double laterals sigmas and weigths).
% In addition the whole content of the TRiP planning folder can be parsed, plotted and saved as 
% *.mat file for further processing.

%% set global path which is valid for all subsections
clc,
clear 
close all
pathTRiP = 'E:\TRiP98DATA_HIT-20131120';

%% extract carbon and proton depth dose profiles

FocusIdx = 1;
Offset   = -2.89; % in mm
visBool  = 1;
Identifier = 'p'; % either p for protons, C for carbons
% parse and save proton ddd's
baseData = matRad_getDDDfromTxt(Identifier,pathTRiP,FocusIdx,Offset,visBool);
save(['..' filesep 'protonBaseDataHIT'],'baseData');
% parse and save carbon ddd's
Identifier = 'C';
baseData = matRad_getDDDfromTxt(Identifier,pathTRiP,FocusIdx,Offset,visBool);
save(['..' filesep 'carbonBaseDataHIT'],'baseData');


%% interpolate double gaussian data froms sparse sigma1, sigma2 and weight matrix
load(['..' filesep 'protonBaseDataHIT.mat']);
%path to sampling points/Stützstellen provided by Katia P.
pathToSparseData = [pathTRiP '\DDD\p\HIT_2D_DB_p_NEW'];
Identifier = 'p';
% if visBool is on then dont forget to press a key to step to the next plot
baseData = matRad_interpLateralBaseData(baseData,pathTRiP,pathToSparseData,Identifier,1);


%% parse single spc files 
pathToSPC = [pathTRiP filesep 'SPC\12C\RF3MM\FLUKA_NEW3_12C.H2O.MeV09000.spc'];
[meta,SPC] = matRad_readSPC(pathToSPC);

%% parse all spc files in a specific folder
pathToSPCfiles = [pathTRiP filesep 'SPC\12C\RF3MM\'];
dirInfo = dir([pathToSPCfiles '*.spc']);

for i = 1:length(dirInfo)
    [MetaSPC,SPC] = matRad_readSPC([pathToSPCfiles dirInfo(i).name]);
    Filename = dirInfo(i).name;
    save([pathToSPCfiles Filename(1:end-4) '.mat'],'MetaSPC','SPC');
end


 %% parse dEdx file
pathdEdx = [pathTRiP filesep '\DEDX\dEdxFLUKAxTRiP.dedx'];
[Meta, dEdx ] = matRad_readdEdx(pathdEdx);
save('dEdx.mat','dEdx');


 %% parse RBEinitial data
pathRBE = [pathTRiP filesep 'RBE'];
[RBE] = matRad_readRBE(pathRBE);
save('RBEinitial.mat','RBE');


%% get dEdx* alpha from CNAO files or get them from the 37 spc files, for now the rapid 
% Scholz algorithm is implemented

% get dose averaged alpha and beta depth curves from the 37spc files, RBE
% inital file and dEdx file using the rapidScholz algorithm
pathToSPCFiles = [pathTRiP filesep 'SPC\12C\RF3MM'];
[ sDataHIT ] = matRad_combineBioSPC(pathToSPCFiles,0);
save('sDataHIT.mat','sDataHIT');


% get dEdx* alpha and dEdx*sqBeta curves from CNAO files. Please keep in
% mind that these curves needs to be divided by the ddd in order to use
% them in matRad. Takes about 10s per alphaBetaRatio. The fields
% sDataCNAO.alphaX and sDataCNAO.betaX have to be added manually for each
% cellline (requirement by matRad).

% TODO
pathCNAO = '\\Mac\Home\Documents\Heidelberg\CNAO_baseData';
[ sDataCNAO ] = matRad_ParseBioDataCNAO(pathCNAO,'C',0);

%% interpolate each entry in the ddd the corresponding depth alpha and depth beta curve
load('sDataHIT.mat');
%load('sDataCNAO.mat');
load(['..' filesep 'carbonBaseDataHIT.mat']);
CNAOisUsed = false;
if CNAOisUsed
    [baseData] = matRad_interpDoseAvgBioData(baseData, sDataCNAO ,CNAOisUsed, 1);
else
    [baseData] = matRad_interpDoseAvgBioData(baseData, sDataHIT ,CNAOisUsed, 1);
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% old code %%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% each *xls file in the provided folder will be converted to *.mat files to use
% them further on for the bio base data generation.
% clc,
% clear 
% close all
% addpath([pwd filesep 'BioDataGeneration']);
% % provide the folder of the *.xls spectra files and the destination path
% matRad_XlsSpectra2Mat('E:\TRiP98DATA_HIT-20131120\SPC\12C\RF3MM',[pwd filesep 'baseDataHIT'])
% 


