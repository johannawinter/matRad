% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad_BaseDataGeneration script
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
% This script can be used to generate machine files from monte carlo
% simulations stored as *.ddd. Each section in this file contains one specific step in the 
% base data generation procedure. First, depth dose distributions for
% protons or carbon ions can be extracted (including single and double laterals sigmas and weigths).
% In addition the content of the TRiP planning folder can be parsed,
% plotted and saved. Biological base data can be generated from RBEinital
% tables

%% set global path which is valid for all subsections
clc,
clear 
close all
pathTRiP = 'E:\TRiP98DATA_HIT-20131120';

%% extract carbon and proton depth dose profiles

FocusIdx = 0;
Offset   = -2.89; % in mm
visBool  = 0;
Identifier = 'C';

metaInformation.SAD             = 6509;   %[mm]
metaInformation.BAMStoIsoDist   = 1226;   %[mm]
metaInformation.minIniBeamSigma = 6;      %[mm]
metaInformation.machineName     = 'HIT';
metaInformation.description     = ['carbon baseData from TRiP98 combined with KatjaP ' ...
                               'lateral double gauss data considering beam widening in air'];

% parse and save carbon ddd's
machine = matRad_getDDDfromTxt(Identifier,pathTRiP,FocusIdx,Offset,metaInformation,visBool);
Name = [machine.meta.radiationMode  '_'  machine.meta.name];
 
%% parse beam widening and inital foki size from LPD.xml files stemming from HIT
% an additional field named iniFocus will be added to the base data set
% which holds for each energy, 4 lookup tables which correspond to four
% different foki. Please note that the initial beam width is already
% considered

load('carbon_HIT.mat');
PathToXMLFile = ['C:\Users\wieserh\Documents\matRad\BaseDataGeneration' filesep 'CarbonLPD_Rifi3mm.xml'];
machine = matRad_readBeamWidthHIT(machine,PathToXMLFile);

%% interpolate double gaussian data from sparse sigma1, sigma2 and weight matrix
% lateral data from katja only describes scattering within the patient

%path to sampling points/Stützstellen provided by Katia P.
pathToSparseData = [pathTRiP '\DDD\12C\HIT_2D_DB_Cwith_KatjaP'];
% if visBool is on then dont forget to press a key to step to the next plot
machine = matRad_interpLateralBaseData(machine,pathTRiP,pathToSparseData,Identifier,FocusIdx,0);


%% parse all spc files in a specific folder
pathToSPCfiles = [pathTRiP filesep 'SPC\12C\RF3MM\'];
dirInfo = dir([pathToSPCfiles '*.spc']);

for i = 1:length(dirInfo)
    [MetaSPC,SPC] = matRad_readSPC([pathToSPCfiles dirInfo(i).name]);
    Filename = dirInfo(i).name;
    save([pathToSPCfiles Filename(1:end-4) '.mat'],'MetaSPC','SPC');
end

 %% parse dEdx file
[Meta, dEdx ] = matRad_readdEdx(pathTRiP);
save('dEdx.mat','dEdx');

 %% parse RBEinitial data
[RBE] = matRad_readRBE(pathTRiP);
save('RBEinitial.mat','RBE');

%% get dEdx* alpha from CNAO files or get them from the 37 spc files,
% for now the rapidScholz algorithm is implemented

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

% % TODO
% pathCNAO = '\\Mac\Home\Documents\Heidelberg\CNAO_baseData';
% [ sDataCNAO ] = matRad_ParseBioDataCNAO(pathCNAO,'C',0);

%% interpolate each entry in the ddd the corresponding depth alpha and depth beta curve

CNAOisUsed = false;
if CNAOisUsed
    [machine] = matRad_interpDoseAvgBioData(machine, sDataCNAO ,CNAOisUsed, 0);
else
    [machine] = matRad_interpDoseAvgBioData(machine, sDataHIT ,CNAOisUsed, 0);
end


%% save file
save(Name,'machine');

