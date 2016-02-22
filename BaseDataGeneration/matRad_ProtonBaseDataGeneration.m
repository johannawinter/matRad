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
% protons can be extracted (including single and double laterals sigmas and weigths).
% In addition, the content of the TRiP planning folder can be parsed,
% plotted and saved. 

%% set global path which is valid for all subsections
clc,
clear 
close all
pathTRiP   = 'E:\TRiP98DATA_HIT-20131120';
matRadRoot = 'C:\Users\wieserh\Documents\matRad';

%% extract carbon and proton depth dose profiles
FocusIdx = 0;
Offset   = -2.89; % in mm
visBool  = 0;
Identifier = 'p';
metaInformation.SAD             = 6509;   %[mm]
metaInformation.BAMStoIsoDist   = 1226;   %[mm]
metaInformation.minIniBeamSigma = 6;      %[mm]
metaInformation.machineName     = 'HIT';
metaInformation.description = ['proton baseData from TRiP98 combined with KatjaP ' ...
                               'lateral double gauss data considering beam widening in air'];

 
% parse and save proton ddd's
machine = matRad_getDDDfromTxt(Identifier,pathTRiP,FocusIdx,Offset,metaInformation,visBool);
Name = [machine.meta.radiationMode  '_'  machine.meta.name];


%% parse beam widening and inital foki size from LPD.xml files stemming from HIT
% an additional field named iniFocus will be added to the base data set
% which holds for each energy, 4 lookup tables which correspond to four
% different foki. Please note that the initial beam width is already
% considered
PathToXMLFile = [matRadRoot filesep 'BaseDataGeneration' filesep 'ProtonLPD.xml'];
machine = matRad_readBeamWidthHIT(machine,PathToXMLFile);


%% interpolate double gaussian data from sparse sigma1, sigma2 and weight matrix
% lateral data from katja only describes scattering within the patient

%path to sampling points/Stützstellen provided by Katia P.
pathToSparseData = [pathTRiP '\DDD\12C\HIT_2D_DB_Cwith_KatjaP'];

% if visBool is on then dont forget to press a key to step to the next plot
machine = matRad_interpLateralBaseData(machine,pathTRiP,pathToSparseData,Identifier,FocusIdx,0);

%% save file
save(Name,'machine');
