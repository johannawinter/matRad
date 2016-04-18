% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad_ProtonBaseDataGeneration script
%
% to run the full script you need:
%   - TRIP treatment planning folder (including DDD and SIS folder)
%   - double gauss data stores as *xml (provided by. Katja.P
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2016, Hans-Peter Wieser
%
% h.wieser@dkfz-heidelberg.de
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
% In addition, the content of the TRiP planning folder can be parsed (e.g.
% sis files).

clc
clear 
close all

%% set global path which is valid for all subsections

pathTRiP   = 'E:\TRiP98DATA_HIT-20131120';
matRadRoot = 'C:\Users\wieserh\Documents\matRad';

%% define some meta information on how the base data set should be generated
SisFocusIdx                 = 0;      % 0-> use no initalBeam, 1...7 determines the FocusIndex which should
                                      % be added to the scattering caused by the patient
BeamOffset                  = -2.89;  %[mm]
visBool                     = 0;
saveToDiskBool              = 0;
Identifier                  = 'p';    % p indicates protons
metaInfo.SAD                = 6509;   %[mm] distance from source to iso-center
metaInfo.BAMStoIsoDist      = 1126;   %[mm] distance beam nozzle to iso-center
metaInfo.LUT_bxWidthminFWHM = [1 2 3 4 6 8 10 20;8 8 8 8 8 8 8 8];      %[mm] 
metaInfo.machine            = 'HIT';
metaInfo.description        = ['proton baseData from TRiP98 combined with KatjaP ' ...
                               'lateral double gauss data considering beam widening in air'];
                           
 
%% parse and save proton ddd's
machine     = matRad_getDDDfromFile(Identifier,pathTRiP,SisFocusIdx,BeamOffset,metaInfo,visBool);
fileName    = [machine.meta.radiationMode  '_'  machine.meta.machine];

%% parse beam widening in air and from LPD.xml file  which are stemming from HIT
% the field machine.data(:).iniFocus will hold for each energy, lookups tables 
% which correspond to different foci. Make sure that the SisFocusIdx is set to 0
% Note that data is only available for foci index 1-4
PathToXMLFile = [matRadRoot filesep 'BaseDataGeneration' filesep 'ProtonLPD.xml'];
machine       = matRad_readBeamWideningAIR(machine,PathToXMLFile);

%% interpolate double gaussian data from sparse sigma1, sigma2 and weight matrix
% lateral data from katja only describes scattering within the patient

%path to sampling points/Stützstellen provided by Katia P.
% if visBool is on then dont forget to press a key to step to the next plot
pathToSparseData = [pathTRiP filesep 'DDD\p\HIT_2D_DB_p_KatjaP'];
machine = matRad_interpLateralBaseData(machine,pathTRiP,pathToSparseData,Identifier,SisFocusIdx,visBool);


%% parse spectra files for protons to enable LET calculations
machine = matRad_getLETfromSPC(machine,Identifier,pathTRiP);

%% save file
save(fileName,'machine');

