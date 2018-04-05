% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad_CarbonBaseDataGeneration script
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

pathTRiP   = 'C:\Users\wieserh\Documents\TRiP98DATA_HIT-20131120';
matRadRoot = 'C:\Users\wieserh\Documents\matlab\matRad';

%% define some meta information on how the base data set should be generated
SisFocusIdx                 = 0;      % 0-> use no initalBeam, 1...7 determines the FocusIndex which should
                                      % be added to the scattering caused by the patient
BeamOffset                  = -2.89;  %[mm]
visBool                     = 0;
saveToDiskBool              = 0;
Identifier                  = 'C';    % p indicates protons
metaInfo.SAD                = 43416.5;   %[mm] distance from source to iso-center  fixedBeamLine 6850  Gantry 43416.5
metaInfo.BAMStoIsoDist      = 993.9;   %[mm] distance beam nozzle to iso-center  fixedBeamLine 1126  Gantry 993.9
metaInfo.LUT_bxWidthminFWHM = [1 2 3 4 5 8 10 20; 6 6 6 6 6 10 10 10];      %[mm] 
metaInfo.machine            = 'HIT';
metaInfo.description        = ['carbon baseData from TRiP98 combined with KatjaP ' ...
                               'lateral double gauss data considering beam widening in air, updated beam widening in air considered'];

%% parse and save carbon ddd's
machine     = matRad_getDDDfromFile(Identifier,pathTRiP,SisFocusIdx,BeamOffset,metaInfo,visBool);
fileName    = [machine.meta.radiationMode  '_'  machine.meta.machine];
 
%% parse beam widening in air and from LPD.xml file  which are stemming from HIT
% the field machine.data(:).iniFocus will hold for each energy, lookups tables 
% which correspond to different foci. Make sure that the SisFocusIdx is set to 0
% Note that data is only available for foci index 1-4

%PathToXMLFile = [matRadRoot filesep 'BaseDataGeneration' filesep 'Basisdaten_HIT' filesep ...
%                 'LPD' filesep 'SingleGaussianBeamWidthSets_vor_12.9.2016' filesep 'CarbonLPD_Rifi3mm'];

PathToXMLFile = [matRadRoot filesep 'BaseDataGeneration' filesep 'Basisdaten_HIT' filesep ...
                 'LPD' filesep 'SingleGaussianBeamWidthSets_nach_12.9.2016' filesep 'CarbonGantryLPD_Rifi3mm_2016'];

machine       = matRad_readBeamWideningAIR(machine,PathToXMLFile);

%% interpolate double gaussian data from sparse sigma1, sigma2 and weight matrix
% lateral data from katja only describes scattering within the patient

% path to sampling points/Stützstellen provided by Katia P.
% if visBool is on then dont forget to press a key to step to the next plot
pathToSparseData = [pathTRiP filesep 'DDD' filesep '12C' filesep 'HIT_2D_DB_Cwith_KatjaP'];
machine = matRad_interpLateralBaseData(machine,pathTRiP,pathToSparseData,Identifier,SisFocusIdx,visBool);

%% parse spectra files for protons to enable LET calculations
machine = matRad_getLETfromSPC(machine,Identifier,pathTRiP);

%% parse dEdx file
[Meta, dEdx ] = matRad_readdEdx(pathTRiP);
save('dEdx.mat','dEdx');

%% parse RBEinitial data
[RBE] = matRad_readRBE(pathTRiP);
save('RBE.mat','RBE');

%% get dose averaged alpha and beta depth curves from the 37spc files, RBE
% inital file and dEdx file using the rapidScholz algorithm
pathToSPCFiles = [pathTRiP filesep 'SPC' filesep '12C' filesep 'RF3MM'];
[ BioDataHIT ] = matRad_getDepthDoseAvgLQM(pathToSPCFiles,[],0);
save('BioDataHIT.mat','BioDataHIT');

%% interpolate each entry in the ddd the corresponding depth alpha and depth beta curve
[machine] = matRad_interpDepthDoseAvgData(machine, BioDataHIT ,false, 0);

%% adjust distances of initial focus manually if you create base data for the old version
for i = 1:length(machine.data)
    machine.data(i).initFocus.dist(:,1) = 0;
    machine.data(i).initFocus.dist(:,2) = 20000;
end


%% save file
save(fileName,'machine');
