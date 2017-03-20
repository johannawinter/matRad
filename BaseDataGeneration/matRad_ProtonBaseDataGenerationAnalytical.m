% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad_ProtonBaseDataGenerationAnalytical script
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
% This script can be used to generate depth dose and depth let profiles from protons

clc
clear 
close all

%% set global path which is valid for all subsections

load protons_Generic.mat

for i = 1:numel(machine.data)
   cla
    [ddd] = matRad_getDDDfromAnalyCalc('p',machine.data(i).range, 1, 1, machine.data(i).depths);
     machine.data(i).LET = ddd.LETd_RT;
     machine.data(i).Z   = ddd.dose;
end

