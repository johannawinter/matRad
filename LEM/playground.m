clc
clear
close all

%% set parameters
addpath(['..' filesep 'baseDataHIT'])
load dEdx
Energy = 10; % MeV/u;
Particle = 'C';
MaterialRho = 1;
ImpactParameter = [0]; % in �m = radial position of hit
RadiusTarget    = [5]; % in �m

[ TotalResponse ] = LEM_singelHIT( ImpactParameter(1), RadiusTarget(1), Energy, dEdx.(Particle));