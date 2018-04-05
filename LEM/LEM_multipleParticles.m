% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LEM_multipleParticles
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2016, Hans-Peter Wieser
%
% h.wieser@dkfz.de
%
% This file is NOT part of the official matRad release. 
% This file has to be used only for internal purposes! 
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This simulates an ensemble of particles. THIS IS SCRIPT IS NOT WORKING
% RIGHT NOW - FURTHER WORK NEEDS TO BE DONE

clc
clear
close all

matRadDir = ['C:\Users\wieserh\Documents\matRad'];

%pathToResource = 'E:\TRiP98DATA_HIT-20131120';
TRiPdir = 'E:\TRiP98DATA_HIT-20131120';
addpath(matRadDir);
%% if biological base data does not exist - create it

% RBEinitial = matRad_readRBE(TRiPdir);
% [metadEdx,dEdx]   = matRad_readdEdx(TRiPdir);
% save('RBE','RBEinitial')
% save('metadEdx','metadEdx');
% save('dEdx','dEdx');


load('RBE');
load('dEdx');
load('metadEdx');

%% define tiusse
tissue.sAlphaX         = 0.1;  %Gy-1
tissue.sBetaX          = 0.05; %Gy-2
tissue.sDcut           = 30;   %Gy
tissue.RadiusTarget_um = 5; % µm

%% set meta parameters
visBool           = 0;
Particle          = 'C';
MaterialRho       = 1;
vNumParticles     = 11:1:20;   % num of one random traversal
% number of energy values for which cell survival should be calculated  
% only use energies greater than 0.5 MeV
vEnergy           = dEdx.(Particle).Energy;  
vEnergy           = vEnergy(vEnergy>=0.5);
vEnergy           = 20;

h =waitbar(0,'Please wait...');

% plot stopping power of the selected particle
if visBool
    figure('Name','Stopping power'),set(gcf,'Color',[1 1 1]); 
        plot(dEdx.(Particle).Energy,dEdx.(Particle).dEdx,'LineWidth',3),grid on, grid minor
       set(gca,'XScale','log'),set(gca,'YScale','log'),ylabel('dEdx'),xlabel('E in MeVcm^2/g per nucleon');
    h = figure;
end


Cnt = 1;
%loop over all energies
for IdxE =  1:length(vEnergy)

    LET_MeVcm2_g = dEdx.(Particle).dEdx(find(dEdx.(Particle).Energy == vEnergy(IdxE)));
    
    for IdxPart = 1:length(vNumParticles)
        
        
        %determine the maximal range of delta electrons = ion track radius
        RadiusTrack_um = LEM_maxElectronRange(vEnergy(IdxE),0);
        AreaTarget_cm2 = (pi*(RadiusTrack_um + tissue.RadiusTarget_um)^2)* 1e-8;
        %calculate current dose in Gy
        Dose_Gy(Cnt)  = LEM_LET2Dose(LET_MeVcm2_g, vNumParticles(IdxPart), AreaTarget_cm2);
        
        N_TE(Cnt)     = LEM_Dose2CntHits(1,LET_MeVcm2_g,AreaTarget_cm2);
                    

        % get impact parameter steps to determine the average damage to cell by one random traversal
        %[vImpactParameter, dr] = LEM_getImpactParameterSteps(tissue.RadiusTarget_um,RadiusTrack_um,0);
        vImpactParameter = 0;
        % assess for each spatial configuration(impact parameter) the biological effect
        for i = 1:length(vImpactParameter)
             vBioEffect(i) = LEM_singelHIT( vImpactParameter(i),...
                                            tissue.RadiusTarget_um, ...
                                            RadiusTrack_um,tissue,...
                                            vEnergy(IdxE), ...
                                            dEdx.(Particle),....
                                            vNumParticles(IdxPart),...
                                            LET_MeVcm2_g,visBool);
        end
        
        % calculate the survival propability
        vSurvivalProp = exp(-vBioEffect);
        
        S_TE_central     = vSurvivalProp(1);
        alpha_TE_central = vBioEffect(1)/Dose_Gy(Cnt);
        
        ExpSurvival(Cnt).NumPart    = vNumParticles(IdxPart);
        ExpSurvival(Cnt).Dose       = Dose_Gy(Cnt);
        ExpSurvival(Cnt).S          = S_TE;
        ExpSurvival(Cnt).Energy     = vEnergy(Cnt);
        ExpSurvival(Cnt).alpha_TE   =  alpha_TE_central;

        waitbar(IdxE/length(vEnergy));
        Cnt= Cnt +1;
    end
    
end  
close(h);

