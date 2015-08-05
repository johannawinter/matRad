clc
clear
close all

%% set parameters
addpath(['..' filesep 'baseDataHIT'])
load dEdx
vEnergy        = [20];
vFluence_cm2   = [1];
Particle       = 'C';
MaterialRho    = 1;
RadiusTarget_um   = 10;    % in µm
NumImpactSteps = 100;
xRayData.sDcut   = 30;     % Gy
xRayData.sAlphaX = 0.18;   % Gy^-1
xRayData.sBetaX  = 0.028;  % Gy^-2
%%


[~,idx] = min(abs(dEdx.(Particle).energy-vEnergy));
LET_MeVcm2_g = dEdx.(Particle).dEdx(idx);
% calculate the dose for Number of particles from 1:50
for i = 1:50
    Dose_Gy(i,1) = LET2Dose(LET_MeVcm2_g, i, ((RadiusTarget_um*1e-4)^2));
end

%% mean number of traversals per unit dose
Dose_GyR = round(Dose_Gy);
[counts,centers] = hist(Dose_GyR,Dose_GyR(end));
N_TE = mean(counts);

for IdxE = 1:length(vEnergy)

    vBioEffect  = zeros(NumImpactSteps,1);
    weight      = zeros(NumImpactSteps,1);
    %determine the maximal range of delta electrons = track radius
    RadiusTrack_um = LEM_maxElectronRange(vEnergy(IdxE),0);
    % get impact parameter steps to determine the average damage to cell by
    % a random traversal
    [ImpactParameter, vDelta] = LEM_getImpactParameterSteps(RadiusTarget_um,RadiusTrack_um,NumImpactSteps);

    for i = 1:length(ImpactParameter)
     vBioEffect(i) = LEM_singelHIT( ImpactParameter(i), RadiusTarget_um, RadiusTrack_um,xRayData,vEnergy(IdxE), dEdx.(Particle),0);
    end

    vCellSurvival      = exp(-vBioEffect);
    TotSurvivalCentral =  exp(-vBioEffect(1));

    %% weightening according to code snipped from S.Greilich
    % only S_r dr needs to be calculated due to normalizing the contribution in
    % LEM_shellIntegration
    for i = 1:length(ImpactParameter)-1
     weight(i) = ImpactParameter(i+1)^2-ImpactParameter(i)^2;
    end
    weight = weight./(sum(weight));
    weight(end) = weight(end-1);
    expectedBioEffect    = sum(vBioEffect.*weight);
    expectedCellSurvival = exp(-expectedBioEffect);

    alpha_TE(IdxE) = (1-expectedCellSurvival)*N_TE;
    alpha_z(IdxE)  = expectedBioEffect/(sLET*1e-3);


end


figure,plot(vEnergy,alpha_TE),grid on, grid minor, set(gca,'xscale','log')
xlabel('Energy in MeV'), ylabel('alpha_{T,E}')

figure,plot(vEnergy,alpha_z),grid on, grid minor, set(gca,'xscale','log')
xlabel('Energy in MeV'), ylabel('alpha_{Z}')
