clc
clear
close all

%% set parameters
addpath(['..' filesep 'baseDataHIT'])
load dEdx
vEnergy        = [1 5 10 20 30 40 50];
vFluence_cm2       = [0.5 1 2 3 5 10];
Particle       = 'C';
MaterialRho    = 1;
RadiusTarget   = 10;    % in µm
NumImpactSteps = 100;
xRay.sDcut   = 30;     % Gy
xRay.sAlphaX = 0.18;   % Gy^-1
xRay.sBetaX  = 0.028;  % Gy^-2
%%


[~,idx] = min(abs(dEdx.(Particle).energy-vEnergy(2)));
sLET = dEdx.(Particle).dEdx(idx);
for i = 1:50
    Dose_Gy(i,1) = LET2Dose(sLET, i, ((RadiusTarget*1e-4)^2));
end

Dose_GyR = round(Dose_Gy);
[counts,centers] = hist(Dose_GyR,Dose_GyR(end));
%% mean number of traversals per unit dose
N_TE = mean(counts);

for Eidx = 1:length(vEnergy)

    vBioEffect  = zeros(NumImpactSteps,1);
    weight      = zeros(NumImpactSteps,1);
    RadiusTrack = LEM_maxElectronRange(vEnergy(Eidx),0);
    [ImpactParameter, vDelta] = LEM_getImpactParameterSteps(RadiusTarget,RadiusTrack,NumImpactSteps);

    for i = 1:length(ImpactParameter)
     vBioEffect(i) = LEM_singelHIT( ImpactParameter(i), RadiusTarget, RadiusTrack,xRay,vEnergy(Eidx), dEdx.(Particle),0);
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

    alpha_TE(Eidx) = (1-expectedCellSurvival)*N_TE;
    alpha_z(Eidx)  = expectedBioEffect/(sLET*1e-3);


end


figure,plot(vEnergy,alpha_TE),grid on, grid minor, set(gca,'xscale','log')
xlabel('Energy in MeV'), ylabel('alpha_{T,E}')

figure,plot(vEnergy,alpha_z),grid on, grid minor, set(gca,'xscale','log')
xlabel('Energy in MeV'), ylabel('alpha_{Z}')
