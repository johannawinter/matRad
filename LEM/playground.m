clc
clear
close all

%% set parameters
addpath(['..' filesep 'baseDataHIT'])
load dEdx
vEnergy        = 10;
Particle       = 'C';
MaterialRho    = 1;
RadiusTarget   = 10;    % in µm
NumImpactSteps = 100;
xRay.sDcut   = 30;     % Gy
xRay.sAlphaX = 0.18;   % Gy^-1
xRay.sBetaX  = 0.028;  % Gy^-2
%%
vBioEffectTrack = zeros(NumImpactSteps,1);
RadiusTrack = LEM_maxElectronRange(vEnergy,0);
[ImpactParameter, vDelta] = LEM_getImpactParameterSteps(RadiusTarget,RadiusTrack,NumImpactSteps);

for i = 1:length(ImpactParameter)
 vBioEffectTrack(i) = LEM_singelHIT( ImpactParameter(i), RadiusTarget, RadiusTrack,xRay,vEnergy, dEdx.(Particle),0);
end


vSurvivalTrack = exp(-vBioEffectTrack);
vTotSurvivalCentral =  exp(-vBioEffectTrack(1));
expS = sum((2*ImpactParameter.*exp(-vBioEffectTrack)'.*vDelta)./((RadiusTrack^2)));

%% Monte Carlo Integration
b = ImpactParameter(end);
a = ImpactParameter(1);
n = 10000;
randIdx = round((length(ImpactParameter)-1).*rand(n,1))+1;
S_TE = 0;
for i = 1:n
    S_TE = S_TE + vSurvivalTrack(randIdx(i))*ImpactParameter(randIdx(i));
end
Res = (b/n)*(2*S_TE/RadiusTrack^2);

%% according to existing code
for i = 1:length(ImpactParameter)-1
 weight(i) = ImpactParameter(i+1)^2-ImpactParameter(i)^2;
end
weight = weight./(sum(weight));
weight(end+1)=weight(end);
weight = weight/sum(weight);

Stot = exp(-sum(vBioEffectTrack.*weight'));


figure,plot(vEnergy,vAlpha_TE),grid on, grid minor, set(gca,'xscale','log')
xlabel('Energy in MeV'), ylabel('alpha_{T,E}')

figure,plot(vEnergy,vAlpha_Z),grid on, grid minor, set(gca,'xscale','log')
xlabel('Energy in MeV'), ylabel('alpha_{Z}')
