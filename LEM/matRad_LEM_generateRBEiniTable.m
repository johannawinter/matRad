% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad_LEM_generateRBEiniTable
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
% This function simulates the LQM radiosensitivity parameters for a
% selected ion type based on the central hit approach. In the end the simulated data
% is stored in *rbe files

clc
clear
close all
%% define paths paths
matRadDir = ['\\Mac\Home\Documents\Heidelberg\matRad'];
TRiPdir   = '\\Mac\Home\Documents\Heidelberg\TRiP98DATA';
addpath(matRadDir);

%% if biological base data does not exist - create it
if ~exist('RBEinitial.mat','file') || ~exist('metadEdx.mat','file') ||  ~exist('dEdx.mat','file')
    RBEinitial       = matRad_readRBE(TRiPdir);
    [metadEdx,dEdx]  = matRad_readdEdx(TRiPdir);
    save('RBEinitial','RBEinitial')
    save('metadEdx','metadEdx');
    save('dEdx','dEdx');
else
    load('RBEinitial');
    load('dEdx');
    load('metadEdx');
end

%% set meta parameters
visBool           = 0;
Particle          = {'H','He','Li','Be','B','C'};  % choose a particle type
MaterialRho       = 1;
vNumParticles     = [1];   % number of one random traversal - multiple travelersals is not yet implemented
h =waitbar(0,'Please wait...');
%% get tiusse parameters - LEM input
[tissue] = LEM_getTissueParameter('MDACC'); % V79 % CHO %MDACC

Cnt                 = 1;
ExpSurvivalCentTrav = struct;

bioDataSimu = struct;

for IdxPart = 1:size(Particle,2)
    
    vEnergy           = dEdx.(Particle{IdxPart}).Energy;   % number of energy values for which cell survival should be calculated  
    vEnergy           = vEnergy(vEnergy>=0.1 & vEnergy<=1000);       % only use energies greater than 0.1 MeV

    alpha_z = zeros(length(vEnergy),1);
    beta_z  = zeros(length(vEnergy),1);

     for IdxE = 1:length(vEnergy)

            [~,idx] = min(abs(dEdx.(Particle{IdxPart}).Energy-vEnergy(IdxE)));
            LET_MeVcm2_g = dEdx.(Particle{IdxPart}).dEdx(idx);

            vNumParticles = 1;
            AreaTarget_cm2 = (pi*(tissue.RadiusTarget_um)^2)* 1e-8;
            Dose_Gy        = LEM_LET2Dose(LET_MeVcm2_g, vNumParticles, AreaTarget_cm2); 

            %determine the maximal range of delta electrons = ion track radius
            RadiusTrack_um = LEM_maxElectronRange(vEnergy(IdxE),0);

            vImpactParameter = 0;

            vBioEffect = LEM_singelHIT(vImpactParameter,...
                                       tissue.RadiusTarget_um, ...
                                       RadiusTrack_um,tissue,...
                                       vEnergy(IdxE), ...
                                       dEdx.(Particle{IdxPart}),....
                                       vNumParticles,...
                                       LET_MeVcm2_g,visBool);

             alpha_z(IdxE) = vBioEffect/Dose_Gy; 
             Smax    = tissue.sAlphaX+(2*tissue.sBetaX)*tissue.sDcut;
             beta_z(IdxE)  = (Smax-alpha_z(IdxE))/(2*tissue.sDcut); 

     end
    
     bioDataSimu.(Particle{IdxPart}).energy = vEnergy;
     bioDataSimu.(Particle{IdxPart}).alpha_z = alpha_z;
    [bioDataSimu.(Particle{IdxPart}).alphaD, bioDataSimu.(Particle{IdxPart}).betaD]  = ...
        matRad_rapidScholz(RBEinitial,dEdx,Particle{IdxPart},tissue,vEnergy,alpha_z',beta_z');   

    waitbar(Cnt/(length(vEnergy)*size(Particle,2)));
    Cnt= Cnt +1;

end  
 
ixPart = 1;
[~,ABratioIdx] = min(abs(([RBEinitial.alpha]./ [RBEinitial.beta]) - tissue.sAlphaBetaRatio));
alpha_z     = RBEinitial(ABratioIdx).(Particle{ixPart}).RBE * RBEinitial(1).alpha;
[alpha_D,~] = matRad_rapidScholz(RBEinitial,dEdx,Particle{ixPart},tissue,RBEinitial(ABratioIdx).(Particle{ixPart}).Energy,alpha_z,[]);

figure, set(gcf,'Color',[1 1 1]); 
plot(RBEinitial(ABratioIdx).(Particle{ixPart}).Energy,alpha_z,'r:','LineWidth',4),hold on
plot(RBEinitial(ABratioIdx).(Particle{ixPart}).Energy,alpha_D,'b:','LineWidth',4),hold on
plot(bioDataSimu.(Particle{ixPart}).energy,bioDataSimu.(Particle{ixPart}).alpha_z  ,'r'),grid on
plot(bioDataSimu.(Particle{ixPart}).energy,bioDataSimu.(Particle{ixPart}).alphaD  ,'b'),grid on
set(gca,'xlim',[0 1000]),set(gca,'xscale','log')
legend({'TRiP $$\alpha_z$$','TRiP $$\alpha_D$$',...
        '$$\alpha_z$$','$$\alpha_D$$'},'Interpreter','Latex')


%save to RBEini file
metaInfo.filetype    = 'RBE';
metaInfo.fileversion = '20040715';
metaInfo.celltype    = 'IhaveNOclue';
metaInfo.mean        = '';
metaInfo.alpha       = ['   ' num2str(tissue.sAlphaX,'%1.4f')];
metaInfo.beta        = ['    ' num2str(tissue.sBetaX,'%1.4f')];
metaInfo.cut         = ['     ' num2str(tissue.sDcut,'%2.3f')];
metaInfo.rnucleus    = num2str(tissue.RadiusTarget_um,'%1.4f');

fNames = fieldnames(metaInfo);
FileName = ['custom.rbe'];
fileID = fopen(FileName,'w');
for j = 1:length(fNames)
    fprintf(fileID,['!' fNames{j,1} '    ' metaInfo.(fNames{j,1}) '\n']);
end

for i = 1:size(Particle,2)
   
    fprintf(fileID,['!projectile  '  [num2str(i*2) Particle{i}] '\n']);
    fprintf(fileID,'!rbe\n');
    fprintf(fileID,'# E/(MeV/u)  RBE\n');

    A = [[bioDataSimu.(Particle{i}).energy]'   [(bioDataSimu.(Particle{i}).alpha_z./tissue.sAlphaX)]]';
    fprintf(fileID,'%1.3E %1.3f \r\n',A);

end
fclose(fileID);
clear metaInfo










