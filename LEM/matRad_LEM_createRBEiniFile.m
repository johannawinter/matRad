% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad_LEM_createRBEiniFile
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2016, Hans-Peter Wieser
%
% h.wieser@dkfz.de
%
%
% This file is NOT part of the official matRad release. 
% This file has to be used only for internal purposes! 
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function simulates the LQM radiosensitivity parameters for a
% selected ion type based on different approaches 
clc
clear
close all
visBool = false;
color = colorspecs;

%% define paths paths
matRadDir = 'C:\Users\wieserh\Documents\matlab\matRad';
TRiPdir   = 'C:\Users\wieserh\Documents\TRiP98DATA_HIT-20131120';
addpath(matRadDir);

load('RBE');
load('dEdx');
load('metadEdx');

cParticles = {'H','He','Li','Be','B','C'};

for i = 1:size(RBE,2)

 
    tissue.sAlphaX = RBE(i).alpha;  %Gy-1
    tissue.sBetaX  = RBE(i).beta; %Gy-2
    tissue.sDcut   = RBE(i).cut;;   %Gy
    tissue.Smax    = tissue.sAlphaX + 2*tissue.sBetaX*tissue.sDcut;
    tissue.RadiusTarget_um = RBE(i).rnucleus; % µm
    tissue.sAlphaBetaRatio = tissue.sAlphaX / tissue.sBetaX;
        
    for j = 1:size(cParticles,2)
        
        Particle          = cParticles{1,j};  % choose a particle type {'H','He','Li','Be','B','C'};
        
        MaterialRho       = 1;
        vNumParticles     = [1];   % number of one random traversal - multiple travelersals is not yet implemented
        % number of energy values for which cell survival should be calculated  
        % only use energies greater than 0.1 MeV
        vEnergy           = dEdx.(Particle).Energy;  
        vEnergy           = vEnergy(vEnergy>=0.1);

        h =waitbar(0,'Please wait...');

        Dose_Gy     = zeros(1,length(vEnergy));
        N_TE        = zeros(1,length(vEnergy));
        ExpSurvival = struct;

        Cnt                 = 1;
        ExpSurvival = struct;
        
         for IdxE = 1:length(vEnergy)

            [~,idx] = min(abs(dEdx.(Particle).Energy-vEnergy(IdxE)));
            LET_MeVcm2_g = dEdx.(Particle).dEdx(idx);

            for IdxPart = 1:length(vNumParticles)

                AreaTarget_cm2 = (pi*(tissue.RadiusTarget_um)^2)* 1e-8;
                Dose_Gy(IdxE)  = LEM_LET2Dose(LET_MeVcm2_g, vNumParticles(IdxPart), AreaTarget_cm2); 

                %determine the maximal range of delta electrons = ion track radius
                RadiusTrack_um = LEM_maxElectronRange(vEnergy(IdxE),0);

                vImpactParameter = 0;

                vBioEffect = LEM_singelHit(vImpactParameter,...
                                           tissue.RadiusTarget_um, ...
                                           RadiusTrack_um,tissue,...
                                           vEnergy(IdxE), ...
                                           dEdx.(Particle),....
                                           vNumParticles(IdxPart),...
                                           LET_MeVcm2_g,visBool);

                % calculate the survival propability
                AvgCellSuvival = exp(-(vBioEffect));

                ExpSurvival(Cnt).NumPart    = vNumParticles(IdxPart);
                ExpSurvival(Cnt).Dose       = Dose_Gy(IdxE);
                ExpSurvival(Cnt).S          = AvgCellSuvival;
                ExpSurvival(Cnt).Energy     = vEnergy(IdxE);
                ExpSurvival(Cnt).alpha_z   = vBioEffect/Dose_Gy(IdxE);
                ExpSurvival(Cnt).beta_z    = (tissue.Smax-ExpSurvival(Cnt).alpha_z)/(2*tissue.sDcut);

                ExpSurvival(Cnt).LET_MeVcm2_g = LET_MeVcm2_g;

            end
            waitbar(IdxE/length(vEnergy));
            Cnt= Cnt +1;

         end    

        close(h);
        
        
        RBEini = interp1([ExpSurvival.Energy],[ExpSurvival.alpha_z]./ RBE(i).alpha,RBE(i).(Particle).Energy);

        RBE(i).(Particle).RBE = RBEini;
        
%         figure(),set(gcf,'Color',[1 1 1],'Units', 'normalized', 'Position', [0.1,0.1,0.8,0.8] ),
%         subplot(211),
%         plot(RBE(i).(Particle).Energy,RBE(i).(Particle).RBE),set(gca,'xScale','log');grid on, hold on
%         plot(RBE(i).(Particle).Energy,RBEini); 
% 
%         [alpha_D,~]     = matRad_rapidScholz(RBE,dEdx,Particle,tissue,[ExpSurvival.Energy],[ExpSurvival.alpha_z],[]);
% 
%         alpha_z         = RBE(i).(Particle).RBE * RBE(i).alpha;
%         [alpha_D_HIT,~] = matRad_rapidScholz(RBE,dEdx,Particle,tissue,RBE(i).(Particle).Energy,alpha_z,[]);
%         
%         subplot(212),
%         plot(RBE(i).(Particle).Energy,alpha_D_HIT,'LineWidth',3,'color',color.dkfzmB),set(gca,'xScale','log');grid on, hold on
%         plot(vEnergy,alpha_D,'--','LineWidth',3,'color',color.dre),
%         title(Particle),xlabel('Energy in MeV','Interpreter','Latex'),ylabel('$\alpha_D$ $[Gy^{-1}]$','Interpreter','Latex'),set(gca,'FontSize',16)
        
    end
    
    i
    
end


save('RBE_LEM1','RBE');




