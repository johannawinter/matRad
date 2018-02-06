% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad_LEM_singleParticles
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
%% define paths paths
matRadDir = ['C:\Users\wieserh\Documents\matlab\matRad'];
TRiPdir = 'C:\Users\wieserh\Documents\TRiP98DATA_HIT-20131120';
addpath(matRadDir);

%% if biological base data does not exist - create it
if ~exist('RBE.mat','file') || ~exist('metadEdx.mat','file') ||  ~exist('dEdx.mat','file')
    RBEinitial       = matRad_readRBE(TRiPdir);
    [metadEdx,dEdx]  = matRad_readdEdx(TRiPdir);
    save('RBE','RBE')
    save('metadEdx','metadEdx');
    save('dEdx','dEdx');
else
    load('RBE');
    load('dEdx');
    load('metadEdx');
end


%% get tiusse parameters - LEM input
[tissue] = LEM_getTissueParameter('V79'); % V79 % CHO

%% set meta parameters
visBool           = 0;
Particle          = 'H';  % choose a particle type {'H','He','Li','Be','B','C'};
MaterialRho       = 1;
vNumParticles     = [1];   % number of one random traversal - multiple travelersals is not yet implemented
% number of energy values for which cell survival should be calculated  
% only use energies greater than 0.1 MeV
vEnergy           = dEdx.(Particle).Energy;  
vEnergy           = vEnergy(vEnergy>=0.1);

h =waitbar(0,'Please wait...');

% plot stopping power of the selected particle
if visBool
    figure('Name','Stopping power'),set(gcf,'Color',[1 1 1]); 
        plot(dEdx.(Particle).Energy,dEdx.(Particle).dEdx,'LineWidth',3),grid on, grid minor
       set(gca,'XScale','log'),set(gca,'YScale','log'),ylabel('dEdx'),xlabel('E in MeVcm^2/g per nucleon');
    h = figure;
end

Cnt         =  1;
Dose_Gy     = zeros(1,length(vEnergy));
N_TE        = zeros(1,length(vEnergy));
ExpSurvival = struct;

%loop over all energies
for IdxE =  1:length(vEnergy)
        
        NumParticles = 1;
        LET_MeVcm2_g = dEdx.(Particle).dEdx((dEdx.(Particle).Energy == vEnergy(IdxE)));   
        %determine the maximal range of delta electrons = ion track radius
        RadiusTrack_um = LEM_maxElectronRange(vEnergy(IdxE),0);
        AreaTarget_cm2 = (pi*(RadiusTrack_um + tissue.RadiusTarget_um)^2)* 1e-8;
        %calculate current dose in Gy
        Dose_Gy(IdxE)  = LEM_LET2Dose(LET_MeVcm2_g, NumParticles, AreaTarget_cm2);       
        N_TE(IdxE)     = LEM_Dose2CntHits(1,LET_MeVcm2_g,AreaTarget_cm2);
               
        % get impact parameter steps to determine the average damage to cell by one random traversal
        [vImpactParameter, dr] = LEM_getImpactParameterSteps(tissue.RadiusTarget_um,RadiusTrack_um,0);

        % assess for each spatial configuration(impact parameter) the biological effect
        vBioEffect  = zeros(1,length(vImpactParameter));
        for i = 1:length(vImpactParameter)
             vBioEffect(i) = LEM_singelHit( vImpactParameter(i),...
                                            tissue.RadiusTarget_um, ...
                                            RadiusTrack_um,tissue,...
                                            vEnergy(IdxE), ...
                                            dEdx.(Particle),....
                                            NumParticles,...
                                            LET_MeVcm2_g,visBool);
        end
        
        if visBool
            figure,plot(vImpactParameter,vBioEffect)
        end
        % calculate the survival propability
        vSurvivalProp    = exp(-vBioEffect);
        
        % expected survival of one randome traversal 
        CellSurvivalRef = (sum(2*pi*vImpactParameter .*vSurvivalProp .* dr));
        Atot            = sum(2*pi.*vImpactParameter .*dr);
        S_TE            = CellSurvivalRef/ Atot;
      
        % expected cell survival according to the paper
        AvgCellSurvival2                = trapz(vImpactParameter,vSurvivalProp)./((RadiusTrack_um)^2/2);
        ExpSurvival(IdxE).alpha_TEorg   =  (1-AvgCellSurvival2)*N_TE(IdxE);

        % this code is only here for validation
        alpha_target                  = interp1(RBE(1).(Particle).Energy,RBE(1).(Particle).RBE * RBE(1).alpha,vEnergy(IdxE));
        ExpSurvival(Cnt).N_TE_target  = alpha_target/(1-S_TE);
        
        ExpSurvival(IdxE).NumPart      = NumParticles;
        ExpSurvival(IdxE).Dose         = Dose_Gy(IdxE);
        ExpSurvival(IdxE).S            = S_TE;
        ExpSurvival(IdxE).Energy       = vEnergy(IdxE);
        ExpSurvival(IdxE).alpha_TE     = (1-S_TE)*N_TE(IdxE);
        ExpSurvival(IdxE).N_TE         = N_TE(IdxE);
        ExpSurvival(IdxE).LET_MeVcm2_g = LET_MeVcm2_g;

        if visBool
                subplot(121),plot(h,vImpactParameter,vSurvivalProp,'LineWidth',3),grid on, grid minor
                xlabel('impact parameter in µm'),ylabel(' S_{r},  cell survival of one particle');hold on
                title('cell suvival of one particle traversal','Interpreter','Latex')
                strPlot1{IdxE} = ['expectecd survival: ' num2str(AvgCellSuvival)];
                strPlot2{IdxE} = ...
                    [ num2str(vEnergy(IdxE)) ' MeV; NumPart: ' num2str(vNumParticles(IdxPart)) '; Dose: ' num2str( Dose_Gy(IdxPart)) ];

                subplot(122),plot(h,vImpactParameter,vBioEffect,'LineWidth',3),grid on, grid minor
                xlabel('impact parameter in µm'),ylabel('biological effect of one particle');hold on
                title('biological effect of one particle traversal','Interpreter','Latex')
        end

        waitbar(IdxE/length(vEnergy));
    
end  
close(h);


%% calculate cell survival just for one central traversal
Cnt                 = 1;
ExpSurvivalCentTrav = struct;

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

        ExpSurvivalCentTrav(Cnt).NumPart    = vNumParticles(IdxPart);
        ExpSurvivalCentTrav(Cnt).Dose       = Dose_Gy(IdxE);
        ExpSurvivalCentTrav(Cnt).S          = AvgCellSuvival;
        ExpSurvivalCentTrav(Cnt).Energy     = vEnergy(IdxE);
        ExpSurvivalCentTrav(Cnt).alpha_z   = vBioEffect/Dose_Gy(IdxE);
        ExpSurvivalCentTrav(Cnt).beta_z    = (tissue.Smax-ExpSurvivalCentTrav(Cnt).alpha_z)/(2*tissue.sDcut);
        
        ExpSurvivalCentTrav(Cnt).LET_MeVcm2_g = LET_MeVcm2_g;
                                             
    end
    waitbar(IdxE/length(vEnergy));
    Cnt= Cnt +1;

 end  
 


%% monte carlo run for single particle

CntRuns = 500;
ExpSurvivalMonteCarlo = struct;
h =waitbar(0,'Monte Carlo run...');
Cnt = 1;

 for IdxE = 1:length(vEnergy)

    [~,idx]      = min(abs(dEdx.(Particle).Energy-vEnergy(IdxE)));
    LET_MeVcm2_g = dEdx.(Particle).dEdx(idx);
    
    for IdxPart = 1:length(vNumParticles)
        
        Dose_Gy(IdxE)  = LEM_LET2Dose(LET_MeVcm2_g, vNumParticles(IdxPart), pi*((tissue.RadiusTarget_um*1e-4)^2)); 
        N_te(IdxE)     = LEM_Dose2CntHits(Dose_Gy(IdxE),LET_MeVcm2_g,pi*((tissue.RadiusTarget_um*1e-4)^2));

        %determine the maximal range of delta electrons = ion track radius
        RadiusTrack_um = LEM_maxElectronRange(vEnergy(IdxE),0);

        % get impact parameter steps to determine the average damage to cell by a random traversal
        UpperLim = tissue.RadiusTarget_um + RadiusTrack_um ;
        LowerLim = 0;
        %How to sample?
        vImpactParametery = ( UpperLim - LowerLim) .* rand(CntRuns,1) + LowerLim;
        vImpactParameterx = ( UpperLim - LowerLim) .* rand(CntRuns,1) + LowerLim;
        vImpactParameter = sqrt(vImpactParameterx.^2 + vImpactParametery.^2);
        %figure,scatter(vImpactParameterx,vImpactParametery);
        
        for i = 1:length(vImpactParameter)
              vBioEffect(i) = LEM_singelHit(vImpactParameter(i),...
                                            tissue.RadiusTarget_um, ...
                                            RadiusTrack_um,tissue,...
                                            vEnergy(IdxE), ...
                                            dEdx.(Particle),....
                                            vNumParticles(IdxPart),...
                                            LET_MeVcm2_g,visBool);
        end

        % calculate the survival propability
        AvgCellSuvival = exp(-mean(vBioEffect));

        ExpSurvivalMonteCarlo(Cnt).NumPart    = vNumParticles(IdxPart);
        ExpSurvivalMonteCarlo(Cnt).Dose       = Dose_Gy(IdxE);
        ExpSurvivalMonteCarlo(Cnt).S          = AvgCellSuvival;
        ExpSurvivalMonteCarlo(Cnt).Energy     = vEnergy(IdxE);
        ExpSurvivalMonteCarlo(Cnt).alpha_TE   =  (1-AvgCellSuvival)*N_TE(IdxE);
        ExpSurvivalMonteCarlo(Cnt).alpha      =  -log(AvgCellSuvival)/Dose_Gy(IdxE);
    end
    waitbar(IdxE/length(vEnergy));
    Cnt= Cnt +1;

 end  
 
 %% plot alpha
[~,ABratioIdx] = min(abs(([RBE.alpha]./ [RBE.beta]) - tissue.sAlphaBetaRatio));

alpha_z     = RBE(ABratioIdx).(Particle).RBE * RBE(1).alpha;
[alpha_D,~] = matRad_rapidScholz(RBE,dEdx,Particle,tissue,RBE(ABratioIdx).(Particle).Energy,alpha_z,[])

figure, set(gcf,'Color',[1 1 1]); 
plot(RBE(ABratioIdx).(Particle).Energy,alpha_z,'r-','LineWidth',4),hold on
plot(RBE(ABratioIdx).(Particle).Energy,alpha_D,'b-','LineWidth',4),hold on
plot([ExpSurvivalCentTrav.Energy],[ExpSurvivalCentTrav.alpha_z],'k-.','LineWidth',3),hold on
plot([ExpSurvival.Energy],[ExpSurvival.alpha_TE],'k--','LineWidth',4),hold on

% self validation - convert alpha_z to alpha_ion
[alpha_ionRS,beta_ionRS] = matRad_rapidScholz( RBE,dEdx,Particle,tissue,[ExpSurvival.Energy],[ExpSurvivalCentTrav.alpha_z],[ExpSurvivalCentTrav.beta_z]);
plot([ExpSurvival.Energy],alpha_ionRS,'go','LineWidth',3),hold on

legend({'TRiP $$\alpha_z$$','TRiP $$\alpha_D$$','$$\alpha_z$$ only central traversal','$$\alpha_D$$ full calculation','$$\alpha_D$$ RapidScholz based on $$\alpha_z$$'},'Interpreter','Latex')
set(gca,'xscale','log'),grid on, grid minor;
set(gca,'FontSize',12);
xlabel('energy [MeV/u]','Interpreter','Latex');
ylabel('$$ \alpha$$ [$$Gy^{-1}$$]','Interpreter','Latex');
title(['particle: ' Particle ', $\alpha_x=$' num2str(tissue.sAlphaX) ', $\beta_x=$' num2str(tissue.sBetaX) ', $r_{nuc}=$' num2str(tissue.RadiusTarget_um) '$\mu m$'],'Interpreter','Latex');
set(gca,'FontSize',26)


RBEval = LEM_calcRBE(tissue.sAlphaX,tissue.sBetaX, alpha_ionRS, beta_ionRS,[ExpSurvival.Dose]);
% plot RBE vs LET
figure, set(gcf,'Color',[1 1 1])
        plot([ExpSurvival.LET_MeVcm2_g],RBEval,'LineWidth',3),hold on, 
        set(gca,'xscale','log'),grid on, grid minor
        title(['particle: ' Particle ', $\alpha_x=$' num2str(tissue.sAlphaX) ', $\beta_x=$' num2str(tissue.sBetaX) ', $r_{nuc}=$' num2str(tissue.RadiusTarget_um) '$\mu m$'],'Interpreter','Latex');
        ylabel('$\alpha_{Ion}$','Interpreter','Latex');
        xlabel('LET [keV/um]','Interpreter','Latex');
        set(gca,'FontSize',26)

%% helper plots
figure, set(gcf,'Color',[1 1 1])
       plot(vEnergy,[ExpSurvival(:).S],'LineWidth',3),hold on, 
       plot(vEnergy,[ExpSurvivalCentTrav(:).S],'LineWidth',3), hold on
      % plot(vEnergy,[ExpSurvivalMonteCarlo(:).S],'--','LineWidth',3), hold on 
       set(gca,'xscale','log'),
       xlabel('energy'),ylabel('cell survival');
legend({'survival full','survival central'}),grid on,


figure, set(gcf,'Color',[1 1 1])
       plot(vEnergy,[ExpSurvival(:).N_TE],'LineWidth',3),hold on, 
       plot(vEnergy,[ExpSurvival(:).N_TE_target],'LineWidth',3), hold on
      % plot(vEnergy,[ExpSurvivalMonteCarlo(:).S],'--','LineWidth',3), hold on 
       set(gca,'xscale','log'), set(gca,'yscale','log'),
       xlabel('energy','Interpreter','Latex'),ylabel('$$N_{TE}$$','Interpreter','Latex');
legend({'$$N_{TE}$$','$$N_{TE target}$$'},'Interpreter','Latex'),grid on,
set(gca,'FontSize',16)


%% nice plot
% figure,set(gcf,'Color',[1 1 1]); 
% 
% load('data.mat');
% 
% data.(Particle).Trip.E     = RBEinitial(1).(Particle).Energy;
% data.(Particle).Trip.alpha = RBEinitial(1).(Particle).RBE * RBEinitial(1).alpha;
% data.(Particle).Own.E     = [ExpSurvivalMonteCarlo.Energy];
% data.(Particle).Own.alpha = [ExpSurvivalMonteCarlo.alpha];
% save('data.mat','data');
% 
% plot(data.('H').Trip.E,data.('H').Trip.alpha,'r','LineWidth',3),hold on
% plot(data.('H').Own.E,data.('H').Own.alpha,'b--','LineWidth',3),hold on
% plot(data.('C').Trip.E,data.('C').Trip.alpha,'r','LineWidth',3),hold on
% plot(data.('C').Own.E,data.('C').Own.alpha,'b--','LineWidth',3),hold on
% 
% 
% legend({' alpha from TRiP98',' recalc. alpha'});
% 
% xlabel('energy [MeV]'),ylabel('alpha in [Gy^-1]'),title('carbon ion alpha from LEM 1')
% set(gca,'FontSize',16)
% set(gca,'xscale','log'),grid on, grid minor;

%% sort data
vX = [0 ExpSurvival.Dose];
vY = [1 ExpSurvival.S];
[val ,ix] = sort(vX);

figure, 
plot(val,vY(ix),'-','LineWidth',3),grid on, grid minor, hold on
set(gca,'Yscale','log'), set(gca,'Ylim',[0 1]);
xlabel('Dose in Gy'); ylabel('Cell Survival'), title([Particle ' ions'])

vX = [0 ExpSurvivalMonteCarlo.Dose];
vY = [1 ExpSurvivalMonteCarlo.S];
[val ,ix] = sort(vX);
plot(val,vY(ix),'-','LineWidth',3),hold on

vX = [0 ExpSurvivalMonteCarlo.Dose];
vY = [1 ExpSurvivalMonteCarlo.S];
[val ,ix] = sort(vX);
plot(val,vY(ix),'-','LineWidth',3),hold on


legend({'analytical','MC','central traversal'})

 


%% test monte carlo run for multiple particles
% CntRuns = 10 * NumImpactSteps;
% h =waitbar(0,'Monte Carlo run...');
% vNumParticles = [2 3];
% Cnt = 1;
%  for IdxE = 10%1:length(vEnergy)
% 
%     [~,idx] = min(abs(dEdx.(Particle).Energy-vEnergy(IdxE)));
%     LET_MeVcm2_g = dEdx.(Particle).dEdx(idx);
%     
%     for IdxPart = 1:length(vNumParticles)
%         
%         Dose_Gy(IdxE)  = LEM_LET2Dose(LET_MeVcm2_g, vNumParticles(IdxPart), pi*((tissue.RadiusTarget_um*1e-4)^2)); 
%         N_te(IdxE)     = LEM_Dose2CntHits(Dose_Gy(IdxE),LET_MeVcm2_g,pi*((tissue.RadiusTarget_um*1e-4)^2));
% 
%         %initialize bio effet vector
%         vBioEffect  = zeros(NumImpactSteps,1);
%         %determine the maximal range of delta electrons = ion track radius
%         RadiusTrack_um = LEM_maxElectronRange(vEnergy(IdxE),0);
% 
%         % get impact parameter steps to determine the average damage to cell by a random traversal
%         UpperLim = tissue.RadiusTarget_um + RadiusTrack_um ;
%         LowerLim = 0;
%         vImpactParameter = [];
%         for j = 1:vNumParticles(IdxPart)
%             vImpactParameter(:,j) = ( UpperLim - LowerLim) .* rand(CntRuns,1) + LowerLim;
%         end
%         % assess for each spatial configuration(impact parameter) the
%         % biological effect
%         for i = 1:length(vImpactParameter)
%          vBioEffect(i) = LEM_multipleHIT( vImpactParameter(i,:), tissue.RadiusTarget_um, ...
%                                         RadiusTrack_um,tissue,vEnergy(IdxE), ...
%                                         dEdx.(Particle),vNumParticles(IdxPart),0);
%         end
% 
%         % calculate the survival propability
%         AvgCellSuvival = exp(-mean(vBioEffect));
% 
%         ExpSurvivalMonteCarlo(Cnt).NumPart = vNumParticles(IdxPart);
%         ExpSurvivalMonteCarlo(Cnt).Dose    = Dose_Gy(IdxE);
%         ExpSurvivalMonteCarlo(Cnt).S       = AvgCellSuvival;
%         ExpSurvivalMonteCarlo(Cnt).Energy  = vEnergy(IdxE);
%         ExpSurvivalMonteCarlo(Cnt).alpha   =  (1 - AvgCellSuvival) * N_te;
%     end
%     waitbar(IdxE/length(vEnergy));
%     Cnt= Cnt +1;
% 
%  end  
