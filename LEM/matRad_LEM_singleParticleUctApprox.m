% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad_LEM_singleParticleUctApprox
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
% This function performs sampling from LEM input parameters propagates
% these realizations through LEM, dose averaging into a 1D C12 SOBP

clc
clear
close all

matRadDir = ['C:\Users\wieserh\Documents\matlab\matRad'];
TRiPdir   = 'C:\Users\wieserh\Documents\TRiP98DATA_HIT-20131120';
addpath(matRadDir);


% %% if biological base data does not exist - create it
% 
% RBEinitial = matRad_readRBE(TRiPdir);
% [metadEdx,dEdx]   = matRad_readdEdx(TRiPdir);
% save('RBE','RBEinitial')
% save('metadEdx','metadEdx');
% save('dEdx','dEdx');

load('RBE');
load('dEdx');
load('metadEdx');

%%
FLAG_PLOT       = false;
FLAG_SAVE       = true; 
foldername      = 'bioUCT';
defaultFontSize = 16;
exportPath      = [pwd filesep 'LEM' filesep 'expFigures' filesep foldername];   
color           = colorspecs();   
addpath(exportPath) ; showInfoFlag = false; addpath('LEM'); addpath(['APM' filesep 'expFigures']); 

extraAxisAptions = ['title style={font=\normalsize},'...
                    'xlabel style={font=\normalsize},'...
                    'ylabel style={font=\normalsize},',...
                    'legend style={font=\normalsize},',...
                    'ticklabel style={font=\small},'];

%% define tiusse
tissue.sAlphaX         = 0.1;  %Gy-1
tissue.sBetaX          = 0.05; %Gy-2
tissue.sDcut           = 30;   %Gy
tissue.RadiusTarget_um = 5; % µm
tissue.sAlphaBetaRatio = tissue.sAlphaX / tissue.sBetaX;

%% define uncertainties
tissue.Realisations = 1;

tissue.sAlphaXvar   = tissue.sAlphaX*0;%;0.1;
tissue.vAlphaUct    = ((randn(tissue.Realisations,1)*tissue.sAlphaXvar) + tissue.sAlphaX);
tissue.sAlphaXnom   = tissue.sAlphaX;
%hist(tissue.vAlphaUct);

tissue.sBetaXvar    = tissue.sBetaX*0;%0.2;
tissue.vBetaUct     = ((randn(tissue.Realisations,1)*tissue.sBetaXvar) + tissue.sBetaX) ;
tissue.sBetaXnom    = tissue.sBetaX;

tissue.sDcutVar     = tissue.sDcut*0;
tissue.vDcutUct     = ((randn(tissue.Realisations,1)*tissue.sDcutVar) + tissue.sDcut);
tissue.sDcutNom     = tissue.sDcut;

tissue.sRnucVar     = tissue.RadiusTarget_um*0;
tissue.vRnucUct     = ((randn(tissue.Realisations,1)*tissue.sRnucVar) + tissue.RadiusTarget_um); 
tissue.sRnucNom     = tissue.RadiusTarget_um;


%% set meta parameters
visBool           = 0;
Particle          = {'H','He','Li','Be','B','C'};
MaterialRho       = 1;
vNumParticles     = [1];   

% initialize some variables
for IdxPart = 1:size(Particle,2)
    vEnergy                           = dEdx.(Particle{IdxPart}).Energy;  
    vEnergy                           = vEnergy(vEnergy>=0.2);
    UctDataAlphaZ.(Particle{IdxPart}) = zeros(tissue.Realisations,length(vEnergy));
    UctDataBetaZ.(Particle{IdxPart})  = zeros(tissue.Realisations,length(vEnergy));
    UctDataAlphaD.(Particle{IdxPart}) = zeros(tissue.Realisations,length(vEnergy));
    UctDataBetaD.(Particle{IdxPart})  = zeros(tissue.Realisations,length(vEnergy));
end

h =waitbar(0,'Please wait...');
% alpha uncertainty - only central traversal 
for IdxPart = 1:size(Particle,2)

vEnergy  = dEdx.(Particle{IdxPart}).Energy;  
vEnergy  = vEnergy(vEnergy>=0.2);

for IdxReal = 1:tissue.Realisations

     tissue.sAlphaX         = tissue.vAlphaUct(IdxReal);
     tissue.sBetaX          = tissue.vBetaUct(IdxReal);
     tissue.sDcut           = tissue.vDcutUct(IdxReal);
     tissue.RadiusTarget_um = tissue.vRnucUct(IdxReal);
            
     for IdxE = 1:length(vEnergy)

        [~,idx] = min(abs(dEdx.(Particle{IdxPart}).Energy-vEnergy(IdxE)));
        LET_MeVcm2_g = dEdx.(Particle{IdxPart}).dEdx(idx);

        vNumParticles = 1;
        AreaTarget_cm2 = (pi*(tissue.RadiusTarget_um)^2)* 1e-8;
        Dose_Gy  = LEM_LET2Dose(LET_MeVcm2_g, vNumParticles, AreaTarget_cm2); 

        %determine the maximal range of delta electrons = ion track radius
        RadiusTrack_um = LEM_maxElectronRange(vEnergy(IdxE),0);

        vImpactParameter = 0;

        vBioEffect = LEM_singelHit(vImpactParameter,...
                                   tissue.RadiusTarget_um, ...
                                   RadiusTrack_um,tissue,...
                                   vEnergy(IdxE), ...
                                   dEdx.(Particle{IdxPart}),....
                                   vNumParticles,...
                                   LET_MeVcm2_g,visBool);
                               
         alpha_z = vBioEffect/Dose_Gy; 
         Smax    = tissue.sAlphaX+(2*tissue.sBetaX)*tissue.sDcut;
         beta_z  = (Smax-alpha_z)/(2*tissue.sDcut); 
         UctDataAlphaZ.(Particle{IdxPart})(IdxReal,IdxE) = alpha_z;      
         UctDataBetaZ.(Particle{IdxPart})(IdxReal,IdxE)  = beta_z;    
     end
     
     % apply rapid scholz algorithm
     [UctDataAlphaD.(Particle{IdxPart})(IdxReal,:),UctDataBetaD.(Particle{IdxPart})(IdxReal,:) ] =...
                matRad_rapidScholz(RBE,dEdx,Particle{IdxPart},tissue,vEnergy,UctDataAlphaZ.(Particle{IdxPart})(IdxReal,:),...
                UctDataBetaZ.(Particle{IdxPart})(IdxReal,:));   
 
     waitbar(IdxReal/(tissue.Realisations*size(Particle,2)));
    
     if IdxReal == tissue.Realisations
       tissue.sAlphaX         = tissue.sAlphaXnom;
       tissue.sBetaX          = tissue.sBetaXnom;
       tissue.sDcut           = tissue.sDcutNom;
       tissue.RadiusTarget_um = tissue.sRnucNom;
     end
     
end

end
close(h)

%% plot a specific particle
if FLAG_PLOT
    CurrentParticle = 'C';
    for j = 1:2
        if j == 1
            UctDataAlphaLoop = UctDataAlphaZ.(CurrentParticle);
        elseif j == 2
            UctDataAlphaLoop = UctDataAlphaD.(CurrentParticle);
        end
        figure('units','normalized','outerposition',[0 0 1 1]),set(gcf,'Color',[1 1 1]) ,hold on
        set(gca,'xscale','log'),
        H(1) = shadedErrorBar(vEnergy,UctDataAlphaLoop,{@mean, @(x) 2*std(x)}, '-r',0);
        H(2) = shadedErrorBar(vEnergy,UctDataAlphaLoop,{@mean, @(x) 1*std(x)}, '-g',0);
        H(3) = shadedErrorBar(vEnergy,UctDataAlphaLoop,{@mean, @(x) 0.5*std(x)}, '-b',0);
        grid on,grid minor
        legend([H(3).mainLine, H.patch], '\mu','2\sigma','\sigma','0.5\sigma');
        xlabel('energy in MeV/u','Interpreter','Latex');
        if j  == 1
             ylabel('$$ \alpha_z$$  in $$Gy^{-1}$$','Interpreter','Latex');
        elseif j == 2
             ylabel('$$ \alpha_D$$  in $$Gy^{-1}$$','Interpreter','Latex');
        end
        title(['particle: ',CurrentParticle ', $\alpha_x=$' num2str(tissue.sAlphaX) ', $\beta_x=$' num2str(tissue.sBetaX) ', $r_{nuc}=$' num2str(tissue.RadiusTarget_um) '$\mu m$'],'Interpreter','Latex');
        set(gca,'FontSize',26)

        newStr = strrep(num2str(tissue.sAlphaXvar),'.','_');
        if j == 1
            filename  = ['alpha_z_uct' newStr '.tex'];
        elseif j == 2
            filename  = ['alpha_D_uct' newStr '.tex'];
        end
        latexPath = [exportPath filesep filename];

        if FLAG_SAVE && exist('matlab2tikz','file') == 2
            cleanfigure;
            matlab2tikz([latexPath],'height', '12cm', 'width', '21cm','showInfo',showInfoFlag,'standalone', true,...
                'extraaxisoptions',extraAxisAptions);         

            currPath = pwd; cd(exportPath);
            if ispc  
                 command = sprintf('pdflatex %s',filename);
            elseif ismac
                 command = sprintf('/Library/Tex/texbin/pdflatex %s',[filename]);
            end
            [status,cmdout] = system(command); cd(currPath);
            if status > 0
                warning(['couldnt compile pdf: ' cmdout]);
            end
        end


    end

end


%% apply dose averaging with spc files
load('carbon_HIT.mat');
[machine,BioDataHIT] = LEM_DoseAvg(TRiPdir,machine,UctDataAlphaD,UctDataBetaD,vEnergy,tissue);

% plot HIT vs. myLEM
machineOrg = load('carbon_HIT.mat');
machineOrg = machineOrg.machine;
ix = 198;
vDepth    = machine.data(ix).depths'./machine.data(ix).peakPos;
figure('units','normalized','outerposition',[0 0 1 1]),set(gcf,'Color',[1 1 1]) ,hold on
plot(vDepth,machine.data(ix).alpha(:,1)),hold on
plot(vDepth,machineOrg.data(ix).alpha(:,1)),hold on

%save('BioDataHIT.mat','BioDataHIT');
%save('carbon_HIT_LEM.mat','machine');

        
%% plot depth depended dose averaged alpha and beta parameters  
ix = 198;
vDepth    = machine.data(ix).depths'./machine.data(ix).peakPos;
figure('units','normalized','outerposition',[0 0 1 1]),set(gcf,'Color',[1 1 1]) ,hold on
subplot(211),plot(vDepth,machine.data(ix).alpha)
xlabel('z [mm]','Interpreter','Latex');grid on, grid minor
ylabel('$Gy^{-1}$','Interpreter','Latex');
title(['depth-dependend dose-averaged $\alpha$-profiles;',' C12, $E_0$=' num2str(machine.data(ix).energy)],'Interpreter','Latex')
set(gca,'FontSize',20),set(gca,'YLim',[0 1.4]);

subplot(212),hold on
alpha_avg = machine.data(ix).alpha';

H(1) = shadedErrorBar(vDepth,alpha_avg,{@mean, @(x) 2*std(x)}, '-r',0);
H(2) = shadedErrorBar(vDepth,alpha_avg,{@mean, @(x) 1*std(x)}, '-g',0);
H(3) = shadedErrorBar(vDepth,alpha_avg,{@mean, @(x) 0.5*std(x)}, '-b',0);
xlabel('z [mm]','Interpreter','Latex'),grid on,grid minor
legend([H(3).mainLine, H.patch], {'\mu','2\sigma','\sigma','0.5\sigma'},'FontSize',16);
xlabel('depth','Interpreter','Latex');
ylabel('depth-depend dose-avg $$ \alpha$$  in $$Gy^{-1}$$','Interpreter','Latex');
 title(['pristine carbon ion beam 350MeV/u, $$ \alpha_x$$=0.1, $$ \beta_x$$=0.05,' ...
    'Realizations = ' num2str(tissue.Realisations) ', $$\sigma({\alpha_x})=' num2str(tissue.sAlphaXvar*100) '\%$$'],'Interpreter','Latex')
set(gca,'FontSize',20),set(gca,'YLim',[0 1.4]);

filename  = ['doseAVG_alpha_profiles_cut' newStr '.tex'];
latexPath = [exportPath filesep filename];

if FLAG_SAVE && exist('matlab2tikz','file') == 2
    cleanfigure;
    matlab2tikz([latexPath],'height', '17cm', 'width', '21cm','showInfo',showInfoFlag,'standalone', true,...
        'extraaxisoptions',extraAxisAptions);         

    currPath = pwd; cd(exportPath);
    if ispc  
         command = sprintf('pdflatex %s',filename);
    elseif ismac
         command = sprintf('/Library/Tex/texbin/pdflatex %s',[filename]);
    end
    [status,cmdout] = system(command);
    delete('*.aux');delete('*.log');
    cd(currPath);
    if status > 0
        warning(['couldnt compile pdf: ' cmdout]);
    end
end


% posBig = get(gca, 'Position');hold on;
% subFig(1) = axes('Position',[0.6 0.25 0.15 0.15]);
% set(subFig(1),'FontSize',12)
figure('units','normalized','outerposition',[0 0 1 1]),set(gcf,'Color',[1 1 1]) ,hold on
plot(vDepth,std(alpha_avg,1,1),'LineWidth',3),grid on, grid minor,hold on,
title('std');


meanStd = mean(std(alpha_avg,1,1));
plot(vDepth,meanStd*ones(numel(vDepth)),'LineWidth',3),grid on, grid minor,hold on,

%% create 1D treatment plan
%load('carbon_HIT_LEM.mat');
VoxelSize                  = 1;   %[mm]
PatientLength              = 200; %[mm]
TargetLength               = 51;  %[mm]
TargetEntry                = 99;
Voxel.Position             = (VoxelSize/2):VoxelSize:PatientLength;
Voxel.NumberOfVoxels       = numel(Voxel.Position);
Voxel.IxNT                 = (Voxel.Position < TargetEntry) | (Voxel.Position > (TargetEntry + TargetLength));
Voxel.IxT                  = ~Voxel.IxNT ;
Voxel.presDose(Voxel.IxNT) = 1;
Voxel.presDose(Voxel.IxT)  = 3;
Voxel.penalty(Voxel.IxNT)  = 10;
Voxel.penalty(Voxel.IxT)   = 100;
Voxel.presEffect = tissue.sAlphaXnom .* Voxel.presDose + tissue.sBetaXnom .* Voxel.presDose.^2; 

availablePeakPos  = [machine.data.peakPos];
availableEnergies = [machine.data.peakPos];
Spot.Position     = availablePeakPos(availablePeakPos > TargetEntry & availablePeakPos < TargetEntry + TargetLength);
Spot.IxEnergy     = find(availablePeakPos > TargetEntry & availablePeakPos < TargetEntry + TargetLength);
Spot.NumOfSpots   = length(Spot.Position);
Spot.Weights      = ones(Spot.NumOfSpots,1);

%% matRad world

[pln,cst] = matRad_getDummy1D_Plan(Voxel,tissue);

dij.totalNumOfBixels = Spot.NumOfSpots;
dij.numOfVoxels      = Voxel.NumberOfVoxels;
dij.dimensions       = [Voxel.NumberOfVoxels 1];
dij.dose  = zeros(Voxel.NumberOfVoxels,Spot.NumOfSpots);
dij.alpha = zeros(Voxel.NumberOfVoxels,Spot.NumOfSpots);
dij.beta  = zeros(Voxel.NumberOfVoxels,Spot.NumOfSpots);

%% calculate mean scenario / nominal scenario
for j = 1:Spot.NumOfSpots
    baseEntry                = machine.data(Spot.IxEnergy(j));
    ix                       = baseEntry.depths(end) > Voxel.Position;
    dij.physicalDose(ix,j)   = interp1(baseEntry.depths,baseEntry.Z,Voxel.Position(ix),'linear');
    dij.mAlphaDose(ix,j)     = (interp1(baseEntry.depths,mean(baseEntry.alpha,2),Voxel.Position(ix),'linear')) .* dij.physicalDose(ix,j)' ;
    dij.mSqrtBetaDose(ix,j)  = sqrt(interp1(baseEntry.depths,mean(baseEntry.beta,2),Voxel.Position(ix),'linear')) .* dij.physicalDose(ix,j)';
end

addpath('C:\Users\wieserh\Documents\matRad')
resultGUI = matRad_fluenceOptimization(dij,cst,pln);

for i = 1:tissue.Realisations
    for j = 1:Spot.NumOfSpots
        baseEntry                   = machine.data(Spot.IxEnergy(j));
        ix                          = baseEntry.depths(end) > Voxel.Position;
        dijUCT.physicalDose(ix,j)   = interp1(baseEntry.depths,baseEntry.Z,Voxel.Position(ix),'linear');
        dijUCT.mAlphaDose(ix,j)     = (interp1(baseEntry.depths,baseEntry.alpha(:,i),Voxel.Position(ix),'linear')) .* dij.physicalDose(ix,j)' ;
        dijUCT.mSqrtBetaDose(ix,j)  = sqrt(interp1(baseEntry.depths,baseEntry.beta(:,i),Voxel.Position(ix),'linear')) .* dij.physicalDose(ix,j)';
    end  
    dijUCT.dimensions       = [Voxel.NumberOfVoxels 1];
    cst{1,5}.alphaX         = tissue.vAlphaUct(i);
    cst{2,5}.alphaX         = tissue.vAlphaUct(i);
    cst{1,5}.betaX          = tissue.vBetaUct(i);
    cst{2,5}.betaX          = tissue.vBetaUct(i);
    
    UCTresultGUI(i)         = matRad_calcCubes(resultGUI.w,dijUCT,cst);
    
end


figure('units','normalized','outerposition',[0 0 1 1]),set(gcf,'Color',[1 1 1]) ,hold on
H(1) = shadedErrorBar(Voxel.Position,[UCTresultGUI.RBExDose]',{@mean, @(x) 2*std(x)}, '-r',0);
H(2) = shadedErrorBar(Voxel.Position,[UCTresultGUI.RBExDose]',{@mean, @(x) 1*std(x)}, '-g',0);
H(3) = shadedErrorBar(Voxel.Position,[UCTresultGUI.RBExDose]',{@mean, @(x) 0.5*std(x)}, '-b',0);
grid on,grid minor
legend([H(3).mainLine, H.patch], '\mu','2\sigma','\sigma','0.5\sigma');
xlabel('depth [mm]','Interpreter','Latex');
ylabel('RBE x dose [Gy(RBE)]','Interpreter','Latex');
title('C12 SOPB with varied $$\alpha_x$$ and $$\beta_x$$','Interpreter','Latex');
set(gca,'FontSize',20)
hold on, plot(Voxel.Position,resultGUI.physicalDose,'k','LineWidth',2)




figure,set(gcf,'Color',[1 1 1]), hold on
Cnt = 1;
for i = 1:8:tissue.Realisations
    plot(Voxel.Position,[UCTresultGUI(i).RBExDose],'LineWidth',3);
    str{Cnt} = [num2str(fix((1-(0.1./tissue.vAlphaUct(i))) * 100)) ' % ,' num2str(fix((1-(0.05./tissue.vBetaUct(i))) * 100)) ' % '];
    Cnt = Cnt +1;
end

legend(str),
grid on, grid minor

    
    