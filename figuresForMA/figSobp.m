% figure for MA
% SOBP
clear, close all
addpath(genpath('submodules'))

% load('protons_HIT_APMgantry')
load('BOXPHANTOM')

% set up proton plan
pln.radiationMode   = 'protons';
pln.machine         = 'HIT_APMgantry';
pln.numOfFractions  = 30;

pln.propStf.bixelWidth      = 5;
pln.propStf.gantryAngles    = [270];
pln.propStf.couchAngles     = [0];
pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);

pln.propOpt.bioOptimization = 'const_RBExD';
pln.propOpt.runDAO          = false;
pln.propOpt.runSequencing   = false;

% calculate proton plan
stf = matRad_generateStf(ct,cst,pln);
dij = matRad_calcParticleDose(ct,stf,pln,cst);
resultGUI = matRad_fluenceOptimization(dij,cst,pln);

% % get proton depth dose curve through the isocenter
% ddPr = resultGUIPr.physicalDose(plnPr.propStf.isoCenter(1)/ct.resolution.x, :, ...
%     plnPr.propStf.isoCenter(3)/ct.resolution.z);
% ddPr = ddPr(40:120);
% depthsPr = ct.resolution.x * [0:80];


%% plot single BP that compose the SOBP
% from weightAnalysis

% load machine
fileName = [pln.radiationMode '_' pln.machine];
load(['C:\Matlab\matrad\' fileName]);
offset = machine.data(1).offset;                % offset to add to depth dose curves due to BAMS

% get all available energies from machine file
availableEnergies = zeros(1,size(machine.data,2));
for j = 1:size(machine.data,2)
    availableEnergies(j) = machine.data(j).energy;
end

% find indices of energies that were used in treatment plan
numEnergiesPerRay = zeros(1,size(stf.ray,2));
for j = 1:size(stf.ray,2)
    numEnergiesPerRay(j) = size(stf.ray(j).energy,2);
end
numEnergiesMax = max(numEnergiesPerRay);

energyIxAllRays = zeros(numEnergiesMax,stf.numOfRays);
for j = 1:numEnergiesMax
    for k = 1:stf.numOfRays
        try
            energyIxAllRays(j,k) = find(availableEnergies == stf.ray(k).energy(j));
        catch
            energyIxAllRays(j,k) = 0;
        end
    end
end

% pool used energy indices in one vector and find corresponding energies
energyIxUsed = unique(energyIxAllRays(:));
energyIxUsed = energyIxUsed(find(energyIxUsed));

energiesUsed = zeros(size(energyIxUsed));
for j = 1:length(energyIxUsed)
    energiesUsed(j) = machine.data(energyIxUsed(j)).energy;
end

% get tabulated depths for used energies
z = zeros(numEnergiesMax, length(machine.data(energyIxUsed(1)).depths));
for j = 1:numEnergiesMax
    z(j,:) = machine.data(energyIxUsed(j)).depths;
end

% separate weights in rays (columns, j) and energies (rows,k)
% The weight vector w is divided into beams which are divided into rays
% which are divided into energies.
weightsAllRaysSeparate = zeros(size(energyIxAllRays));
weightIx = 0;
for j = 1:size(weightsAllRaysSeparate,2)
    for k = 1:find(energyIxAllRays(:,j),1,'last')
        weightIx = weightIx + 1;
        weightsAllRaysSeparate(k,j) = resultGUI.w(weightIx);
    end
end


% plot energy-separated dose profiles summed up over all rays
% sum up weights over all rays for each energy
weightsAllRays = zeros(length(energyIxUsed),1);
for j = 1:length(energyIxUsed)
    weightsAllRays(j) = sum(weightsAllRaysSeparate(energyIxAllRays==energyIxUsed(j)));
end

% calculate weighted dose profiles
unweightedDose = zeros(numEnergiesMax, size(z,2));
weightedDoseAllRays = zeros(numEnergiesMax,size(z,2));
for j = 1:numEnergiesMax
    unweightedDose(j,:) = machine.data(energyIxUsed(j)).Z.doseAPM;
    weightedDoseAllRays(j,:) = weightsAllRays(j)*unweightedDose(j,:);
end

% get physical dose integrated over all rays
doseSumAllRays = sum(sum(resultGUI.RBExDose(:,:,:),1),3);

coords_water = [3:3:3*160]-120;
startTarget = 208.5-120;
endTarget = 271.5-120;

% plot energy-separated dose profiles summed up over all rays
weightedDoseProfilesAllRaysFig = figure;
hold on
for j = 1:size(weightedDoseAllRays,1)
    plot(z(j,:) + offset, weightedDoseAllRays(j,:)/max(weightedDoseAllRays(:)));
end
doseSumAllRaysPlot = plot(coords_water,doseSumAllRays/max(doseSumAllRays(:)),'k--');
targetStartPlot = plot([startTarget startTarget],[0 1],'k-','linewidth',2);
targetEndPlot = plot([endTarget endTarget],[0 1],'k-','linewidth',2);
xlabel('depth [mm]')
ylabel('relative dose')
axis([0 endTarget+20 0 1])
legend([doseSumAllRaysPlot,targetStartPlot], ...
    'sum of dose profiles = SOBP', 'target boundaries', 'location','northwest')


%% save
savefig(myFig,'X:\Masterarbeit\figures\sobp.fig')
matlab2tikz('X:\Masterarbeit\figures\sobp.tex','width','\fwidth')
