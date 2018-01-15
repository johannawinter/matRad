% Analysis of the relation between modulation and falloff / DeltaD95
% therefore analysis of the weights from optimization

clear
close all

breastThickness = 30;       % [mm]
targetThickness = 40;       % [mm]
lungGeoThickness = [2 7 20 30 40 50 60 70 80 90 100]; % [mm]
% breastThickness = 30;
% targetThickness = 80;
% lungGeoThickness = [5 10 17 30 40 50 60 70 80 90]; % 100]; % [mm]
% breastThickness = 70;       % [mm]
% targetThickness = 40;       % [mm]
% lungGeoThickness = [5 10 17 30 40 50 60 70 80 90 100]; % [mm]
% breastThickness = 70;
% targetThickness = 80;
% lungGeoThickness = [5 10 20 31 40]; % 50 60 70 80 90 100]; % [mm]

rho = .306;                 % [g/cm^3] lung density

for i = 1:length(lungGeoThickness)
    close all
    
    % load results from optimized phantom treatment plan
    load(['C:\Matlab\Analysis phantom degradation\fallOff_D95_accordingToSigmaAnalysis\breast' ...
        num2str(breastThickness) '_target' num2str(targetThickness) ...
        '\results_breastThickness_' num2str(breastThickness) ...
        '_targetThickness_' num2str(targetThickness) '_lungThickness_' ...
        num2str(lungGeoThickness(i)) '.mat'])
    
    
    %% find energies for each ray and separate weight vector into rays and energies
    % load machine
    fileName = [pln.radiationMode '_' pln.machine];
    load(['C:\Matlab\matrad\' fileName]);
    
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
    
    
    %% plot weight over energy for central ray
    centralRayNum = round(size(stf.ray,2)/2);
    weightsCentralRay = weightsAllRaysSeparate(:,centralRayNum);
    
    % plot weight over energy in MeV
    weightsCentralRayFig = figure;
    title(['Weights for central ray (# ' num2str(centralRayNum) '), ' ...
        num2str(breastThickness) ' mm breast, ' ...
        num2str(targetThickness) ' mm target, ' num2str(lungGeoThickness(i)) ...
        ' mm geometrical lung'])
    hold on
    plot(stf.ray(centralRayNum).energy,weightsCentralRay,'x-')
    xlabel('energy [MeV]')
    ylabel('weight = # particles*10^6')
    
    
    %% plot sum of weight of all rays over energy
    % sum up weights over all rays for each energy
    weightsAllRays = zeros(length(energyIxUsed),1);
    for j = 1:length(energyIxUsed)
        weightsAllRays(j) = sum(weightsAllRaysSeparate(energyIxAllRays==energyIxUsed(j)));
    end
    
    % plot sum of weight of all rays over energy in MeV
    weightsAllRaysFig = figure;
    title(['Sum of weights over all rays, ' ...
        num2str(breastThickness) ' mm breast, ' ...
        num2str(targetThickness) ' mm target, ' num2str(lungGeoThickness(i)) ...
        ' mm geometrical lung'])
    hold on
    plot(energiesUsed,weightsAllRays,'x-')
    xlabel('energy [MeV]')
    ylabel('sum of weights = # particles*10^6')
    
    
    %% plot energy-separated dose profiles for central ray
    % calculate weighted dose profiles
    unweightedDose = zeros(numEnergiesMax, size(z,2));
    weightedDoseCentralRay = zeros(numEnergiesMax, size(z,2));
    for j = 1:numEnergiesMax
        unweightedDose(j,:) = machine.data(energyIxUsed(j)).Z.doseAPM;
        weightedDoseCentralRay(j,:) = weightsAllRaysSeparate(j,centralRayNum)*unweightedDose(j,:);
    end
    
    % get physical dose, that is sum of dose profiles, on central ray
    centralRay.x = round(pln.isoCenter(2)/2);
    centralRay.z = round(pln.isoCenter(3)/2);
    doseSumCentralRay = resultGUI.physicalDose_noHeterogeneity(centralRay.x, :, centralRay.z);
    coords_matRad = 1:1:250;            % [mm*2]
    
    % calculate water coordinates for doseSum in mm (different from coords_matRad in mm*2!)
    % breast starts at z = 2 --> coords_water starts there --> all coords -2
    % resolution: 2 mm
    if mod(lungGeoThickness(i),2) == 1
        coords_water = [0, coords_matRad(2:breastThickness/2+1)*2-2, ...
            (breastThickness+rho*2): rho*2 :breastThickness + lungGeoThickness(i)*rho, ...
            breastThickness + lungGeoThickness(i)*rho, ...
            breastThickness+lungGeoThickness(i)*rho+2: 2 :(coords_matRad(end-1)*2-lungGeoThickness(i)*(1-rho))];
    else
        coords_water = [0, coords_matRad(2:breastThickness/2+1)*2-2, ...
            (breastThickness+rho*2): rho*2 :breastThickness + lungGeoThickness(i)*rho, ...
            breastThickness+lungGeoThickness(i)*rho+2: 2 :(coords_matRad(end-1)*2-lungGeoThickness(i)*(1-rho))];
    end
    
    % get start and end points of target
    startTarget = breastThickness + rho*lungGeoThickness(i);
    endTarget = startTarget + targetThickness;
    
    
    % plot separate dose profiles, sum of them and target boundaries
    weightedDoseProfilesCentralRayFig = figure;
    title(['Weighted dose profiles for central ray (# ' ...
        num2str(centralRayNum) '), ' num2str(breastThickness) ...
        ' mm breast, ' num2str(targetThickness) ' mm target, ' ...
        num2str(lungGeoThickness(i)) ' mm geometrical lung'])
    hold on
    for j = 1:size(weightedDoseCentralRay,1)
        plot(z(j,:),weightedDoseCentralRay(j,:)/max(weightedDoseCentralRay(:)))
    end
    doseSumPlot = plot(coords_water,doseSumCentralRay/max(doseSumCentralRay(:)),'+k--');
    targetStartPlot = plot([startTarget startTarget],[0 2.5],'k-','linewidth',2);
    targetEndPlot = plot([endTarget endTarget],[0 2.5],'k-','linewidth',2);
    xlabel('z in water [mm]')
    ylabel('relative dose')
    axis([0 endTarget+20 0 1.2])
    legend([doseSumPlot,targetStartPlot], ...
        'sum of dose profiles = SOBP', 'target boundaries','location','northwest')
    
    
    %% plot energy-separated dose profiles summed up over all rays
    % calculate weighted dose profiles
    weightedDoseAllRays = zeros(numEnergiesMax,size(z,2));
    for j = 1:numEnergiesMax
        weightedDoseAllRays(j,:) = weightsAllRays(j)*unweightedDose(j,:);
    end
    
    % get physical dose integrated over all rays
    doseSumAllRays = sum(sum(resultGUI.physicalDose_noHeterogeneity(:,:,:),1),3);
    
    % plot energy-separated dose profiles summed up over all rays
    weightedDoseProfilesAllRaysFig = figure;
    title(['Energy-separated weighted dose profiles summed up over all rays, ' ...
        num2str(breastThickness) ' mm breast, ' num2str(targetThickness) ...
        ' mm target, ' num2str(lungGeoThickness(i)) ' mm geometrical lung'])
    hold on
    for j = 1:size(weightedDoseAllRays,1)
        plot(z(j,:),weightedDoseAllRays(j,:)/max(weightedDoseAllRays(:)));
    end
    doseSumAllRaysPlot = plot(coords_water,doseSumAllRays/max(doseSumAllRays(:)),'+k--');
    targetStartPlot = plot([startTarget startTarget],[0 2.5],'k-','linewidth',2);
    targetEndPlot = plot([endTarget endTarget],[0 2.5],'k-','linewidth',2);
    xlabel('z in water [mm]')
    ylabel('relative dose')
    axis([0 endTarget+20 0 1.2])
    legend([doseSumAllRaysPlot,targetStartPlot], ...
        'sum of dose profiles = SOBP', 'target boundaries', 'location','northwest')
    
    
    %% save results and plots
    save(['C:\Matlab\Analysis phantom degradation\weight_analysis\breast' ...
        num2str(breastThickness) '_target' num2str(targetThickness) ...
        '\results_lung' num2str(lungGeoThickness(i))],...
        'breastThickness','targetThickness','lungGeoThickness','i','rho',...
        'stf','centralRayNum','startTarget','endTarget','z',...
        'energyIxAllRays','energyIxUsed','energiesUsed',...
        'weightsAllRaysSeparate','weightsCentralRay','weightsAllRays',...
        'weightedDoseCentralRay','doseSumCentralRay','coords_water',...
        'weightedDoseAllRays','doseSumAllRays',...
        '-v7.3')
    
    savefig([weightsCentralRayFig,weightsAllRaysFig,weightedDoseProfilesCentralRayFig,weightedDoseProfilesAllRaysFig], ...
        ['C:\Matlab\Analysis phantom degradation\weight_analysis\breast' ...
        num2str(breastThickness) '_target' num2str(targetThickness) ...
        '\figs_lung' num2str(lungGeoThickness(i)) '.fig'])
    
end

