

load('protons_HIT_APM')

plotFits = 1;       % 1: true / 0: false
saveFigs = 0;       % 1: true / 0: false

PmodMin = 50;      % modulation power
PmodPhan = 256;
PmodMax = 750;

breastThickness = [70 110 150]; % 70 / 110 / 150 / 200    % [mm]
lungThickness = 100;                                % [mm]
geoThickness = breastThickness + lungThickness;

rho = .297; % formerly: .306         % relative electron density of lung phantom

% Gaussian error function: erf(x) =  -0.5 .* erf( (x-mu)/(sqrt(2)*sigma) ) + 0.5;
% coeffErrorFun(1) = mu; coeffErrorFun(2) = sigma;
gaussErrorFitFunction = @(coeffErrorFun,x)...
    -0.5 .* erf( (x-coeffErrorFun(1))/(sqrt(2)*coeffErrorFun(2)) ) + 0.5;


%% calculate sigma range straggeling by fitting an error function to base data
% lung part
for i = 1:length(breastThickness)
   zGeoTotal(i,1:(geoThickness(i)-breastThickness(i))/5+1) = linspace(breastThickness(i), geoThickness(i), lungThickness/5+1); 
end
zGeoLung = zGeoTotal(1,:) - breastThickness(1);	% same for all breast thicknesses as they cancel
for i = 1:length(breastThickness)
    wetTotal(i,:) = (zGeoTotal(i,:) - breastThickness(i)) * rho + breastThickness(i);
end

%%%%%% ....


%% theroetical approach using [Mev], [cm]!!!
% according to Bortfeld_1997_Med.Phys.24_12: An analytical approximation of 
% the Bragg curve for therapeutic proton beams (eq. 17 / B5)

% breast (water) part
alphaWater = 2.2e-3;        % [cm MeV^(-p)] for p+ in water, ~ sqrt(Aeff) [14], ~ 1/mass density
p = 1.77;                   % for protons between 10 and 250 MeV, p = 1.8 according to sources [14],[15] 
alphaPrimeWater = 0.087;    % [MeV^2/cm]
% lung part
alphaLung = alphaWater / 0.317;             % alpha ~ 1/mass density, mass density phantom = 0.317 g/cm^3
alphaPrimeLung = alphaPrimeWater * rho;     % alphaPrime ~ electron density


ixBeamEnergy = NaN(size(wetTotal));
for i = 1:length(breastThickness)
    for j = 1:nnz(wetTotal(i,:)>0)
        ixBeamEnergyTemp = zeros(size(machine.data,2),1);
        for k = 1:size(machine.data,2)
            ixBeamEnergyTemp(k) = abs(machine.data(k).peakPos + machine.data(k).offset - wetTotal(i,j));
        end
        [~,ixBeamEnergy(i,j)] = min(ixBeamEnergyTemp);
    end
end

% set beam energies
E0 = NaN(size(ixBeamEnergy));
for i = 1:length(breastThickness)
    for j = 1:nnz(wetTotal(i,:)>0)
        E0(i,j) = machine.data(ixBeamEnergy(:)).energy;
    end
end

% determine range in lung (geom)
R0Lung = NaN(size(wetTotal)); 
for i = 1:length(breastThickness)
    for j = 1:nnz(wetTotal(i,:)>0)
        R0Lung(i,j) = alphaLung/alphaWater * ...
            ( (alphaWater*E0(i,j)^p)^(1/p) - (alphaWater*E0(i,j).^p - breastThickness(i)).^(1/p) )^p;
    end
end

