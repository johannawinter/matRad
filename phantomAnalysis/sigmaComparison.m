% Comparison of sigma due to range straggeling and due to heterogeneity over z_geo
% with range straggling derived from base data and from Bortfeld's theory; 
% theoretical straggling sigma used for plotting; 
% fit theoretical range straggling sigma; 
% get intersections between straggling and heterogeneity sigma;
% test fit to base data for water by plotting line dose profiles

clear

load('protons_HIT_APM')

plotFits = 0;       % 1: true / 0: false
saveFigs = 1;       % 1: true / 0: false

PmodMin = 150;      % modulation power
PmodPhan = 256;
PmodMax = 750;

breastThickness = [70 110 150]; % 70 / 110 / 150 / 200    % [mm]
lungThickness = 100;                                % [mm]
geoThickness = breastThickness + lungThickness;

rho = .306;         % relative electron density of lung phantom

% Gaussian error function: erf(x) =  -0.5 .* erf( (x-mu)/(sqrt(2)*sigma) ) + 0.5;
% coeffErrorFun(1) = mu; coeffErrorFun(2) = sigma;
gaussErrorFitFunction = @(coeffErrorFun,x)...
    -0.5 .* erf( (x-coeffErrorFun(1))/(sqrt(2)*coeffErrorFun(2)) ) + 0.5;


%% calculate sigma lung for Pmod1 and Pmod2
zGeoLungHetero = linspace(0, lungThickness, lungThickness+1);
wetLungHetero = zGeoLungHetero*rho;

sigmaLungMin = sqrt(matRad_getHeterogeneityCorrSigmaSq(wetLungHetero,PmodMin));
sigmaLungPhan = sqrt(matRad_getHeterogeneityCorrSigmaSq(wetLungHetero,PmodPhan));
sigmaLungMax = sqrt(matRad_getHeterogeneityCorrSigmaSq(wetLungHetero,PmodMax));


%% calculate sigma range straggeling by fitting an error function to base data
% lung part
for i = 1:length(breastThickness)
   zGeoTotal(i,1:(geoThickness(i)-breastThickness(i))/5+1) = linspace(breastThickness(i), geoThickness(i), (geoThickness(i)-breastThickness(i))/5+1); 
end
zGeoLung = zGeoTotal(1,:) - breastThickness(1);	% same for all breast thicknesses as they cancel
for i = 1:length(breastThickness)
    wetTotal(i,:) = (zGeoTotal(i,:) - breastThickness(i)) * rho + breastThickness(i);
end

sigmaRSLung = zeros(size(wetTotal));
for i = 1:length(breastThickness)
    
    if plotFits && exist('fitFig')
        close(fitFig)
    end
    
    for j = 1:nnz(wetTotal(i,:)>0)
        ixBeamEnergyTemp = zeros(size(machine.data,2),1);
        for k = 1:size(machine.data,2)
            ixBeamEnergyTemp(k) = abs(machine.data(k).peakPos - wetTotal(i,j));
        end
        [~,ixBeamEnergy] = min(ixBeamEnergyTemp);
        
        xLinear2 = machine.data(ixBeamEnergy).depths;
        braggCurveLung = machine.data(ixBeamEnergy).Z.doseORG;
        [maxBraggCurveLung,ixMaxBraggCurveLung] = max(braggCurveLung);
        
        % start values for fit
        muStart = machine.data(ixBeamEnergy).peakPos;
        sigmaStart = machine.data(ixBeamEnergy).Z.width(end);
        
        % fit error function
        coeffErrorFun = lsqcurvefit(gaussErrorFitFunction,[muStart,sigmaStart],xLinear2,braggCurveLung./maxBraggCurveLung);
        muFit = coeffErrorFun(1);
        sigmaFit = coeffErrorFun(2);
        
        % show error function fit
        if plotFits
            gaussFit = -0.5 .* erf( (xLinear2-muFit)/(sqrt(2)*sigmaFit) ) + 0.5;
            fitFig(j) = figure;
            hold on
            plot(xLinear2,braggCurveLung)
            plot(xLinear2,gaussFit.*maxBraggCurveLung)
            legend('base data','error function fit','location','west')
            xlabel('WET [mm]')
            ylabel('dose [Gy]')
            hold off
        end
        
        sigmaRSLung(i,j) = sigmaFit;  
    end
    
    if plotFits && saveFigs
        savefig(fitFig,['C:\Matlab\Analysis phantom degradation\sigma_analysis\fitFigures_breast' num2str(breastThickness(i)) 'mm.fig'])
    end
    
end


%% add sigma range shifter if range shifter used
sigmaRangeShifter = 0;
if exist('stf','var')
    for j = 1:size(stf.ray,2)
        for k = 1:size(stf.ray(j).rangeShifter,2)
            if any(stf.ray(j).rangeShifter(k).eqThickness)
                sigmaRangeShifter = 0;
                warning('Sigma for range shifter not yet implemented.')
            end
        end
    end
end

sigmaRSLungTotal = sqrt(sigmaRSLung.^2 + sigmaRangeShifter.^2);



%% theroetical approach using [Mev], [cm]!!!
% according to Bortfeld_1997_Med.Phys.24_12: An analytical approximation of 
% the Bragg curve for therapeutic proton beams (eq. 17 / B5)

% breast (water) part
alphaWater = 2.2e-3;        % [cm MeV^(-p)] for p+ in water, ~ sqrt(Aeff) [14], ~ 1/mass density
p = 1.77;                   % for protons between 10 and 250 MeV, p = 1.8 according to sources [14],[15] 
alphaPrimeWater = 0.087;    % [MeV^2/cm]

for i = 1:length(breastThickness)
   zGeoWater(i,1:floor((breastThickness(i)-21)/3+1)) = linspace(21, breastThickness(i), (breastThickness(i)-21)/3+1);	% [mm]
end

R0Breast = zeros(size(zGeoWater));
for i = 1:length(breastThickness)
    for j = 1:nnz(zGeoWater(i,:)>0)
        ixBeamEnergyBreastTemp = zeros(size(machine.data,2),1);
        for k = 1:size(machine.data,2)
            ixBeamEnergyBreastTemp(k) = abs(machine.data(k).peakPos - zGeoWater(i,j));
        end
        [~,ixBeamEnergyBreast] = min(ixBeamEnergyBreastTemp);
        
        % mean range = alpha*E0^p [cm]
        R0Breast(i,j) = alphaWater * machine.data(ixBeamEnergyBreast).energy ^p;
    end
end
    
sigmaTheoWaterSq = alphaPrimeWater * (p^3*alphaWater^(2/p))/(3*p-2) * R0Breast.^(3-2/p);
sigmaTheoWater = sqrt(sigmaTheoWaterSq) *10;        % [mm]

%%% test
% for water, R0 and sigma in cm, eq. 18
sigmaTheoWaterSimple = 0.012*R0Breast.^0.935 *10;     % ex. E = 150 MeV, R0 = 15.64 cm, sigma = 0.16 cm;
%%%


%% lung part (theoretical)
alphaLung = alphaWater / 0.317;             % alpha ~ 1/mass density, mass density phantom = 0.317 g/cm^3
alphaPrimeLung = alphaPrimeWater * rho;     % alphaPrime ~ electron density

R0LungWithoutBreast = zeros(size(wetTotal));
for i = 1:length(breastThickness)
    for j = 1:nnz(wetTotal(i,:)>0)
        ixBeamEnergyLungTemp = zeros(size(machine.data,2),1);
        for k = 1:size(machine.data,2)
            ixBeamEnergyLungTemp(k) = abs(machine.data(k).peakPos - wetTotal(i,j));
        end
        [~,ixBeamEnergyLung] = min(ixBeamEnergyLungTemp);
        
        % mean range = alpha*E0^p [cm];   R0 = alphaLung * E0^p - alphaLung/alphaWater * waterThickness [cm]
        R0LungWithoutBreast(i,j) = alphaLung * machine.data(ixBeamEnergyLung).energy ^p - ...
            alphaLung/alphaWater * breastThickness(i) / 10;
    end
    
    R0Lung(i,:) = R0LungWithoutBreast(i,:) + breastThickness(i) / 10;	% [cm]
end

for i = 1:size(R0LungWithoutBreast,1)
    for j = 1:size(R0LungWithoutBreast,2)
        if R0LungWithoutBreast(i,j) <= 0
            R0LungWithoutBreast(i,j) = 0.0001;
        end
    end
end

%%% correct?
% sigmaTheoRSLungSq = alphaPrimeLung * (p^3*alphaLung^(2/p))/(3*p-2) * R0Lung.^(3-2/p);
% sigmaTheoRSLungSq = alphaPrimeLung * (p^3*alphaLung^(2/p))/(3*p-2) * R0LungWithoutBreast.^(3-2/p) + sigmaTheoWaterSq(end);
% sigmaTheoRSLung = sqrt(sigmaTheoRSLungSq) *10;

sigmaTheoRSLungSq = alphaPrimeWater * (p^3*alphaWater^(2/p))/(3*p-2) * R0Lung.^(3-2/p);
sigmaTheoRSLungWet = sqrt(sigmaTheoRSLungSq) *10;                       % [mm]
for i = 1:length(breastThickness)
    sigmaTheoRSLung(i,:) = sigmaTheoRSLungWet(i,:)*rho - ...
        sigmaTheoWater(i,find(sigmaTheoWater(i,:),1,'last'))*rho + ...
        sigmaTheoWater(i,find(sigmaTheoWater(i,:),1,'last'));           % [mm]
end
%%%


%% fit theoretical range straggling sigma
powerFitFun = @(x,xdata)x(1)*xdata.^x(2)+x(3);   % ydata = x(1)*xdata^.x(2) + x(3)
% for i = 1:length(breastThickness)
%     coeffFitSigmaRSLung(i,:) = lsqcurvefit(powerFitFun,[0.0001, 1, 0],zGeoLung,sigmaRSLung(i,:))
% end

for i = 1:length(breastThickness)
    coeffFitSigmaTheoRSLung(i,:) = lsqcurvefit(powerFitFun,[0.0001, 1, 0],zGeoLung,sigmaTheoRSLung(i,:))
end


%% get intersections
for i = 1:length(breastThickness)
    try
        [xIntersectMin(i),yIntersectMin(i)] = intersections(zGeoLung,sigmaTheoRSLung(i,:), zGeoLungHetero,sigmaLungMin,1);
    catch
        xIntersectMin(i) = []; yIntersectMin(i) = [];
    end
    [xIntersectMax(i),yIntersectMax(i)] = intersections(zGeoLung,sigmaTheoRSLung(i,:), zGeoLungHetero,sigmaLungMax,1);
end


%% combination of heterogeneity effect and range straggeling
sigmaLungMinCombi = sqrt(matRad_getHeterogeneityCorrSigmaSq(zGeoLung*rho,PmodMin));
sigmaLungMaxCombi = sqrt(matRad_getHeterogeneityCorrSigmaSq(zGeoLung*rho,PmodMax));

sigmaCombiMin = sqrt(sigmaLungMinCombi.^2 + sigmaTheoRSLung(2,:).^2);
sigmaCombiMax = sqrt(sigmaLungMaxCombi.^2 + sigmaTheoRSLung(2,:).^2);


%% plot sigmas - with breast in front of lung
sigmaFig = figure('Name','sigma comparison');
hold on
title('\sigma comparison - heterogeneity effect and range straggeling (RS) in lung')
xlabel('z_g_e_o lung [mm]')
ylabel('\sigma [mm] in water')
grid on, grid minor
box on
plot(zGeoLungHetero, sigmaLungMin, 'b--')
plot(zGeoLungHetero, sigmaLungPhan, 'b-')
plot(zGeoLungHetero, sigmaLungMax, 'b--')

plot(zGeoLung, sigmaRSLung(1,:), 'o','color',[.7,0,1])     % purple
plot(zGeoLung, sigmaRSLung(2,:), 'ro')
plot(zGeoLung, sigmaRSLung(3,:), 'o','color',[1,.7,0])     % orange
% plot(zGeoLung, sigmaRSLung(4,:), 'o','color',[.5,.5,0])	% olive

plot(zGeoLung, sigmaTheoRSLung(1,:),'*','color',[.7,0,1])	% purple
plot(zGeoLung, sigmaTheoRSLung(2,:),'r*')
plot(zGeoLung, sigmaTheoRSLung(3,:),'*','color',[1,.7,0])	% orange
% plot(zGeoLung, sigmaTheoRSLung(4,:),'*','color',[.5,.5,0])	% olive

plot(zGeoLung, powerFitFun(coeffFitSigmaTheoRSLung(1,:),zGeoLung),'color',[.7,0,1])     % purple
plot(zGeoLung, powerFitFun(coeffFitSigmaTheoRSLung(2,:),zGeoLung),'r')
plot(zGeoLung, powerFitFun(coeffFitSigmaTheoRSLung(3,:),zGeoLung),'color',[1,.7,0])     % orange
% plot(zGeoLung, powerFitFun(coeffFitSigmaTheoRSLung(4,:),zGeoLung),'color',[.5,.5,0])	% olive

plot(zGeoLung, sigmaCombiMin ,'k:')
plot(zGeoLung, sigmaCombiMax ,'k:')

plot(xIntersectMin(1),yIntersectMin(1),'kx','Linewidth',2,'MarkerSize',16)
plot(xIntersectMax(1),yIntersectMax(1),'kx','Linewidth',2,'MarkerSize',16)
plot(xIntersectMin(2),yIntersectMin(2),'kx','Linewidth',2,'MarkerSize',16)
plot(xIntersectMax(2),yIntersectMax(2),'kx','Linewidth',2,'MarkerSize',16)
plot(xIntersectMin(3),yIntersectMin(3),'kx','Linewidth',2,'MarkerSize',16)
plot(xIntersectMax(3),yIntersectMax(3),'kx','Linewidth',2,'MarkerSize',16)
% plot(xIntersectMin(4),yIntersectMin(4),'kx','Linewidth',2,'MarkerSize',16)
% plot(xIntersectMax(4),yIntersectMax(4),'kx','Linewidth',2,'MarkerSize',16)

ax1 = gca;
set(ax1,'XMinorTick','on')
ax2 = axes('Position',[ax1.Position(1),ax1.Position(2)-.07,ax1.Position(3),1e-20],...
    'color','none','xlim',[0 max(wetLungHetero)],'XMinorTick','on');
xlabel(ax2,'WET lung [mm]')

legend(ax1,'heterogeneity - P_m_o_d_,_m_i_n = 150 µm',...
    'heterogeneity - P_m_o_d_,_p_h_a_n_t_o_m = 256 µm','heterogeneity - P_m_o_d_,_m_a_x = 750 µm',...
    ['RS - water thickness ' num2str(breastThickness(1)) ' mm'],...
    ['RS - water thickness ' num2str(breastThickness(2)) ' mm'],...
    ['RS - water thickness ' num2str(breastThickness(3)) ' mm'],... %     ['RS - water thickness ' num2str(breastThickness(4)) ' mm'],...
    ['RS (theoretical)  - water thickness ' num2str(breastThickness(1)) ' mm'],...
    ['RS (theoretical)  - water thickness ' num2str(breastThickness(2)) ' mm'],...
    ['RS (theoretical)  - water thickness ' num2str(breastThickness(3)) ' mm'],... %     ['RS (theoretical)  - water thickness ' num2str(breastThickness(4)) ' mm'],...
    ['power fit with a = ' num2str(coeffFitSigmaTheoRSLung(1,1),2) ', b = ' num2str(coeffFitSigmaTheoRSLung(1,2),3)],...
    ['power fit with a = ' num2str(coeffFitSigmaTheoRSLung(2,1),2) ', b = ' num2str(coeffFitSigmaTheoRSLung(2,2),3)],...
    ['power fit with a = ' num2str(coeffFitSigmaTheoRSLung(3,1),2) ', b = ' num2str(coeffFitSigmaTheoRSLung(3,2),3)],... %     ['power fit with a = ' num2str(coeffFitSigmaTheoRSLung(4,1),2) ', b = ' num2str(coeffFitSigmaTheoRSLung(4,2),3)],...
    'combination RS (70 mm) with P_m_o_d_,_m_i_n','combination RS (70 mm) with P_m_o_d_,_m_a_x',...
    'location','northwest')

if saveFigs
    savefig(sigmaFig,['C:\Matlab\Analysis phantom degradation\sigma_analysis\sigmaComparison_' num2str(length(breastThickness)) 'breastThicknesses.fig'])
end

%% test for water
%% calculate sigma range straggeling by fitting an error function to base data

% breast part - 21 mm corresponds to lowest proton energy
% for i = 1:length(breastThickness)
%    zGeoWater(i,1:floor((breastThickness(i)-21)/3+1)) = linspace(21, breastThickness(i), (breastThickness(i)-21)/3+1);	% [mm]
% end

sigmaWater = zeros(size(zGeoWater));

for i = 1:length(breastThickness)
    
    for j = 1:nnz(zGeoWater(i,:)>0)
        ixBeamEnergyTemp = zeros(size(machine.data,2),1);
        for k = 1:size(machine.data,2)
            ixBeamEnergyTemp(k) = abs(machine.data(k).peakPos - zGeoWater(i,j));
        end
        [~,ixBeamEnergyTest] = min(ixBeamEnergyTemp);
        ixBeamEnergy = ixBeamEnergyTest;
        
        xWater1 = machine.data(ixBeamEnergy).depths;
        braggCurveWater = machine.data(ixBeamEnergy).Z.doseORG;
        [maxBraggCurveWater,ixMaxBraggCurveWater] = max(braggCurveWater);
        
        % start values for fit
        muStart = machine.data(ixBeamEnergy).peakPos;
        sigmaStart = machine.data(ixBeamEnergy).Z.width(end);
        
        % fit error function
        coeffErrorFun = lsqcurvefit(gaussErrorFitFunction,[muStart,sigmaStart],xWater1,braggCurveWater./maxBraggCurveWater);
        muFit = coeffErrorFun(1);
        sigmaFit = coeffErrorFun(2);
        
        % show error function fit
        if plotFits
            gaussFit = -0.5 .* erf( (xWater1-muFit)/(sqrt(2)*sigmaFit) ) + 0.5;
            fitFigWater(j) = figure;
            hold on
            plot(xWater1,braggCurveWater)
            plot(xWater1,gaussFit.*maxBraggCurveWater)
            legend('base data','error function fit','location','west')
            xlabel('z_g_e_o [mm]')
            ylabel('dose [Gy]')
            hold off
        end
        
        sigmaWater(i,j) = sigmaFit;
    end
    
    if plotFits && i == 3 && saveFigs
        savefig(fitFigWater,['C:\Matlab\Analysis phantom degradation\sigma_analysis\fitFiguresWaterTest_breast' num2str(breastThickness(i)) 'mm.fig'])
    end
end

% add sigma from range shifter if necessary
sigmaRSWaterTotal = sqrt(sigmaWater.^2 + sigmaRangeShifter.^2);


% fit
coeffFitSigmaTheoWater = zeros(length(breastThickness),3);
for i = 1:length(breastThickness)
    coeffFitSigmaTheoWater(i,:) = lsqcurvefit(powerFitFun,[0.1, 1, 0],...
        zGeoWater(i,find(sigmaTheoWater(i,:),1,'last')),...
        sigmaTheoWater(i,find(sigmaTheoWater(i,:),1,'last')))
end

coeffFitSigmaWater = zeros(length(breastThickness),3);
for i = 1:length(breastThickness)
    coeffFitSigmaWater(i,:) = lsqcurvefit(powerFitFun,[0.1, 1, 0],...
        zGeoWater(i,find(sigmaWater(i,:),1,'last')),...
        sigmaWater(i,find(sigmaWater(i,:),1,'last')))
end

%% plot water
waterTest = figure;
hold on
xlabel('z_g_e_o water [mm]')
ylabel('\sigma [mm] in water')
ylim([0,5])
grid on, grid minor
box on
plot(zGeoWater(1,find(zGeoWater(1,:),1,'last')),sigmaWater(1,find(sigmaTheoWater(1,:),1,'last')), 'o','color',[.7,0,1]) % purple
plot(zGeoWater(2,find(zGeoWater(2,:),1,'last')),sigmaWater(2,find(sigmaTheoWater(2,:),1,'last')), 'ro')
plot(zGeoWater(3,find(zGeoWater(3,:),1,'last')),sigmaWater(3,find(sigmaTheoWater(3,:),1,'last')), 'o','color',[1,.7,0]) % orange
plot(zGeoWater(4,end),sigmaWater(4,end), 'o','color',[.5,.5,0]) % olive

plot(zGeoWater(1,find(zGeoWater(1,:),1,'last')),sigmaTheoWater(1,find(sigmaTheoWater(1,:),1,'last')), '*','color',[.7,0,1])
plot(zGeoWater(2,find(zGeoWater(2,:),1,'last')),sigmaTheoWater(2,find(sigmaTheoWater(2,:),1,'last')), 'r*')
plot(zGeoWater(3,find(zGeoWater(3,:),1,'last')),sigmaTheoWater(3,find(sigmaTheoWater(3,:),1,'last')), '*','color',[1,.7,0])
plot(zGeoWater(4,end),sigmaTheoWater(4,end), '*','color',[.5,.5,0])
% plot(zGeoLinear1,sigmaTheoWaterSimple,'gx')

plot(zGeoWater(1,1:find(zGeoWater(1,:),1,'last')), powerFitFun(coeffFitSigmaTheoWater(1,:),zGeoWater(1,1:find(sigmaTheoWater(1,:),1,'last'))),'color',[.7,0,1])
plot(zGeoWater(2,1:find(zGeoWater(2,:),1,'last')), powerFitFun(coeffFitSigmaTheoWater(2,:),zGeoWater(2,1:find(sigmaTheoWater(2,:),1,'last'))),'r')
plot(zGeoWater(3,1:find(zGeoWater(3,:),1,'last')), powerFitFun(coeffFitSigmaTheoWater(3,:),zGeoWater(3,1:find(sigmaTheoWater(3,:),1,'last'))),'color',[1,.7,0])
plot(zGeoWater(4,1:end), powerFitFun(coeffFitSigmaTheoWater(4,:),zGeoWater(4,1:end)),'color',[.5,.5,0])

legend('base data','base data','base data','base data',...
    'theoretical','theoretical','theoretical','theoretical',...
    'location','northwest')

if saveFigs
    savefig(waterTest,['C:\Matlab\Analysis phantom degradation\sigma_analysis\waterTest_' num2str(length(breastThickness)) 'breastThicknesses.fig'])
end
