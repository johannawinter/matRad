% Analyse falloff and Delta D95 of DD for a specific compbination of 
% breast wall thickness, target size and lung thickness

clear
close all
load PHANTOM_for_falloffs.mat

breastThickness = 70;   % [mm]
targetThickness = 80;   % [mm]
lungGeoThickness = 60; %[5 17 40 60 80 100];	% [mm]
% Pmod = 256;             % [µm]

plotDD = 1;           % true / false

for h = 1:length(lungGeoThickness)
close all

%% adjust target size, water phantom, breast wall and lung (w/o heterogeneity)
% breast wall
cst{3,1} = 2;
cst{3,2} = 'BreastWall';
cst{3,3} = 'OAR';
cst{3,5} = cst{1,5};
cst{3,5}.visibleColor = [.7 .7 0];

A = zeros(250,250,250);
A(:,2:breastThickness/2+1,:) = 1;
ix = find(A > 0);
cst{3,4}{1} = ix;
cst{3,7} = [];

ct.cube{1}(:) = 0;

ct.cube{1}(cst{3,4}{1}) = 1;

% lung
cst{4,1} = 3;
cst{4,2} = 'Lung';
cst{4,3} = 'OAR';
cst{4,5} = cst{1,5};
cst{4,5}.visibleColor = [.5 .5 .5];

A = zeros(250,250,250);
A(:, (breastThickness/2+2) : round(breastThickness/2+2 + lungGeoThickness(h)/2-1), :) = 1;
ix = find(A > 0);
cst{4,4}{1} = ix;

ct.cube{1}(cst{4,4}{1}) = .306;

% target size
x1target = 125 - round(targetThickness/4) + 1;
x2target = 125 + round(targetThickness/4);
y1target = round(breastThickness/2+2 + lungGeoThickness(h)/2-1 + 1);
y2target = round(breastThickness/2+2 + lungGeoThickness(h)/2-1 + 1 + targetThickness/2 - 1);

A = zeros(250,250,250);
A(x1target:x2target, y1target:y2target, x1target:x2target) = 1;
ix = find(A > 0);
cst{2,4}{1} = ix;
cst{2,7} = [];

% water phantom
A = zeros(250,250,250);
A(:, y1target:end, :) = 1;
ix = find(A > 0);
cst{1,4}{1} = ix;
cst{1,7} = [];

ct.cube{1}(cst{1,4}{1}) = 1;


%% optimization without lung heterogeneity
matRad;

resultGUI.physicalDose_noHeterogeneity = resultGUI.physicalDose;

% create depth dose curves (DD)
coords_matRad = 1:1:250;       % [mm*2]
coords_spline = .05:.05:250;       % [mm*2]

dd_0 = resultGUI.physicalDose_noHeterogeneity(round(pln.isoCenter(2)/2), :, round(pln.isoCenter(3)/2));
dd_0_spline = spline(coords_matRad,dd_0,coords_spline);

% calculate falloff 80%-20%
[~,ix_peak] = max(dd_0_spline);

R80_0 = 2 * .8;                         % nominal dose = 2 Gy
[~,ix_R80_behind] = min(abs(dd_0_spline(ix_peak:end)-R80_0));
ix_R80 = ix_R80_behind + ix_peak - 1;
coord_R80 = coords_spline(ix_R80);

R20_0 = 2 * .2;
[~,ix_R20_behind] = min(abs(dd_0_spline(ix_peak:end)-R20_0));
ix_R20 = ix_R20_behind + ix_peak - 1;
coord_R20 = coords_spline(ix_R20);

z8020(1,1) = 0.0001;                    % thickness of lung tissue [mm]
z8020(1,2) = (coord_R20-coord_R80)*2;	% falloff [mm]

% calculate DeltaD95 [mm]
D95_0 = 2 * .95;
[~,ix_D95_0_behind] = min(abs(dd_0_spline(ix_peak:end)-D95_0));
ix_D95_0 = ix_D95_0_behind + ix_peak - 1;
coord_D95_0 = coords_spline(ix_D95_0);

DeltaD95(1,1) = 0.0001;                 % thickness of lung tissue [mm]
DeltaD95(1,2) = 0;                      % difference to D95 without lung material [mm]


%% add heterogeneity
cst{4,5}.HeterogeneityCorrection = 'Lung';

resultGUI_lung = matRad_calcDoseDirect(ct,stf,pln,cst,resultGUI.w);
resultGUI.physicalDose_Lung = resultGUI_lung.physicalDose;

dd_lung = resultGUI.physicalDose_Lung(round(pln.isoCenter(2)/2), :, round(pln.isoCenter(3)/2));
dd_lung_spline = spline(coords_matRad,dd_lung,coords_spline);

if plotDD
    dd = figure;
    title(['DD: p+ on ' num2str(breastThickness) ' mm breast wall and ' ...
        num2str(lungGeoThickness(h)) ' mm lung, target size ' ...
        num2str(targetThickness) ' mm'])   
    hold on
    plot(coords_matRad*2,dd_0,'ob')
    plot(coords_spline*2,dd_0_spline,'b')
    plot(coords_matRad*2,dd_lung,'or')
    plot(coords_spline*2,dd_lung_spline,'r')
    plot([y1target*2,y1target*2],[0,2.3],'k')
    plot([y2target*2,y2target*2],[0,2.3],'k')
    legend('without heterogeneity','spline','with lung heterogeneity','spline','target boundary','location','northeast')
%     legend('without heterogeneity','with lung heterogeneity','target boundary','location','north')
    xlabel('depth in water [mm]')
    ylabel('dose [Gy]')
    axis([0 500 0 max(dd_0_spline)+.2])
    grid on, grid minor
    box on
    
%     savefig(dd,['C:\Matlab\Analysis phantom degradation\fallOff_D95_accordingToSigmaAnalysis\DD_breastThickness' num2str(breastThickness) '_targetThickness_' num2str(targetThickness) '_lungThickness_' num2str(lungGeoThickness(h)) '.fig']) 
     savefig(dd,['DD_breastThickness' num2str(breastThickness) '_targetThickness_' num2str(targetThickness) '_lungThickness_' num2str(lungGeoThickness(h)) '.fig']) 
end

% calculate falloff 80%-20% [mm]
[~,ix_peak] = max(dd_lung_spline);

R80_lung = 2 * .8;
[~,ix_R80_behind] = min(abs(dd_lung_spline(ix_peak:end)-R80_lung));
ix_R80 = ix_R80_behind + ix_peak - 1;
coord_R80 = coords_spline(ix_R80);

R20_lung = 2 * .2;
[~,ix_R20_behind] = min(abs(dd_lung_spline(ix_peak:end)-R20_lung));
ix_R20 = ix_R20_behind + ix_peak - 1;
coord_R20 = coords_spline(ix_R20);

z8020(2,1) = lungGeoThickness(h);
z8020(2,2) = (coord_R20-coord_R80)*2;


% calculate DeltaD95 [mm]
D95_lung = 2 * .95;
[~,ix_D95_lung_behind] = min(abs(dd_lung_spline(ix_peak:end)-D95_lung));
ix_D95_lung = ix_D95_lung_behind + ix_peak - 1;
coord_D95_lung = coords_spline(ix_D95_lung);

DeltaD95_lung = coord_D95_0 - coord_D95_lung;
DeltaD95(2,1) = lungGeoThickness(h);
DeltaD95(2,2) = DeltaD95_lung*2;


%% show absolute dose differences in matRadGUI
absDiffCube = resultGUI.physicalDose_Lung - resultGUI.physicalDose_noHeterogeneity;
resultGUI.physicalDose_absDiffHeteroHomo = absDiffCube;


%% save results
% save(['C:\Matlab\Analysis phantom degradation\fallOff_D95_accordingToSigmaAnalysis\results_breastThickness_' ...
%     num2str(breastThickness) '_targetThickness_' num2str(targetThickness) '_lungThickness_' num2str(lungGeoThickness(h))],...
%     'ct','cst','pln','resultGUI','stf','z8020','DeltaD95','-v7.3')
save(['results_breastThickness_' ...
    num2str(breastThickness) '_targetThickness_' num2str(targetThickness) '_lungThickness_' num2str(lungGeoThickness(h))],...
    'ct','cst','pln','resultGUI','stf','z8020','DeltaD95','-v7.3')


%% include DVH comparison no lung vs. lung
dvh_0 = matRad_calcDVH(cst,resultGUI.physicalDose_noHeterogeneity,'cum');
qi_0  = matRad_calcQualityIndicators(cst,pln,resultGUI.physicalDose_noHeterogeneity);

dvh_lung = matRad_calcDVH(cst,resultGUI.physicalDose_Lung,'cum');
qi_lung = matRad_calcQualityIndicators(cst,pln,resultGUI.physicalDose_Lung);


dvhTitle = 'DVH comparison - solid line: no heterogeneity, dotted line: heterogeneous lung';
dvhFig = figure('Name','DVH comparison','Color',[0.5 0.5 0.5],'Position',([300 300 800 600]));
hold on
subplot(311)
matRad_showDVH(dvh_0,cst,pln,1,dvhTitle)
matRad_showDVH(dvh_lung,cst,pln,2)

subplot(312)
matRad_showQualityIndicators(qi_0);

subplot(313)
matRad_showQualityIndicators(qi_lung);

% savefig(dvhFig,['C:\Matlab\Analysis phantom degradation\fallOff_D95_accordingToSigmaAnalysis\DVH_breastThickness_' num2str(breastThickness) '_targetThickness_' num2str(targetThickness) '_lungThickness_' num2str(lungGeoThickness(h)) '.fig'])
savefig(dvhFig,['DVH_breastThickness_' num2str(breastThickness) '_targetThickness_' num2str(targetThickness) '_lungThickness_' num2str(lungGeoThickness(h)) '.fig'])

end



%% compare falloffs and DeltaD95
% clear
% breastThickness = 30;   % [mm]
% targetThickness = 40;   % [mm]
% lungGeoThickness = [2 7 20 40 60 80 100];	% [mm]
% 
% % breastThickness = 70;   % [mm]
% % targetThickness = 40;   % [mm]
% % lungGeoThickness = [5 17 40 60 80 100];	% [mm]
% 
% for h = 1:length(lungGeoThickness)
%     results(h) = load(['C:\Matlab\Analysis phantom degradation\fallOff_D95_accordingToSigmaAnalysis\breast' ...
%         num2str(breastThickness) '_target' num2str(targetThickness) '\results_breastThickness_' ...
%         num2str(breastThickness) '_targetThickness_' num2str(targetThickness) '_lungThickness_' ...
%         num2str(lungGeoThickness(h)) '.mat']);
% end
% 
% dc = figure;
% hold on
% title(['DeltaD95 for p+ on ' num2str(breastThickness) ' mm breast wall and ' num2str(targetThickness) ' mm target size'])
% for h = 1:length(lungGeoThickness)
%     plot(results(h).DeltaD95(:,1), results(h).DeltaD95(:,2),'x')
% end
% xlabel('z_g_e_o lung [mm]')
% ylabel('DeltaD95 [mm]')
% % legend(' ','location','northwest')
% grid on
% box on
% 
% fc = figure;
% hold on
% title(['falloff widening by swtiching on heterogeneity in the lung for p+ on ' num2str(breastThickness) ' mm breast wall and ' num2str(targetThickness) ' mm target size'])
% for h = 1:length(lungGeoThickness)
%     plot(results(h).z8020(:,1), results(h).z8020(:,2)-results(h).z8020(1,2),'x')
% end
% xlabel('z_g_e_o lung [mm]')
% ylabel('80% - 20% [mm]')
% % legend(' ','location','northwest')
% grid on
% box on
