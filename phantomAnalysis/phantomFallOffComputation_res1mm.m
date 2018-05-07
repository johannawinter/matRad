% Analyse falloff and Delta D95 of DD curves for a specific compbination of 
% breast wall thickness, target size and lung thickness

clear
close all
load PHANTOM_for_falloffs.mat

ct = struct('cube',cell(1,1),'resolution',struct('x',1,'y',1,'z',1),...
    'cubeDim',[500,500,500],'numOfCtScen',1);
ct.cube{1,1} = zeros(500,500,500);

breastThickness = 30;   % [mm]
targetThickness = 40;   % [mm]
lungGeoThickness = 50;  % [mm]
% Pmod = 256;         % [µm]

plotDD = 1;           % true / false

for h = 1:length(lungGeoThickness)
close all

%% adjust target size, water phantom, breast wall and lung (w/o heterogeneity)
% breast wall
cst{3,1} = 2;
cst{3,2} = 'ChestWall';
cst{3,3} = 'OAR';
cst{3,5} = cst{1,5};
cst{3,5}.visibleColor = [.7 .7 0];	% olive

A = zeros(500,500,500);
A(:,2:breastThickness+1,:) = 1;
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
cst{4,5}.visibleColor = [.7,0,1];	% purple

A = zeros(500,500,500);
A(:, (breastThickness+2) : round(breastThickness+2 + lungGeoThickness(h)-1), :) = 1;
ix = find(A > 0);
cst{4,4}{1} = ix;

ct.cube{1}(cst{4,4}{1}) = .297;     % formerly .306

% target size
x1target = 250 - round(targetThickness/2) + 1;
x2target = 250 + round(targetThickness/2);
y1target = round(breastThickness+2 + lungGeoThickness(h)-1 + 1);
y2target = round(breastThickness+2 + lungGeoThickness(h)-1 + 1 + targetThickness-1);

A = zeros(500,500,500);
A(x1target:x2target, y1target:y2target, x1target:x2target) = 1;
ix = find(A > 0);
cst{2,4}{1} = ix;
cst{2,7} = [];

% water phantom
A = zeros(500,500,500);
A(:, y1target:end, :) = 1;
ix = find(A > 0);
cst{1,4}{1} = ix;
cst{1,7} = [];

ct.cube{1}(cst{1,4}{1}) = 1;

% add margin around target for optimization
margin = 10;    % [mm]
vMargin.x = margin;
vMargin.y = margin;
vMargin.z = margin;
% assign ones to target voxels
target = cst{2,4}{:};
targetMask = zeros(ct.cubeDim);
targetMask(target) = 1;
% add margin
targetEnlargedVoi = matRad_addMargin(targetMask,cst,ct.resolution,vMargin,true);
targetEnlarged = find(targetEnlargedVoi>0);
% assign enlarged target to cst
cst{5,1} = 4;
cst{5,2} = ['TargetMargin' num2str(margin) 'mm'];
cst{5,3} = 'OAR';
cst{5,4}{1} = targetEnlarged;
cst{5,5} = cst{1,5};      % same as WaterPhantom
cst{5,5}.visibleColor = [.5 .5 .5];
cst{5,6} = cst{1,6};
cst{5,6}.dose = 50;
cst{5,6}.penalty = 40;


%% optimization without lung heterogeneity
matRad
resultGUI = matRad_calcDoseDirect(ct,stf,pln,cst,resultGUI.w);

resultGUI.RBExDose_homo = resultGUI.RBExDose;

% create depth dose curves (DD)
coords_matRad = 1:1:500;        % [mm]
coords_spline = .05:.001:500;	% [mm]

dd_0 = resultGUI.RBExDose_homo(round(pln.propStf.isoCenter(2)), :, round(pln.propStf.isoCenter(3)));
dd_0_spline = spline(coords_matRad,dd_0,coords_spline);

% calculate falloff 80%-20%
[~,ix_peak] = max(dd_0_spline);
R80 = 2 * .8;                 % nominal dose = 2 Gy
R20 = 2 * .2;

coord_R80 = matRad_interp1(dd_0_spline(end:-1:ix_peak)',coords_spline(end:-1:ix_peak)',R80);
coord_R20 = matRad_interp1(dd_0_spline(end:-1:ix_peak)',coords_spline(end:-1:ix_peak)',R20);

% [~,ix_peak_direct] = max(dd_0);
% coords_R80_direct = matRad_interp1(dd_0(end:-1:ix_peak_direct)',coords_matRad(end:-1:ix_peak_direct)',R80_0);

z8020(1,1) = 0.0001;                    % thickness of lung tissue [mm]
z8020(1,2) = (coord_R20-coord_R80)*2;	% falloff [mm]

% calculate DeltaD95 [mm]
D95 = 2 * .95;
coord_D95_0 = matRad_interp1(dd_0_spline(end:-1:ix_peak)',coords_spline(end:-1:ix_peak)',D95);

DeltaD95(1,1) = 0.0001;                 % thickness of lung tissue [mm]
DeltaD95(1,2) = 0;                      % difference to D95 without lung material [mm]


%% add heterogeneity
cst{4,5}.HeterogeneityCorrection = 'Lung';

resultGUI_hetero = matRad_calcDoseDirect(ct,stf,pln,cst,resultGUI.w);
resultGUI.RBExDose_hetero = resultGUI_hetero.RBExDose;

dd_lung = resultGUI.RBExDose_hetero(round(pln.propStf.isoCenter(2)), :, round(pln.propStf.isoCenter(3)));
dd_lung_spline = spline(coords_matRad,dd_lung,coords_spline);

if plotDD
    ddFig = figure;
    title(['DD: p+ on ' num2str(breastThickness) ' mm breast wall and ' ...
        num2str(lungGeoThickness(h)) ' mm lung, target size ' ...
        num2str(targetThickness) ' mm'])   
    hold on
    plot(coords_matRad+2,dd_0,'ob')
    plot(coords_spline+2,dd_0_spline,'b')
    plot(coords_matRad+2,dd_lung,'or')
    plot(coords_spline+2,dd_lung_spline,'r')
    plot([y1target+2,y1target+2],[0,2.3],'k')
    plot([y2target+2,y2target+2],[0,2.3],'k')
    legend('without heterogeneity','spline','with lung heterogeneity','spline','target boundary',...
        'location','northeast')
%     legend('without heterogeneity','with lung heterogeneity','target boundary','location','north')
    xlabel('depth in water [mm]')
    ylabel('dose [Gy]')
    axis([0 500 0 max(dd_0_spline)+.2])
    grid on, grid minor
    box on
    
    savefig(ddFig,...
        ['D:\analyzed matRad data\Analysis phantom degradation\fallOff_D95_bugfix\'...
        'resolution1mm\DD_breastThickness' ...
        num2str(breastThickness) '_targetThickness_' num2str(targetThickness) ...
        '_lungThickness_' num2str(lungGeoThickness(h)) '.fig']) 
%     savefig(dd,['DD_breastThickness' num2str(breastThickness) '_targetThickness_' ...
%     num2str(targetThickness) '_lungThickness_' num2str(lungGeoThickness(h)) '.fig']) 
end

% calculate falloff 80%-20% [mm]
[~,ix_peak] = max(dd_lung_spline);
coord_R80 = matRad_interp1(dd_lung_spline(end:-1:ix_peak)',coords_spline(end:-1:ix_peak)',R80);
coord_R20 = matRad_interp1(dd_lung_spline(end:-1:ix_peak)',coords_spline(end:-1:ix_peak)',R20);

z8020(2,1) = lungGeoThickness(h);
z8020(2,2) = (coord_R20-coord_R80)*2;


% calculate DeltaD95 [mm]
coord_D95_lung = matRad_interp1(dd_lung_spline(end:-1:ix_peak)',coords_spline(end:-1:ix_peak)',D95);
DeltaD95_lung = coord_D95_0 - coord_D95_lung;

DeltaD95(2,1) = lungGeoThickness(h);
DeltaD95(2,2) = DeltaD95_lung*2;


%% show absolute dose differences in matRadGUI
absDiffCube = resultGUI.RBExDose_hetero - resultGUI.RBExDose_homo;
resultGUI.RBExDose_diffHeteroHomo = absDiffCube;


%% save results
save(['D:\analyzed matRad data\Analysis phantom degradation\fallOff_D95_bugfix\'...
    'resolution1mm\results_breastThickness_' ...
    num2str(breastThickness) '_targetThickness_' num2str(targetThickness) ...
    '_lungThickness_' num2str(lungGeoThickness(h))],...
    'ct','cst','dij','pln','resultGUI','stf',...
    'breastThickness','targetThickness','lungGeoThickness','h','coords_matRad','coords_spline',...
    'z8020','DeltaD95','-v7.3')
% save(['results_breastThickness_' ...
%     num2str(breastThickness) '_targetThickness_' num2str(targetThickness) ...
%     '_lungThickness_' num2str(lungGeoThickness(h))],...
%     'ct','cst','pln','resultGUI','stf','z8020','DeltaD95','-v7.3')


%% include DVH comparison no lung vs. lung
dvh_0 = matRad_calcDVH(cst,resultGUI.RBExDose_homo,'cum');
qi_0  = matRad_calcQualityIndicators(cst,pln,resultGUI.RBExDose_homo);

dvh_lung = matRad_calcDVH(cst,resultGUI.RBExDose_hetero,'cum');
qi_lung = matRad_calcQualityIndicators(cst,pln,resultGUI.RBExDose_hetero);


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

savefig(dvhFig,...
    ['D:\analyzed matRad data\Analysis phantom degradation\fallOff_D95_bugfix\'...
    'resolution1mm\DVH_breastThickness_' ...
    num2str(breastThickness) '_targetThickness_' num2str(targetThickness) ...
    '_lungThickness_' num2str(lungGeoThickness(h)) '.fig'])
% savefig(dvhFig,...
%     ['DVH_breastThickness_' num2str(breastThickness) '_targetThickness_' ...
%     num2str(targetThickness) '_lungThickness_' num2str(lungGeoThickness(h)) '.fig'])

end


%% compare falloffs and DeltaD95
% clear
% close all
% 
% % breastThickness = 30;   % [mm]
% % targetThickness = 40;   % [mm]
% % lungGeoThickness = [2 7 20 30 40 50 60 70 80 90 100];	% [mm]
% 
% % breastThickness = 30;   % [mm]
% % targetThickness = 80;   % [mm]
% % lungGeoThickness = [5 10 17 30 40 50 60 70 80 90]; % 100];	% [mm]
% 
% % breastThickness = 70;   % [mm]
% % targetThickness = 40;   % [mm]
% % lungGeoThickness = [5 10 17 30 40 50 60 70 80 90 100];	% [mm]
% 
% breastThickness = 70;   % [mm]
% targetThickness = 80;   % [mm]
% lungGeoThickness = [5 10 20 31 40]; % 50 60 70 80 90 100];	% [mm]
% 
% 
% for h = 1:length(lungGeoThickness)
%     results(h) = load(['C:\Matlab\Analysis phantom degradation\implementationComparison\breast' ...
%         num2str(breastThickness) '_target' num2str(targetThickness) '\results_breastThickness_' ...
%         num2str(breastThickness) '_targetThickness_' num2str(targetThickness) '_lungThickness_' ...
%         num2str(lungGeoThickness(h)) '.mat']);
% end
% 
% dc = figure;
% hold on
% title(['DeltaD95 for p+ on ' num2str(breastThickness) ' mm breast wall and ' ...
%     num2str(targetThickness) ' mm target size'])
% for h = 1:length(lungGeoThickness)
%     plot(results(h).DeltaD95(:,1), results(h).DeltaD95(:,2),'bx')
% end
% xlabel('z_g_e_o lung [mm]')
% ylabel('DeltaD95 [mm]')
% % legend(' ','location','northwest')
% grid on
% box on
% 
% fdc = figure;
% hold on
% title(['falloff widening by swtiching on heterogeneity in the lung for p+ on ' ...
%     num2str(breastThickness) ' mm breast wall and ' num2str(targetThickness) ...
%     ' mm target size'])
% for h = 1:length(lungGeoThickness)
%     plot(results(h).z8020(:,1), results(h).z8020(:,2)-results(h).z8020(1,2),'bx')
% end
% xlabel('z_g_e_o lung [mm]')
% ylabel('80% - 20% [mm]')
% % legend(' ','location','northwest')
% grid on
% box on
% 
% fc = figure;
% hold on
% title(['falloff for p+ on ' num2str(breastThickness) ' mm breast wall and ' ...
%     num2str(targetThickness) ' mm target size'])
% for h = 1:length(lungGeoThickness)
%     plot(results(h).z8020(:,1), results(h).z8020(:,2),'bx')
% end
% xlabel('z_g_e_o lung [mm]')
% ylabel('80% - 20% [mm]')
% % legend(' ','location','northwest')
% grid on
% box on
% 
% 
% savefig(dc,...
%     ['C:\Matlab\Analysis phantom degradation\implementationComparison\DeltaD95Comparison_breastThickness_' ...
%     num2str(breastThickness) '_targetThickness_' num2str(targetThickness) '.fig'])
% savefig(fdc,...
%     ['C:\Matlab\Analysis phantom degradation\implementationComparison\falloffDifferenceComparison_breastThickness_' ...
%     num2str(breastThickness) '_targetThickness_' num2str(targetThickness) '.fig'])
% savefig(fc,...
%     ['C:\Matlab\Analysis phantom degradation\implementationComparison\falloffComparison_breastThickness_' ...
%     num2str(breastThickness) '_targetThickness_' num2str(targetThickness) '.fig'])

