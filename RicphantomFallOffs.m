% Analyse falloff and Delta D95 of DD for different lung thicknesses
% Adjust y-value of target acoordingly to the depth of the BP / beam energy
clear
close all
load RICPHANTOM_none.mat

% create target of 20x20x20 mm^3 - its depth (y) refers to energies
% use (20 mm,) (50 mm,) 60mm, 100 mm, 150 mm, 200 mm, 250 mm depth, water tank starts at y = 151
A = zeros(500,500,500);
A(241:260,391:410,241:260) = 1; % middle of water tank y = 291:310 (150 mm depth)
ix = find(A > 0);
cst{2,4}{1} = ix;
cst{2,5}.visibleColor = [0 1 1];

% geom. lung thickness in mm, dependent on target depth
lungThickness = [2 5:5:140];        % not longer than 140! default: [2 5:5:depth of target (middle)]??


matRad;
targetDepth = pln.isoCenter(1) - 150.5;
resultGUI.physicalDose_0mm = resultGUI.physicalDose;

% create depth dose curves (DD)
coords_matRad = .5:1:350;       % [mm]
coords_spline = .5:.1:350;      % [mm]

dd_0 = resultGUI.physicalDose(251,151:end,251);
dd_0_spline = spline(coords_matRad,dd_0,coords_spline);

% save depth dose curve to cell array
dd.lungThickness(1) = 0;
dd.dose(1,:) = dd_0;
dd.spline(1,:) = dd_0_spline;

h(1) = figure;
hold on
plot(coords_matRad,dd_0,'ob')
plot(coords_spline,dd_0_spline,'b')
title(['DD: p+ target size 2x2x2 cm^3 - 0 mm lung - target depth ' num2str(targetDepth) ' mm'])
xlabel('depth in water [mm]')
ylabel('dose [Gy]')
axis([0 350 0 2.5])
grid on
box on


% calculate falloff 80%-20%
[~,ix_peak] = max(dd_0_spline);

R80_0 = max(dd_0_spline)*.8;
[~,ix_R80_behind] = min(abs(dd_0_spline(ix_peak:end)-R80_0));
ix_R80 = ix_R80_behind + ix_peak - 1;
coord_R80 = coords_spline(ix_R80);

R20_0 = max(dd_0_spline)*.2;
[~,ix_R20_behind] = min(abs(dd_0_spline(ix_peak:end)-R20_0));
ix_R20 = ix_R20_behind + ix_peak - 1;
coord_R20 = coords_spline(ix_R20);

z8020(1,1) = 0.001;                 % thickness of lung tissue [mm]
z8020(1,2) = coord_R20-coord_R80;   % falloff [mm]

% calculate DeltaD95 [mm]
D95_0 = max(dd_0_spline)*.95;
[~,ix_D95_0_behind] = min(abs(dd_0_spline(ix_peak:end)-D95_0));
ix_D95_0 = ix_D95_0_behind + ix_peak - 1;
coord_D95_0 = coords_spline(ix_D95_0);

DeltaD95(1,1) = 0.001;   % thickness of lung tissue [mm]
DeltaD95(1,2) = 0;      % difference to D95 without lung material [mm]


%% add heterogeneous sample
rho = 0.306;    % relative electron density of lung phantom

cst{3,1} = 2;
cst{3,2} = 'Sample';
cst{3,3} = 'OAR';
cst{3,5} = cst{1,5};
cst{3,5}.Priority = 1;
cst{3,5}.visibleColor = [.5 .5 .5];
cst{3,5}.HeterogeneityCorrection = 'Lung';


for i = 1:size(lungThickness,2)
    
    A = zeros(500,500,500);
    A(:,10:lungThickness(i)+10,:) = 1;
    ix = find(A > 0);
    cst{3,4}{1} = ix;
    ct.cube{1}(cst{3,4}{1}) = rho;
    
    resultGUI_lung = matRad_calcDoseDirect(ct,stf,pln,cst,resultGUI.w);
    
    dd_lung = resultGUI_lung.physicalDose(251,151:end,251);
    dd_lung_spline = spline(coords_matRad,dd_lung,coords_spline);
    
    % save depth dose curve to cell array
    dd.lungThickness(i+1,:) = lungThickness(i);
    dd.dose(i+1,:) = dd_lung;
    dd.spline(i+1,:) = dd_lung_spline;
    
    h(i+1) = figure;
    hold on
    plot(coords_matRad,dd_lung,'ob')
    plot(coords_spline,dd_lung_spline,'b')
    title(['DD: p+ target size 2x2x2 cm^3 - ' num2str(lungThickness(i)) ' mm lung - target depth ' num2str(targetDepth) ' mm'])
    xlabel('depth in water [mm]')
    ylabel('dose [Gy]')
    axis([0 350 0 2.5])
    grid on
    box on
    
    
    % calculate falloff 80%-20% [mm]
    [~,ix_peak] = max(dd_lung_spline);
    
    R80_lung = max(dd_lung_spline)*.8;
    [~,ix_R80_behind] = min(abs(dd_lung_spline(ix_peak:end)-R80_lung));
    ix_R80 = ix_R80_behind + ix_peak - 1;
    coord_R80 = coords_spline(ix_R80);
    
    R20_lung = max(dd_lung_spline)*.2;
    [~,ix_R20_behind] = min(abs(dd_lung_spline(ix_peak:end)-R20_lung));
    ix_R20 = ix_R20_behind + ix_peak - 1;
    coord_R20 = coords_spline(ix_R20);
    
    z8020(i+1,1) = lungThickness(i);
    z8020(i+1,2) = coord_R20-coord_R80;
    
    
    % calculate DeltaD95 [mm]
    D95_lung = max(dd_lung_spline)*.95;
    [~,ix_D95_lung_behind] = min(abs(dd_lung_spline(ix_peak:end)-D95_lung));
    ix_D95_lung = ix_D95_lung_behind + ix_peak - 1;
    coord_D95_lung = coords_spline(ix_D95_lung);
    
    WetSample = lungThickness(i)*rho;         % WET of sample that shifts the Bragg curve to the front
    
    DeltaD95_lung = coord_D95_0 - coord_D95_lung - WetSample;
    DeltaD95(i+1,1) = lungThickness(i);
    DeltaD95(i+1,2) = DeltaD95_lung;
    
end

save(['C:\Matlab\Analysis Ricphantom\fallOff_D95\P750\dd_targetDepth_' num2str(targetDepth)],'dd','-v7.3');
savefig(h,['C:\Matlab\Analysis Ricphantom\fallOff_D95\P750\dd_targetDepth_' num2str(targetDepth) '.fig']);
save(['C:\Matlab\Analysis Ricphantom\fallOff_D95\P750\falloff_targetDepth_' num2str(targetDepth)],'z8020','-v7.3');
save(['C:\Matlab\Analysis Ricphantom\fallOff_D95\P750\DeltaD95_targetDepth_' num2str(targetDepth)],'DeltaD95','-v7.3');


%% plot falloff for various lung thicknesses

% fit falloff dependence to different models
fitFalloffLinear = fit(z8020(:,1),z8020(:,2),'poly1');
fitFalloffSqrt_temp = [sqrt(z8020(:,1)) ones(size(z8020(:,2)))] \ z8020(:,2);
fitFalloffSqrt = fitFalloffSqrt_temp(1)*z8020(:,1).^.5 + fitFalloffSqrt_temp(2); % y = K+x^.5+c
fitFalloffPower = fit(z8020(:,1),z8020(:,2),'power2');      % y = a*x^b+c
abc = coeffvalues(fitFalloffPower);
a = abc(1); b = abc(2); c = abc(3);

% plot z8020
f = figure;
hold on
title(['falloff dependence on thickness of traversed lung material - target depth ' num2str(targetDepth) ' mm'])
plot(z8020(:,1),z8020(:,2),'bo')
plot(fitFalloffLinear,'g')
plot(z8020(:,1),fitFalloffSqrt,'c')
plot(fitFalloffPower,'r')
xlabel('z_g_e_o lung [mm]')
ylabel('80% - 20% [mm]')
xticks(0:10:max(lungThickness))
grid on
box on
legend('simulation data','linear fit','square root fit',['power fit with a = ' num2str(a) ', b = ' num2str(b) ', c = ' num2str(c)],'location','northwest')


% fit falloff dependence to different models
fitDeltaD95Linear = fit(DeltaD95(:,1),DeltaD95(:,2),'poly1');
fitDeltaD95Sqrt_temp = [sqrt(DeltaD95(:,1)) ones(size(DeltaD95(:,2)))] \ DeltaD95(:,2);
fitDeltaD95Sqrt = fitFalloffSqrt_temp(1)*DeltaD95(:,1).^.5 + fitDeltaD95Sqrt_temp(2); % y = K+x^.5+c
fitDeltaD95Power = fit(DeltaD95(:,1),DeltaD95(:,2),'power2');      % y = a*x^b+c
abc = coeffvalues(fitDeltaD95Power);
a = abc(1); b = abc(2); c = abc(3);

% plot DeltaD95
d = figure;
hold on
title(['Delta D95 dependence on thickness of traversed lung material - target depth ' num2str(targetDepth) ' mm'])
plot(DeltaD95(:,1),DeltaD95(:,2),'bo')
plot(fitDeltaD95Linear,'g')
plot(DeltaD95(:,1),fitDeltaD95Sqrt,'c')
plot(fitDeltaD95Power,'r')
xlabel('z_g_e_o lung [mm]')
ylabel('Delta D95 [mm]')
xticks(0:10:max(lungThickness))
grid on
box on
legend('simulation data','linear fit','square root fit',['power fit with a = ' num2str(a) ', b = ' num2str(b) ', c = ' num2str(c)],'location','northwest')


savefig(f,['C:\Matlab\Analysis Ricphantom\fallOff_D95\P750\falloff_targetDepth_' num2str(targetDepth) '.fig'])
savefig(d,['C:\Matlab\Analysis Ricphantom\fallOff_D95\P750\deltaD95_targetDepth_' num2str(targetDepth) '.fig'])

