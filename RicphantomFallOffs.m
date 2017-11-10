% Analyse falloff and Delta D95 of DD for different lung thicknesses
% Adjust y-value of target acoordingly to the depth of the BP / beam energy
clear
close all
load RICPHANTOM_none_res2mm.mat

% create target of 60x60x60 mm^3 - its depth (y) refers to energies
% use (20 mm,) (50 mm,) 60mm, 100 mm, 150 mm, 200 mm, 250 mm depth, water tank starts at y = 152

beginTarget = 136;
endTarget = 165;

A = zeros(250,250,250);
A(111:140,beginTarget:endTarget,111:140) = 1; % middle of water tank y = 135:165 (150 mm depth)
ix = find(A > 0);
cst{2,4}{1} = ix;

% geom. lung thickness in mm*2, dependent on target depth
lungThickness = [1:2:21 30:10:70];          % not longer than 70!

matRad;
targetDepth = pln.isoCenter(1) - 151;
resultGUI.physicalDose_0mm = resultGUI.physicalDose;

% create depth dose curves (DD)
coords_matRad = 1:1:175;       % [mm*2]
coords_spline = 1:.05:175;       % [mm*2]

dd_0 = resultGUI.physicalDose_0mm(126,76:end,126);
dd_0_spline = spline(coords_matRad,dd_0,coords_spline);

% save depth dose curve to cell array
dd.lungThickness(1) = 0; 
dd.dose(1,:) = dd_0;
dd.spline(1,:) = dd_0_spline;

h(1) = figure;
hold on
plot(coords_matRad*2,dd_0,'ob')
plot(coords_spline*2,dd_0_spline,'b')
plot([beginTarget*2-150,beginTarget*2-150],[0,2.3],'k')
plot([endTarget*2-150,endTarget*2-150],[0,2.3],'k')
title(['DD: p+ target size 6x6x6 cm^3 - 0 mm lung - target depth ' num2str(targetDepth) ' mm'])
xlabel('depth in water [mm]')
ylabel('dose [Gy]')
axis([0 350 0 2.5])
grid on
box on


% calculate falloff 80%-20%
[~,ix_peak] = max(dd_0_spline);

% R80_0 = max(dd_0_spline)*.8;
R80_0 = 2 * .8;                       % nominal dose = 2 Gy
[~,ix_R80_behind] = min(abs(dd_0_spline(ix_peak:end)-R80_0));
ix_R80 = ix_R80_behind + ix_peak - 11;
coord_R80 = coords_spline(ix_R80);

% R20_0 = max(dd_0_spline)*.2;
R20_0 = 2 * .2;
[~,ix_R20_behind] = min(abs(dd_0_spline(ix_peak:end)-R20_0));
ix_R20 = ix_R20_behind + ix_peak - 11;
coord_R20 = coords_spline(ix_R20);

z8020(1,1) = 0.0001;                 % thickness of lung tissue [mm]
z8020(1,2) = (coord_R20-coord_R80)*2;   % falloff [mm]

% calculate DeltaD95 [mm]
% D95_0 = max(dd_0_spline)*.95;
D95_0 = 2 * .95;
[~,ix_D95_0_behind] = min(abs(dd_0_spline(ix_peak:end)-D95_0));
ix_D95_0 = ix_D95_0_behind + ix_peak - 11;
coord_D95_0 = coords_spline(ix_D95_0);

DeltaD95(1,1) = 0.0001;   % thickness of lung tissue [mm]
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

resultGUI.physicalDoseLung = cell(size(lungThickness,2),1);

for i = 1:size(lungThickness,2)
    WetSample = lungThickness(i)*rho;         % WET of sample that shifts the Bragg curve to the front

    A = zeros(250,250,250);
    A(:,5:lungThickness(i)+5,:) = 1;
    ix = find(A > 0);
    cst{3,4}{1} = ix;
    ct.cube{1}(cst{3,4}{1}) = rho;
    
    resultGUI_lung = matRad_calcDoseDirect(ct,stf,pln,cst,resultGUI.w);
    resultGUI.physicalDoseLung{i} = resultGUI_lung.physicalDose;
    
    dd_lung = resultGUI.physicalDoseLung{i}(126,76:end,126);
    dd_lung_spline = spline(coords_matRad,dd_lung,coords_spline);
    
    % save depth dose curve to cell array
    dd.lungThickness(i+1,:) = lungThickness(i);
    dd.dose(i+1,:) = dd_lung;
    dd.spline(i+1,:) = dd_lung_spline;
    
    h(i+1) = figure;
    hold on
    plot(coords_matRad*2,dd_lung,'ob')
    plot(coords_spline*2,dd_lung_spline,'b')
    plot([beginTarget*2-150-WetSample*2,beginTarget*2-150-WetSample*2],[0,2.3],'k')
    plot([endTarget*2-150-WetSample*2,endTarget*2-150-WetSample*2],[0,2.3],'k')
    title(['DD: p+ target size 6x6x6 cm^3 - ' num2str(lungThickness(i)*2) ' mm lung - target depth ' num2str(targetDepth) ' mm'])
    xlabel('depth in water [mm]')
    ylabel('dose [Gy]')
    axis([0 350 0 2.5])
    grid on
    box on
    
    
    % calculate falloff 80%-20% [mm]
    [~,ix_peak] = max(dd_lung_spline);
    
%     R80_lung = max(dd_lung_spline)*.8;
    R80_lung = 2 * .8;
    [~,ix_R80_behind] = min(abs(dd_lung_spline(ix_peak:end)-R80_lung));
    ix_R80 = ix_R80_behind + ix_peak - 11;
    coord_R80 = coords_spline(ix_R80);
    
%     R20_lung = max(dd_lung_spline)*.2;
    R20_lung = 2 * .2;
    [~,ix_R20_behind] = min(abs(dd_lung_spline(ix_peak:end)-R20_lung));
    ix_R20 = ix_R20_behind + ix_peak - 11;
    coord_R20 = coords_spline(ix_R20);
    
    z8020(i+1,1) = lungThickness(i)*2;
    z8020(i+1,2) = (coord_R20-coord_R80)*2;
    
    
    % calculate DeltaD95 [mm]
%     D95_lung = max(dd_lung_spline)*.95;
    D95_lung = 2 * .95;
    [~,ix_D95_lung_behind] = min(abs(dd_lung_spline(ix_peak:end)-D95_lung));
    ix_D95_lung = ix_D95_lung_behind + ix_peak - 11;
    coord_D95_lung = coords_spline(ix_D95_lung);
        
    DeltaD95_lung = coord_D95_0 - coord_D95_lung - WetSample;
    DeltaD95(i+1,1) = lungThickness(i)*2;
    DeltaD95(i+1,2) = DeltaD95_lung*2;
    
end


% comparison DD and IDD with 0 mm and 140 mm lung 
idd_0 = sum(sum(resultGUI.physicalDose_0mm(:,76:end,:),3),1);
idd_0_spline = spline(coords_matRad,idd_0,coords_spline);
idd_lung = sum(sum(resultGUI.physicalDoseLung{i}(:,76:end,:),3),1);
idd_lung_spline = spline(coords_matRad,idd_lung,coords_spline);

idd(1) = figure;
hold on; 
p1 = plot(coords_matRad*2,dd_0,'ob');
p2 = plot(coords_spline*2,dd_0_spline,'b');
p3 = plot(coords_matRad*2,dd_lung,'or');
p4 = plot(coords_spline*2,dd_lung_spline,'r');
p5 = plot(coords_matRad*2 + rho*lungThickness(i)*2,dd_lung,'oc');
p6 = plot(coords_spline*2 + rho*lungThickness(i)*2,dd_lung_spline,'c');
legend([p1,p3,p5],'0 mm lung',[num2str(lungThickness(i)*2) ' mm lung'],[num2str(lungThickness(i)*2) ' mm lung shifted by WET_l_u_n_g'],'location','southwest')
title(['DD: p+ target size 6x6x6 cm^3 - target depth ' num2str(targetDepth) ' mm'])
xlabel('depth in water [mm]')
ylabel('dose [Gy]')
axis([0 350 0 2.5])
grid on
box on

idd(2) = figure;
hold on; 
p1 = plot(coords_matRad*2,idd_0,'ob');
p2 = plot(coords_spline*2,idd_0_spline,'b');
p3 = plot(coords_matRad*2,idd_lung,'or');
p4 = plot(coords_spline*2,idd_lung_spline,'r');
p5 = plot(coords_matRad*2 + rho*lungThickness(i)*2,idd_lung,'oc');
p6 = plot(coords_spline*2 + rho*lungThickness(i)*2,idd_lung_spline,'c');
legend([p1,p3,p5], '0 mm lung',[num2str(lungThickness(i)*2) ' mm lung'],[num2str(lungThickness(i)*2) ' mm lung shifted by WET_l_u_n_g'],'location','southwest')
title(['IDD: p+ target size 6x6x6 cm^3 - target depth ' num2str(targetDepth) ' mm'])
xlabel('depth in water [mm]')
ylabel('dose [Gy]')
axis([0 350 0 3500])
grid on
box on


save(['C:\Matlab\Analysis phantom degradation\fallOff_D95\resultsFallOffD95_targetDepth_' num2str(targetDepth)],'cst','ct','pln','resultGUI','stf','-v7.3');
save(['C:\Matlab\Analysis phantom degradation\fallOff_D95\dd_targetDepth_' num2str(targetDepth)],'dd','-v7.3');
savefig(h,['C:\Matlab\Analysis phantom degradation\fallOff_D95\dd_targetDepth_' num2str(targetDepth) '.fig']);
savefig(idd,['C:\Matlab\Analysis phantom degradation\fallOff_D95\idd_vs_dd_targetDepth_' num2str(targetDepth) '.fig']);


%% plot falloff for various lung thicknesses

% fit falloff dependence to different models
% fitFalloffLinear = fit(z8020(:,1),z8020(:,2),'poly1');
% fitFalloffSqrt_temp = [sqrt(z8020(:,1)) ones(size(z8020(:,2)))] \ z8020(:,2);
% fitFalloffSqrt = fitFalloffSqrt_temp(1)*z8020(:,1).^.5 + fitFalloffSqrt_temp(2); % y = K+x^.5+c
fitFalloffPower = fit(z8020(:,1),z8020(:,2),'power2');      % y = a*x^b+c
abc = coeffvalues(fitFalloffPower);
a = abc(1); b = abc(2); c = abc(3);

% plot z8020
f = figure;
hold on
title(['falloff dependence on thickness of traversed lung material - target depth ' num2str(targetDepth) ' mm'])
plot(z8020(:,1),z8020(:,2),'bo')
% plot(fitFalloffLinear,'g')
% plot(z8020(:,1),fitFalloffSqrt,'c')
plot(fitFalloffPower,'r')
xlabel('z_g_e_o lung [mm]')
ylabel('80% - 20% [mm]')
xticks(0:10:max(lungThickness)*2)
grid on
box on
% legend('simulation data','linear fit','square root fit',['power fit with a = ' num2str(a) ', b = ' num2str(b) ', c = ' num2str(c)],'location','northwest')
legend('simulation data',['power fit with a = ' num2str(a) ', b = ' num2str(b) ', c = ' num2str(c)],'location','northwest')


% fit DeltaD95 dependence to different models
% fitDeltaD95Linear = fit(DeltaD95(:,1),DeltaD95(:,2),'poly1');
% fitDeltaD95Sqrt_temp = [sqrt(DeltaD95(:,1)) ones(size(DeltaD95(:,2)))] \ DeltaD95(:,2);
% fitDeltaD95Sqrt = fitFalloffSqrt_temp(1)*DeltaD95(:,1).^.5 + fitDeltaD95Sqrt_temp(2); % y = K+x^.5+c
fitDeltaD95Power = fit(DeltaD95(:,1),DeltaD95(:,2),'power2');      % y = a*x^b+c
abc = coeffvalues(fitDeltaD95Power);
a = abc(1); b = abc(2); c = abc(3);

% plot DeltaD95
d = figure;
hold on
title(['Delta D95 dependence on thickness of traversed lung material - target depth ' num2str(targetDepth) ' mm'])
plot(DeltaD95(:,1),DeltaD95(:,2),'bo')
% plot(fitDeltaD95Linear,'g')
% plot(DeltaD95(:,1),fitDeltaD95Sqrt,'c')
plot(fitDeltaD95Power,'r')
xlabel('z_g_e_o lung [mm]')
ylabel('Delta D95 [mm]')
xticks(0:10:max(lungThickness)*2)
grid on
box on
% legend('simulation data','linear fit','square root fit',['power fit with a = ' num2str(a) ', b = ' num2str(b) ', c = ' num2str(c)],'location','northwest')
legend('simulation data',['power fit with a = ' num2str(a) ', b = ' num2str(b) ', c = ' num2str(c)],'location','northwest')


save(['C:\Matlab\Analysis phantom degradation\fallOff_D95\falloff_targetDepth_' num2str(targetDepth)],'z8020','fitFalloffPower','-v7.3');
save(['C:\Matlab\Analysis phantom degradation\fallOff_D95\DeltaD95_targetDepth_' num2str(targetDepth)],'DeltaD95','fitDeltaD95Power','-v7.3');
savefig(f,['C:\Matlab\Analysis phantom degradation\fallOff_D95\falloff_targetDepth_' num2str(targetDepth) '.fig'])
savefig(d,['C:\Matlab\Analysis phantom degradation\fallOff_D95\deltaD95_targetDepth_' num2str(targetDepth) '.fig'])


%% include DVH comparison no lung vs. 38 mm lung
ixForDvh = 10;      % different behaviour / changes for differen lung thicknesses!
lungThicknessForDvh = lungThickness(ixForDvh); % length(lungThickness) = 16

resultGUI_0mm = resultGUI;
resultGUI_0mm.physicalDose = resultGUI.physicalDose_0mm;
cst_0ForDvh = matRad_indicatorWrapper(cst,pln,resultGUI_0mm);

% for lung dose distribution, target must be shifted by WetSample
WetSampleForDvh = lungThicknessForDvh*rho;
cst_lungForDvh = cst;

A = zeros(250,250,250);
A(111:140,beginTarget-round(WetSampleForDvh):endTarget-round(WetSampleForDvh),111:140) = 1;
ix = find(A > 0);
cst_lungForDvh{2,4}{1} = ix;

resultGUI_lungForDvh = resultGUI;
resultGUI_lungForDvh.physicalDose = resultGUI.physicalDoseLung{ixForDvh};
cst_lungForDvh = matRad_indicatorWrapper(cst_lungForDvh,pln,resultGUI_lungForDvh);

for i = 1:size(cst,1)
    cst{i,8}{1} = cst_0ForDvh{i,8}{1};
    cst{i,8}{2} = cst_lungForDvh{i,8}{1};
    
    cst{i,9}{1} = cst_0ForDvh{i,9}{1};
    cst{i,9}{2} = cst_lungForDvh{i,9}{1};
end

dvhTitle = ['DVH comparison - solid line: 0 mm lung, dotted line: ' num2str(lungThicknessForDvh*2) ' mm lung'];
dvh = figure('Name','DVH comparison','Color',[0.5 0.5 0.5],'Position',([300 300 800 600]));
hold on
for scenIx = [1,2]
    matRad_showDVH(cst,pln,scenIx,scenIx,dvhTitle)   % scenIx = 1: 0 mm, scenIx = 2: lung
end

savefig(dvh,['C:\Matlab\Analysis phantom degradation\fallOff_D95\DVH_targetDepth_' num2str(targetDepth) '_lung_' num2str(lungThicknessForDvh*2) '.fig'])


% %% falloff and DeltaD95 comparison P256 vs. P750
% % clear
% % targetDepth = 250;
% 
% % falloff
% z8020_P256 = load(['C:\Matlab\Analysis phantom degradation\fallOff_D95\P256\falloff_targetDepth_' num2str(targetDepth)]);
% z8020_P750 = load(['C:\Matlab\Analysis phantom degradation\fallOff_D95\P750\falloff_targetDepth_' num2str(targetDepth)]);
% 
% abc_P256 = coeffvalues(z8020_P256.fitFalloffPower);
% a_P256 = abc_P256(1); b_P256 = abc_P256(2); c_P256 = abc_P256(3);
% abc_P750 = coeffvalues(z8020_P750.fitFalloffPower);
% a_P750 = abc_P750(1); b_P750 = abc_P750(2); c_P750 = abc_P750(3);
% 
% 
% fc = figure;
% hold on
% title(['falloff comparison of modulation powers P - target depth ' num2str(targetDepth) ' mm'])
% plot(z8020_P256.z8020(:,1),z8020_P256.z8020(:,2),'bo')
% plot(z8020_P256.fitFalloffPower,'b')
% plot(z8020_P750.z8020(:,1),z8020_P750.z8020(:,2),'rx')
% plot(z8020_P750.fitFalloffPower,'r')
% xlabel('z_g_e_o lung [mm]')
% ylabel('80% - 20% [mm]')
% legend('P_m_o_d = 256 µm',['power fit with a = ' num2str(a_P256) ', b = ' num2str(b_P256) ', c = ' num2str(c_P256)],...
%     'P_m_o_d = 750 µm',['power fit with a = ' num2str(a_P750) ', b = ' num2str(b_P750) ', c = ' num2str(c_P750)],...
%     'location','northwest')
% grid on
% box on
% 
% 
% % DeltaD95
% DeltaD95_P256 = load(['C:\Matlab\Analysis phantom degradation\fallOff_D95\P256\DeltaD95_targetDepth_' num2str(targetDepth)]);
% DeltaD95_P750 = load(['C:\Matlab\Analysis phantom degradation\fallOff_D95\P750\DeltaD95_targetDepth_' num2str(targetDepth)]);
% 
% abc_P256 = coeffvalues(DeltaD95_P256.fitDeltaD95Power);
% a_P256 = abc_P256(1); b_P256 = abc_P256(2); c_P256 = abc_P256(3);
% abc_P750 = coeffvalues(DeltaD95_P750.fitDeltaD95Power);
% a_P750 = abc_P750(1); b_P750 = abc_P750(2); c_P750 = abc_P750(3);
% 
% dc = figure;
% hold on
% title(['Delta D95 comparison of modulation powers P - target depth ' num2str(targetDepth) ' mm'])
% plot(DeltaD95_P256.DeltaD95(:,1),DeltaD95_P256.DeltaD95(:,2),'bo')
% plot(DeltaD95_P256.fitDeltaD95Power,'b')
% plot(DeltaD95_P750.DeltaD95(:,1),DeltaD95_P750.DeltaD95(:,2),'rx')
% plot(DeltaD95_P750.fitDeltaD95Power,'r')
% xlabel('z_g_e_o lung [mm]')
% ylabel('Delta D95 [mm]')
% legend('P_m_o_d = 256 µm',['power fit with a = ' num2str(a_P256) ', b = ' num2str(b_P256) ', c = ' num2str(c_P256)],...
%     'P_m_o_d = 750 µm',['power fit with a = ' num2str(a_P750) ', b = ' num2str(b_P750) ', c = ' num2str(c_P750)],...
%     'location','northwest')
% grid on
% box on
% 
% 
% savefig(fc,['C:\Matlab\Analysis phantom degradation\fallOff_D95\falloff_P_comparison_targetDepth_' num2str(targetDepth) '.fig'])
% savefig(dc,['C:\Matlab\Analysis phantom degradation\fallOff_D95\DeltaD95_P_comparison_targetDepth_' num2str(targetDepth) '.fig'])
% 
% 
