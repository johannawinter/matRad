% Evaluate changes in Ricphantom with target volume / several hundred PB
 
% load RICPHANTOM_none.mat;
% A = zeros(500,500,500);
% A(241:260,291:310,241:260) = 1;
% ix = find(A > 0);
% cst{2,4}{1} = ix;
% cst{2,5}.visibleColor = [0 1 1];
% 
% matRad;
% 
% coords_matRad = .5:1:350;       % [mm]
% coords_spline = .5:.1:350;      % [mm]
% 
% matRad_dd_0 = resultGUI.physicalDose(251,151:end,251);
% matRad_dd_0_spline = spline(coords_matRad,matRad_dd_0,coords_spline);
% 
% matRad_idd_0 = sum(sum(resultGUI.physicalDose(:,151:end,:),3),1);
% matRad_idd_0_spline = spline(coords_matRad,matRad_idd_0,coords_spline);


%% include phantom A but use same weights
load RICPHANTOM_target.mat;
cst{3,5}.visibleColor = [0 1 1];

% remove heterogeneity
cst{2,5}.HeterogeneityCorrection = [];
matRad;

% add heterogeneity
cst{2,5}.HeterogeneityCorrection = 'Lung';
resultGUI_A = matRad_calcDoseDirect(ct,stf,pln,cst,resultGUI.w);
resultGUI.physicalDose_A = resultGUI_A.physicalDose;

% cst_A = matRad_indicatorWrapper(cst,pln,resultGUI_A);
% matRad_showDVH(cst_A,pln)

resultGUI.physicalDose_diff = resultGUI.physicalDose_A - resultGUI.physicalDose;


% create depth dose curves (DD) and integral depth dose curves (IDD)
coords_matRad = .5:1:350;       % [mm]
coords_spline = .5:.1:350;      % [mm]

matRad_dd_0 = resultGUI.physicalDose(251,151:end,251);
matRad_dd_0_spline = spline(coords_matRad,matRad_dd_0,coords_spline);

matRad_idd_0 = sum(sum(resultGUI.physicalDose(:,151:end,:),3),1);
matRad_idd_0_spline = spline(coords_matRad,matRad_idd_0,coords_spline);

matRad_dd_A_same_w = resultGUI.physicalDose_A(251,151:end,251);
matRad_dd_A_same_w_spline = spline(coords_matRad,matRad_dd_A_same_w,coords_spline);

matRad_idd_A_same_w = sum(sum(resultGUI.physicalDose_A(:,151:end,:),3),1);
matRad_idd_A_same_w_spline = spline(coords_matRad,matRad_idd_A_same_w,coords_spline);

%%
figure
hold on
title(['DD: p+ - same weights - target size 2x2x2 cm^3'])
plot(coords_matRad,matRad_dd_0,'bo')
plot(coords_spline,matRad_dd_0_spline,'b')
plot(coords_matRad,matRad_dd_A_same_w,'ro')
plot(coords_spline,matRad_dd_A_same_w_spline,'r')
plot(coords_matRad + 9.2,matRad_dd_A_same_w,'co')
plot(coords_spline + 9.2,matRad_dd_A_same_w_spline,'c')
legend('no phantom','spline','phantom A (P_m_o_d = 256 µm)','spline','phantom A shifted','spline shifted','location','northwest')
axis([0 200 0 2])
xlabel('depth in water [mm]')
ylabel('dose [Gy]')
grid on
box on

%%
figure
hold on
title(['IDD: p+ - same weights - target size 2x2x2 cm^3'])
plot(coords_matRad,matRad_idd_0,'bo')
plot(coords_spline,matRad_idd_0_spline,'b')
plot(coords_matRad,matRad_idd_A_same_w,'ro')
plot(coords_spline,matRad_idd_A_same_w_spline,'r')
plot(coords_matRad + 9.2,matRad_idd_A_same_w,'co')
plot(coords_spline + 9.2,matRad_idd_A_same_w_spline,'c')
legend('no phantom','spline','phantom  A (P_m_o_d = 256 µm)','spline','phantom A shifted','spline shifted','location','northwest')
axis([0 200 0 2100])
xlabel('depth in water [mm]')
ylabel('dose [Gy]')
grid on
box on


%% new optimization
load RICPHANTOM_target.mat;
cst{3,5}.visibleColor = [0 1 1];

matRad;

matRad_dd_A_new_w = resultGUI.physicalDose(251,151:end,251);
matRad_dd_A_new_w_spline = spline(coords_matRad,matRad_dd_A_new_w,coords_spline);

matRad_idd_A_new_w = sum(sum(resultGUI.physicalDose(:,151:end,:),3),1);
matRad_idd_A_new_w_spline = spline(coords_matRad,matRad_idd_A_new_w,coords_spline);


figure
hold on
title(['DD: p+ - new optimization - target size 2x2x2 cm^3'])
plot(coords_matRad,matRad_dd_0,'bo')
plot(coords_spline,matRad_dd_0_spline,'b')
plot(coords_matRad,matRad_dd_A_new_w,'ro')
plot(coords_spline,matRad_dd_A_new_w_spline,'r')
legend('no phantom','spline','phantom ','spline','location','northwest')
axis([0 200 0 2])
xlabel('depth in water [mm]')
ylabel('dose [Gy]')
grid on
box on


figure
hold on
title(['IDD: p+ - new optimization - target size 2x2x2 cm^3'])
plot(coords_matRad,matRad_idd_0,'bo')
plot(coords_spline,matRad_idd_0_spline,'b')
plot(coords_matRad,matRad_idd_A_new_w,'ro')
plot(coords_spline,matRad_idd_A_new_w_spline,'r')
legend('no phantom','spline','phantom A','spline','location','northwest')
axis([0 200 0 2100])
xlabel('depth in water [mm]')
ylabel('dose [Gy]')
grid on
box on
