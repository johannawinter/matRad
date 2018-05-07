% testRicDataProton_new
% Calculate plan without and with phantom A and compare them with
% measurements and MC simulations.

clear
close all

addpath(genpath('phantomAnalysis'))

%% Calculate plan without a phantom
clear
energyStep = 112;   % 58/70/112

load('RicData.mat');
pEnergyStep0exp = eval(['pE' num2str(energyStep) '0exp']);
% All simulation data is in the format:
% start_bin_z - end_bin_z - dose - error_dose 
% where the bins in z are in cm for simulated data and in mm for experimental data;
% the dose is the absolute dose for 1 primary for the simulation and 
% the relative dose IonizationChamber1/IonizationChamber2 for the experiment.

% The files with the index "0" refer to the measurement without any phantom.
% Phantom A (Gammex lung) has the following properties:
Lz_A  = 30;     % [mm] physical length of Phantom A
% WET_A = 8.9;    % [mm] with protons

load RICPHANTOM_none.mat;
setUpPlanAndStf

% restrict matRad coordinates and idd to water phantom (and behind)
coords_matRad = .5:1:350;
% spline interpolation
coords_spline_complete = .5:.1:350;
coords_spline = coords_spline_complete(446:1396);       % restrict to reasonable depths where we have measurement data: 45-140 mm

% compute matRad idd (phyical dose!) and interpolate with spline
matRad_idd_pEnergyStep0        = sum(sum(resultGUI.physicalDose(:,151:end,:),3),1);
matRad_idd_pEnergyStep0_spline = spline(coords_matRad,matRad_idd_pEnergyStep0,coords_spline);

if energyStep == 70
    % compute middle of bins for MC simulation in [mm]
    coords_pE700sim             = 10*mean(pE700sim(:,[1 2]),2);
    % interpolate exp and sim data with spline, preparation for evaluating R80
    pE700expStephan_splinePre	= spline(pE700expStephan(:,1),pE700expStephan(:,2),coords_spline);
    pE700sim_splinePre          = spline(coords_pE700sim,pE700sim(:,3),coords_spline);
end
pE700exp_splinePre              = spline(pEnergyStep0exp(:,1),pEnergyStep0exp(:,2),coords_spline);


% calculate offsets to R80 of matRad
% find R80 value of matRad idd, both experimental and simulation data
[maxMatradIddSpline, maxMatradIddSplineIx] = max(matRad_idd_pEnergyStep0_spline);
R80matrad = matRad_interp1(matRad_idd_pEnergyStep0_spline(end:-1:maxMatradIddSplineIx)', ...
    coords_spline(end:-1:maxMatradIddSplineIx)', maxMatradIddSpline*.8);

[maxExpIddSpline, maxExpIddSplineIx] = max(pE700exp_splinePre);
R80exp = matRad_interp1(pE700exp_splinePre(end:-1:maxExpIddSplineIx)', ...
    coords_spline(end:-1:maxExpIddSplineIx)', maxExpIddSpline*.8);

if energyStep == 70
    [maxExpStephanIddSpline, maxExpStephanIddSplineIx] = max(pE700expStephan_splinePre);
    R80expStephan = matRad_interp1(pE700expStephan_splinePre(end:-1:maxExpStephanIddSplineIx)', ...
        coords_spline(end:-1:maxExpStephanIddSplineIx)', maxExpStephanIddSpline*.8);
    
    [maxSimIddSpline, maxSimIddSplineIx] = max(pE700sim_splinePre);
    R80sim = matRad_interp1(pE700sim_splinePre(end:-1:maxSimIddSplineIx)', ...
        coords_spline(end:-1:maxSimIddSplineIx)', maxSimIddSpline*.8);
end

offset_expMatrad = R80matrad-R80exp;
if energyStep == 70
offset_expStephanMatrad = R80matrad-R80expStephan;
offset_simMatrad = R80matrad-R80sim;
% check BAMS offset
offset_BAMS = machine.data(energyStep).offset;
end

% re-interpolate exp and sim data with correct offset
pEnergyStep0exp_spline      = spline(pEnergyStep0exp(:,1) + offset_expMatrad, pEnergyStep0exp(:,2), coords_spline);
if energyStep == 70
    pE700expStephan_spline  = spline(pE700expStephan(:,1) + offset_expStephanMatrad, pE700expStephan(:,2), coords_spline);
    pE700sim_spline         = spline(coords_pE700sim + offset_simMatrad, pE700sim(:,3), coords_spline);
end


% calculate trapezoidal integral over 45-95/135 mm (x-dir.)
if energyStep == 58 || energyStep == 70
    trapz_matRad_pE700 = trapz(coords_spline(1:501),matRad_idd_pEnergyStep0_spline(1:501));
    trapz_pE700exp = trapz(coords_spline(1:501) + offset_expMatrad, pEnergyStep0exp_spline(1:501));
elseif energyStep == 112
    trapz_matRad_pE700 = trapz(coords_spline(1:901),matRad_idd_pEnergyStep0_spline(1:901));
    trapz_pE700exp = trapz(coords_spline(1:901) + offset_expMatrad, pEnergyStep0exp_spline(1:901));
end
% calculate normalization factors to matRad
normFactorExp = trapz_matRad_pE700/trapz_pE700exp;
if energyStep == 70
    trapz_pE700expStephan = trapz(coords_spline(1:501) + offset_expStephanMatrad, pE700expStephan_spline(1:501));
    trapz_pE700sim = trapz(coords_spline(1:501) + offset_simMatrad, pE700sim_spline(1:501));
    normFactorExpStephan = trapz_matRad_pE700/trapz_pE700expStephan;
    normFactorSim = trapz_matRad_pE700/trapz_pE700sim;
end


% plot all curves without phantom
validationFig = figure;
hold on
if energyStep == 58 || energyStep == 70
    title(['IDD: p+ Energy ' num2str(energyStep) ' - solid: no sample, dashed: phantom A - normalization to R_{80} (x) and ' ...
        'by trapezoidal integration over x = 45-95 mm (y) - Pmod = 256 µm'])
elseif energyStep == 112
    title(['IDD: p+ Energy ' num2str(energyStep) ' - solid: no sample, dashed: phantom A - normalization to R_{80} (x) and ' ...
        'by trapezoidal integration over x = 45-135 mm (y) - Pmod = 256 µm'])
end
plot(coords_matRad,                             matRad_idd_pEnergyStep0, 'xr')
plot(coords_spline,                             matRad_idd_pEnergyStep0_spline, '-r')
plot(pEnergyStep0exp(:,1) + offset_expMatrad,   pEnergyStep0exp(:,2) * normFactorExp, 'xb')
plot(coords_spline,                             pEnergyStep0exp_spline * normFactorExp,'-b')
legend('matRad','spline','measurement (Riccardo)','spline','Location','northwest')
% legend('matRad','measurement (Riccardo)','Location','northwest')
axis([max(45,round(machine.data(energyStep).peakPos,-1)-45) round(machine.data(energyStep).peakPos,-1)+10 0 max(matRad_idd_pEnergyStep0_spline)*1.1])
xlabel('depth in water [mm]')
ylabel('dose [a.u.]')
grid on
box on

if energyStep == 70
    plot(pE700expStephan(:,1) + offset_expStephanMatrad,    pE700expStephan(:,2) * normFactorExpStephan,'xk')
    plot(coords_spline,                                     pE700expStephan_spline * normFactorExpStephan, '-k')
    plot(coords_pE700sim + offset_simMatrad,                pE700sim(:,3) * normFactorSim,'xg')
    plot(coords_spline,                                     pE700sim_spline * normFactorSim, '-g')
    legend('matRad','spline','measurement (Riccardo)','spline',...
        'Stephan''s measurement','spline',...
        'MC simulation (Riccardo)','spline','Location','northwest')
%     legend('matRad','measurement (Riccardo)','Stephan''s measurement',...
%         'MC simulation (Riccardo)','Location','northwest')
end


%% Calculate plan with phantom A
% create phantom in cst
A = zeros(500,500,500);
A(:,21:50,:) = 1;
ix = find(A>0);

cst{3,1} = 2;
cst{3,2} = 'Lung phantom';
cst{3,3} = 'OAR';
cst{3,4}{1} = ix;
cst{3,5} = cst{1,5};
cst{3,5}.visibleColor = [1 1 0.3333];
cst{3,5}.HeterogeneityCorrection = 'Lung';
cst{3,6} = cst{1,6};

ct.cube{1}(cst{3,4}{1}) = 0.297;        % so far: rho_lung = 0.306;

% compute dose distribution
setUpPlanAndStf

% compute matRad idd with phantom and interpolate with spline
matRad_idd_pEnergyStepA        = sum(sum(resultGUI.physicalDose(:,151:end,:),3),1);
matRad_idd_pEnergyStepA_spline = spline(coords_matRad,matRad_idd_pEnergyStepA,coords_spline);


% find correct measurement data phantom A
pEnergyStepAexp = eval(['pE' num2str(energyStep) 'Aexp']); % any other solutions without the eval function?!

if energyStep == 70
    % compute middle of bins for MC simulation in [mm]
    coords_pE70Asim             = 10*mean(pE70Asim(:,[1 2]),2);
    % interpolate exp and sim data with offset from above
    pE70Asim_spline             = spline(coords_pE70Asim + offset_simMatrad - Lz_A, pE70Asim(:,3), coords_spline);
end
pEnergyStepAexp_spline          = spline(pEnergyStepAexp(:,1) + offset_expMatrad, pEnergyStepAexp(:,2), coords_spline);

% calculate trapezoidal integral over 45-95 mm  (x-dir.) for normalization of simulation
if energyStep == 58 || energyStep == 70
    trapz_matRad_pEnergyStepA = trapz(coords_spline(1:501),matRad_idd_pEnergyStepA_spline(1:501));
    trapz_pEnergyStepAexp = trapz(coords_spline(1:501),pEnergyStepAexp_spline(1:501));
elseif energyStep == 112
    trapz_matRad_pEnergyStepA = trapz(coords_spline(1:901),matRad_idd_pEnergyStepA_spline(1:901));
    trapz_pEnergyStepAexp = trapz(coords_spline(1:901),pEnergyStepAexp_spline(1:901));
end
% calculate normalization factors to matRad
normFactorExpA = trapz_matRad_pEnergyStepA/trapz_pEnergyStepAexp;
if energyStep == 70
    trapz_pE70Asim = trapz(coords_spline(1:501),pE70Asim_spline(1:501));
    normFactorSimA = trapz_matRad_pEnergyStepA/trapz_pE70Asim;
end

% add phantom curves to plot
plot(coords_matRad,                             matRad_idd_pEnergyStepA, 'xr')
plot(coords_spline,                             matRad_idd_pEnergyStepA_spline, '--r')
plot(pEnergyStepAexp(:,1) + offset_expMatrad,   pEnergyStepAexp(:,2) * normFactorExpA, 'xb')
plot(coords_spline,                             pEnergyStepAexp_spline * normFactorExpA,'--b')
legend('matRad','spline','measurement (Riccardo)','spline','Location','northwest')
% legend('matRad','measurement (Riccardo)','Location','northwest')
if energyStep == 70
    plot(coords_pE70Asim + offset_simMatrad - Lz_A, pE70Asim(:,3) * normFactorSimA, 'xg')
    plot(coords_spline,                         pE70Asim_spline * normFactorSimA, '--g')
    legend('matRad','spline','measurement (Riccardo)','spline',...
        'Stephan''s measurement','spline',...
        'MC simulation (Riccardo)','spline','Location','northwest')
%     legend('matRad','measurement (Riccardo)','Stephan''s measurement','MC simulation (Riccardo)', ...
%         'Location','northwest')
end

%%
%%%%%%%%%%%%%%%%%%%%%%%% test
figure; 
hold on
title(['IDD: p+ Energy ' num2str(energyStep) ' - solid: no sample, dashed: phantom A - normalization to R_{80} (x) and ' ...
'maximum - Pmod = 256 µm'])
plot(coords_matRad, matRad_idd_pEnergyStep0 / max(matRad_idd_pEnergyStep0_spline), 'xr')
plot(coords_spline, matRad_idd_pEnergyStep0_spline / max(matRad_idd_pEnergyStep0_spline), '-r')
plot(pEnergyStep0exp(:,1) + offset_expMatrad, pEnergyStep0exp(:,2) / max(pEnergyStep0exp_spline), 'xb')
plot(coords_spline, pEnergyStep0exp_spline / max(pEnergyStep0exp_spline) ,'-b')
plot(coords_pE700sim + offset_simMatrad, pE700sim(:,3) / max(pE700sim_spline),'xg')
plot(coords_spline, pE700sim_spline / max(pE700sim_spline), '-g')

axis([70 92 0 1.05])
xlabel('depth in water [mm]')
ylabel('dose [a.u.]')
grid on
box on

plot(coords_matRad, matRad_idd_pEnergyStepA / max(matRad_idd_pEnergyStepA_spline), 'xr')
plot(coords_spline, matRad_idd_pEnergyStepA_spline / max(matRad_idd_pEnergyStepA_spline), '--r')
plot(pEnergyStepAexp(:,1) + offset_expMatrad, pEnergyStepAexp(:,2) / max(pEnergyStepAexp_spline), 'xb')
plot(coords_spline, pEnergyStepAexp_spline / max(pEnergyStepAexp_spline), '--b')
plot(coords_pE70Asim + offset_simMatrad - Lz_A, pE70Asim(:,3) / max(pE70Asim_spline), 'xg')
plot(coords_spline, pE70Asim_spline / max(pE70Asim_spline), '--g')
%%%%%%%%%%%%%%%%%%%%%%%% end test


%% save results for voxelwise convolution
if energyStep == 70
    save(['C:\Matlab\Analysis RicData bugfix\results_EIx' num2str(energyStep)],...
        'cst','ct','pln','resultGUI','stf',...
        'energyStep','Lz_A',...
        'coords_matRad','coords_spline','matRad_idd_pEnergyStep0','matRad_idd_pEnergyStep0_spline',...
        'matRad_idd_pEnergyStepA','matRad_idd_pEnergyStepA_spline',...
        'coords_pE700sim','coords_pE70Asim',...
        'pEnergyStep0exp','pEnergyStep0exp_spline','pEnergyStepAexp','pEnergyStepAexp_spline',...
        'pE700expStephan','pE700expStephan_spline','pE700sim','pE700sim_spline','pE70Asim','pE70Asim_spline',...
        'normFactorExp','normFactorExpA','normFactorExpStephan','normFactorSim','normFactorSimA',...
        'offset_expMatrad','offset_expStephanMatrad','offset_simMatrad','offset_BAMS','-v7.3')
else
    save(['C:\Matlab\Analysis RicData bugfix\results_EIx' num2str(energyStep)],...
        'cst','ct','pln','resultGUI','stf',...
        'energyStep','Lz_A',...
        'coords_matRad','coords_spline','matRad_idd_pEnergyStep0','matRad_idd_pEnergyStep0_spline',...
        'matRad_idd_pEnergyStepA','matRad_idd_pEnergyStepA_spline',...
        'pEnergyStep0exp','pEnergyStep0exp_spline','pEnergyStepAexp','pEnergyStepAexp_spline',...
        'normFactorExp','normFactorExpA','offset_expMatrad','-v7.3')
end

savefig(validationFig, ['C:\Matlab\Analysis RicData bugfix\validation_rho0,297_P256_EIx'...
    num2str(energyStep) '_complete.fig'])


%% compare different implementations
clear
close all

energyStep = 70;

complete = load(['C:\Matlab\Analysis RicData bugfix\completeConvolution\results_EIx' num2str(energyStep)]);
depthBased = load(['C:\Matlab\Analysis RicData bugfix\depthBasedConvolution\results_EIx' num2str(energyStep)]);
voxelwise = load(['C:\Matlab\Analysis RicData bugfix\voxelwiseConvolution\results_EIx' num2str(energyStep)]);

diff_depthBased_complete_0 = depthBased.matRad_idd_pEnergyStep0 - complete.matRad_idd_pEnergyStep0;
diff_voxelwise_complete_0  = voxelwise.matRad_idd_pEnergyStep0 - complete.matRad_idd_pEnergyStep0;
if sum(diff_depthBased_complete_0) ~= 0 || sum (diff_voxelwise_complete_0) ~= 0
    warning('matRad data without phantom is not the same in all implementations.')
end

diff_depthBased_complete_A  = depthBased.matRad_idd_pEnergyStepA - complete.matRad_idd_pEnergyStepA;
maxDiff_depthBased_complete = max(abs(diff_depthBased_complete_A));
fprintf('Maximum absolute difference between depth based convolution and complete convolution is %.2e.\n', ...
    maxDiff_depthBased_complete)

diff_voxelwise_complete_A  = voxelwise.matRad_idd_pEnergyStepA - complete.matRad_idd_pEnergyStepA;
maxDiff_voxelwise_complete = max(abs(diff_voxelwise_complete_A));
fprintf('Maximum absolute difference between voxelwise convolution and complete convolution is %.2e.\n', ...
    maxDiff_voxelwise_complete)

diff_depthBased_voxelwise_A  = depthBased.matRad_idd_pEnergyStepA - voxelwise.matRad_idd_pEnergyStepA;
maxDiff_depthBased_voxelwise = max(abs(diff_depthBased_voxelwise_A));
fprintf('Maximum absolute difference between depth based convolution and voxelwise convolution is %.2e.\n', ...
    maxDiff_depthBased_voxelwise)


figure
hold on
plot(complete.coords_matRad+151, diff_depthBased_complete_A, 'xb')
plot(complete.coords_matRad+151, diff_voxelwise_complete_A, 'xr')
plot(complete.coords_matRad+151, diff_depthBased_voxelwise_A, 'xg')
legend('depthBased - complete','voxelwise - complete','depthBased - voxelwise')
xlim([150 300])

plot(complete.coords_matRad+151, complete.matRad_idd_pEnergyStepA, 'k')
legend('depthBased - complete','voxelwise - complete','depthBased - voxelwise','complete')
