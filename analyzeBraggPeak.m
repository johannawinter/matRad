% Analyze Bragg peak regarding peak position [mm], R80 position [mm], falloff
% 80\%-20\% (z8020) [mm], width R80-R80 (before peak - behind peak) [mm] and 
% gradient at R80 behind the peak [mm/dose].

%% Choose energy step of HIT_APM data

energyStep = 58;

%% Generate plan
load('RicData.mat');
load RICPHANTOM_extension.mat;
setUpPlanAndStf;

% Lz_A  = 30;      % [mm]  physical length of Phantom A
offset_1 = 0.833; % [mm]  offset between Riccardo's experimental data and his MC simulation [MA p. 53]
offset_2 = 1.1;   % [mm]  offset to correct for discrepancies between Riccardo's measurements and MC simulations and Stephan Brons' data
% offset_3 = 2.89; % [mm]  offset to correct for BAMS and air from nozzle to isocenter, which is considered in matRad, 
                        % value taken from protons_HITfixedBL.data.offset

coords_matRad = .5:1:350;

% % compute middle of bins for MC simulation in [mm]
% coords_pE70Asim = 10*mean(pE70Asim(:,[1 2]),2);

% find correct measurement data
pEnergyStepAexp = eval(['pE' num2str(energyStep) 'Aexp']); % any other solutions without the eval function?!

% compute matRad idd with phantom
matRad_idd_pEnergyStepA = sum(sum(resultGUI.physicalDose(:,151:end,:),3),1);

% spline interpolation
coords_spline_complete = .5:.1:350;                     % step size 0.1 mm
coords_spline = coords_spline_complete(446:1396);       % restrict to 45-140 mm where we have measurement data

pEnergyStepAexp_spline         = spline(pEnergyStepAexp(:,1) - offset_1 - offset_2,pEnergyStepAexp(:,2),coords_spline);
% pE70Asim_spline               = spline(coords_pE70Asim - Lz_A - offset_2,pE70Asim(:,3),coords_spline);
matRad_idd_pEnergyStepA_spline = spline(coords_matRad,matRad_idd_pEnergyStepA,coords_spline);


%% Find peak position
[~,ix_peak_exp] = max(pEnergyStepAexp_spline);
coord_peak_exp = coords_spline(ix_peak_exp);

% [~,ix_peak_sim] = max(pE70Asim_spline);
% coord_peak_sim = coords_spline(ix_peak_sim);

[~,ix_peak_matRad] = max(matRad_idd_pEnergyStepA_spline);
coord_peak_matRad = coords_spline(ix_peak_matRad);

peakPos(1) = coord_peak_exp;
peakPos(2) = coord_peak_matRad;
% peakPos(3) = coord_peak_sim;

format compact
disp(['Energy step E ' num2str(energyStep) ':'])
disp(['Peak position for Riccardo''s measurement (1), matRad simulation (2) and if applicable MC simulation (3).']);
disp(peakPos);


%% Find range R80
R80_exp = max(pEnergyStepAexp_spline)*.8;                                          % determine R80
[~,ix_R80_exp_behind] = min(abs(pEnergyStepAexp_spline(ix_peak_exp:end)-R80_exp)); % find x-value closest to R80 (behind peak)
ix_R80_exp = ix_R80_exp_behind + ix_peak_exp - 1;
coord_R80_exp = coords_spline(ix_R80_exp);                                  % find coordinate

% R80_sim = max(pE70Asim_spline)*.8;
% [~,ix_R80_sim_behind] = min(abs(pE70Asim_spline(ix_peak_sim:end)-R80_sim));
% ix_R80_sim = ix_R80_sim_behind + ix_peak_sim - 1;
% coord_R80_sim = coords_spline(ix_R80_sim);

% [~,x_peak_matRad] = max(matRad_idd_pE70A_spline);
R80_matRad = max(matRad_idd_pEnergyStepA_spline)*.8;
[~,ix_R80_matRad_behind] = min(abs(matRad_idd_pEnergyStepA_spline(ix_peak_matRad:end)-R80_matRad));
ix_R80_matRad = ix_R80_matRad_behind + ix_peak_matRad - 1;
coord_R80_matRad = coords_spline(ix_R80_matRad);

R80Pos(1) = coord_R80_exp;
R80Pos(2) = coord_R80_matRad;
% R80Pos(3) = coord_R80_sim;

disp('Range (R80) for Riccardo''s measurement (1), matRad simulation (2) and if applicable MC simulation (3).');
disp(R80Pos);


%% Calculate falloff 80%-20%
R20_exp = max(pEnergyStepAexp_spline)*.2;                                  % determine 20% of peak
[~,ix_R20_exp_behind] = min(abs(pEnergyStepAexp_spline(ix_peak_exp:end)-R20_exp));  % find x-value closest to 20% (behind peak)
ix_R20_exp = ix_R20_exp_behind + ix_peak_exp - 1;
coord_R20_exp = coords_spline(ix_R20_exp);                          % find coordinate for 20%
z8020exp = coord_R20_exp-coord_R80_exp;

% R20_sim = max(pE70Asim_spline)*.2;
% [~,ix_R20_sim_behind] = min(abs(pE70Asim_spline(ix_peak_sim:end)-R20_sim));
% ix_R20_sim = ix_R20_sim_behind + ix_peak_sim - 1;
% coord_R20_sim = coords_spline(ix_R20_sim);
% z8020sim = coord_R20_sim-coord_R80_sim;

R20_matRad = max(matRad_idd_pEnergyStepA_spline)*.2;
[~,ix_R20_matRad_behind] = min(abs(matRad_idd_pEnergyStepA_spline(ix_peak_matRad:end)-R20_matRad));
ix_R20_matRad = ix_R20_matRad_behind + ix_peak_matRad - 1;
coord_R20_matRad = coords_spline(ix_R20_matRad);
z8020matRad = coord_R20_matRad-coord_R80_matRad;

z8020(1) = z8020exp;
z8020(2) = z8020matRad;
% z8020(3) = z8020sim;

disp('Distal falloff widths z_(80-20) for Riccardo''s measurement (1), matRad simulation (2) and if applicable MC simulation (3).');
disp(z8020);


%% Calculate width of Bragg peak with z difference of R80 in front of and behind peak
[~,ix_R80_exp_front] = min(abs(pEnergyStepAexp_spline(1:ix_peak_exp)-R80_exp));  % find x-value closest to R80 (in front of peak)
coord_R80_exp_front = coords_spline(ix_R80_exp_front);                 % find coordinate
width8080exp = coord_R80_exp-coord_R80_exp_front;

% [~,ix_R80_sim_front] = min(abs(pE70Asim_spline(1:ix_peak_sim)-R80_sim));
% coord_R80_sim_front = coords_spline(ix_R80_sim_front);
% width8080sim = coord_R80_sim-coord_R80_sim_front;

[~,ix_R80_matRad_front] = min(abs(matRad_idd_pEnergyStepA_spline(1:ix_peak_matRad)-R80_matRad));
coord_R80_matRad_front = coords_spline(ix_R80_matRad_front);
width8080matRad = coord_R80_matRad-coord_R80_matRad_front;

width8080(1) = width8080exp;
width8080(2) = width8080matRad;
% width8080(3) = width8080sim;

disp('Bragg peak width from R80 in front of peak to R80 behind peak for Riccardo''s measurement (1), matRad simulation (2) and if applicable MC simulation (3).');
disp(width8080);


%% Calculate gradients at R80 (behind peak), normalized to max(pEnergyStepAexp_spline)
grad_exp = gradient(pEnergyStepAexp_spline./max(pEnergyStepAexp_spline))*10;  % *10 as pE70Aexp_spline (or coords_spline) is given in [0.1 mm] 
gradR80exp = grad_exp(ix_R80_exp);

% grad_sim = gradient(pE70Asim_spline./max(pE70Asim_spline))*10;
% gradR80sim = grad_sim(ix_R80_sim);

grad_matRad = gradient(matRad_idd_pEnergyStepA_spline./max(matRad_idd_pEnergyStepA_spline))*10;
gradR80matRad = grad_matRad(ix_R80_matRad);

gradR80(1) = gradR80exp;
gradR80(2) = gradR80matRad;
% gradR80(3) = gradR80sim;

disp('Gradient at R80 behind peak for Riccardo''s measurement (1), matRad simulation (2) and if applicable MC simulation (3).');
disp(gradR80);

