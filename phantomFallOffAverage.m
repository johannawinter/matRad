% calculation of DeltaD95 and falloff and averaging over several depth dose
% curves around central ray

% define phantom setup
breastThickness = 30;
targetThickness = 40;
lungGeoThickness = [2 7 20 30 40 50 60 70 80 90 100];

% load results
result = struct(length(lungGeoThickness),1);
for h = 1:length(lungGeoThickness)
    result(h) = load(['C:\Matlab\Analysis phantom degradation\fallOff_D95_accordingToSigmaAnalysis\breast'...
        num2str(breastThickness) '_target' num2str(targetThickness) ...
        '\results_breastThickness_' num2str(breastThickness) ...
        '_targetThickness_' num2str(targetThickness) ...
        '_lungThickness_' num2str(lungGeoThickness(h)) '.mat']);
end

%% compute 30 depth dose curves around central ray
% define coordinates
coords_matRad = 1:1:250;       % [mm*2]
coords_spline = .05:.0005:250;       % [mm*2]

% no heterogeneity: calculate depth dose curves
% get central ray
% calculate 30 DD around central ray
for i = 30
    tempIx = -15+i;
    
    dd_0_spline.i = 
    
    
    
    % calculate DeltsD95 each
    [~,ix_peak] = max(dd_0_spline);
    
    D95_0 = 2 * .95;
    [~,ix_D95_0_behind] = min(abs(dd_0_spline(ix_peak:end)-D95_0));
    ix_D95_0 = ix_D95_0_behind + ix_peak - 1;
    coord_D95_0 = coords_spline(ix_D95_0);
    
    DeltaD95_0(1,1) = 0.0001;
    DeltaD95_0(1,2) = 0;
    
    % falloff
    
end

% average DeltaD95 and falloff



% calculate for different lung thicknesses
h = 1;

% get central ray
centralRay.x = round(result(h).pln.isoCenter(2)/2);
centralRay.y = round(result(h).pln.isoCenter(3)/2);

doseHomo = result(h).resultGUI.physicalDose_noHeterogeneity;
doseLung = result(h).resultGUI.physicalDose_Lung;
% calculate 30 DD around central ray
for i = 1:30
    tempIx = -15+i;
    dd.i = doseLung(round(result(h).pln.isoCenter(2)/2)+tempIx, :, round(result(h).pln.isoCenter(3)/2)+tempIx);
    dd_spline.i = spline(coords_matRad,dd.i,coords_spline);
    
    % comupte Dalta D95 each
    [~,ix_peak] = max(dd_spline.i);

    D95 = 2 * .95;
    [~,ix_D95_behind] = min(abs(dd_spline.i(ix_peak:end)-D95));
    ix_D95 = ix_D95_behind + ix_peak - 1;
    coord_D95 = coords_spline(ix_D95);

    DeltaD95 = coord_D95_0 - coord_D95;
    DeltaD95.i(2,1) = lungGeoThickness(h);
    DeltaD95.i(2,2) = DeltaD95*2;
end



