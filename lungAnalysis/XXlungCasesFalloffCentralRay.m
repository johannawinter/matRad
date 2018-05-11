% Plot depth dose curves for lung cancer cases - first only 1 field. 
% Calculate (DeltaD95 and) falloff z80%-20% with and without
% heterogeneity correction.

clear
close all

testPlots = 1;


%% prepare patient data
% load treatment case with all dose distributions, only 1 field
load('C:\Matlab\HIT-Lung\H03368\1_field\results_1fields_P256')
% load('C:\Matlab\HIT-Lung\S00003\2_fields\results_2fields_P256')
% load('C:\Matlab\HIT-Lung\S00002\results_3fields_P256')
% load('C:\Matlab\HIT-Lung\S00001\results_3fields_P256')

patientID = 'H03368';
% patientID = 'S00003';
% patientID = 'S00002';
% patientID = 'S00001';
beam = 1;                       % choose which field to be analyzed
completeDoseDistribution = 1;   % chose if falloff of complete dose distribution should be analyzed (1) or not (0)

% get ct cube
ctCube = ct.cube{1};

% get lung cube
lungCube = zeros(ct.cubeDim);
for i = 1:size(cst,1)
   isLung = contains(cst{i,2},'lung','IgnoreCase',true);
   if isLung
       lungCube(cst{i,4}{1}) = 1;
       fprintf(['Structure "' cst{i,2} '" added to binary lung cube.\n']);
   end
end

% get dose cubes
if completeDoseDistribution
    doseCube{1} = resultGUI.RBExDose;           % Original
    doseCube{2} = resultGUI.matRadRecalc;       % Recalculated
    doseCube{3} = resultGUI.matRadHeteroRecalc; % Heterogeneity Correction
else
    doseCube{1} = resultGUI.RBExDose_BEAM_4;     % Original for patient S00002
    warning('Falloff analysis for separate beams not implemented yet.')
end

if testPlots
    % test if lung cube is in the right position
    addpath(genpath('tools'))
    addpath(genpath('plotting'))
    slice = round(stf(1).isoCenter(3)/ct.resolution.z);
    
    figure
    matRad_plotSliceWrapper(gca,ct,cst,1,doseCube{1},3,slice,[],[],colorcube);
    
    figure
    matRad_plotSliceWrapper(gca,ct,cst,1,lungCube,3,slice,[],[],colorcube,[],[-.01 1.01]);
end


%% get depth dose curves along central ray
% find ray corresponding to isocenter
cnt = 1;
for i = 1:length(stf(beam).ray)
    if abs(stf(beam).ray(i).rayPos_bev) <= [.1 .1 .1]
        isocenterRay(cnt) = i;
        cnt = cnt + 1;
    end
end

% check that only one ray is going through isocenter
if length(isocenterRay) ~= 1
    warning(['There is more than one ray going through the isocenter. ' ...
        'Only first ray is considered for falloff analysis.'])
    isocenterRay = isocenterRay(1);
end

% perform Siddon ray tracing on central axis
[~,l,rho,~,ix] = matRad_siddonRayTracer(stf(beam).isoCenter, ...
                                   ct.resolution, ...
                                   stf(beam).sourcePoint, ...
                                   stf(beam).ray(isocenterRay).targetPoint, ...
                                   [{ctCube},{lungCube},{doseCube{1}},{doseCube{2}},{doseCube{3}}]);

% Calculate geometrical, radiological and lung depth
geomDepth = cumsum(l) - l/2;

d = [0 l.*rho{1}];
radDepth = cumsum(d(1:end-1));

dLung = [0 l.*rho{2}];
lungDepth = cumsum(dLung(1:end-1));

% calculate depth dose curves for all dose cubes
dd{1} = rho{3};     % Original
dd{2} = rho{4};     % Recalc
dd{3} = rho{5};     % Hetero


% create spline of dd curves
coordsSplineGeom = geomDepth(1):.1:geomDepth(end);
coordsSplineRad = radDepth(1):.1:radDepth(end);

ddSplineGeom = cell(1,3);
for i = 1:length(ddSplineGeom)
    ddSplineGeom{i} = spline(geomDepth,dd{i},coordsSplineGeom);
end

ddSplineRad = cell(1,3);
for i = 1:length(ddSplineRad)
    ddSplineRad{i} = spline(radDepth,dd{i},coordsSplineRad);
end


% create dd curve for linearly interpolated geometrical coordinates
coordsInterpGeom = linspace(geomDepth(1),geomDepth(end),5000);
ddInterpGeom = cell(1,3);
for i = 1:length(ddInterpGeom)
    ddInterpGeom{i} = matRad_interp1(geomDepth', dd{i}', coordsInterpGeom');
end


%% geometrical depth as reference
% plot dose vs. geometrical depth
ddGeomFig = figure;
title(['DD of case ' patientID ', beam ' num2str(beam) ', along central ray ' num2str(isocenterRay)])
hold on
plot(geomDepth, dd{1},'-xb');
plot(geomDepth, dd{2},'-xg');
plot(geomDepth, dd{3},'-xr');
legend('original dose','matRad recalculated dose',...
    'matRad dose with heterogeneity correction')
xlabel('geometrical depth [mm]')
ylabel('dose [Gy]')
axis([0 max(geomDepth) 0 max(dd{1})+.2])
grid on
grid minor
box on


% asses target and CTV coordinates
tmpPrior = intmax;
tmpSize = 0;
for i = 1:size(cst,1)
    if strcmp(cst{i,3},'TARGET') && tmpPrior >= cst{i,5}.Priority && tmpSize<numel(cst{i,4}{1})
        linIdxTarget = unique(cst{i,4}{1});
        tmpPrior=cst{i,5}.Priority;
        tmpSize=numel(cst{i,4}{1});
        VOI = cst{i,2};
    end
    
    if strcmp(cst{i,2},'CTV')
        linIdxCtv = cst{i,4}{1};
    end
%     if strfind(cst{i,2},'CTV') > 0
%         linIdxCtv = cst{i,4}{1};
%     end
end

% plot target and CTV boundaries
mTargetCube = zeros(ct.cubeDim);
mTargetCube(linIdxTarget) = 1;
vProfileTarget = mTargetCube(ix);

geomTargetEntry = geomDepth(find(vProfileTarget,1));
geomTargetExit  = geomDepth(find(vProfileTarget,1,'last'));

mCtvCube = zeros(ct.cubeDim);
mCtvCube(linIdxCtv) = 1;
vProfileCtv = mCtvCube(ix);

geomCtvEntry = geomDepth(find(vProfileCtv,1));
geomCtvExit  = geomDepth(find(vProfileCtv,1,'last'));

plot([geomCtvEntry geomCtvEntry],       [0 max(dd{1}(:))+.2],'--k')
plot([geomTargetEntry geomTargetEntry], [0 max(dd{1}(:))+.2],'--','color',[0 .6 0])
plot([geomCtvExit geomCtvExit],         [0 max(dd{1}(:))+.2],'--k')
plot([geomTargetExit geomTargetExit],   [0 max(dd{1}(:))+.2],'--','color',[0 .6 0])
legend('original dose','matRad recalculated dose',...
    'matRad dose with heterogeneity correction',...
    'CTV boundaries', ['target boundaries: ' VOI])


% plot lung depth vs. geom depth
lungDepthGeomFig = figure;
title(['Lung depth vs. geometrical depth along central ray ' num2str(isocenterRay)])
hold on
plot(geomDepth, lungDepth);
xlabel('geometrical depth [mm]')
ylabel('geometrical lung depth [mm]')
grid on
grid minor
box on


%% radiological depth as reference
% plot dose vs. water equivalent depth
ddRadFig = figure;
title(['DD of case ' patientID ', beam ' num2str(beam) ', along central ray ' num2str(isocenterRay)])
hold on
plot(radDepth, dd{1},'-xb');                 % Original
plot(radDepth, dd{2},'-xg');                 % Recalc
plot(radDepth, dd{3},'-xr');                 % Hetero
legend('original dose','matRad recalculated dose',...
    'matRad dose with heterogeneity correction')
xlabel('radiological depth [mm]')
ylabel('dose [Gy]')
axis([0 max(radDepth) 0 max(dd{1})+.2])
grid on
grid minor
box on

% plot target and CTV boundaries
radTargetEntry = radDepth(find(vProfileTarget,1));
radTargetExit  = radDepth(find(vProfileTarget,1,'last'));

radCtvEntry = radDepth(find(vProfileCtv,1));
radCtvExit  = radDepth(find(vProfileCtv,1,'last'));

plot([radCtvEntry radCtvEntry],         [0 max(dd{1}(:))+.2], '--k')
plot([radTargetEntry radTargetEntry],   [0 max(dd{1}(:))+.2], '--','color',[0 .6 0])
plot([radCtvExit radCtvExit],           [0 max(dd{1}(:))+.2], '--k')
plot([radTargetExit radTargetExit],     [0 max(dd{1}(:))+.2], '--','color',[0 .6 0])
legend('original dose','matRad recalculated dose',...
    'matRad dose with heterogeneity correction',...
    'CTV boundaries', ['target boundaries: ' VOI])


% plot lung depth vs. rad depth
lungDepthRadFig = figure;
title(['Lung depth vs. radiological depth along central ray ' num2str(isocenterRay)])
hold on
plot(radDepth, lungDepth);
xlabel('radiological depth [mm]')
ylabel('geometrical lung depth [mm]')
grid on
grid minor
box on


%% find distal falloff z80%-20%
% use spline coordinates
for i = 1:length(ddSplineGeom)
    % get peak value
    peakValue(i) = max(dd{i});
    R80(i) = peakValue(i) * .8;
    R20(i) = peakValue(i) * .2;
    
    % get spline coordinates for peak, R80, R20
    [~,ixPeakSpline(i)] = min(abs(ddSplineGeom{i} - peakValue(i)));
    
    [~,ixR80behind] = min(abs(ddSplineGeom{i}(ixPeakSpline(i):end)-R80(i)));
    ixR80 = ixR80behind + ixPeakSpline(i) - 1;
    coordR80Spline(i) = coordsSplineGeom(ixR80);
    
    [~,ixR20behind] = min(abs(ddSplineGeom{i}(ixPeakSpline(i):end)-R20(i)));
    ixR20 = ixR20behind + ixPeakSpline(i) - 1;
    coordR20Spline(i) = coordsSplineGeom(ixR20);
    
    % calculate z8020
    z8020Spline(i) = (coordR20Spline(i)-coordR80Spline(i));
end

% use linear interpolation coordinates
for i = 1:length(ddInterpGeom)
    % get interpolated coordinates for peak, R80, R20
    [~,ixPeakInterp(i)] = min(abs(ddInterpGeom{i} - peakValue(i)));
    
    [~,ixR80behind] = min(abs(ddInterpGeom{i}(ixPeakInterp(i):end)-R80(i)));
    ixR80 = ixR80behind + ixPeakInterp(i) - 1;
    coordR80Interp(i) = coordsInterpGeom(ixR80);
    
    [~,ixR20behind] = min(abs(ddInterpGeom{i}(ixPeakInterp(i):end)-R20(i)));
    ixR20 = ixR20behind + ixPeakInterp(i) - 1;
    coordR20Interp(i) = coordsInterpGeom(ixR20);
    
    % calculate z8020
    z8020Interp(i) = (coordR20Interp(i)-coordR80Interp(i));
end

% plot z8020 into plots of all dose curves - from spline coordinates
falloffTestFig = figure;

subplot(311)
title(['DD of case ' patientID ', beam ' num2str(beam) ...
    ', along central ray ' num2str(isocenterRay) ', original dose'])
hold on
ddLine = plot(geomDepth, dd{1},'-xb');
plot(coordR80Spline(1),R80(1),'xr');
plot(coordR20Spline(1),R20(1),'xr');
z8020Line = plot([coordR80Spline(1) coordR80Spline(1)+z8020Spline(1)], [R80(1) R20(1)],'r');
legend([ddLine, z8020Line], ...
    'dose values',['z8020 = ' num2str(z8020Spline(1)) ' mm from spline'])
xlabel('geometrical depth [mm]')
ylabel('dose [Gy]')
axis([0 max(geomDepth) 0 max(ddSplineGeom{1})+.2])
grid on
grid minor
box on

subplot(312)
title('matRad recalculated dose')
hold on
ddLine = plot(geomDepth, dd{2},'-xb');
plot(coordR80Spline(2),R80(2),'xr');
plot(coordR20Spline(2),R20(2),'xr');
z8020Line = plot([coordR80Spline(2) coordR80Spline(2)+z8020Spline(2)], [R80(2) R20(2)],'r');
legend([ddLine, z8020Line], ...
    'dose values',['z8020 = ' num2str(z8020Spline(2)) ' mm from spline'])
xlabel('geometrical depth [mm]')
ylabel('dose [Gy]')
axis([0 max(geomDepth) 0 max(ddSplineGeom{1})+.2])
grid on
grid minor
box on

subplot(313)
title('matRad recalculated dose with heterogeneity correction')
hold on
ddLine = plot(geomDepth, dd{3},'-xb');
plot(coordR80Spline(3),R80(3),'xr');
plot(coordR20Spline(3),R20(3),'xr');
z8020Line = plot([coordR80Spline(3) coordR80Spline(3)+z8020Spline(3)], [R80(3) R20(3)],'r');
legend([ddLine, z8020Line], ...
    'dose values',['z8020 = ' num2str(z8020Spline(3)) ' mm from spline'])
xlabel('geometrical depth [mm]')
ylabel('dose [Gy]')
axis([0 max(geomDepth) 0 max(ddSplineGeom{1})+.2])
grid on
grid minor
box on


%% save results
save(['C:\Matlab\HIT-Lung_falloff\' patientID '\results_falloff'],...
    'cst','ct','stf','lungCube','doseCube','isocenterRay',...
    'l','rho','ix','geomDepth','radDepth','lungDepth',...
    'coordsSplineGeom','coordsSplineRad','coordsInterpGeom',...
    'dd','ddSplineGeom','ddSplineRad','ddInterpGeom',...
    'geomCtvEntry','geomCtvExit','radCtvEntry','radCtvExit',...
    'geomTargetEntry','geomTargetExit','radTargetEntry','radTargetExit','VOI',...
    'peakValue','ixPeakSpline','ixPeakInterp',...
    'R80','R20','coordR80Spline','coordR80Interp','coordR20Spline','coordR20Interp',...
    'z8020Spline','z8020Interp')

save(['C:\Matlab\HIT-Lung_falloff\' patientID '\z8020'],...
    'z8020Spline','z8020Interp')

savefig([ddGeomFig,lungDepthGeomFig,ddRadFig,falloffTestFig], ...
    ['C:\Matlab\HIT-Lung_falloff\' patientID '\dds.fig'])

