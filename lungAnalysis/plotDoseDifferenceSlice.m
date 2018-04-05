% plot slice of S00002 for DGMP abstract
clear

load('C:\Matlab\HIT-Lung\S00002\results_3fields_P256.mat')
plane = 3;
slice = round(pln.propStf.isoCenter(1,3)./ct.resolution.z);
borderDoseWindow = min(abs(min(min(resultGUI.diff_matRadHetero_matRadRecalc(:,:,slice)))), max(max(resultGUI.diff_matRadHetero_matRadRecalc(:,:,slice))));
doseWindow = [-borderDoseWindow borderDoseWindow];
doseIsoLevels = [-95 -70 -30 -10 10 30 70 95]/100 * max(abs(min(doseWindow)),max(doseWindow));
voiSelection = [1 0 1 1 1 1 0 0 0 0 1 1 0 0 0 0 1 1];
thresh = 0.05;

doseFig = figure;
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.diff_matRadHetero_matRadRecalc,plane,slice,thresh,1,...
[],redblue,doseWindow,doseIsoLevels,voiSelection,'Gy(RBE)',1,'Linewidth',2);
matRad_plotIsoCenterMarker(gca,pln,ct,plane,slice);
axis([100 225 50 175])
title(' ')

str = {'Dose difference slice, negative values show a lower dose with heterogeneity correction than without. ', ...
    'Target dose was 11.07 Gy(RBE), beam directions: 10°, 65°, 110°.'};
title(str)

% dim = [0.1 0.0 0.5 0.3];
% t = annotation('textbox',dim,'String',str,'FitBoxToText','on');
% t.FontSize = 12;

savefig(doseFig,'C:\Matlab\HIT-Lung\S00002\doseDifferenceSlice')
