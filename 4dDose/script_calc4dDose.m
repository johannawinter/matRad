% 4D dose calculation workflow

fileName = 'phantom30402';

% plan treatment or load ct, cst, dij, pln, stf, resultGUI

% post processing
% This step is necessary to remove beam spots with too few particles that
% could not be delivered
resultGUI = matRad_postprocessing(resultGUI,dij,pln,cst,stf);

% export plan
% Number of beams files called PBP_0X_Plan01.xml are generated the order 
% of spot positions defines the delivery path - you can choose between 
% 'stfMode'- order of spots as in stf file - line wise 'backforth' - first 
% row is delivered from left to right, next right to left and so on 'TSP' 
% (attention, very slow) the shortest path between all spots in each 
% energy slice is calculated.
matRad_export_HITXMLPlan_modified(fileName,pln,stf,resultGUI,'stfMode')

% make Lmd out
% Need to call external program makeLmdout (installed on sievert22) to 
% calculate the time structure of the delivery 
% makeLmdout -p PBP_0X_Plan01.xml -o D_0X_Plan01 -y v2015 for each beam X
% The name of the output file should be the same as for the plan file 
% (for matRad_calc4dDose).

% calculate 4D dose
% Make sure that the correct pln, dij and stf are loeaded in the workspace.
[resultGUI, delivery] = matRad_calc4dDose(ct,pln,dij,stf,cst,resultGUI,fileName);

% Plot the result in comparison to the static dose
slice = round(pln.isoCenter(1,3)./ct.resolution.z);
figure
subplot(2,2,1)
imagesc(resultGUI.RBExDose(:,:,slice)),colorbar, colormap(jet);
title('Static dose distribution')
subplot(2,2,2)
imagesc(resultGUI.accRBExDose(:,:,slice)),colorbar, colormap(jet);
title('4D dose distribution')
subplot(2,2,3)
imagesc(resultGUI.RBExDose(:,:,slice) - resultGUI.accRBExDose(:,:,slice)) ,colorbar, colormap(jet);
title('Difference')
