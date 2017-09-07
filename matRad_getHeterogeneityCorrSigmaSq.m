function sigmaSq = matRad_getHeterogeneityCorrSigmaSq(WET)
% WET [mm] corresponds to the water equivalent thickness of penetrated lung
% tissue

% sigmaSq [mm^2]
% a = 1.43 [mm/Sqr(cm)] from Riccardos thesis, independent from proton energy
% (smaller a tightens the curve (smaller degradation), peak stays at same
% position)
% a = 1.60 [mm/Sqr(cm)] fits the measurements better

aSq = (1.6/sqrt(10))^2; % [mm]

sigmaSq = aSq * WET;