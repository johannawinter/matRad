function sigmaSq = matRad_getHeterogeneityCorrSigmaSq(WET)
% WET [mm] corresponds to the water equivalent thickness of penetrated lung
% tissue
% sigmaSq [mm^2]

% a = 1.43 [mm/sqrt(cm)] from Riccardos thesis, independent from proton energy
% (smaller a tightens the curve (smaller degradation), R80 stays at same position)
% a = 1.60 [mm/sqrt(cm)] = 1.6/sqrt(10) [sqrt(mm)] fits the measurements better

% Pmod = a^2; 150-750 micrometer for swine lung from Witt et al.
% Pmod = 256 is equivalent to a = 1.6

Pmod = 256; % [micrometer]

sigmaSq = Pmod/1000 * WET;