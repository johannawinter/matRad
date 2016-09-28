function ObjVal = matRad_1D_bioEffectObjFunc(vWeight, dij, Voxel)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LEM_LET2Dose
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2016, Hans-Peter Wieser
%
% h.wieser@dkfz.de
%
% This file is NOT part of the official matRad release. 
% This file has to be used only for internal purposes! 
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% objective function calculation for the biological effect

ObjVal = (dij.alpha .* dij.dose) * vWeight + (sqrt(dij.beta) .* dij.dose * vWeight).^2  -  Voxel.presEffect';

ObjVal(ObjVal < 0 & Voxel.IxNT) = 0;
ObjVal = (Voxel.Penalty'.* ObjVal)' * ObjVal;