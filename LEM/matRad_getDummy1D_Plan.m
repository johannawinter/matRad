function [ pln,cst ] = matRad_getDummy1D_Plan(Voxel,tissue)
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
% creates a dummy pln and cst for the matRad world

cst{1,1} = 0; cst{2,1} = 1;
cst{1,2} = 'Body'; cst{2,2} = 'Target';
cst{1,3} = 'OAR';  cst{2,3} = 'TARGET';
cst{1,4}{1} = find(Voxel.Position);  cst{2,4}{1} = find(Voxel.IxT);
cst{1,5}.alphaX       = tissue.sAlphaXnom;
cst{1,5}.betaX        = tissue.sBetaXnom;
cst{1,5}.Priority     = 2;
cst{2,5}.TissueClass  = 1;
cst{2,5}.alphaX       = tissue.sAlphaXnom;
cst{2,5}.betaX        = tissue.sBetaXnom;
cst{2,5}.Priority     = 1;
cst{2,5}.TissueClass  = 1;

cst{1,6}.penalty      = 10;
cst{1,6}.dose         = 1;
cst{1,6}.type         = 'square overdosing';
cst{1,6}.EUD          = NaN;
cst{1,6}.volume       = NaN;
cst{1,6}.robustness   = 'none';

cst{2,6}.penalty      = 500;
cst{2,6}.dose         = 3;
cst{2,6}.type         = 'square deviation';
cst{2,6}.EUD          = NaN;
cst{2,6}.volume       = NaN;
cst{2,6}.robustness   = 'none';

%%
pln.radiationMode   = 'carbon';
pln.bioOptimization = 'LEMIV_effect';
pln.numOfFractions  = 1;
pln.runDAO          = false;


end

