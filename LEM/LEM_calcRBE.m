function [ RBE ] = LEM_calcRBE( sAlphaX,sBetaX,sAlphaI,sBetaI,physDose )
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LEM_calcRBE script
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
% This script calculates the RBE for given radiosensitivity parameters (photons and ions)
% as well as a given physical dose


RBE = (-sAlphaX + sqrt(sAlphaX.^2 + 4.*sBetaX.*(sAlphaI.*physDose + sBetaI.*physDose.^2)))./(2.*sBetaX.*physDose);

end

