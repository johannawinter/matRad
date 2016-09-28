function [ tissue ] = LEM_getTissueParameter( CellLineName )
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LEM_getTissueParameter
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
% This function represents a simple biological cell line library 


switch CellLineName
    case{'V79'}
        tissue.sAlphaX = 0.1;  %Gy-1
        tissue.sBetaX  = 0.05; %Gy-2
        tissue.sDcut   = 30;   %Gy
        tissue.Smax    = tissue.sAlphaX + 2*tissue.sBetaX*tissue.sDcut;
        tissue.RadiusTarget_um = 5; % µm
        tissue.sAlphaBetaRatio = tissue.sAlphaX / tissue.sBetaX;
    case {'CHO'}
        tissue.sAlphaX = 0.228;  %Gy-1
        tissue.sBetaX  = 0.02;   %Gy-2
        tissue.sDcut   = 30;     %Gy
        tissue.Smax    = tissue.sAlphaX + 2*tissue.sBetaX*tissue.sDcut;
        tissue.RadiusTarget_um = 5.5; % µm
        tissue.sAlphaBetaRatio = tissue.sAlphaX / tissue.sBetaX;
    case {'MDACC'}
        tissue.sAlphaX = 0.101;  %Gy-1
        tissue.sBetaX  = 0.055;  %Gy-2
        tissue.sDcut   = 30;     %Gy
        tissue.Smax    = tissue.sAlphaX + 2*tissue.sBetaX*tissue.sDcut;
        tissue.RadiusTarget_um = 8; % µm
        tissue.sAlphaBetaRatio = tissue.sAlphaX / tissue.sBetaX;
        
end

end

