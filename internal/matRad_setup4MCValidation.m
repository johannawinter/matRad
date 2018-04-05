function [ct,cst,pln] = matRad_setup4MCValidation(MCcube)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% call
%   [ct,cst,pln] = matRad_setup4MCValiadtion(MCcube)
%
% input
%   MCcube:         Monte Carlo dose cube
%
% output
%   ct:             matRad ct struct
%   cst:            matRad cst struct
%   pln:            matRad pln struct only containing iso center
%
% References
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015, Daniel Björkman, on behalf of the matRad development team
%
% d.bjoerkman@dkfz.de
%
% This file is not part of the official matRad release
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ct assuming water
ct.cube{1}     = ones(size(MCcube.cube));
ct.resolution  = MCcube.resolution;
ct.cubeDim     = size(ct.cube{1});
ct.numOfCtScen = 1;
% pln
pln.isoCenter = [MCcube.isoCenter];

% cst
cst       = cell(1,6);
cst{1,1}  = 0;
cst{1,2}  = 'BODY';
cst{1,3}  = 'OAR';
cst{1,4}{1} = (1:numel(ct.cube{1}))';
