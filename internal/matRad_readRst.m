function [stf, pln, w] = matRad_readRst(ct,pln,filename)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad read RST file and generates stearing file
%
% call
%   [stf, pln] = matRad_readRst(pln,filename)
%
% input
%   ct:         matRad's ct struct
%   pln:        plan struct
%   filename:   name of RST file
%
% output
%   stf: steering struct
%   pln: pln struct
%   w:   weight vector
%
% References
%
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

% read RST into cell array
RST = strsplit(fileread(filename),char(10))';
%RST = strsplit(fileread(filename),'\n')';

% Parser only applicable to single gantry angles! The expected gantry
% angles have to specified in the pln sturct, abort of more than one angle
% specified
numOfGantryAnglesInRST = sum(strncmp('gantryangle',RST,10));
if numel(pln.gantryAngles) > 1 || numOfGantryAnglesInRST > 1
    error('parser not validated for more than one beam direction\n');
end

% prepare structures necessary for particles
fileName = [pln.radiationMode '_' pln.machine];
try
   load([pwd filesep fileName]);
   pln.SAD = machine.meta.SAD;
catch
   error(['Could not find the following machine file: ' fileName ]); 
end

% Define steering file like struct. Prellocating for speed.
stf = struct;
w   = [];

% get number of machines in RST
machineStartIx = [find(strncmp('machine#',RST,7)) numel(RST)];
numOfMachinesInRST = numel(machineStartIx)-1;

load([pln.radiationMode '_' pln.machine])

% generate voi surface cube for SSD calculation
stf.SAD = machine.meta.SAD;


% loop over all machines
for i = 1:numOfMachinesInRST
    
    RSTi = RST(machineStartIx(i):machineStartIx(i+1));
    stf(i).isoCenter     = pln.isoCenter;

    % check if couch and gantry angles are consitent
    gantryAngleIx = find(strncmp('gantryangle',RSTi,10));
    couchAngleIx = find(strncmp('couchangle',RSTi,9));
    gantryAngleRST = cell2mat(textscan(RSTi{gantryAngleIx},'gantryangle %f'));
    couchAngleRST = cell2mat(textscan(RSTi{couchAngleIx},'couchangle %f'));
    if gantryAngleRST == pln.gantryAngles(i)
        stf(i).gantryAngle  = pln.gantryAngles(i);
    else
        warning(['Gantry angle mismatch. A Gantry angle of ' num2str(gantryAngleRST) ' degrees used from the RST file.']);
        stf(i).gantryAngle  = gantryAngleRST;
        pln.gantryAngles(i) = gantryAngleRST;
    end
    if couchAngleRST == pln.couchAngles(i)
        stf(i).couchAngle  = pln.couchAngles(i);
    else
        warning(['Couch angle mismatch. A Couch angle of ' num2str(couchAngleRST) ' degrees used from the RST file.']);
        stf(i).couchAngle  = couchAngleRST;
        pln.couchAngles(i) = couchAngleRST;
    end
    
    % check lateral spot spacing i.e. bixel width
    bixelWidthIx = find(strncmp('stepsize',RSTi,8));
    if numel(unique(RSTi(bixelWidthIx))) > 1
        warning(['Different bixel width used for different energies. Setting pln.bixelWidth = NaN.\n']);
        pln.bixelWidth = 3;
        stf.bixelWidth = 3;
        warning('setting bixel width to 3');
    elseif numel(unique(cell2mat(textscan(RSTi{bixelWidthIx(1)},'stepsize %f %f')))) > 1
        warning(['Different bixel width used in x and y direction. Setting pln.bixelWidth = NaN.\n']);
        pln.bixelWidth = NaN;
        stf.bixelWidth = NaN;
    elseif unique(cell2mat(textscan(RSTi{bixelWidthIx(1)},'stepsize %f %f'))) ~= pln.bixelWidth
        warning(['Different bixel width used in RST file. Setting pln.bixelWidth = ' ...
                  num2str(unique(cell2mat(textscan(RSTi{bixelWidthIx(1)},'stepsize %f %f')))) '.\n']);
        pln.bixelWidth = unique(cell2mat(textscan(RSTi{bixelWidthIx(1)},'stepsize %f %f')));
        stf.bixelWidth = unique(cell2mat(textscan(RSTi{bixelWidthIx(1)},'stepsize %f %f')));
    else
        stf.bixelWidth = unique(cell2mat(textscan(RSTi{bixelWidthIx(1)},'stepsize %f %f')));        
    end
    
    % set radiation modality
    projectile = textscan(RSTi{2,:},'projectile %s');
    if strcmp(pln.radiationMode,'protons') && strcmp(projectile{1} , '1H')
        stf(i).radiationMode = pln.radiationMode;
    elseif strcmp(pln.radiationMode,'carbon') && (strcmp(projectile{1}, 'C') || strcmp(projectile{1} , '12C'))
        stf(i).radiationMode = pln.radiationMode;
    elseif strcmp(projectile{1} , '12C') || strcmp(projectile{1}, 'C')
        warning(['Radiationmode mismatch. Carbon used as radiation mode.']);
        stf(i).radiationMode = 'carbon';
        pln.radiationMode = 'carbon';
    elseif strcmp(projectile{1} , '1H')
        warning(['Radiationmode mismatch. Protons used as radiation mode.']);
        stf(i).radiationMode = 'protons';
        pln.radiationMode = 'protons';
    else
        error(['Unknown radiationmode.']);
    end
    
    % get the number and positions of submachines
    submachineIx = find(strncmp('submachine#',RSTi,11));
    numOfSubmachines = numel(submachineIx);
    
    % total number of bixels
    numOfBixelsPerSubmachine = cell2mat(cellfun(@(x)textscan(x,'#points %f'),RSTi(submachineIx+3)));
    totalNumOfBixels = sum(numOfBixelsPerSubmachine);
    vector = zeros(totalNumOfBixels,4);
    
    % read in positions
    counter = 1;
    for j = 1:numOfSubmachines
        submachineAttributes = textscan(RSTi{submachineIx(j)},'submachine# %f %f %f %f');
        data = cellfun(@(x)textscan(x,'%f %f %f'),RSTi(submachineIx(j)+4:submachineIx(j)+3+numOfBixelsPerSubmachine(j)),'UniformOutput',false);
        for k = 1:numOfBixelsPerSubmachine(j)
            vector(counter,:) = [data{k}{1} data{k}{2} submachineAttributes{2} data{k}{3}/1000000];
            counter = counter + 1;
        end     
    end    
    
    % finds all unique rays and saves them in to the stf
    [RayPosTmp, ~, ic] = unique(vector(:,1:2), 'rows');
    [rows, ~] = size(RayPosTmp);
    for j = 1: rows
        stf(i).ray(j).rayPos_bev = [RayPosTmp(j,1) 0 RayPosTmp(j,2)];
        stf(i).ray(j).energy = [];
        stf(i).ray(j).weight = [];
    end
    
    % saves all energies and weights to their corresponding ray
    [rows, ~] = size(vector);
    for j = 1: rows
        k = ic(j);
        stf(i).ray(k).energy = [stf(i).ray(k).energy vector(j,3)];
        stf(i).ray(k).weight = [stf(i).ray(k).weight vector(j,4)];
    end
   
    
    % save the number of rays
    stf(i).numOfRays = numel(stf.ray);
    
    % Save ray and target position in beam eye´s view (bev)
    for j = 1:stf(i).numOfRays
        stf(i).ray(j).targetPoint_bev = [2*stf(i).ray(j).rayPos_bev(1) ...
                                         stf.SAD ...
                                         2*stf(i).ray(j).rayPos_bev(3)];
    end
    
    % source position in bev
    stf.sourcePoint_bev = [0 -stf.SAD 0];
  
    % compute coordinates in lps coordinate system, i.e. rotate beam
    % geometry around fixed patient
    
    % Rotation around Z axis (gantry)
    rotMx_XY_rotated = [ cosd(pln.gantryAngles(i)) sind(pln.gantryAngles(i)) 0;
                        -sind(pln.gantryAngles(i)) cosd(pln.gantryAngles(i)) 0;
                         0                         0                         1];
    
    % Rotation around Y axis (couch)
    rotMx_XZ_rotated = [ cosd(pln.couchAngles(i)) 0 -sind(pln.couchAngles(i));
                         0                        1                        0;
                         sind(pln.couchAngles(i)) 0 cosd(pln.couchAngles(i))];
    
    % Rotated Source point, first needs to be rotated around gantry, and then
    % couch.
    stf(i).sourcePoint =  stf.sourcePoint_bev * rotMx_XY_rotated * rotMx_XZ_rotated;
   
    % Save ray and target position in lps system.
    for j = 1:stf(i).numOfRays
        
        stf(i).ray(j).rayPos      = stf(i).ray(j).rayPos_bev      * rotMx_XY_rotated * rotMx_XZ_rotated;
        stf(i).ray(j).targetPoint = stf(i).ray(j).targetPoint_bev * rotMx_XY_rotated * rotMx_XZ_rotated;
        
        
        % ray tracing necessary to determine depth of the target
        [alpha,~,rho,~,~] = matRad_siddonRayTracer(stf(i).isoCenter, ...
                             ct.resolution, ...
                             stf(i).sourcePoint, ...
                             stf(i).ray(j).targetPoint, ...
                             [ct.cube]);
                         
                            
        DensityThresholdSSD = 0.05;              
        ixSSD = find(rho{1} > DensityThresholdSSD,1,'first');
                           
        if isempty(ixSSD)== 1
            warning('Surface for SSD calculation starts directly in first voxel of CT\n');
        end

        % calculate SSD
        stf(i).ray(j).SSD = 2 * stf(i).SAD * alpha(ixSSD);

        % book keeping & calculate focus index
        stf(i).numOfBixelsPerRay(j) = numel([stf(i).ray(j).energy]);
        currentMinimumFWHM = matRad_interp1(machine.meta.LUT_bxWidthminFWHM(1,:),...
                                     machine.meta.LUT_bxWidthminFWHM(2,:),...
                                     pln.bixelWidth);
        focusIx  =  ones(stf(i).numOfBixelsPerRay(j),1);
        [~, vEnergyIx] = min(abs(bsxfun(@minus,[machine.data.energy]',...
                        repmat(stf(i).ray(j).energy,length([machine.data]),1))));

        % get for each spot the focus index
        for k = 1:stf(i).numOfBixelsPerRay(j)                    
            focusIx(k) = find(machine.data(vEnergyIx(k)).initFocus.SisFWHMAtIso > currentMinimumFWHM,1,'first');
        end

        stf(i).ray(j).focusIx = focusIx';

    end
        
    % save total number of bixels
    numOfBixels = 0;
    for j = 1:numel(stf.ray)
        numOfBixels = numOfBixels + numel(stf.ray(j).energy);
        stf.numOfBixelsPerRay(j) = numel(stf.ray(j).energy);
        w = [w stf(i).ray(j).weight];
    end
    
    stf(i).totalNumOfBixels = numOfBixels;
    
end

    