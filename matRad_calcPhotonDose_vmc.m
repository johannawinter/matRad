function dij = matRad_calcPhotonDose_vmc(ct,stf,pln,cst)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad vmc++ photon dose calculation wrapper
% 
% call
%   dij = matRad_calcPhotonDose(ct,stf,pln,cst,visBool)
%
% input
%   ct:         ct cube
%   stf:        matRad steering information struct
%   pln:        matRad plan meta information struct
%   cst:        matRad cst struct
%
% output
%   dij:        matRad dij struct
%
% References
%   
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015, Mark Bangert, on behalf of the matRad development team
%
% m.bangert@dkfz.de
%
% This file is part of matRad.
%
% matrad is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free
% Software Foundation, either version 3 of the License, or (at your option)
% any later version.
%
% matRad is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License in the
% file license.txt along with matRad. If not, see
% <http://www.gnu.org/licenses/>.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if visBool not set toogle off visualizationpat
if nargin < 5
    visBool = 0;
end

% initialize waitbar
figureWait=waitbar(0,'photon dij-calculation..');
% meta information for dij
dij.numOfBeams         = pln.numOfBeams;
dij.numOfVoxels        = pln.numOfVoxels;
dij.resolution         = ct.resolution;
dij.numOfRaysPerBeam   = [stf(:).numOfRays];
dij.totalNumOfRays     = sum(dij.numOfRaysPerBeam);
dij.totalNumOfBixels   = sum([stf(:).totalNumOfBixels]);
dij.dimensions         = pln.voxelDimensions;

%% (A)
% set environment variables for vmc++
VMCPath     = fullfile(pwd , 'vmc++', '');
runsPath    = fullfile(VMCPath, 'runs', '');
phantomPath = fullfile(VMCPath, 'phantoms');

setenv('vmc_home',VMCPath);
setenv('vmc_dir',runsPath);
setenv('xvmc_dir',VMCPath);

% export CT cube as binary file for vmc++
matRad_export_CT_vmc(ct, fullfile(phantomPath, 'matRad_CT.ct'));

% set general vmc++ parameters
% 1 source
VMC_options.beamlet_source.my_name       = 'source 1';                          % name of source
VMC_options.beamlet_source.monitor_units = 1;                                   % ?
VMC_options.beamlet_source.spectrum      = 'var_6MV.spectrum';                  % energy spectrum source (only used if no mono-Energy given)
VMC_options.beamlet_source.charge        = 0;                                   % charge (-1,0,1)
% 2 transport parameter
VMC_options.MC_parameter.automatic_parameter = 'yes';                           % if yes, automatic transport parameters are used
% 3 MC control
VMC_options.MC_control.ncase     = 5000;                                        % number of histories
VMC_options.MC_control.nbatch    = 10;                                          % ?
VMC_options.MC_control.rng_seeds = [9722,14369];                                % initialization of pseudo random number
% 4 variance reduction
VMC_options.variance_reduction.repeat_history      = 0.251;                     % 
VMC_options.variance_reduction.split_photons       = 'yes';                     % 
VMC_options.variance_reduction.photon_split_factor = -40;                       %
% 5 quasi random numbers
VMC_options.quasi.base      = 2;                                                %   
VMC_options.quasi.dimension = 60;                                               %
VMC_options.quasi.skip      = 1;                                                %
% 6 geometry
VMC_options.geometry.XYZ_geometry.method_of_input = 'CT-PHANTOM';                 % input method ('CT-PHANTOM', 'individual', 'groups') 
VMC_options.geometry.XYZ_geometry.CT              = 'CT';                         % name of geometry
VMC_options.geometry.XYZ_geometry.CT_file         = './phantoms/matRad_CT.ct';    % path of density matrix (only needed if input method is 'CT-PHANTOM')
% 7 scoring manager
VMC_options.scoring_options.start_in_geometry                = 'CT';    % geometry in which partciles start their transport
VMC_options.scoring_options.dose_options.score_in_geometries = 'CT'; % geometry in which dose is recorded
VMC_options.scoring_options.dose_options.score_dose_to_water = 'yes';   % if yes output is dose to water
VMC_options.scoring_options.output_options.name              = 'CT'; % geometry for which dose output is created (geometry has to be scored)
VMC_options.scoring_options.output_options.dump_dose         = 2;       % output format (1: format=float, Dose + deltaDose; 2: format=short int, Dose)
%% (A)

% set up arrays for book keeping
dij.bixelNum = NaN*ones(dij.totalNumOfRays,1);
dij.rayNum   = NaN*ones(dij.totalNumOfRays,1);
dij.beamNum  = NaN*ones(dij.totalNumOfRays,1);

% Allocate space for dij.physicalDose sparse matrix
dij.physicalDose = spalloc(numel(ct.cube),dij.totalNumOfBixels,1);

% Allocate memory for dose_temp cell array
numOfBixelsContainer = ceil(dij.totalNumOfBixels/10);
doseTmpContainer = cell(numOfBixelsContainer,1);

% take only voxels inside patient
V = unique([cell2mat(cst(:,4))]);

% Convert CT subscripts to linear indices.
[yCoordsV, xCoordsV, zCoordsV] = ind2sub(size(ct.cube),V);

xCoordsV = xCoordsV(:)*ct.resolution.x-pln.isoCenter(1);
yCoordsV = yCoordsV(:)*ct.resolution.y-pln.isoCenter(2);
zCoordsV = zCoordsV(:)*ct.resolution.z-pln.isoCenter(3);
coordsV  = [xCoordsV yCoordsV zCoordsV];

% set lateral cutoff value
lateralCutoff = 20; % [mm]

% toggle custom primary fluence on/off. if 0 we assume a homogeneous
% primary fluence, if 1 we use measured radially symmetric data
% useCustomPrimFluenceBool = 0;
% 
% %% kernel convolution
% % load polynomial fits for kernels ppKernel1, ppKernel2, ppKernel3
% load photonPencilBeamKernels_6MV.mat;
% 
% % Make a 2D grid extending +/-100mm with 0.1 mm resolution
% convLimits = 100; % [mm]
% convResolution = .5; % [mm]
% [X,Z] = meshgrid(-convLimits:convResolution:convLimits);
%                           
% % Evaluate piecewise polynomial kernels
% kernel1Mx = ppval(ppKernel1,sqrt(X.^2+Z.^2));
% kernel2Mx = ppval(ppKernel2,sqrt(X.^2+Z.^2));
% kernel3Mx = ppval(ppKernel3,sqrt(X.^2+Z.^2));
% 
% % Create zero matrix for the Fluence
% F = zeros(size(X));
% 
% % set bixel opening to one
% F(abs(X)<=pln.bixelWidth/2 & abs(Z)<=pln.bixelWidth/2) = 1;
% 
% % gaussian convolution of field to model penumbra
% sigmaGauss = 2.1/convResolution; % [mm] / see diploma thesis siggel 4.1.2
% gaussFilter =  convResolution^2/(2*pi*sigmaGauss^2) * exp( -(X.^2+Z.^2)/(2*sigmaGauss^2) );
% F = real(fftshift(ifft2(fft2( ifftshift(F) ).*fft2( ifftshift(gaussFilter) ))));
% 
% if ~useCustomPrimFluenceBool % pre-compute konvolution matrices for idealized homogeneous primary fluence
%     
%     % Display console message.
%     fprintf('matRad: Uniform primary photon fluence -> pre-compute kernel convolution... \n');    
% 
%     % 2D convolution of Fluence and Kernels in fourier domain
%     convMx1 = real(fftshift(ifft2(fft2( ifftshift(F) ).*fft2( ifftshift(kernel1Mx) ) )));
%     convMx2 = real(fftshift(ifft2(fft2( ifftshift(F) ).*fft2( ifftshift(kernel2Mx) ) )));
%     convMx3 = real(fftshift(ifft2(fft2( ifftshift(F) ).*fft2( ifftshift(kernel3Mx) ) )));
% 
%     % Creates an interpolant for kernes from vectors position X and Z
%     if exist('griddedInterpolant','class') % use griddedInterpoland class when available 
%         Interp_kernel1 = griddedInterpolant(X',Z',convMx1','linear');
%         Interp_kernel2 = griddedInterpolant(X',Z',convMx2','linear');
%         Interp_kernel3 = griddedInterpolant(X',Z',convMx3','linear');
%     else
%         Interp_kernel1 = @(x,y)interp2(X(1,:),Z(:,1),convMx1,x,y,'linear');
%         Interp_kernel2 = @(x,y)interp2(X(1,:),Z(:,1),convMx2,x,y,'linear');
%         Interp_kernel3 = @(x,y)interp2(X(1,:),Z(:,1),convMx3,x,y,'linear');
%     end
% 
% end
% 
% % define source position for beam eye view.
% sourcePoint_bev = [0 -pln.SAD 0];
% 
counter = 0;

fprintf('matRad: Photon dose calculation... ');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:dij.numOfBeams; % loop over all beams
    
    % gantry and couch roation matrices according to IEC 61217 standard
    % instead of moving the beam around the patient, we perform an inverse
    % rotation of the patient, i.e. we consider a beam's eye view
    % coordinate system
    
%     % rotation around Z axis (gantry)
%     rotMx_XY = [cosd(pln.gantryAngles(i)) -sind(pln.gantryAngles(i)) 0;
%                 sind(pln.gantryAngles(i))  cosd(pln.gantryAngles(i)) 0;
%                                         0                          0 1];
%     
%     % rotation around Y axis (couch)
%     rotMx_XZ = [ cosd(pln.couchAngles(i)) 0 sind(pln.couchAngles(i));
%                                         0 1                         0;
%                 -sind(pln.couchAngles(i)) 0 cosd(pln.couchAngles(i))];
%     
%     % rotate target coordinates around Y axis and then around Z axis
%     % i.e. 1st couch, 2nd gantry; matrix multiplication not cummutative
%     rot_coordsV = coordsV*rotMx_XZ*rotMx_XY;
%     
%     rot_coordsV(:,1) = rot_coordsV(:,1)-sourcePoint_bev(1);
%     rot_coordsV(:,2) = rot_coordsV(:,2)-sourcePoint_bev(2);
%     rot_coordsV(:,3) = rot_coordsV(:,3)-sourcePoint_bev(3);
%     
    for j = 1:stf(i).numOfRays % loop over all rays / for photons we only have one bixel per ray!
        
        counter = counter + 1;
        
%         if useCustomPrimFluenceBool % use custom primary fluence if specifried
%             
%             r     = sqrt( (X-stf(i).ray(j).rayPos(1)).^2 + (Z-stf(i).ray(j).rayPos(3)).^2 );
%             Psi   = interp1(primaryFluence(:,1),primaryFluence(:,2),r);
%             FxPsi = F .* Psi;
%         
%             % 2D convolution of Fluence and Kernels in fourier domain
%             convMx1 = real(fftshift(ifft2(fft2( ifftshift(FxPsi) ).*fft2( ifftshift(kernel1Mx) ) )));
%             convMx2 = real(fftshift(ifft2(fft2( ifftshift(FxPsi) ).*fft2( ifftshift(kernel2Mx) ) )));
%             convMx3 = real(fftshift(ifft2(fft2( ifftshift(FxPsi) ).*fft2( ifftshift(kernel3Mx) ) )));
% 
%             % Creates an interpolant for kernes from vectors position X and Z
%             if exist('griddedInterpolant','class') % use griddedInterpoland class when available 
%                 Interp_kernel1 = griddedInterpolant(X',Z',convMx1','linear');
%                 Interp_kernel2 = griddedInterpolant(X',Z',convMx2','linear');
%                 Interp_kernel3 = griddedInterpolant(X',Z',convMx3','linear');
%             else
%                 Interp_kernel1 = @(x,y)interp2(X(1,:),Z(:,1),convMx1,x,y,'linear');
%                 Interp_kernel2 = @(x,y)interp2(X(1,:),Z(:,1),convMx2,x,y,'linear');
%                 Interp_kernel3 = @(x,y)interp2(X(1,:),Z(:,1),convMx3,x,y,'linear');
%             end
% 
%         end

        % Display progress
        matRad_progress(counter,dij.totalNumOfBixels);
        % update waitbar only 100 times
        if mod(counter,round(dij.totalNumOfBixels/100)) == 0
            waitbar(counter/dij.totalNumOfBixels);
        end
        
        % remember beam and bixel number
        dij.beamNum(counter)  = i;
        dij.rayNum(counter)   = j;
        dij.bixelNum(counter) = j;
        
%         % Ray tracing for beam i and bixel j
%         [ix,radDepths,geoDists,latDistsX,latDistsZ] = matRad_calcRadGeoDists(ct.cube, ...
%                                                         V, ...
%                                                         pln.isoCenter, ...
%                                                         rot_coordsV, ...
%                                                         ct.resolution, ...
%                                                         stf(i).sourcePoint, ...
%                                                         stf(i).ray(j).targetPoint, ...
%                                                         sourcePoint_bev, ...
%                                                         stf(i).ray(j).targetPoint_bev, ...
%                                                         coordsV, ...
%                                                         lateralCutoff, ...
%                                                         visBool);
%         
%         % calculate photon dose for beam i and bixel j
%         bixelDose = matRad_calcPhotonDoseBixel(pln.SAD,m,betas, ...
%                                                Interp_kernel1,...
%                                                Interp_kernel2,...
%                                                Interp_kernel3,...
%                                                radDepths,...
%                                                geoDists,...
%                                                latDistsX,...
%                                                latDistsZ);

        %% (B)
        % set ray specific vmc++ parameters
        % a) change coordinate system (Isocenter cs-> physical cs) and units mm -> cm
        ray_corner_1 = (stf(i).ray(j).rayCorners_SCD(1,:) + pln.isoCenter)/10;              
        ray_corner_2 = (stf(i).ray(j).rayCorners_SCD(2,:) + pln.isoCenter)/10;
        ray_corner_3 = (stf(i).ray(j).rayCorners_SCD(3,:) + pln.isoCenter)/10; %vmc needs only three corners (counter-clockwise)
        beam_source  = (stf(i).sourcePoint + pln.isoCenter)/10;
        
        % b) swap x and y (CT-standard = [y,x,z])
        ray_corner_1 = ray_corner_1([2,1,3]);              
        ray_corner_2 = ray_corner_2([2,1,3]);
        ray_corner_3 = ray_corner_3([2,1,3]);
        beam_source  = beam_source([2,1,3]);
        
        % c) set vmc++ parameters
        VMC_options.beamlet_source.mono_energy                   = stf(i).ray(j).energy;                        % photon energy
        VMC_options.beamlet_source.beamlet_edges                 = [ray_corner_1,ray_corner_2,ray_corner_3];    % counter-clockwise beamlet edges
        VMC_options.beamlet_source.virtual_point_source_position = beam_source;                                 % virtual beam source position
        
        % create inputfile with vmc++ parameters
        outfile = 'MCpencilbeam_temp';
        matRad_create_VMC_input(VMC_options,fullfile(runsPath, [outfile,'.vmc']));
        
        % perform vmc++ simulation
        current = pwd;
        cd(VMCPath);
        dos(['start /NORMAL /B /WAIT ' fullfile('.', 'bin', 'vmc_Windows.exe') ' ' outfile '']); % (D)
        cd(current);
        
        % import calculated dose
        [bixelDose,~] = matRad_read_dose_vmc(fullfile(VMCPath, 'runs',...
                                             [outfile, '_', VMC_options.scoring_options.output_options.name, '.dos']));
        
        % Save dose for every bixel in cell array
        doseTmpContainer{mod(counter-1,numOfBixelsContainer)+1,1} = sparse(V,1,bixelDose(V),numel(ct.cube),1);
        %% (B)
        
        % Save dose for every bixel in cell array
        % doseTmpContainer{mod(counter-1,numOfBixelsContainer)+1,1} = sparse(V(ix),1,bixelDose,numel(ct.cube),1);
                
        % save computation time and memory by sequentially filling the 
        % sparse matrix dose.dij from the cell array
        if mod(counter,numOfBixelsContainer) == 0 || counter == dij.totalNumOfBixels
            dij.physicalDose(:,(ceil(counter/numOfBixelsContainer)-1)*numOfBixelsContainer+1:counter) = [doseTmpContainer{1:mod(counter-1,numOfBixelsContainer)+1,1}];
        end
        
    end
end

%% (C)
% delete phantom and run files
delete(fullfile(phantomPath, 'matRad_CT.ct')); % phantom file
delete(fullfile(runsPath, [outfile,'.vmc']));  % vmc input file
delete(fullfile(runsPath, [outfile,'_',VMC_options.scoring_options.dose_options.score_in_geometries,'.dos'])); % vmc outputfile

%% (C)

close(figureWait);