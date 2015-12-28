function dij = matRad_calcPhotonDose_vmc(ct,stf,pln,cst)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad vmc++ photon dose calculation wrapper
% 
% call
%   dij = matRad_calcPhotonDose(ct,stf,pln,cst)
%
% input
%   ct:         matRad ct struct
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

% initialize waitbar
figureWait = waitbar(0,'VMC++ photon dij-calculation..');
% meta information for dij
dij.numOfBeams         = pln.numOfBeams;
dij.numOfVoxels        = pln.numOfVoxels;
dij.resolution         = ct.resolution;
dij.numOfRaysPerBeam   = [stf(:).numOfRays];
dij.totalNumOfRays     = sum(dij.numOfRaysPerBeam);
dij.totalNumOfBixels   = sum([stf(:).totalNumOfBixels]);
dij.dimensions         = pln.voxelDimensions;

% set environment variables for vmc++
VMCPath     = fullfile(pwd , 'vmc++');
runsPath    = fullfile(VMCPath, 'runs');
phantomPath = fullfile(VMCPath, 'phantoms');

setenv('vmc_home',VMCPath);
setenv('vmc_dir',runsPath);
setenv('xvmc_dir',VMCPath);

% export CT cube as binary file for vmc++
matRad_export_CT_vmc(ct, fullfile(phantomPath, 'matRad_CT.ct'));

% set general vmc++ parameters
% 1 source
VMC_options.beamlet_source.my_name       = 'source 1';                                               % name of source
VMC_options.beamlet_source.monitor_units = 1;                                                        % ?
VMC_options.beamlet_source.spectrum      = './spectra/Artiste_Geant4_8mm_phsp_tacke.spectrum';       % energy spectrum source (only used if no mono-Energy given)
VMC_options.beamlet_source.charge        = 0;                                                        % charge (-1,0,1)
% 2 transport parameter
VMC_options.MC_parameter.automatic_parameter = 'yes';                          % if yes, automatic transport parameters are used
% 3 MC control
VMC_options.MC_control.ncase     = 5000;                                       % number of histories
VMC_options.MC_control.nbatch    = 10;                                         % ?
VMC_options.MC_control.rng_seeds = [9722,14369];                               % initialization of pseudo random number
% 4 variance reduction
VMC_options.variance_reduction.repeat_history      = 0.251;                    % 
VMC_options.variance_reduction.split_photons       = 'yes';                    % 
VMC_options.variance_reduction.photon_split_factor = -40;                      %
% 5 quasi random numbers
VMC_options.quasi.base      = 2;                                               %   
VMC_options.quasi.dimension = 60;                                              %
VMC_options.quasi.skip      = 1;                                               %
% 6 geometry
VMC_options.geometry.XYZ_geometry.method_of_input = 'CT-PHANTOM';              % input method ('CT-PHANTOM', 'individual', 'groups') 
VMC_options.geometry.XYZ_geometry.CT              = 'CT';                      % name of geometry
VMC_options.geometry.XYZ_geometry.CT_file         = './phantoms/matRad_CT.ct'; % path of density matrix (only needed if input method is 'CT-PHANTOM')
% 7 scoring manager
VMC_options.scoring_options.start_in_geometry                = 'CT';           % geometry in which partciles start their transport
VMC_options.scoring_options.dose_options.score_in_geometries = 'CT';           % geometry in which dose is recorded
VMC_options.scoring_options.dose_options.score_dose_to_water = 'yes';          % if yes output is dose to water
VMC_options.scoring_options.output_options.name              = 'CT';           % geometry for which dose output is created (geometry has to be scored)
VMC_options.scoring_options.output_options.dump_dose         = 2;              % output format (1: format=float, Dose + deltaDose; 2: format=short int, Dose)

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

counter = 0;

fprintf('matRad: VMC++ photon dose calculation... ');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:dij.numOfBeams; % loop over all beams
       
    for j = 1:stf(i).numOfRays % loop over all rays / for photons we only have one bixel per ray!
        
        counter = counter + 1;

        % Display progress
        fprintf(['finished ' num2str(counter/dij.totalNumOfBixels*100) '%% \n']);
        % update waitbar only 100 times
        if mod(counter,round(dij.totalNumOfBixels/100)) == 0
            waitbar(counter/dij.totalNumOfBixels);
        end
        
        % remember beam and bixel number
        dij.beamNum(counter)  = i;
        dij.rayNum(counter)   = j;
        dij.bixelNum(counter) = j;
        
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
        %VMC_options.beamlet_source.mono_energy                   = stf(i).ray(j).energy;                       % photon energy
        VMC_options.beamlet_source.mono_energy                   = []                  ;                        % use photon spectrum
        VMC_options.beamlet_source.beamlet_edges                 = [ray_corner_1,ray_corner_2,ray_corner_3];    % counter-clockwise beamlet edges
        VMC_options.beamlet_source.virtual_point_source_position = beam_source;                                 % virtual beam source position
        
        % create inputfile with vmc++ parameters
        outfile = 'MCpencilbeam_temp';
        matRad_create_VMC_input(VMC_options,fullfile(runsPath, [outfile,'.vmc']));
        
        % perform vmc++ simulation
        current = pwd;
        cd(VMCPath);
        dos(['start /NORMAL /B /WAIT ' fullfile('.', 'bin', 'vmc_Windows.exe') ' ' outfile '']);
        cd(current);
        
        % import calculated dose
        [bixelDose,~] = matRad_read_dose_vmc(fullfile(VMCPath, 'runs',...
                                             [outfile, '_', VMC_options.scoring_options.output_options.name, '.dos']));
                                         
        % apply relative dose cutoff
        Dose_cutoff                        = 10^(-4)*max(bixelDose);
        bixelDose(bixelDose < Dose_cutoff) = 0;
        
        % apply conversion factor (enables comparability of dose calculations)
        bixelDose = bixelDose*91.876665940287400;
                                         
        % Save dose for every bixel in cell array
        doseTmpContainer{mod(counter-1,numOfBixelsContainer)+1,1} = sparse(V,1,bixelDose(V),numel(ct.cube),1);

        if mod(counter,numOfBixelsContainer) == 0 || counter == dij.totalNumOfBixels
            dij.physicalDose(:,(ceil(counter/numOfBixelsContainer)-1)*numOfBixelsContainer+1:counter) = [doseTmpContainer{1:mod(counter-1,numOfBixelsContainer)+1,1}];
        end
        
    end
end

% delete phantom and run files
delete(fullfile(phantomPath, 'matRad_CT.ct')); % phantom file
delete(fullfile(runsPath, [outfile,'.vmc']));  % vmc input file
delete(fullfile(runsPath, [outfile,'_',VMC_options.scoring_options.dose_options.score_in_geometries,'.dos'])); % vmc outputfile

close(figureWait);