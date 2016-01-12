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

% set random seed (enables reproducibility)
rng(0);

% set number of parallel simulations
max_parallel_simulations = 8;

% set relative dose cutoff
rel_Dose_cutoff = 10^(-3);

% set absolute calibration factor
% CALCULATION
% absolute_calibration_factor = 1/D(depth = 100,5mm) -> D(depth = 100,5mm) = 1Gy
% SETUP
% SAD = 1000mm, SCD = 500mm, bixelWidth = 5mm, IC = [240mm,240mm,240mm]
% fieldsize@IC = 105mm x 105mm, phantomsize = 81 x 81 x 81 = 243mm x 243mm x 243mm
% rel_Dose_cutoff = 10^(-3), ncase = 500000/bixel
absolute_calibration_factor_vmc = 99.818252282632300;

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

% export CT cube as binary file for vmc++
matRad_export_CT_vmc(ct, fullfile(phantomPath, 'matRad_CT.ct'));

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

counter                     = 0;
counter2                    = 0;
max_no_executed_simulations = 0;

fprintf('matRad: VMC++ photon dose calculation... ');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:dij.numOfBeams; % loop over all beams
       
    for j = 1:stf(i).numOfRays % loop over all rays / for photons we only have one bixel per ray!
        
        counter = counter + 1;
        
        % create different seeds for every bixel
        VMC_options.MC_control.rng_seeds = [randi(30000),randi(30000)];

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
        outfile = ['MCpencilbeam_temp_',num2str(mod(counter-1,max_parallel_simulations)+1)];
        matRad_create_VMC_input(VMC_options,fullfile(runsPath, [outfile,'.vmc']));
        
        if mod(counter,max_parallel_simulations) == 0 || counter == dij.totalNumOfBixels
            % create batch file (enables parallel processes)
            if counter == dij.totalNumOfBixels && mod(counter,max_parallel_simulations) ~= 0
                parallel_simulations = mod(counter,max_parallel_simulations);
            else
                parallel_simulations = max_parallel_simulations;
            end
            matRad_create_batch_file(parallel_simulations,fullfile(VMCPath,'run_parallel_simulations.bat'));
            
            % save max number of executed parallel simulations
            if parallel_simulations > max_no_executed_simulations 
                max_no_executed_simulations = parallel_simulations;
            end
            
            % perform vmc++ simulation
            current = pwd;
            cd(VMCPath);
            dos('run_parallel_simulations.bat')
            cd(current);
            
            for k = 1:parallel_simulations
                counter2 = counter2 + 1;
                
                % import calculated dose
                idx = regexp(outfile,'_');
                [bixelDose,~] = matRad_read_dose_vmc(fullfile(VMCPath, 'runs',...
                                                     [outfile(1:idx(2)),num2str(k), '_', VMC_options.scoring_options.output_options.name, '.dos']));

                % apply relative dose cutoff
                Dose_cutoff                        = rel_Dose_cutoff*max(bixelDose);
                bixelDose(bixelDose < Dose_cutoff) = 0;

                % apply absolute calibration factor
                bixelDose = bixelDose*absolute_calibration_factor_vmc;

                % Save dose for every bixel in cell array
                doseTmpContainer{mod(counter2-1,numOfBixelsContainer)+1,1} = sparse(V,1,bixelDose(V),numel(ct.cube),1);

                if mod(counter2,numOfBixelsContainer) == 0 || counter2 == dij.totalNumOfBixels
                dij.physicalDose(:,(ceil(counter2/numOfBixelsContainer)-1)*numOfBixelsContainer+1:counter2) = [doseTmpContainer{1:mod(counter2-1,numOfBixelsContainer)+1,1}];
                end
            end
        end
        
    end
end

% delete phantom and run files
delete(fullfile(VMCPath, 'run_parallel_simulations.bat')); % batch file
delete(fullfile(phantomPath, 'matRad_CT.ct')); % phantom file
for j=1:max_no_executed_simulations
    delete(fullfile(runsPath, ['MCpencilbeam_temp_',num2str(mod(j-1,max_parallel_simulations)+1),'.vmc'])); % vmc input file
    delete(fullfile(runsPath, ['MCpencilbeam_temp_',num2str(mod(j-1,max_parallel_simulations)+1),'_',VMC_options.scoring_options.dose_options.score_in_geometries,'.dos'])); % vmc outputfile
end

close(figureWait);