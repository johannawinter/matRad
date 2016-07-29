function dAcc = matRad_doseAcc(d,ct,accMethod)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad dose accumulation function
% 
% call
%   dAcc = matRad_doseAcc(d,dvf)
%
% input
%   d:                  4D dose cube, size(d) = [dimX dimY dimZ nPhases]
%   ct:                 matRad ct struct inclduing 4d ct, deformation
%                       vector fields, and meta information
%
%   +++ Attention +++   the deformation vector fields are in [mm]
%
% output
%   dAcc:               accumulated dose, size(dAcc) = [dimX dimY dimZ]kes place)
%
% References
%   [1] http://iopscience.iop.org/0031-9155/59/21/6401/
%   [2] http://iopscience.iop.org/0031-9155/59/1/173/
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015, Mark Bangert
%
% m.bangert@dkfz.de
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


nPhases = size(d.physicalDose,4);
dimensions = [size(d.physicalDose,1) size(d.physicalDose,2) size(d.physicalDose,3)];

xGridVec = 1:dimensions(1);
yGridVec = 1:dimensions(2);
zGridVec = 1:dimensions(3);

dAcc = zeros(dimensions);

% DDM:  direct dose mapping
% DDMM: divergent dose mapping
% EMT:  energy mass transfer algorithm
if nargin < 3 % set default accumulation method
    accMethod = 'DDM';
end

    
if strcmp(accMethod,'DDM')
    
    if ~strcmp(ct.dvfType,'pull');
        error('dose accumulation via direct dose mapping (DDM) requires pull dvfs');
    end
    
    [Y,X,Z] = meshgrid(xGridVec,yGridVec,zGridVec);
 
    for i = 1:nPhases
        
        dvf_x_i = squeeze(ct.dvf(1,:,:,:,i));
        dvf_y_i = squeeze(ct.dvf(2,:,:,:,i));
        dvf_z_i = squeeze(ct.dvf(3,:,:,:,i));
        
        ix  = d.physicalDose(:,:,:,i) > 0;
                
        d_ref = interp3(yGridVec,xGridVec',zGridVec,d.physicalDose(:,:,:,i), ...
                        Y(ix) + dvf_y_i(ix)/ct.resolution.y, ...
                        X(ix) + dvf_x_i(ix)/ct.resolution.x, ...
                        Z(ix) + dvf_z_i(ix)/ct.resolution.z, ...
                        'linear',0);
                    
        dAcc(ix) = dAcc(ix) + d_ref;
        
    end
    
elseif strcmp(accMethod,'EMT')
   
    if ~strcmp(ct.dvfType,'push');
        error('dose accumulation via interpolation requires push dvfs');
    end

    [X,Y,Z] = ndgrid(xGridVec,yGridVec,zGridVec);

    for i = 1:nPhases
        
        m_i     = ct.cube(:,:,:,i);
        e_i     = d.physicalDose(:,:,:,i).*m_i;
        dvf_x_i = squeeze(ct.dvf(1,:,:,:,i))/ct.resolution.x;
        dvf_y_i = squeeze(ct.dvf(2,:,:,:,i))/ct.resolution.y;
        dvf_z_i = squeeze(ct.dvf(3,:,:,:,i))/ct.resolution.z;
        
        ix = e_i>0;
        
        m_ref = zeros(dimensions);
        e_ref = zeros(dimensions);

        X_i = X(ix) + dvf_x_i(ix);
        Y_i = Y(ix) + dvf_y_i(ix);
        Z_i = Z(ix) + dvf_z_i(ix);
        
        ix_i = sub2ind(dimensions,floor(X_i)  ,floor(Y_i)+1,floor(Z_i)  );
        overlap = (floor(X_i)-X_i+1) .* (Y_i-floor(Y_i))   .* (floor(Z_i)-Z_i+1);
        
        m_ref(ix_i) = m_ref(ix_i) + overlap .* m_i(ix);
        e_ref(ix_i) = e_ref(ix_i) + overlap .* e_i(ix);
        
        ix_i = sub2ind(dimensions,floor(X_i)  ,floor(Y_i)+1,floor(Z_i)+1);
        overlap = (floor(X_i)-X_i+1) .* (Y_i-floor(Y_i))   .* (Z_i-floor(Z_i));
        
        m_ref(ix_i) = m_ref(ix_i) + overlap .* m_i(ix);
        e_ref(ix_i) = e_ref(ix_i) + overlap .* e_i(ix);
        
        ix_i = sub2ind(dimensions,floor(X_i)  ,floor(Y_i)  ,floor(Z_i)+1);
        overlap = (floor(X_i)-X_i+1) .* (floor(Y_i)-Y_i+1) .* (Z_i-floor(Z_i));
        
        m_ref(ix_i) = m_ref(ix_i) + overlap .* m_i(ix);
        e_ref(ix_i) = e_ref(ix_i) + overlap .* e_i(ix);
        
        ix_i = sub2ind(dimensions,floor(X_i)  ,floor(Y_i)  ,floor(Z_i)  );
        overlap = (floor(X_i)-X_i+1) .* (floor(Y_i)-Y_i+1) .* (floor(Z_i)-Z_i+1);
        
        m_ref(ix_i) = m_ref(ix_i) + overlap .* m_i(ix);
        e_ref(ix_i) = e_ref(ix_i) + overlap .* e_i(ix);
        
        ix_i = sub2ind(dimensions,floor(X_i)+1,floor(Y_i)+1,floor(Z_i)  );
        overlap = (X_i-floor(X_i))   .* (Y_i-floor(Y_i))   .* (floor(Z_i)-Z_i+1);
        
        m_ref(ix_i) = m_ref(ix_i) + overlap .* m_i(ix);
        e_ref(ix_i) = e_ref(ix_i) + overlap .* e_i(ix);
        
        ix_i = sub2ind(dimensions,floor(X_i)+1,floor(Y_i)+1,floor(Z_i)+1);
        overlap = (X_i-floor(X_i))   .* (Y_i-floor(Y_i))   .* (Z_i-floor(Z_i));
        
        m_ref(ix_i) = m_ref(ix_i) + overlap .* m_i(ix);
        e_ref(ix_i) = e_ref(ix_i) + overlap .* e_i(ix);
        
        ix_i = sub2ind(dimensions,floor(X_i)+1,floor(Y_i)  ,floor(Z_i)+1);
        overlap = (X_i-floor(X_i))   .* (floor(Y_i)-Y_i+1) .* (Z_i-floor(Z_i));
        
        m_ref(ix_i) = m_ref(ix_i) + overlap .* m_i(ix);
        e_ref(ix_i) = e_ref(ix_i) + overlap .* e_i(ix);
        
        ix_i = sub2ind(dimensions,floor(X_i)+1,floor(Y_i)  ,floor(Z_i)  );
        overlap = (X_i-floor(X_i))   .* (floor(Y_i)-Y_i+1) .* (floor(Z_i)-Z_i+1);
                
        m_ref(ix_i) = m_ref(ix_i) + overlap .* m_i(ix);
        e_ref(ix_i) = e_ref(ix_i) + overlap .* e_i(ix);
        
        dAcc = dAcc + e_ref./m_ref;
                
    end
    
elseif strcmp(accMethod,'DDMM')
    
    % this implementation is experimental
    warning('The implementation of divergent dose mapping (DDMM) has never been debugged. Use with utmost care...');
    
    if ~strcmp(ct.dvfType,'pull');
        error('dose accumulation via divergent dose mapping requires pull dvfs');
    end

      
    massWeightedDose = d.physicalDose .* ct.cube;
    
    l_x = 2; % number of sub samples along x dir
    l_y = 3; % number of sub samples along y dir
    l_z = 4; % number of sub samples along z dir
    l   = l_x*l_y*l_z;
        
    x_steps = round(linspace(0,dimensions(1),l_x+1));
    y_steps = round(linspace(0,dimensions(2),l_y+1));
    z_steps = round(linspace(0,dimensions(3),l_z+1));
    
    counter = 0;
    
    for i = 1:nPhases
        
        curr_dm = massWeightedDose(:,:,:,i);
        curr_ct = ct.cube(:,:,:,i);
        
        dm_i = zeros(dimensions);
        m_i  = zeros(dimensions);
        
        
        for x = 1:l_x
            
            xStart = xGridVec(x_steps(x)+1) - 1/2 + 1/2/l_x;
            xEnd   = xGridVec(x_steps(x+1)) + 1/2 - 1/2/l_x;
            
            highRes_xGridVec = xStart:1/l_x:xEnd;
            
            for y = 1:l_y

                yStart = yGridVec(y_steps(y)+1) - 1/2 + 1/2/l_y;
                yEnd   = yGridVec(y_steps(y+1)) + 1/2 - 1/2/l_y;
    
                highRes_yGridVec = yStart:1/l_y:yEnd;

                for z = 1:l_z
                    
                    zStart = zGridVec(z_steps(z)+1) - 1/2 + 1/2/l_z;
                    zEnd   = zGridVec(z_steps(z+1)) + 1/2 - 1/2/l_z;
    
                    highRes_zGridVec = zStart:1/l_z:zEnd;
                    
                    % interpolation of sub cube of dose to high res
                    %highRes_dose_jkl = interp3(xGridVec,yGridVec',zGridVec, ...
                    %                           d.physicalDose(:,:,:,i), ...
                    %                           xGridVec_j,yGridVec_k',zGridVec_l,'linear',0);
                    
                    % interpolation of sub cube of deformation vector fields to high res
                    highRes_dvf_x    = interp3(yGridVec,xGridVec',zGridVec, ...
                                               squeeze( ct.dvf(1,:,:,:,i) ), ...
                                               highRes_yGridVec,highRes_xGridVec',highRes_zGridVec,'linear',0);
                    
                    highRes_dvf_y    = interp3(yGridVec,xGridVec',zGridVec, ...
                                               squeeze( ct.dvf(2,:,:,:,i) ), ...
                                               highRes_yGridVec,highRes_xGridVec',highRes_zGridVec,'linear',0);
                    
                    highRes_dvf_z    = interp3(yGridVec,xGridVec',zGridVec, ...
                                               squeeze( ct.dvf(3,:,:,:,i) ), ...
                                               highRes_yGridVec,highRes_xGridVec',highRes_zGridVec,'linear',0);
                    
                    % conversion to indices
                    [highRes_X,highRes_Y,highRes_Z] = ndgrid(highRes_xGridVec,highRes_yGridVec',highRes_zGridVec);
                    
                    origin_ix = sub2ind(dimensions,round(highRes_X(:)),round(highRes_Y(:)),round(highRes_Z(:)));
                    
                    highRes_dvf_x = round(highRes_X+highRes_dvf_x/ct.resolution.x);
                    highRes_dvf_y = round(highRes_Y+highRes_dvf_y/ct.resolution.y);
                    highRes_dvf_z = round(highRes_Z+highRes_dvf_z/ct.resolution.z);
                    
                    target_ix = sub2ind(dimensions,highRes_dvf_x(:),highRes_dvf_y(:),highRes_dvf_z(:));

                    
                    for j = 1:numel(target_ix)
                        
                        dm_i(target_ix(j)) = dm_i(target_ix(j)) + curr_dm(origin_ix(j));
                        m_i (target_ix(j)) = m_i (target_ix(j)) + curr_ct(origin_ix(j));
                        
                    end
                    
                    counter = counter + 1;
                    
                end
            end
        end

    end
        
end