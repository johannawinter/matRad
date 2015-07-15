function [ rGrid ] = LEM_getOrginalIntegrationSteps(r_min, r_max )

    rGrid = 0;
	log_step_factor		= 1.005; % equals 0.5 %
	lin_step_width  	= 0.005; % equals 5nm
	r_trans         	= 0.50;  % transition between log and lin steps
                                 % currently at 0.5 um
	n_log_steps			= 0;
	n_lin_steps			= 0;

	if r_min < r_trans 
		n_log_steps = floor(log(r_trans/r_min)/log(log_step_factor)) + 1;
    end
    
    if r_max >= r_trans
       n_lin_steps = floor((r_max-r_trans) / lin_step_width); 
    end
    
    n_steps	= n_log_steps + n_lin_steps;
    
    if n_steps == 0
        return
    end
    
    if n_log_steps == 0
        n_lin_steps = n_lin_steps-1;
    end
    
    if n_lin_steps == 0
        n_log_steps = n_log_steps-1;
    end
    
    if isnan(n_steps) || n_steps == 0 || isinf(n_steps)
       st = 2; 
    end
    rGrid = zeros(n_steps,1);
    rGrid(1) = r_min;
    
    %% generate grid, first make logaritmic steps then >r_trans make linear steps
    if n_log_steps > 0 
        for i = 1:n_log_steps
            rGrid(i+1) = rGrid(i)*log_step_factor;
        end
    end
    
    if n_lin_steps > 0
       for i = 1:n_lin_steps-1
           rGrid(i+1+n_log_steps) = rGrid(i+n_log_steps) + lin_step_width;
       end
    end
    
    rGrid(end) = r_max;
    
end

