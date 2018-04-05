function [ rGrid ] = LEM_getOrginalIntegrationSteps(r_min, r_max,visBool )
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
% This function determines a grid for the shell integration


    rGrid = 0;
	log_step_factor		= 1.003; % equals 0.5 %
	lin_step_width  	= 0.004; % equals 5nm
	r_trans         	= 0.35;  % transition between log and lin steps
                                 % currently at 0.5 um
                                 
    linGrid    = [];
    logGrid    = [];
    logGrid(1) = r_min;
    
    %
   if r_min < r_trans && r_trans < r_max
         n_log_steps = floor(log(r_trans/r_min)/log(log_step_factor));
   elseif r_min < r_trans && r_trans > r_max  
         n_log_steps = floor(log(r_max/r_min)/log(log_step_factor));
   else 
         n_log_steps = 0;
   end
   % create log grid
   for j = 1 : n_log_steps
       logGrid(j+1) =  logGrid(j) * log_step_factor;
   end
    
    
    if r_max >= r_trans && r_trans > r_min
       n_lin_steps = floor((r_max-r_trans) / lin_step_width); 
       linGrid = linspace(r_trans,r_max,n_lin_steps);
    elseif r_max >= r_trans && r_trans < r_min
       n_lin_steps = floor((r_max-r_min) / lin_step_width); 
       linGrid = linspace(r_min,r_max,n_lin_steps);
    end
    
    rGrid = [logGrid linGrid];
    

    % consistency checks
    
    if min(rGrid) < r_min || max(rGrid) > r_max
        warning('shell integration steps exceed boarders')
    end
    
    if visBool
        figHandles = get(0,'Children');
        NewFigure = true;
        for j = 1 : length(figHandles)
            if ~isempty(findstr(figHandles(j).Name,'ShellIntegrationSteps'))
                NewFigure = false;
            end
        end
    
        if NewFigure
            Num = numel(rGrid)-numel(linGrid);
            figure('Name','ShellIntegrationSteps'),set(gcf,'Color',[1 1 1]); 
            plot(rGrid,'LineWidth',3),grid on, grid minor, hold on
            plot(logGrid,'--','LineWidth',2)
            plot(Num:1:numel(rGrid)-1,linGrid,'--','LineWidth',2)
            xlabel('number of grid points','Interpreter','Latex');
            ylabel('shell integration grid points','Interpreter','Latex');
            title(['Total number of points: ' num2str(numel(rGrid))],'Interpreter','Latex')
            if isempty(linGrid)
                legend('total grid','log grid')
            else
                legend('total grid','log grid','lin grid')
            end
            
        end
    end
    
end

