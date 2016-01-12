function matRad_create_batch_file(parallel_simulations, filepath)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad batchfile creation
%
% call
%   matRad_create_batch_file(parallel_simulations, filepath)
%
% input
%   parallel_simulations:   no of parallel simulations
%   filepath:               path where batchfile is created (this has to be the 
%                           path of the vmc++ folder)
%
%
% References
%   
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parallel_processes = cell(1,parallel_simulations);
for i=1:parallel_simulations
    parallel_processes{1,i} = ['start "" 9>"%lock%',num2str(i),'" .\bin\vmc_Windows.exe MCpencilbeam_temp_',num2str(i)];
end

batch_file = {...
    ['@echo off'],...
    ['setlocal'],...
    ['set "lock=%temp%\wait%random%.lock"'],...
    [''],...
    parallel_processes{:},...
    [''],...
    [':Wait for both processes to finish (wait until lock files are no longer locked)'],...
    ['1>nul 2>nul ping /n 2 ::1'],...
    ['for %%N in (',strjoin(arrayfun(@(x) num2str(x),(1:parallel_simulations),'UniformOutput',false),' '),') do ('],...
    ['  (call ) 9>"%lock%%%N" || goto :Wait'],...
    [') 2>nul'],...
    [''],...
    ['del "%lock%*"'],...
    [''],...
    ['echo Done - ready to continue processing']
    };

% write batch file
fid = fopen(filepath,'wt');
for i = 1 : length(batch_file)
  fprintf(fid,'%s\n',batch_file{i});
end
fclose(fid);

end
