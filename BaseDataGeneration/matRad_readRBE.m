function [RBE] = matRad_readRBE(path)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad_readdEdx script to parse RBE files from TRiP
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015, Hans-Peter Wieser
%
% h.wieser@dkfz.de
%
% This file is NOT part of the official matRad release. 
% This file has to be used only for internal purposes! 
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script can be used to parse the *.rbe files
path = [path filesep 'RBE'];
folderInfo = dir([path filesep '*.rbe']);

for i = 1:length(folderInfo)

    fileID = fopen([path filesep folderInfo(i).name]);
    currentline = fgetl(fileID);
    % parse meta
    while ischar(currentline)

        if strfind(currentline,'projectile')
            break;
        end

        if strcmp(currentline(1),'!')
           vTmp = strsplit(currentline(2:end));
           if length(vTmp)>1
             if sum(strcmp(vTmp{1,1},{'alpha','beta','cut','rnucleus'}))>0
                vTmp{1,2} = str2num(vTmp{1,2});
             end
             RBE(i).(vTmp{1,1}) = vTmp{1,2};
           end
        end

        currentline = fgetl(fileID);
    end


    ParticleCnt = 1;
    % parse data block
    while ischar(currentline)

        if strfind(currentline,'projectile')>0
            Cnt = 1;
            InnerCnt = 1;
            vTmp = strsplit(currentline(2:end));
            % delete number from string - otherwise it cannot be used as
            % fieldname
            projectile = vTmp{1,2};
            currentParticle = (regexprep(projectile,'[^a-zA-Z]',''));
            RBE(i).particle{1,ParticleCnt} = currentParticle;
            ParticleCnt = ParticleCnt + 1;
        end

        if strcmp(currentline(1),'#')
            unit = strsplit(currentline(3:end));
            % save units to meta
            for j=1:length(unit)
               RBE(i).(['unit' num2str(j)]) = unit{j}; 
            end
        end

        if Cnt > 3
            row = strsplit(currentline(3:end));
            RBE(i).(currentParticle).Energy(InnerCnt) = str2double(row{1});
            RBE(i).(currentParticle).RBE(InnerCnt)    = str2double(row{2});
            InnerCnt = InnerCnt +1;
        end

        currentline = fgetl(fileID);
        Cnt = Cnt+1;
    end

end

end

