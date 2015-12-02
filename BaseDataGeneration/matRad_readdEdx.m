function [Meta, dEdx ] = matRad_readdEdx(path)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse TRiP's dEdx file
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
% This script can be used to parse the stopping power ratios from dEdx
% files

path = [path filesep 'DEDX' filesep];
folderInfo = dir([path '*.dedx']);

if isempty(folderInfo)
    error(' cant find a dedx file');
elseif length(folderInfo)>1
    warning('more than one dedx exists');
    fileID = fopen(folderInfo(1).name);
else
    fileID = fopen([path folderInfo(1).name]);
end

currentline = fgetl(fileID);

% parse Meta
while ischar(currentline)

    if strfind(currentline,'projectile')
        break;
    end
    
    if strcmp(currentline(1),'!')
       vTmp = strsplit(currentline(2:end));
       Meta.(vTmp{1,1}) = vTmp{1,2};
    end
    
    currentline = fgetl(fileID);
end


% parse data block
while ischar(currentline)

    if strfind(currentline,'projectile')>0
        Cnt = 1;
        InnerCnt = 1;
        currentParticle = '';
        vTmp = strsplit(currentline(2:end));
        % delete number from string - otherwise it cannot be used as
        % fieldname
        projectile = vTmp{1,2};
        currentParticle = (regexprep(projectile,'[^a-zA-Z]',''));
    end
    
    if strcmp(currentline(1),'#')
        unit = strsplit(currentline(3:end));
        % save units to Meta
        for i=1:length(unit)
           Meta.(['unit' num2str(i)]) = unit{i}; 
        end
    end
    
    if Cnt > 3
        row = strsplit(currentline(3:end));
        dEdx.(currentParticle).Energy(InnerCnt)    = str2double(row{1});
        dEdx.(currentParticle).dEdx(InnerCnt) = str2double(row{2});
        dEdx.(currentParticle).range(InnerCnt)= str2double(row{3});
        InnerCnt = InnerCnt +1;
    end
    
    currentline = fgetl(fileID);
    Cnt = Cnt+1;
end

data = textscan(fileID,'%s');

data = data{1,1};
data = data(28:end,:);
Cnt = 1;
InnerCnt=1;
CntParticle = 1;
for i = 1:length(data)
   if Cnt>length(data)
        break; 
   end
   
   if strcmp(data(Cnt,1),'!dedx')
        Cnt=Cnt+1;
        
        while ~strcmp(data(Cnt,1),'!projectile')
            dEdx.(dEdxarticles{CntParticle}).energy{InnerCnt} = str2double(data(Cnt,1));
            dEdx.(dEdxarticles{CntParticle}).dEdx{InnerCnt} = str2double(data(Cnt+1,1));
            dEdx.(dEdxarticles{CntParticle}).range{InnerCnt} = str2double(data(Cnt+2,1));
            InnerCnt = InnerCnt+1;
            Cnt = Cnt+3;
            if Cnt>length(data)
               break; 
            end
        end
        dEdx.(dEdxarticles{CntParticle}).energy =(cell2mat(dEdx.(dEdxarticles{CntParticle}).energy))';
        dEdx.(dEdxarticles{CntParticle}).dEdx =(cell2mat(dEdx.(dEdxarticles{CntParticle}).dEdx))';
        dEdx.(dEdxarticles{CntParticle}).range =(cell2mat(dEdx.(dEdxarticles{CntParticle}).range))';
        CntParticle = CntParticle+1;
        InnerCnt = 1;
   else
       Cnt = Cnt+1;
   end
end



end

