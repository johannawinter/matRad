function [ sData ] = matRad_ParseBioData(PathToBaseData,ParticlePrefix,visBool)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad fucntion to parse CNAO dEdx*alpha and dEdx*sqrtBeta curves from
% monte carlo simulations
%
% call
%   sData = matRad_ParseBioData(PathToBaseData,ParticlePrefix,visBool)
% example
%   sData = matRad_ParseBioData('C:\Users\User\baseData','C',visBool)  
%
% input
%   PathToBaseData:       path to folder containing dEdx*alpha and
%                         dEdx*sqrtbeta files for each energy, whereat the
%                         folder
%   visBool:              toggle on/off visualization (optional)
%
% output
%   sData:                 TODO
%
% References
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
% %%

% parse directories - each folder represents data for one specific
% alpha/beta ratio
Directories = dir([PathToBaseData filesep '*AB*']);
%% TODO sort Directories according to alpha-beta ratio

CntCellLines = length(Directories);
sData=cell(1,CntCellLines);

% loop over all folder 
for i=1:CntCellLines
    Cnt = 1;
    Headers = cell(0);
    if Directories(i).isdir
        
        % process files in the i-th folder
        Path = [PathToBaseData filesep Directories(i).name];
        Files = dir([Path filesep '*' ParticlePrefix '_*.txt']);
        
        % loop over all files in folder i
        for j=1:length(Files)
            
            % construct complete prefix for dEdx*alpha and for
            % dEdx*sqrtBeta
            CompletePrefix = [ParticlePrefix '_E' num2str(sscanf(Files(j).name, [ParticlePrefix '_E%d'])) '_'];
            
            % find file pair for specific energy and add meta data to
            % header info
            if isempty(find(ismember(Headers,CompletePrefix)))
                %TODO sort FoundFilePair according to Energy
                FoundFilePair = dir([Path filesep CompletePrefix '*']);
                Headers{Cnt}=CompletePrefix;
            else
                FoundFilePair = [];
            end
                  
            % the length is supposed to be two as two files files per
            % energy shall be available
            for l=1:length(FoundFilePair)

                if  ~isempty(strfind(FoundFilePair(l).name, 'alpha')) && ~isempty(strfind(FoundFilePair(l).name,  Headers{Cnt}))
                    Import =importdata([Path filesep  FoundFilePair(l).name]);
                    
                    sData{1,i}{Cnt}.energy(:,1) = sscanf(char(Import.textdata(4,1)),'!Energy: %f');
                    sData{1,i}{Cnt}.depths(:,1) = Import.data(:,1)*10;
                    sData{1,i}{Cnt}.dEdxA(:,i) = Import.data(:,2);                    
                    sData{1,i}{Cnt}.alphaBetaRatio=Directories(i).name;
                
                end

                if ~isempty(strfind(FoundFilePair(l).name, 'beta')) && ~isempty(strfind(FoundFilePair(l).name,  Headers{Cnt}))
                     Import =importdata([Path filesep  FoundFilePair(l).name]);
                     sData{1,i}{Cnt}.dEdxB(:,i) = Import.data(:,2);
                end

                if l==length(FoundFilePair)
                    Cnt=Cnt+1;
                end
                
            end

            
        end
    end
end
 

%% convert cell to matrix...
A= (cell2mat(sData{1,1}));
vEnergies=[A(:).energy];
[val idx] = sort(vEnergies);
 
sData= (cell2mat(sData(idx)));

if visBool
    figure,hold on,xlabel('depth in cm'), ylabel('Gy-1')
    for i = 1:5:length(sData)
         str =['dE/dx * alpha with energy: '  num2str(sData(i).energy)];

        plot(sData(i).depths,sData(i).dEdxA), title(str,'FontSize',18),grid on
        waitforbuttonpress
    end


     figure,hold on,xlabel('depth in cm'), ylabel('Gy-1')
    for i = 1:5:length(sData)
         str =['dE/dx * sqrt(beta) with energy: '  num2str(sData(i).Energy)];

        plot(sData(i).Depth,sData(i).dEdxB), title(str,'FontSize',18),grid on
        waitforbuttonpress
    end
end


end
 
