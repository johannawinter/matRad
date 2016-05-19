function [ct, cst, resultGUI] = matRad_convertTPSdata2matRad(PathTo_tps_data)

load(PathTo_tps_data);
addpath([ pwd filesep 'dicomImport'])

sNames = fieldnames(tps_data.structures);

for i = 1:numel(sNames)
    
    cst{i,1} = i-1;
    cst{i,2} = sNames{i,1};
    
    if ~isempty(regexpi(cst{i,2},'tv')) || ...
       ~isempty(regexpi(cst{i,2},'target')) || ...
       ~isempty(regexpi(cst{i,2},'gtv')) || ...
       ~isempty(regexpi(cst{i,2},'ctv')) || ...
       ~isempty(regexpi(cst{i,2},'ptv')) || ...
       ~isempty(regexpi(cst{i,2},'boost')) || ...
       ~isempty(regexpi(cst{i,2},'tumor'))
        cst{i,3} = 'TARGET';
        % default objectives for targets
        cst{i,6}(1).dose = 60; %  
        cst{i,6}(1).penalty = 800; %  
        cst{i,6}(1).type = 'square deviation';
        cst{i,5}.Priority = 1;
    else
        cst{i,3} = 'OAR';
        cst{i,6} = []; % define no OAR dummy objcetives   
        cst{i,5}.Priority = 2;
    end

    cst{i,4}{1} = find(tps_data.structures.(sNames{i,1}).indicator_mask);
    
    % set default parameter for biological planning
    cst{i,5}.TissueClass = 1; 
    cst{i,5}.alphaX = 0.1;
    cst{i,5}.betaX = 0.05;
   
end

CT_Ref         = tps_data.ct;
ct.cube{1}     = CT_Ref.cube;
ct.dicomInfo   = CT_Ref.dicomInfo;
ct.resolution  = CT_Ref.resolution;
ct.cubeDim     = size(ct.cube{1});
ct.numOfCtScen = 1;
ct             = matRad_calcWaterEqD(ct);
resultGUI.physicalDose = tps_data.dose;
