function [cst,dvhFig,qiFig] = ringContoursDVH(patientID,cst,ct,pln,stf,resultGUI,margins,boolSave)
% add ring-shaped VOIs around PTV: inner (small), outer, large ring

% define margins if not specified
if ~exist('margins','var') || isempty(margins)
    margins{1} = 10;    % [mm]
    margins{2} = 20;
end

% find PTV and External in cst
for i = 1:size(cst,1)
    if strcmp(cst{i,2},'PTV')
        PTVix = i;
    elseif strcmp(cst{i,2},'External') || strcmp(cst{i,2},'Aussenkontur')
        ExternalIx = i;
    end
end

% assign ones to PTV voxels
try
    PTVvoxels = cst{PTVix,4}{:};
catch
    % if no PTV could be found, use ITV instead
    for i = 1:size(cst,1)
        if strcmp(cst{i,2},'ITV')
            PTVix = i;
        elseif strcmp(cst{i,2},'External') || strcmp(cst{i,2},'Aussenkontur')
            ExternalIx = i;
        end
    end
    
    warning('ITV is used instead of PTV.')
    PTVvoxels = cst{PTVix,4}{:};
end
PTVmask = zeros(ct.cubeDim);
PTVmask(PTVvoxels) = 1;

% add rings around PTV in margin thicknesses
vMargin         = cell(1,length(margins));
PTVenlargedVoi  = cell(size(vMargin));
PTVenlarged     = cell(size(vMargin));
ringPTV         = cell(size(vMargin));
for j = 1:length(margins)
    vMargin{j}.x = margins{j};
    vMargin{j}.y = margins{j};
    vMargin{j}.z = margins{j};
    
    PTVenlargedVoi{j} = matRad_addMargin(PTVmask,cst,ct.resolution,vMargin{j},true);
    PTVenlarged{j} = find(PTVenlargedVoi{j} > 0);
    ringPTV{j} = setdiff(PTVenlarged{j},PTVvoxels);
    
    % assign ring around PTV to cst
    cst{i+j,1} = i+j-1;
    cst{i+j,2} = ['ring PTV ' num2str(margins{j}) 'mm'];
    cst{i+j,3} = 'OAR';
    cst{i+j,4}{1} = ringPTV{j};
    cst{i+j,5} = cst{ExternalIx,5}; % as external
    cst{i+j,5}.Visible = 1;
    cst{i+j,6} = [];                % no objectives
end

% add outer ring(s) around PTV
if length(margins) > 1
    ringPTVOuter = cell(1,length(margins)-1);
    
    for k = 1:length(margins)-1
        % outer ring contour minus inner ring contour
        ringPTVOuter{k} = setdiff(ringPTV{k+1},ringPTV{k});
        
        % assign outer ring to cst
        cst{i+j+k,1} = i+j+k-1;
        cst{i+j+k,2} = ['outerRing PTV ' num2str(margins{k+1}) '-' num2str(margins{k}) 'mm'];
        cst{i+j+k,3} = 'OAR';
        cst{i+j+k,4}{1} = ringPTVOuter{k};
        cst{i+j+k,5} = cst{ExternalIx,5}; % as external
        cst{i+j+k,5}.Visible = 1;
        cst{i+j+k,6} = [];                % no objectives
    end
end


%% calculate DVH and QI
dvh_homo = matRad_calcDVH(cst,resultGUI.matRadRecalc_RBExDose,'cum');
qi_homo  = matRad_calcQualityIndicators(cst,pln,resultGUI.matRadRecalc_RBExDose);

dvh_hetero = matRad_calcDVH(cst,resultGUI.matRadHetero_RBExDose,'cum');
qi_hetero = matRad_calcQualityIndicators(cst,pln,resultGUI.matRadHetero_RBExDose);

% set several structures to invisible to clearer DVH
switch patientID
    case {'H03368_1','H03368_2'}
        invis = [1 2 3 4 5 6 7 12 13 14 15 16];
    case 'H04889'
        invis = [1 2 3 4 5 6 7 8 10 11 14 15 16 17 18 19 20 21 22 23 24 25 26 27];
    case 'S00001'
        invis = [1 2 3 4 5 6 7 8 9 10 11 15 16 17 18 19];
    case 'S00002'
        invis = [1 2 3 4 5 6 7 8 9 10 11 12 14 15];
    case {'S00003_2','S00003_3'}
        invis = [1 2 3 4 5 6 7 8 9 14 15 16 17 18];
    case 'S00004'
        invis = [1 2 3 4 5 6 7 8 9 14 15 16 17 18];
    case 'S000005'
        invis = [1 2 3 4 5 6 7 8 9 15 16 17 18 19 20];
    case 'S00006'
        invis = [1 2 3 4 5 6 7 8 9 10 11 15 16 17 18 19 20];
end

for i = invis
    cst{i,5}.Visible = 0;
end

% plot DVH comparison
dvhTitle = 'DVH comparison - phantom Pmod - solid: homogeneous lung, dotted: heterogeneous lung';
dvhFig = figure('Name','DVH comparison','Color',[0.5 0.5 0.5],'Position',([300 300 800 600]));
hold on
matRad_showDVH(dvh_homo,cst,pln,1,dvhTitle)
legend('AutoUpdate','off')
matRad_showDVH(dvh_hetero,cst,pln,2)
hold off

% show QI comparison
qiTitle = 'Copmarison quality indicators - phantom Pmod - top: homogeneous lung, bottom: heterogeneous lung';
qiFig = figure('Name','QI comparison','Color',[0.5 0.5 0.5],'Position',([300 300 800 600]));
hold on
subplot(211)
matRad_showQualityIndicators(qi_homo)
title(qiTitle)
subplot(212)
matRad_showQualityIndicators(qi_hetero)
hold off


%% save results and figures
if boolSave
    save(['D:\analyzed matRad data\HIT-Lung\' patientID '\voxelwiseConv\results_' ...
        num2str(size(stf,2)) 'fields_P256_rings'],...
        'cst','ct','patientID','pln','resultGUI','stf')
    
    savefig(dvhFig, ['D:\analyzed matRad data\HIT-Lung\' patientID '\voxelwiseConv\dvh_rings.fig'])
    savefig(qiFig, ['D:\analyzed matRad data\HIT-Lung\' patientID '\voxelwiseConv\qi_rings.fig'])
end
