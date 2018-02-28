% Comparison of differences in falloff values of single rays in all
% patients

% information about z8020 values:
% - negative z8020 values were deleted
% - first column: original imported dose distribution
%   second column: matRad recalculated dose distribution
%   third column: matRad recalculated with heterogeneity correction
% - rows: rays

% patients: 1 - H03368 1 field; 2 - H03368 2 fields; 3 - H04889; 4 - S00001; 
% 5 - S00002; 6 - S00003 2 fields; 7 - S00003 3 fields;


%% load all falloffs from patient data
clear
patientID{1} = 'H03368';
patientID{2} = 'H03368';
patientID{3} = 'H04889';
patientID{4} = 'S00001';
patientID{5} = 'S00002';
patientID{6} = 'S00003';
patientID{7} = 'S00003';


z8020{1,1} = load('C:\Matlab\HIT-Lung_falloff\H03368\1field\z8020.mat');
z8020{2,1} = load('C:\Matlab\HIT-Lung_falloff\H03368\2fields\beam1\z8020.mat');
z8020{2,2} = load('C:\Matlab\HIT-Lung_falloff\H03368\2fields\beam2\z8020.mat');

z8020{3,1} = load('C:\Matlab\HIT-Lung_falloff\H04889\beam1\z8020.mat');
z8020{3,2} = load('C:\Matlab\HIT-Lung_falloff\H04889\beam2\z8020.mat');

z8020{4,1} = load('C:\Matlab\HIT-Lung_falloff\S00001\beam1\z8020.mat');
z8020{4,2} = load('C:\Matlab\HIT-Lung_falloff\S00001\beam2\z8020.mat');
z8020{4,3} = load('C:\Matlab\HIT-Lung_falloff\S00001\beam3\z8020.mat');

z8020{5,1} = load('C:\Matlab\HIT-Lung_falloff\S00002\beam1\z8020.mat');
z8020{5,2} = load('C:\Matlab\HIT-Lung_falloff\S00002\beam2\z8020.mat');
z8020{5,3} = load('C:\Matlab\HIT-Lung_falloff\S00002\beam3\z8020.mat');

z8020{6,1} = load('C:\Matlab\HIT-Lung_falloff\S00003\2fields\beam1\z8020.mat');
z8020{6,2} = load('C:\Matlab\HIT-Lung_falloff\S00003\2fields\beam2\z8020.mat');

z8020{7,1} = load('C:\Matlab\HIT-Lung_falloff\S00003\3fields\beam1\z8020.mat');
z8020{7,2} = load('C:\Matlab\HIT-Lung_falloff\S00003\3fields\beam2\z8020.mat');
z8020{7,3} = load('C:\Matlab\HIT-Lung_falloff\S00003\3fields\beam3\z8020.mat');


%% calculate differences in z8020 for all rays in each patient
deltaZ8020 = cell(size(z8020));
for i = 1:size(z8020,1)         % loop over all patients / plans
    for j = 1:size(z8020,2)     % loop over all beams
        if ~isempty(z8020{i,j})
            tmp = cat(3, z8020{i,j}.z8020(:,3), -z8020{i,j}.z8020(:,2));        % hetero - matRad
            deltaZ8020{i,j} = sum(tmp,3,'omitnan');
            
            % remove 0-values as they come from NaNs in z8020
            deltaZ8020{i,j}(deltaZ8020{i,j} == 0) = NaN;
        end
    end
end


%% plot differences in z8020 vs. z8020
% deltaZ8020Fig = figure;
for i = 1:size(z8020,1)
    deltaZ8020Fig(i) = figure;
    for j = 1:size(z8020,2)
        if ~isempty(deltaZ8020{i,j})
            scatter(z8020{i,j}.z8020(:,2),deltaZ8020{i,j},20,'x')
        end
    end

title(['Falloff changes, all rays, ' patientID{i}])
xlabel('z8020 in water for homogeneous lung [mm]')
ylabel('Delta z8020 in water [mm]')
end


%% save figure results and figures
save('C:\Matlab\HIT-Lung_falloff\differences\results',...
    'patientID','z8020','deltaZ8020')

savefig(deltaZ8020Fig, 'C:\Matlab\HIT-Lung_falloff\differences\deltaZ8020.fig')

