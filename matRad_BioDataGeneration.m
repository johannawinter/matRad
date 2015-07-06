clc,
clear 
close all

addpath([pwd filesep 'BioDataGeneration']);
matRad_XlsSpectra2Mat('E:\TRiP98DATA_HIT-20131120\SPC\12C\RF3MM',[pwd filesep 'baseDataHIT'])


%% determine which foci shall be used
clc,
clear 
close all
protonBaseDataHIT=matRad_getDDDfromTxt('p','E:\TRiP98DATA_HIT-20131120',1);
carbonBaseDataHIT=matRad_getDDDfromTxt('C','E:\TRiP98DATA_HIT-20131120',1);

%% get dEdx* alpha from CNAO files or get them from the 37 spc files
pathCNAO = '\\psf\Home\Documents\Heidelberg\CNAO_baseData';
pathHIT = 'C:\Users\wieserh\Documents\matRad\baseDataHIT'; 

[ sDataHIT ] = matRad_ParseBioDataHIT(pathHIT,0);
%[ sDataCNAO ] = matRad_ParseBioDataCNAO(pathCNAO,'C',0);

%% merge ddd and dose averaged depth alpha curves

 load('C:\Users\wieserh\Documents\matRad\sDataHIT.mat');
 load('C:\Users\wieserh\Documents\matRad\carbonBaseDataHIT.mat');

[ baseData ] = matRad_interpDoseAvgBioData(baseData, sDataHIT ,0);

%% interpolate to deeper depths
clc
clear
close all

% load('\\psf\Home\Documents\Heidelberg\matRad\baseDataHIT2\sDataHIT.mat');
% figure,hold on
% CellLine = 9;
% for i = 1:10
%    plot(sDataHIT{1,CellLine}(1,i).depths,sDataHIT{1,CellLine}(1,i).alpha);
% end
% figure,hold on
% for i = 1:10
%    plot(sDataHITraw{1,CellLine}(1,i).depths,sDataHITraw{1,CellLine}(1,i).alpha);
% end
% 
% 

load('\\psf\Home\Documents\Heidelberg\matRad\carbonBaseDataHITBio.mat');
[ baseData ] = extrapDeeper( baseData,0 );



