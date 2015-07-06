clc,
clear 
close all
%% determine which foci shall be used
protonBaseDataHIT=matRad_getDDDfromTxt('p','\\psf\Home\Documents\Heidelberg\TRiP98DATA');
%carbonBaseDataHIT=matRad_getDDDfromTxt('C','\\psf\Home\Documents\Heidelberg\TRiP98DATA');

%% get dEdx* alpha from CNAO files or get them from the 37 spc files
pathCNAO = '\\psf\Home\Documents\Heidelberg\CNAO_baseData';
pathHIT = '\\psf\Home\Documents\Heidelberg\matRad\baseDataHIT2'; 

[ sDataHIT ] = matRad_ParseBioDataHIT(pathHIT,0);
%[ sDataCNAO ] = matRad_ParseBioDataCNAO(pathCNAO,'C',0);

%% merge ddd and dose averaged depth alpha curves

 load('\\psf\Home\Documents\Heidelberg\matRad\baseDataHIT2\sDataHIT.mat');
 load('\\psf\Home\Documents\Heidelberg\matRad\baseDataHIT2\carbonBaseDataHIT.mat');

[ baseData ] = matRad_interpDoseAvgBioData(carbonBaseDataHIT, sDataHIT ,0);

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



