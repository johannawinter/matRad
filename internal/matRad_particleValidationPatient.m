% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015, Mark Bangert, on behalf of the matRad development team
%
% m.bangert@dkfz.de
%
% This file is not part of matRad.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
clc
rootPath = ['E:\matRad validation\protons\patientData'];
% Syngo dose cube was imported using matRad's dicom import
load([rootPath filesep 'H02333_111mm_CstNew.mat']);
%matRadGUI

ctTmp = ct.cube{1};
ctTmp(isnan(ctTmp)) = 0;
ct.cube{1} = ctTmp;

SyngoDoseCube     = resultGUI.physicalDose;
%matRadGUI

% load Fluka Monte Carlo Cube
MCfilename   = [rootPath filesep 'H02333_01T180_dosePhys' '.txt'];
MCcube       = matRad_readMCdataPatient(MCfilename);
pln.runSequencing   = false; 
pln.runDAO          = false; 


counter = 0;
for i = 1:pln.numOfBeams
    for j = 1:stf(i).numOfRays
        for k = 1:stf(i).numOfBixelsPerRay(j)
            counter = counter + 1;
            w(counter) = stf(i).ray(j).weight(k);
        end
    end
end
   

%% dose calculation
%matRadDoseCubeSupineDirect = matRad_calcDoseDirect(ct,stf,pln,cst,w);  
%save('matRadDoseCubeSupineDirect','matRadDoseCubeSupineDirect');
load([rootPath filesep 'matRadDoseCubeSupineDirect']);
matRadDoseCubeSupine = matRadDoseCubeSupineDirect.physicalDose;
matRadDoseCube       = matRadDoseCubeSupine;

%% flip ct, dose cube and ct to prone


%% Plot all cubes
cellName = {'Syngo','FLUKA','matRad'};
slice       = round(pln.isoCenter(3)/ct.resolution.z);
plane       = 3;
DoseCutOff  = 0;
defFontSize = 20;

% % integral dose
relIntDoseDif = (1-sum(matRadDoseCube(:))/sum(SyngoDoseCube(:)))*100;
fprintf(['Relative difference in integral dose: ' num2str(relIntDoseDif) '%%\n']);

vDim = size(SyngoDoseCube);
defaultLineWidth = 1.5;

vLevels = [0.2 0.4 0.5 0.6 0.70 0.80 0.90 0.95 1 1.1];
MaxDoseSyngo   = max(max(SyngoDoseCube(:)));
MaxDoseMC      = max(max(MCcube.cube(:)));
MaxDoseMatRad  = max(max(matRadDoseCube(:)));

if plane == 1
    mSyngo  = squeeze(SyngoDoseCube(slice,:,:));
    mMC     = squeeze(MCcube.cube(slice,:,:));
    mMatRad = squeeze(matRadDoseCube(slice,:,:));
    ctImg   = squeeze(ct.cube{1}(slice,:,:));
elseif plane == 2
    mSyngo  = squeeze(SyngoDoseCube(:,slice,:));
    mMC     = squeeze(MCcube.cube(:,slice,:));
    mMatRad = squeeze(matRadDoseCube(:,slice,:));
    ctImg   = squeeze(ct.cube{1}(:,slice,:));
elseif plane == 3
    mSyngo  = squeeze(SyngoDoseCube(:,:,slice));
    %mMC     = squeeze(MCcube.cube(:,:,slice));
    mMatRad = squeeze(matRadDoseCube(:,:,slice));
    ctImg   = squeeze(ct.cube{1}(:,:,slice));
end
vLevelsDose = vLevels * MaxDoseSyngo;
%% compare matRad and Syngo dose cube
m=100; s = 4;
cm_magma=magma(m);cm_inferno=inferno(m);cm_plasma=plasma(m);cm_viridis=viridis(m);
CM = cm_viridis;%jet;
%downsample cubes

maxDose = max([MaxDoseSyngo MaxDoseMatRad]);
figure,set(gcf,'Color',[1 1 1]);
CropLim = [120 300 45 265];
ct_rgb = ind2rgb(uint8(100*ctImg/max(ct.cube{1}(:))),bone(100));

ct_rgbCrop  = ct_rgb;%(CropLim(1):CropLim(2),CropLim(3):CropLim(4));
mMatRadCrop = mMatRad;%(CropLim(1):CropLim(2),CropLim(3):CropLim(4));
mSyngoCrop  = mSyngo;%(CropLim(1):CropLim(2),CropLim(3):CropLim(4));
% transversal slices
subplot(131),imagesc(ct_rgbCrop),hold on;colormap(CM),
h1=imagesc(mMatRadCrop);
set(h1,'AlphaData', .6*double(mMatRadCrop>DoseCutOff));
contour(mMatRadCrop,vLevelsDose,'LevelListMode','manual','LineWidth',1.5);
title('matRad','Interpreter','Latex','FontSize',defFontSize); cBarHandel = colorbar(gca);
set(get(cBarHandel,'ylabel'),'String', 'dose [Gy]','FontSize',defFontSize,'Interpreter','Latex');cBarHandel.FontSize = 15;
set(gca,'XTickLabel','');set(gca,'YTickLabel',''), axis equal
mask = zeros(size(ct.cube{1})); 
mask(unique(cst{s,4}{1})) = 1;
contour(squeeze(mask(:,:,slice)),.5*[1 1],'Color',[0 0 0],'LineWidth',1.5,'DisplayName',cst{s,2});
set(gca,'ylim',[CropLim(1) CropLim(2)]),set(gca,'xlim',[CropLim(3) CropLim(4)])

subplot(132),imagesc(ct_rgbCrop),hold on;
h1=imagesc(mMatRadCrop-mSyngoCrop);polarmap(100,1);cBarHandel = colorbar(gca) ;set(gca,'XTickLabel','');set(gca,'YTickLabel',''), axis square
set(get(cBarHandel,'ylabel'),'String', 'dose [Gy]','FontSize',defFontSize,'Interpreter','Latex');cBarHandel.FontSize = 15;
title('matRad - Syngo','Interpreter','Latex','FontSize',defFontSize); , axis equal
mask = zeros(size(ct.cube{1})); 
mask(unique(cst{s,4}{1})) = 1;
contour(squeeze(mask(:,:,slice)),.5*[1 1],'Color',[0 0 0],'LineWidth',1.5,'DisplayName',cst{s,2});
set(gca,'ylim',[CropLim(1) CropLim(2)]),set(gca,'xlim',[CropLim(3) CropLim(4)]),

[gammaCube,costumMap,gammaPassRate] = matRad_gammaIndex(matRadDoseCube,SyngoDoseCube,[1 1 1]);
subplot(133),imagesc(ct_rgbCrop),hold on;
h1=imagesc(gammaCube(:,:,slice),[0 2]);
set(h1,'AlphaData', .6*double(mMatRadCrop>DoseCutOff));
set(gca,'XTickLabel','');set(gca,'YTickLabel','');
colormap(gca,costumMap);cBarHandel = colorbar(gca); , axis equal
set(get(cBarHandel,'ylabel'),'String','passed voxels $<$ 1','fontsize',defFontSize,'Interpreter','Latex');
cBarHandel.FontSize = 15;
title(['$\gamma$ index, pass rate: ' num2str(gammaPassRate,5) '$\%$'],'FontSize',defFontSize,'Interpreter','Latex');
mask = zeros(size(ct.cube{1})); 
mask(unique(cst{s,4}{1})) = 1;
contour(squeeze(mask(:,:,slice)),.5*[1 1],'Color',[0 0 0],'LineWidth',1.5,'DisplayName',cst{s,2});
set(gca,'ylim',[CropLim(1) CropLim(2)]),set(gca,'xlim',[CropLim(3) CropLim(4)]),

export_fig('testt','-eps', '-transparent')
export_fig('testt','-dpdf',' -transparent' )

%%
profileSlice = round(pln.isoCenter(2));49; 
cube1CentralRayProf = (mMatRad(profileSlice,:));
cube2CentralRayProf = (mSyngo(profileSlice,:));

subplot(224),hold on
plot(1:1:size(mMatRad,2),cube1CentralRayProf,'r','LineWidth',defaultLineWidth)
plot(1:1:size(mMatRad,2),cube2CentralRayProf,'b--','LineWidth',defaultLineWidth)
box on,grid on,title(['depth profiles at slice: ' num2str(profileSlice)]),legend({'matRad','Syngo'})

figure,
for k = profileSlice-15:profileSlice+15 
cube1CentralRayProf = (mMatRad(k,:));
cube2CentralRayProf = (mSyngo(k,:));hold on
plot(1:1:size(mMatRad,2),cube1CentralRayProf,'r','LineWidth',defaultLineWidth),hold on
plot(1:1:size(mMatRad,2),cube2CentralRayProf,'b--','LineWidth',defaultLineWidth)
box on,grid on,title(['depth profiles at slice: ' num2str(k)]),legend({'matRad','Syngo'})

waitforbuttonpress;
cla
end
%% plot the result in coronal, sagittal and axial slice
figure,set(gcf,'Color',[1 1 1]);
sliceCoronal = 90;
subplot(131),
ct_rgb = ind2rgb(uint8(63*squeeze(ct.cube{1}(sliceCoronal,:,:))/max(ct.cube{1}(:))),bone);
imagesc(ct_rgb),hold on;
colormap(jet),
h1=imagesc(squeeze(matRadDoseCube(sliceCoronal,:,:)));
set(h1,'AlphaData', .4*double(squeeze(matRadDoseCube(sliceCoronal,:,:))>DoseCutOff));
[C,myContour] = contour(squeeze(matRadDoseCube(sliceCoronal,:,:)),vLevelsDose,'LevelListMode','manual','LineWidth',1.5);
title(['coronal ' cellName{1,3} ': slice ' num2str(sliceCoronal)]);
cBarHandel = colorbar(gca);
set(get(cBarHandel,'ylabel'),'String', 'dose [Gy]','fontsize',15,'Interpreter','Latex');
set(cBarHandel,'YLim',[0 max(max(matRadDoseCube(sliceCoronal,:,:)))]);
% plot structure contours
colors = colorcube; 
mask = zeros(size(ct.cube{1})); 
hold on,
for s = 1:size(cst,1)
    if ~strcmp(cst{s,3},'IGNORED') && ~isempty(cst{s,4}) && sum(strcmp(cst{s,2},{'H_PTV'})>0)
        mask(:) = 0;
        mask(unique(cst{s,4})) = 1;
        if sum(sum(mask(sliceCoronal,:,:))) > 0 && var(var(mask(sliceCoronal,:,:))) > 0
            contour(squeeze(mask(sliceCoronal,:,:)),.5*[1 1],'Color',colors(s,:),'LineWidth',1.5,'DisplayName',cst{s,2});
        end
    end
end
axis square


subplot(132),
sliceSagittal = 90;
ct_rgb = ind2rgb(uint8(63*squeeze(ct.cube{1}(:,sliceSagittal,:))/max(ct.cube{1}(:))),bone);
imagesc(ct_rgb),hold on;
colormap(jet),
h1=imagesc(squeeze(matRadDoseCube(:,sliceSagittal,:)));
set(h1,'AlphaData', .4*double(squeeze(matRadDoseCube(:,sliceSagittal,:))>DoseCutOff));
[C,myContour] = contour(squeeze(matRadDoseCube(:,sliceSagittal,:)),vLevelsDose,'LevelListMode','manual','LineWidth',1.5);
title(['sagittal ' cellName{1,3} ': slice ' num2str(sliceSagittal)]);
cBarHandel = colorbar(gca);
set(get(cBarHandel,'ylabel'),'String', 'dose [Gy]','fontsize',15,'Interpreter','Latex');
axis square
% plot structure contours
colors = colorcube; 
mask = zeros(size(ct.cube{1})); 
hold on,
for s = 1:size(cst,1)
    if ~strcmp(cst{s,3},'IGNORED') && ~isempty(cst{s,4}) && sum(strcmp(cst{s,2},{'H_PTV'})>0)
        mask(:) = 0;
        mask(unique(cst{s,4})) = 1;
        if sum(sum(mask(:,sliceSagittal,:))) > 0 && var(var(mask(:,sliceSagittal,:))) > 0
            contour(squeeze(mask(:,sliceSagittal,:)),.5*[1 1],'Color',colors(s,:),'LineWidth',1.5,'DisplayName',cst{s,2});
        end
    end
end


subplot(133),
sliceAxial = 29;
ct_rgb = ind2rgb(uint8(63*squeeze(ct.cube{1}(:,:,sliceAxial))/max(ct.cube{1}(:))),bone);
imagesc(ct_rgb),hold on;
colormap(jet),
h1=imagesc(squeeze(matRadDoseCube(:,:,sliceAxial)));
set(h1,'AlphaData', .4*double(matRadDoseCube(:,:,sliceAxial)>DoseCutOff));
[C,myContour] = contour(squeeze(matRadDoseCube(:,:,sliceAxial)),vLevelsDose,'LevelListMode','manual','LineWidth',1.5);
title(['axial ' cellName{1,3} ': slice ' num2str(sliceAxial)]);
cBarHandel = colorbar(gca);
set(get(cBarHandel,'ylabel'),'String', 'dose [Gy]','fontsize',15,'Interpreter','Latex');
axis square
% plot structure contours
colors = colorcube; 
mask = zeros(size(ct.cube{1})); 
hold on,
for s = 1:size(cst,1)
    if ~strcmp(cst{s,3},'IGNORED') && ~isempty(cst{s,4}) && sum(strcmp(cst{s,2},{'H_PTV'})>0)
        mask(:) = 0;
        mask(unique(cst{s,4})) = 1;
        if sum(sum(mask(:,:,sliceAxial))) > 1 && var(var(mask(:,:,sliceAxial))) > 0
            contour(squeeze(mask(:,:,sliceAxial)),.5*[1 1],'Color',colors(s,:),'LineWidth',1.5,'DisplayName',cst{s,2});
        end
    end
end


%% second plot 3x3
figure
set(gcf,'Color',[1 1 1]);
ct_rgb = ind2rgb(uint8(63*ctImg/max(ct.cube{1}(:))),bone);
% transversal slices
subplot(3,3,1)
imagesc(ct_rgb),hold on;
colormap(jet),
h1=imagesc(mSyngo);
set(h1,'AlphaData', .6*double(mSyngo>DoseCutOff));
%[C,myContour] = contour(SyngoDoseCube(:,:,slice),vLevelsDose,'LevelListMode','manual','LineWidth',1.5);
title([ cellName{1,1} ': slice ' num2str(slice)]);
cBarHandel = colorbar(gca);
set(get(cBarHandel,'ylabel'),'String', 'dose [Gy]','fontsize',15,'Interpreter','Latex');
set(cBarHandel,'YLim',[0 MaxDoseSyngo]);

% plot structure contours
colors = colorcube; 
InnerCnt = 1;
mask = zeros(size(ct.cube{1})); 
hold on,
for s = 1:size(cst,1)
    if ~strcmp(cst{s,3},'IGNORED') && ~isempty(cst{s,4}) && sum(strcmp(cst{s,2},{'PTV'})>0)
        mask(:) = 0;
        mask(unique(cst{s,4})) = 1;
        if plane == 1 && sum(sum(mask(slice,:,:))) > 0
            contour(squeeze(mask(slice,:,:)),.5*[1 1],'Color',colors(s,:),'LineWidth',1.5,'DisplayName',cst{s,2});
            strLegend{InnerCnt,1} = cst{s,2};
            InnerCnt = InnerCnt +1;
        elseif plane == 2 && sum(sum(mask(:,slice,:))) > 0
            contour(squeeze(mask(:,slice,:)),.5*[1 1],'Color',colors(s,:),'LineWidth',1.5,'DisplayName',cst{s,2});
            strLegend{InnerCnt,1} = cst{s,2};
            InnerCnt = InnerCnt +1;
        elseif plane == 3 && sum(sum(mask(:,:,slice))) > 1 && var(var(mask(:,:,slice))) > 0
            contour(squeeze(mask(:,:,slice)),.5*[1 1],'Color',colors(s,:),'LineWidth',1.5,'DisplayName',cst{s,2});
            strLegend{InnerCnt,1} = cst{s,2};
            InnerCnt = InnerCnt +1;
        end
    end
end
axis equal
myLegend = legend(strLegend,'location','NorthEast');
set(myLegend,'FontSize',7);
set(myLegend,'color','none');
set(myLegend,'TextColor', [1 1 1]);
legend boxoff


subplot(3,3,2)
imagesc(ct_rgb),hold on;
colormap(jet),
h1=imagesc(mMC);
set(h1,'AlphaData', .6*double(mMC>DoseCutOff));
%[C,myContour] = contour(SyngoDoseCube(:,:,slice),vLevelsDose,'LevelListMode','manual','LineWidth',1.5);
title([ cellName{1,2} ': slice ' num2str(slice)]);
cBarHandel = colorbar(gca);
set(get(cBarHandel,'ylabel'),'String', 'dose [Gy]','fontsize',15,'Interpreter','Latex');
set(cBarHandel,'YLim',[0 MaxDoseSyngo]);

% plot structure contours
colors = colorcube; 
InnerCnt = 1;
mask = zeros(size(ct.cube{1})); 
hold on,
for s = 1:size(cst,1)
    if ~strcmp(cst{s,3},'IGNORED') && ~isempty(cst{s,4}) && sum(strcmp(cst{s,2},{'PTV'})>0)
        mask(:) = 0;
        mask(unique(cst{s,4})) = 1;
        if plane == 1 && sum(sum(mask(slice,:,:))) > 0
            contour(squeeze(mask(slice,:,:)),.5*[1 1],'Color',colors(s,:),'LineWidth',1.5,'DisplayName',cst{s,2});
            strLegend{InnerCnt,1} = cst{s,2};
            InnerCnt = InnerCnt +1;
        elseif plane == 2 && sum(sum(mask(:,slice,:))) > 0
            contour(squeeze(mask(:,slice,:)),.5*[1 1],'Color',colors(s,:),'LineWidth',1.5,'DisplayName',cst{s,2});
            strLegend{InnerCnt,1} = cst{s,2};
            InnerCnt = InnerCnt +1;
        elseif plane == 3 && sum(sum(mask(:,:,slice))) > 1 && var(var(mask(:,:,slice))) > 0
            contour(squeeze(mask(:,:,slice)),.5*[1 1],'Color',colors(s,:),'LineWidth',1.5,'DisplayName',cst{s,2});
            strLegend{InnerCnt,1} = cst{s,2};
            InnerCnt = InnerCnt +1;
        end
    end
end
axis equal
myLegend = legend(strLegend,'location','NorthEast');
set(myLegend,'FontSize',7);
set(myLegend,'color','none');
set(myLegend,'TextColor', [1 1 1]);
legend boxoff


subplot(3,3,3)
imagesc(ct_rgb),hold on;
colormap(jet),
h1=imagesc(mMatRad);
set(h1,'AlphaData', .6*double(mMatRad>DoseCutOff));
%[C,myContour] = contour(SyngoDoseCube(:,:,slice),vLevelsDose,'LevelListMode','manual','LineWidth',1.5);
title([ cellName{1,3} ': slice ' num2str(slice)]);
cBarHandel = colorbar(gca);
set(get(cBarHandel,'ylabel'),'String', 'dose [Gy]','fontsize',15,'Interpreter','Latex');
set(cBarHandel,'YLim',[0 MaxDoseSyngo]);

% plot structure contours
colors = colorcube; 
InnerCnt = 1;
mask = zeros(size(ct.cube{1})); 
hold on,
for s = 1:size(cst,1)
    if ~strcmp(cst{s,3},'IGNORED') && ~isempty(cst{s,4}) && sum(strcmp(cst{s,2},{'PTV'})>0)
        mask(:) = 0;
        mask(unique(cst{s,4})) = 1;
        if plane == 1 && sum(sum(mask(slice,:,:))) > 0
            contour(squeeze(mask(slice,:,:)),.5*[1 1],'Color',colors(s,:),'LineWidth',1.5,'DisplayName',cst{s,2});
            strLegend{InnerCnt,1} = cst{s,2};
            InnerCnt = InnerCnt +1;
        elseif plane == 2 && sum(sum(mask(:,slice,:))) > 0
            contour(squeeze(mask(:,slice,:)),.5*[1 1],'Color',colors(s,:),'LineWidth',1.5,'DisplayName',cst{s,2});
            strLegend{InnerCnt,1} = cst{s,2};
            InnerCnt = InnerCnt +1;
        elseif plane == 3 && sum(sum(mask(:,:,slice))) > 1 && var(var(mask(:,:,slice))) > 0
            contour(squeeze(mask(:,:,slice)),.5*[1 1],'Color',colors(s,:),'LineWidth',1.5,'DisplayName',cst{s,2});
            strLegend{InnerCnt,1} = cst{s,2};
            InnerCnt = InnerCnt +1;
        end
    end
end
axis equal
myLegend = legend(strLegend,'location','NorthEast');
set(myLegend,'FontSize',7);
set(myLegend,'color','none');
set(myLegend,'TextColor', [1 1 1]);
legend boxoff



normMC      = max(max(mMC));
normMatRad = max(max(mMatRad));

ax3 = subplot(3,3,4);
imagesc(100*(mSyngo-mMC)./normMC);
myMap = matRad_getCostumColorbarDiff(SyngoDoseCube,MCcube.cube,slice,plane);
colormap(ax3,myMap); colorbar;
title(['rel diff [%] (' cellName{1,1} ' - ' cellName{1,2} ')/ ' cellName{1,2} ' : slice ' num2str(slice)])

normRad = max(max(mMatRad));
ax3 = subplot(3,3,5);
imagesc(100*(mSyngo-mMatRad)./normRad)
myMap = matRad_getCostumColorbarDiff(SyngoDoseCube,matRadDoseCube,slice,plane);
colormap(ax3,myMap); colorbar;
title(['rel diff [%] (' cellName{1,1} ' - ' cellName{1,3} ')/ ' cellName{1,3} '  : slice ' num2str(slice)])


ax3 = subplot(3,3,6);
imagesc(100*(mMatRad-mMC)./normMC)
myMap = matRad_getCostumColorbarDiff(matRadDoseCube,MCcube.cube,slice,plane);
colormap(ax3,myMap); colorbar;
title(['rel diff [%] (' cellName{1,3} ' - ' cellName{1,2} ')/ ' cellName{1,2} ' : slice ' num2str(slice)])

%% profile along depth
profileSlice = 90;
cube1CentralRayProf = squeeze(SyngoDoseCube(profileSlice,:,slice));
cube2CentralRayProf = squeeze(MCcube.cube(profileSlice,:,slice));
cube3CentralRayProf = squeeze(matRadDoseCube(profileSlice,:,slice));

subplot(3,3,7)
hold on
plot(1:1:vDim(2),cube1CentralRayProf,'r','LineWidth',defaultLineWidth)
plot(1:1:vDim(2),cube2CentralRayProf,'b--','LineWidth',defaultLineWidth)
plot(1:1:vDim(2),cube3CentralRayProf,'k:','LineWidth',defaultLineWidth)
box on
grid on
title(['depth profiles at slice: ' num2str(profileSlice)])
legend(cellName)


%% idds

x = 1:1:vDim(2);
cube1IDD = sum(sum(SyngoDoseCube,3),1);
cube2IDD = sum(sum(MCcube.cube,3),1);
cube3IDD = sum(sum(matRadDoseCube,3),1);


subplot(3,3,8)
hold on
plot(x,cube1IDD,'r','LineWidth',defaultLineWidth)
plot(x,cube2IDD,'b--','LineWidth',defaultLineWidth)
plot(x,cube3IDD,'K:','LineWidth',defaultLineWidth)
title(['IDDs'])
box on
grid on
legend(cellName)

%% profile along lateral direction entrance
LatprofileSlice = 80;
cube1LatProfileEnt = squeeze(SyngoDoseCube(:,LatprofileSlice,slice));
cube2LatProfileEnt = squeeze(MCcube.cube(:,LatprofileSlice,slice));
cube3LatProfileEnt = squeeze(matRadDoseCube(:,LatprofileSlice,slice));

subplot(3,3,9)
hold on
plot(1:1:vDim(1),cube1LatProfileEnt,'r','LineWidth',defaultLineWidth)
plot(1:1:vDim(1),cube2LatProfileEnt,'b--','LineWidth',defaultLineWidth)
plot(1:1:vDim(1),cube3LatProfileEnt,'k:','LineWidth',defaultLineWidth)
box on
grid on
title(['lateral profile at: ' num2str(LatprofileSlice)])
legend(cellName)  

%%
[gammaCube,costumMap] = matRad_gammaIndex(MCcube.cubeOrg,matRadDoseCubeOrg,[0.9990 0.9990 3],32);
figure,subplot(131),
set(gcf,'Color',[1 1 1]);
ct_rgb = ind2rgb(uint8(63*ctImg/max(ct.cube{1}(:))),bone);
imagesc(ct_rgb),hold on;
colormap(jet),
h1=imagesc(matRadDoseCubeOrg(:,:,32));
set(h1,'AlphaData', .6*double(matRadDoseCubeOrg(:,:,32)>0));
%[C,myContour] = contour(SyngoDoseCube(:,:,slice),vLevelsDose,'LevelListMode','manual','LineWidth',1.5);
title(['a) matRad'],'FontSize',15);
set(gca,'XTickLabel','');
set(gca,'YTickLabel','');
cBarHandel = colorbar(gca);
set(get(cBarHandel,'ylabel'),'String', 'dose [Gy]','fontsize',16);
cBarHandel.FontSize = 14;
sMax = max(max(matRadDoseCubeOrg(:,:,32)));
set(cBarHandel,'YLim',[0 sMax]);
mask = zeros(size(ct.cube{1})); 
hold on,
s = 10;
mask(unique(cst{s,4})) = 1;
contour(squeeze(mask(:,:,slice)),.5*[1 1],'Color',[0 0 0],'LineWidth',1.5,'DisplayName',cst{s,2});
myLegend = legend({'PTV'},'location','NorthEast');
set(myLegend,'FontSize',15);
set(myLegend,'color','none');
set(myLegend,'TextColor', [1 1 1],'FontWeight','bold');
legend boxoff
set(gca,'FontSize',defFontSize);


vNorm = max(max(MCcube.cubeOrg(:)));
ax3 = subplot(132);
imagesc(ct_rgb),hold on;
colormap(jet),
mI = 100*((matRadDoseCubeOrg(:,:,32)./vNorm)-(MCcube.cubeOrg(:,:,32)./vNorm));
h1=imagesc(mI);
myMap = matRad_getCostumColorbarDiff(matRadDoseCubeOrg,MCcube.cubeOrg,32,3);
colormap(ax3,myMap); 
cBarHandel = colorbar;
cBarHandel.FontSize = 15;
set(get(cBarHandel,'ylabel'),'String', 'rel. diff [%]','fontsize',16);
set(h1,'AlphaData', .6*double(matRadDoseCubeOrg(:,:,32)>0));
mask = zeros(size(ct.cube{1})); 
hold on,
s = 10;
mask(unique(cst{s,4})) = 1;
contour(squeeze(mask(:,:,slice)),.5*[1 1],'Color',[0 0 0],'LineWidth',1.5,'DisplayName',cst{s,2});
title(['b) rel. diff of matRad - Monte Carlo'],'FontSize',15);
set(gca,'XTickLabel','');
set(gca,'YTickLabel','');
myLegend = legend({'PTV'},'location','NorthEast');
set(myLegend,'FontSize',14);
set(myLegend,'color','none');
set(myLegend,'TextColor', [1 1 1],'FontWeight','bold');
legend boxoff
set(gca,'FontSize',defFontSize);


subplot(133),
imagesc(ct_rgb),hold on;
colormap(jet),
h1=imagesc(gammaCube(:,:,32),[0 2]);
set(h1,'AlphaData', .6*double(matRadDoseCubeOrg(:,:,32)>0));
set(gca,'XTickLabel','');
set(gca,'YTickLabel','');
colormap(gca,costumMap);
cBarHandel = colorbar(gca);
absDoseThreshold = 3/100 * max(matRadDoseCubeOrg(:));
doseIx          = MCcube.cubeOrg > absDoseThreshold;
numOfPassGamma  = sum(gammaCube(doseIx) < 1);
gammaPassRate   = 100 * numOfPassGamma / sum(doseIx(:));

set(get(cBarHandel,'ylabel'),'String','passed voxels < 1','fontsize',16);
cBarHandel.FontSize = 15;
title(['c) \gamma index, pass rate: ' num2str(gammaPassRate,5) '%'],'FontSize',15);
mask = zeros(size(ct.cube{1})); 
mask(unique(cst{s,4})) = 1;
contour(squeeze(mask(:,:,slice)),.5*[1 1],'Color',[0 0 0],'LineWidth',1.5,'DisplayName',cst{s,2});

%axis equal
myLegend = legend({'PTV'},'location','NorthEast');
set(myLegend,'FontSize',14);
set(myLegend,'color','none');
set(myLegend,'TextColor', [1 1 1],'FontWeight','bold');
legend boxoff
set(gca,'FontSize',defFontSize);












