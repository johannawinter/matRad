function matRad_plotTwoDoseCubes(ct,cst,pln,doseCube1,doseCube2,cNames,varargin)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots a dose cube using matRad varibales 
% 
% call
%    matRad_plotDoseCube(ct,cst,pln,doseCube,slice,plane,filename)
%
% input
%   ct:          matRads ct struct
%   cst:         matRads cst struct
%   pln:         matRads pln struct
%   doseCube:    3D dose cube having the same size dube dimensions as the ct
%   slice:       scalar to determine the slice
%   plane:       scalar to determine the plane
%               (1=coronal,2=sagital,3=axial);
%   filename     if a filename is provided the figure will be saved as pdf
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2016, Hans-Peter Wieser
%
% This file is NOT part of the official matRad release. 
% This file has to be used only for internal purposes! 
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[folder, ~, ~] = fileparts(mfilename('fullpath'));
addpath(genpath([folder filesep 'colormaps']));
% select a awesome colormap
%cm = jet; cm = inferno(100);
cm       = viridis(100);
FlagSave = false;
FlagCst  = false; 
FlagBBox = true;

% handle variable number of inputs 
if nargin >= 7
    if isnumeric(varargin{1})
        if varargin{1}>3
            slice = varargin{1};
        else
            plane = varargin{1};
        end
    elseif ischar(varargin{1})
        FlagSave = true;
        filename = varargin{1};
    end
end

if nargin >= 8
     if isnumeric(varargin{2})
        if varargin{2}>3
            slice = varargin{2};
        else
            plane = varargin{2};
        end
    elseif ischar(varargin{2})
        FlagSave = true;
        filename = varargin{2};
     end
end
if nargin == 9
    if ischar(varargin{3})
        FlagSave = true;
        filename = varargin{3};
    end
end

if (sum(size(ct.cube{1}) ~= size(doseCube1))>0) || (sum(size(doseCube1)~=size(doseCube2))>1)
    error('inconsistent cube dimensions');
end
    

if ~exist('plane','var')
    plane = 3;
end
if ~exist('slice','var')
    slice = round(pln.isoCenter(plane)/ct.resolution.z);
end

defaultFontSize = 16;
DoseCutOff      = 0;   % relative number between 0 and 1
maskCst         = zeros(ct.cubeDim); V = [cst{:,4}]; V = unique(vertcat(V{:}));
maskCst(V)      = 1; 

if plane == 1
    maskSlice  = squeeze(maskCst(slice,:,:));
    bBox       = [find(sum(maskSlice,1)>0,1,'first') find(sum(maskSlice,2)>0,1,'first')  find(sum(maskSlice,1)>0,1,'last') find(sum(maskSlice,2)>0,1,'last')];
    ctSlice    = squeeze(ct.cube{1}(slice,:,:)); doseSlice1  = squeeze(doseCube1(slice,:,:));  doseSlice2   = squeeze(doseCube2(slice,:,:));     
elseif plane == 2
    maskSlice  = squeeze(maskCst(:,slice,:));
    bBox       = [find(sum(maskSlice,1)>0,1,'first') find(sum(maskSlice,2)>0,1,'first')  find(sum(maskSlice,1)>0,1,'last') find(sum(maskSlice,2)>0,1,'last')];
    ctSlice    = squeeze(ct.cube{1}(:,slice,:)); doseSlice1 = squeeze(doseCube1(:,slice,:)); doseSlice2    = squeeze(doseCube2(:,slice,:));
elseif plane == 3
    maskSlice  = squeeze(maskCst(:,:,slice));
    bBox       = [find(sum(maskSlice,1)>0,1,'first') find(sum(maskSlice,2)>0,1,'first')  find(sum(maskSlice,1)>0,1,'last') find(sum(maskSlice,2)>0,1,'last')];
    ctSlice    = squeeze(ct.cube{1}(:,:,slice));  doseSlice1  = squeeze(doseCube1(:,:,slice)); doseSlice2    = squeeze(doseCube2(:,:,slice));
end

    if FlagBBox
        ctSlice     = ctSlice(bBox(2):bBox(4),bBox(1):bBox(3));
        doseSlice1  = doseSlice1(bBox(2):bBox(4),bBox(1):bBox(3));
        doseSlice2  = doseSlice2(bBox(2):bBox(4),bBox(1):bBox(3));
    end

vLevels     = [ 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.75 0.8 0.9 0.95 1 1.1 ];
maxDose     =  max([max(doseSlice1(:))  max(doseSlice2(:))]);
vLevelDose  = maxDose  * vLevels;
ct_rgb      = ind2rgb(uint8(63*squeeze(ctSlice)/max(ct.cube{1}(:))),bone);

figure('units','normalized','outerposition',[0 0 1 1]),set(gcf,'Color',[1 1 1]);
subplot(121),
imagesc(ct_rgb),hold on,colormap(cm),h1 =imagesc(doseSlice1);
set(gca,'XTickLabel','');set(gca,'YTickLabel','');
set(h1,'AlphaData', .8*(double(doseSlice1)>DoseCutOff*maxDose));
title([cNames{1} ' slice: ' num2str(slice)  ' plane: ' num2str(plane)],'Interpreter','Latex','FontSize',defaultFontSize)   
[~,~] = contour(doseSlice1,vLevelDose,'LevelListMode','manual','LineWidth',1.5);
cBarHandel = colorbar(gca,'FontSize',defaultFontSize);
caxis(gca,[0, maxDose]); set(get(cBarHandel,'ylabel'),'String', 'dose [Gy]','fontsize',defaultFontSize,'Interpreter','Latex');
set(cBarHandel,'YLim',[0 maxDose]); if plane ==3;axis equal;end

colors = colorcube; 
ColorCnt = 1; CntLegend = 1; cHandle = [];
for i = 1:size(cst,1)
         maskCst(:) = 0;
         maskCst(cst{i,4}{1}) = 1;
         if plane == 1 && sum(sum(maskCst(slice,:,:))) > 0
            tmpSlice = squeeze(maskCst(slice,:,:));
         elseif plane == 2 && sum(sum(maskCst(:,slice,:))) > 0
            tmpSlice = squeeze(maskCst(:,slice,:));
        elseif plane == 3 && sum(sum(maskCst(:,:,slice))) > 0
            tmpSlice = squeeze(maskCst(:,:,slice));
         else
             tmpSlice = 0;
         end
         
         if sum(tmpSlice(:)) > 0 && FlagCst
             [c, h] = contour(gca,squeeze(maskCst(:,:,slice)),.5*[1 1],'LineWidth',2,'Color',colors(ColorCnt,:),'DisplayName',cst{i,2});
               if ~isempty(c)
                 cHandle = [cHandle h]; ColorCnt = ColorCnt + 1; Name = regexprep(cst{i,2},'_','.'); 
                 sLegend{CntLegend} = Name; CntLegend = CntLegend +1;
               end
         end
end
if FlagCst
    legend(cHandle,sLegend,'FontSize',defFontSize)
end

subplot(122),imagesc(ct_rgb),hold on,colormap(cm),h1 =imagesc(doseSlice2);
set(gca,'XTickLabel','');set(gca,'YTickLabel','');
set(h1,'AlphaData', .8*(double(doseSlice2)>DoseCutOff*maxDose));
title([cNames{2} ' slice: ' num2str(slice)  ' plane: ' num2str(plane)],'Interpreter','Latex','FontSize',defaultFontSize)   
[~,~] = contour(doseSlice2,vLevelDose,'LevelListMode','manual','LineWidth',1.5);
cBarHandel = colorbar(gca,'FontSize',defaultFontSize);
caxis(gca,[0, maxDose]); set(get(cBarHandel,'ylabel'),'String', 'dose [Gy]','fontsize',defaultFontSize,'Interpreter','Latex');
set(cBarHandel,'YLim',[0 maxDose]); if plane ==3;axis equal;end


%% save first plot
folderPath = [folder filesep 'exports'];
if exist('FolderPath','file')
    status = mkdir(folderPath);
end

if FlagSave
    FullFileName = [folderPath filesep filename '_' datestr(now, 'dd_mmm_yyyy_HH_MM')];
    set(gcf,'PaperOrientation','landscape','PaperUnits','normalized','PaperPosition', [0 0 1 1]);
    print('-dpdf','-r300',FullFileName);
end


figure('units','normalized','outerposition',[0 0 1 1]),set(gcf,'Color',[1 1 1]);
subplot(121),imagesc(ct_rgb),hold on;
h1=imagesc(doseSlice1-doseSlice2); polarmap(100,1); 
caxis([-1,1]*max(abs(caxis))); 
cBarHandel = colorbar(gca,'Location','eastoutside');  %southoutside
set(gca,'XTickLabel','');set(gca,'YTickLabel',''), 
set(h1,'AlphaData', .6*double(doseSlice1>DoseCutOff));
set(get(cBarHandel,'ylabel'),'String', 'dose [Gy]','FontSize',defaultFontSize,'Interpreter','Latex');cBarHandel.FontSize = defaultFontSize; 
if plane ==3;axis equal; axis tight; end 
title([cNames{1} ' - ' cNames{2}],'Interpreter','Latex','Fontsize',defaultFontSize);

criteria = [2 2];
[gammaCube,costumMap,gammaPassRate] = matRad_gammaIndex(doseCube1,doseCube2,[ct.resolution.x ct.resolution.y ct.resolution.z],criteria);

subplot(122),imagesc(ct_rgb),hold on;
if plane     == 1
    gammaSlice = squeeze(gammaCube(slice,:,:));  
elseif plane == 2
    gammaSlice = squeeze(gammaCube(:,slice,:));
elseif plane == 3
    gammaSlice = squeeze(gammaCube(:,:,slice));
end
    
if FlagBBox; 
    gammaSlice = gammaSlice(bBox(2):bBox(4),bBox(1):bBox(3));
end   
h1=imagesc(gammaSlice,[0 2]);
set(h1,'AlphaData', .6*double(doseSlice1>DoseCutOff));
set(gca,'XTickLabel','');set(gca,'YTickLabel','');
colormap(costumMap);cBarHandel = colorbar(gca,'Location','eastoutside');  %southoutside
if plane ==3;axis equal; axis tight; end 
set(get(cBarHandel,'ylabel'),'String','passed voxels $<$ 1','fontsize',defaultFontSize,'Interpreter','Latex');
cBarHandel.FontSize = defaultFontSize;
title(['$\gamma$ index (' num2str(criteria(1)) '$\%$,'  num2str(criteria(2)) 'mm) = ' num2str(gammaPassRate,5) '$\%$'],'FontSize',defaultFontSize,'Interpreter','Latex');

if FlagSave
    set(gcf,'PaperOrientation','landscape','PaperUnits','normalized','PaperPosition', [0 0 1 1]);
    print('-dpdf','-r300',[FullFileName '_2']);
end















