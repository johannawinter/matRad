function matRad_plotDoseCube(ct,cst,pln,doseCube,varargin)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots a dose cube using matRad varibale 
% 
% call
%    matRad_plotDoseCube(ct,cst,pln,doseCube,slice,plane,filename)
%
% input
%   ct:          matRads ct struct
%   cst:         matRads cst struct
%   pln:         matRads pln struct
%   doseCube:    3D dose cube having the same size as the ct
%   slice:       scalar to determine the slice
%   plane:       scalar to determine the plane
%               (1=coronal,2=sagital,3=axial);
%   filename     if a filename is provided the figure will be saved as pdf
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015, Hans-Peter Wieser
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

% handle variable number of inputs 
if nargin >= 5
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

if nargin == 7
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
elseif nargin == 8
    if ischar(varargin{3})
        FlagSave = true;
        filename = varargin{1};
    end
end


DoseCutOff  = 0; % relative number between 0 and 1
defFontSize = 16;

if ~exist('plane','var')
    plane = 3;
end
if ~exist('slice','var')
    slice = round(pln.isoCenter(plane)/ct.resolution.z);
end

defaultFontSize = 16;

if plane == 1
    ctSlice      = squeeze(ct.cube{1}(slice,:,:));
    doseSlice    = squeeze(doseCube(slice,:,:));   
elseif plane == 2
    ctSlice      = squeeze(ct.cube{1}(:,slice,:));        
    doseSlice    = squeeze(doseCube(:,slice,:));

elseif plane == 3
    ctSlice      = squeeze(ct.cube{1}(:,:,slice));
    doseSlice    = squeeze(doseCube(:,:,slice));
end

vLevels     = [ 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.75 0.8 0.9 0.95 1 1.1 ];
maxDose    =  max(doseSlice(:));
vLevelDose  = maxDose  * vLevels;

figure('units','normalized','outerposition',[0 0 1 1]),set(gcf,'Color',[1 1 1]);

ct_rgb = ind2rgb(uint8(63*squeeze(ctSlice)/max(ct.cube{1}(:))),bone);
imagesc(ct_rgb),hold on,colormap(cm),
h1=imagesc(doseSlice);
set(gca,'XTickLabel','');set(gca,'YTickLabel','');
set(h1,'AlphaData', .8*(double(doseSlice)>DoseCutOff*maxDose));
title(['dose slice at: ' num2str(slice)  ' in plane: ' num2str(plane)],'Interpreter','Latex','FontSize',defaultFontSize)   
[~,~] = contour(doseSlice,vLevelDose,'LevelListMode','manual','LineWidth',1.5);
cBarHandel = colorbar(gca,'FontSize',defaultFontSize);
caxis(gca,[0, maxDose]);
set(get(cBarHandel,'ylabel'),'String', 'dose [Gy]','fontsize',defaultFontSize,'Interpreter','Latex');
set(cBarHandel,'YLim',[0 maxDose]);
axis equal

colors = colorcube; maskCst = zeros(ct.cubeDim);
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
         end
         
         if sum(tmpSlice(:)) > 0
             [c, h] = contour(gca,squeeze(maskCst(:,:,slice)),.5*[1 1],'LineWidth',2,'Color',colors(ColorCnt,:),'DisplayName',cst{i,2});
             cHandle = [cHandle h];
             ColorCnt = ColorCnt + 1;
             sLegend{CntLegend} = cst{i,2}; CntLegend = CntLegend +1;
         end
end
legend(cHandle,sLegend,'FontSize',defFontSize)

folderPath = [folder filesep 'exports'];
if exist('FolderPath','file')
    status = mkdir(folderPath);
end



if FlagSave
    FullFileName = [folderPath filesep filename '_' datestr(now, 'dd_mmm_yyyy_HH_MM')];
    set(gcf,'PaperOrientation','landscape','PaperUnits','normalized','PaperPosition', [0 0 1 1]);
    print('-dpdf','-r300',FullFileName);
    % consider downloading the export figure picture 
    % https://github.com/altmany/export_fig
    % export_fig([FullFileName '_' datestr(now, 'dd_mmm_yyyy_HH_MM')],'-png');
    
    % consider downloading the matlab2tiks package to create latex based plots using pgfplots 
    % https://github.com/matlab2tikz/matlab2tikz
    % FullFileName = [FullFileName '.tex']
    % matlab2tikz( 'FullFileName, 'height', '\fheight', 'width', '\fwidth'
    % --> goto https://github.com/matlab2tikz/matlab2tikz/wiki for more
    % information
end



