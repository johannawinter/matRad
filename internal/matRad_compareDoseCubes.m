function matRad_compareDoseCubes(pln,energyIx,cube1,cube2,resolution,cellName,isSOPB,slice,filename)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% comparison of 3D cubes
% 
% call
%   matRad_compareDoseCubes(pln,energyIx,cube1,cube2,resolution,cellName,slice,filename)
%
% input
%   pln:        matRads pln struct
%   energyIx:   only necessary for pristine carbon beam evalution
%   cube1:      first 3D array
%   cube2:      second 3D array
%   resolution: resolution
%   filename:   optional: name of ps file for output
%
% output
%
%
% References
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015, Mark Bangert, on behalf of the matRad development team
%
% m.bangert@dkfz.de
%
% This file is not part of matRad.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if ~exist('slice','var')
   slice = 100; 
end

    
defaultLineWidth = 1.5;

% cut pristin carbon beam cube at last known depth position
if isequal(pln.radiationMode,'carbon') && ~isSOPB
    load([pln.radiationMode '_' pln.machine]);
    maxDepth = round(max(machine.data(energyIx).depths));
    if maxDepth< size(cube1,2)
     cube1 = cube1(:,end-maxDepth:end,:);
        cube2 = cube2(:,end-maxDepth:end,:);
    end
end

% necessary for sobp evaluation
if isSOPB
    cube2 = permute(cube2,[3 1 2]);
end

if ~isequal(size(cube1),size(cube2))
    error('cube dimensions inconsistent');
end

%% integral dose
relIntDoseDif = (1-sum(cube1(:))/sum(cube2(:)))*100;

fprintf(['Relative difference in integral dose: ' num2str(relIntDoseDif) '%%\n']);

figure
set(gcf,'Color',[1 1 1]);
%% transversal slices
subplot(3,2,1)
imagesc(cube1(:,:,slice))
colorbar
title([ cellName{1,1} ': slice ' num2str(slice)])

subplot(3,2,2)
imagesc(cube2(:,:,slice))
colorbar
title([ cellName{1,2} ': slice ' num2str(slice)])

norm = max(max(cube1(:,:,slice)));

ax3 = subplot(3,2,3);
imagesc(100*(cube1(:,:,slice)-cube2(:,:,slice))./norm)
myMap = matRad_getCostumColorbarDiff(cube1,cube2,slice,3);
colormap(ax3,myMap); colorbar;
title(['rel diff [%] ' cellName{1,1} ' - ' cellName{1,2} ' : slice ' num2str(slice)])

%% idds

x = resolution.x*[1/2:1:size(cube1,2)-1/2];

cube1IDD = sum(sum(cube1,3),1);
cube2IDD = sum(sum(cube2,3),1);

subplot(3,2,4)
hold on
plot(x,cube1IDD,'r','LineWidth',defaultLineWidth)
plot(x,cube2IDD,'b--','LineWidth',defaultLineWidth)
title(['IDD: rel diff int. dose ' num2str(relIntDoseDif)])
box on
grid on
legend(cellName,'Location','northeast')

%% profile along depth
slice = 100;
cube1CentralRayProf = squeeze(cube1(slice,:,slice));
cube2CentralRayProf = squeeze(cube2(slice,:,slice));

subplot(3,2,5)
hold on
plot(x,cube1CentralRayProf,'r','LineWidth',defaultLineWidth)
plot(x,cube2CentralRayProf,'b--','LineWidth',defaultLineWidth)
box on
grid on
title('depth profile')
legend(cellName,'Location','southeast')

%% profile along lateral direction entrance

y = resolution.y*[1/2:1:size(cube1,1)-1/2];

cube1LatProfileEnt = squeeze(cube1(:,end,slice));
cube2LatProfileEnt = squeeze(cube2(:,end,slice));

NormFac = 1;%;max(cube1LatProfileEnt);
cube1LatProfileEnt = cube1LatProfileEnt/NormFac;
cube2LatProfileEnt = cube2LatProfileEnt/NormFac;
 
subplot(3,2,6)
hold on
plot(y,cube1LatProfileEnt,'r','LineWidth',defaultLineWidth)
plot(y,cube2LatProfileEnt,'b--','LineWidth',defaultLineWidth)
box on
grid on
title('lateral entrance profile')
legend(cellName,'Location','northeast')

    
if nargin > 8
    h=gcf;
    set(h,'PaperPositionMode','auto');         
    set(h,'PaperOrientation','landscape');
    set(h,'Position',[50 50 1200 800]);
    print(gcf, '-dpdf', [filename '1.pdf'])
end


%% show either SOPB or pristine peak plot
if ~isSOPB
    figure,set(gcf,'Color',[1 1 1]);
    norm = max(max(cube1(:,:,slice)));
    ax1 = subplot(121);
    imagesc(100*(squeeze(cube1(:,:,slice))-squeeze(cube2(:,:,slice)))./norm);
    myMap = matRad_getCostumColorbarDiff(cube1,cube2,slice,3);
    colormap(ax1,myMap); colorbar;
    xlabel('y [mm]'),ylabel('x [mm]'),set(gca,'FontSize',14),title(['axial slice'])

    [~,slice] = max(cube1IDD);
    norm = max(max(cube1(:,slice,:)));
    ax1 = subplot(122);
    imagesc(100*(squeeze(cube1(:,slice,:))-squeeze(cube2(:,slice,:)))./norm);
    myMap = matRad_getCostumColorbarDiff(cube1,cube2,slice,2);
    colormap(ax1,myMap); colorbar;
    title(['sagittal slice at bragg peak'])
    xlabel('x [mm]'),ylabel('z [mm]'),set(gca,'FontSize',16)

    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    text(0.5, 1,'\bf pristine proton beam 156.35 MeV/u; rel. diff. [%] matRad - MC','HorizontalAlignment','center','VerticalAlignment', 'top')
else
    defFontSize = 16;
    figure,set(gcf,'Color',[1 1 1]);
    ax1=subplot(131),
    imagesc(flip(100*(cube1(:,:,slice)-cube2(:,:,slice))./norm,2))
    myMap = matRad_getCostumColorbarDiff(cube1,cube2,slice,3);
    colormap(ax1,myMap); colorbar('FontSize',defFontSize);set(gca,'xlim',[0 150])
    title('a) rel. diff. matRad - Monte Carlo','FontSize',defFontSize)
    rectangle('Position',[60 45 41 150],'Curvature',0.1,'LineWidth',3)
    xlabel('depth [mm]'),ylabel('x [mm]');
    set(gca,'FontWeight','bold','FontSize',defFontSize); 

    MCcentral     = cube2(100,:,100);
    MatradCentral = cube1(100,:,100);
    sMaxVal       = max(MCcentral);
    subplot(132),plot(flip(MatradCentral)/sMaxVal,'r','LineWidth',2),hold on
                 plot(flip(MCcentral)/sMaxVal,'b--','LineWidth',2)
                 xlabel('depth [mm]','LineWidth',defFontSize,'FontWeight','bold')
                 ylabel('rel. dose','LineWidth',defFontSize,'FontWeight','bold')
    set(gca,'xlim',[0 150]),
    l = legend({'matRad','Monte Carlo'},'LineWidth',30,'FontWeight','bold','Location','southwest');
    legend boxoff
    set(gca,'FontWeight','bold');   
    title('b) central profiles of SOBP')

    pSOBPmatRadLatProfile = squeeze(cube1(:,160,100));
    pSOBP_MCLatProfile    = squeeze(cube2(:,160,100));
    maxLat = max(pSOBPmatRadLatProfile);
    set(gca,'FontSize',defFontSize)
    grid on;
    subplot(133),plot(pSOBPmatRadLatProfile/maxLat,'r','LineWidth',2),hold on
                 plot(pSOBP_MCLatProfile/maxLat,'b--','LineWidth',2)
                 xlabel('x [mm]','LineWidth',defFontSize,'FontWeight','bold')
                 ylabel('rel. dose','LineWidth',defFontSize,'FontWeight','bold')
    l = legend({'matRad','Monte Carlo'},'LineWidth',30,'FontWeight','bold','Location','south');
    legend boxoff
    set(gca,'FontWeight','bold');  
    set(gca,'FontSize',defFontSize),
    grid on;
    title('c) lateral profiles of SOBP')  
    set(gca,'xlim',[0 240]),
end

if nargin > 8
    h=gcf;
    set(h,'PaperPositionMode','auto');         
    set(h,'PaperOrientation','landscape');
    set(h,'Position',[50 50 1200 800]);
    print(gcf, '-dpdf', [filename '2.pdf'])
end


%% calculate gamma pass rate
gammaCube = matRad_gammaIndex(cube1,cube2,[resolution.x resolution.y resolution.z] ,slice);
if nargin > 8
    h=gcf;
    set(h,'PaperPositionMode','auto');         
    set(h,'PaperOrientation','landscape');
    set(h,'Position',[50 50 1200 800]);
    print(gcf, '-dpdf', [filename '3.pdf'])
end


%% relative differences
figure,set(gcf,'Color',[1 1 1]);
subplot(2,2,1)
plot(x,100*(cube1IDD-cube2IDD)/max(cube2IDD),'r')
box on
grid on
title('rel diff [%] idd')

subplot(2,2,2)
plot(x,100*(cube1CentralRayProf-cube2CentralRayProf)/max(cube2CentralRayProf),'r')
box on
grid on
title('rel diff [%] central profile')

subplot(2,2,3)
plot(y,100*(cube1LatProfileEnt-cube2LatProfileEnt)/max(cube2LatProfileEnt),'r')
box on
grid on
title('rel diff [%] lateral profile @ entrance')

dim = size(cube1);
depthIx = round(dim(2)/2);
cube1LatProfileIx = squeeze(cube1(:,depthIx,slice));

dim = size(cube2);
depthIx = round(dim(2)/2);
cube2LatProfileIx = squeeze(cube2(:,depthIx,slice));

subplot(2,2,4)
plot(y,100*(cube1LatProfileIx-cube2LatProfileIx)/max(cube2LatProfileIx),'r')
box on
grid on
title('rel diff [%] lateral profile before peak')

if nargin > 8
    h=gcf;
    set(h,'PaperPositionMode','auto');         
    set(h,'PaperOrientation','landscape');
    set(h,'Position',[50 50 1200 800]);
    print(gcf, '-dpdf', [filename '4.pdf'])
end




end