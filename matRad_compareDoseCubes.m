function matRad_compareDoseCubes(cube1,cube2,resolution,filename)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% comparison of 3D cubes
% 
% call
%   matRad_compareDoseCubes(cube1,cube2,resolution,filename)
%
% input
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

if ~isequal(size(cube1),size(cube2))
    error('cube dimensions inconsistent');
end

%% integral dose
relIntDoseDif = (1-sum(cube1(:))/sum(cube2(:)))*100;

fprintf(['Relative difference in integral dose: ' num2str(relIntDoseDif) '%%\n']);

figure

%% transversal slices
slice = 100;
subplot(3,2,1)
imagesc(cube1(:,:,slice))
colorbar
title(['cube 1: slice ' num2str(slice)])

subplot(3,2,2)
imagesc(cube2(:,:,slice))
colorbar
title(['cube 2: slice ' num2str(slice)])

norm = max(max(cube1(:,:,slice)));

subplot(3,2,3)
imagesc(100*(cube1(:,:,slice)-cube2(:,:,slice))./norm)
colorbar
title(['rel diff [%] cube 1 - cube 2: slice ' num2str(slice)])

%% idds
x = resolution.x*[1/2:1:size(cube1,2)-1/2];

cube1IDD = sum(sum(cube1,3),1);
cube2IDD = sum(sum(cube2,3),1);

subplot(3,2,4)
hold on
plot(x,cube1IDD,'r')
plot(x,cube2IDD,'b--')
title('IDD')
box on
grid on
legend({'cube 1','cube 2'})

%% profile along depth

cube1CentralRayProf = squeeze(cube1(slice,:,slice));
cube2CentralRayProf = squeeze(cube2(slice,:,slice));

subplot(3,2,5)
hold on
plot(x,cube1CentralRayProf,'r')
plot(x,cube2CentralRayProf,'b--')
box on
grid on
title('depth profile')
legend({'cube 1','cube 2'})

%% profile along lateral direction entrance

y = resolution.y*[1/2:1:size(cube1,1)-1/2];

cube1LatProfileEnt = squeeze(cube1(:,end,slice));
cube2LatProfileEnt = squeeze(cube2(:,end,slice));

subplot(3,2,6)
hold on
plot(y,cube1LatProfileEnt,'r')
plot(y,cube2LatProfileEnt,'b--')
box on
grid on
title('lateral entrance profile')
legend({'cube 1','cube 2'})

if nargin > 3
    annotation('textbox', [0 0.9 1 0.1],'String', filename, ...
    'EdgeColor', 'none','HorizontalAlignment', 'center')
    set(gcf,'PaperOrientation','landscape','PaperUnits','normalized','PaperPosition', [0 0 1 1]);
    print('-dpsc','-r300','-append',filename)
end
    
%% relative differences
figure
subplot(3,2,1)
plot(x,100*(cube1IDD-cube2IDD)/max(cube2IDD),'r')
box on
grid on
title('rel diff [%] idd')

subplot(3,2,2)
plot(x,100*(cube1CentralRayProf-cube2CentralRayProf)/max(cube2CentralRayProf),'r')
box on
grid on
title('rel diff [%] central profile')

subplot(3,2,3)
plot(y,100*(cube1LatProfileEnt-cube2LatProfileEnt)/max(cube2LatProfileEnt),'r')
box on
grid on
title('rel diff [%] lateral profile @ entrance')

depthIx = 230;

cube1LatProfileIx = squeeze(cube1(:,depthIx,slice));
cube2LatProfileIx = squeeze(cube2(:,depthIx,slice));

subplot(3,2,4)
plot(y,100*(cube1LatProfileIx-cube2LatProfileIx)/max(cube2LatProfileIx),'r')
box on
grid on
title('rel diff [%] lateral profile before peak')

if nargin > 3
    annotation('textbox', [0 0.9 1 0.1],'String', filename, ...
    'EdgeColor', 'none','HorizontalAlignment', 'center')
    set(gcf,'PaperOrientation','landscape','PaperUnits','normalized','PaperPosition', [0 0 1 1]);
    print('-dpsc','-r300','-append',filename)
end
