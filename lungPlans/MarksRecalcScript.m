% set const RBE mode for protons
pln.bioOptimization='const_RBExD';
% recalc!
tic
zzz = matRad_calcDoseDirect(ct,stf,pln,cst);
toc
clear zzz
% copy to resultGUI struct
resultGUI.matRadRecalc = xxx.RBExDose;

%%
% Let's plot single profiles that are perpendicular to the beam direction
% and a gamma index map
ixProfileX = round(pln.isoCenter(1,1)./ct.resolution.x);
ixProfileY = round(pln.isoCenter(1,2)./ct.resolution.y);
ixProfileZ = round(pln.isoCenter(1,3)./ct.resolution.z);
profileOrginalX = squeeze(resultGUI.RBExDose(:,ixProfileX,ixProfileZ));
profilematRadX  = squeeze(resultGUI.matRadRecalc(:,ixProfileX,ixProfileZ));
profileOrginalY = squeeze(resultGUI.RBExDose(ixProfileY,:,ixProfileZ));
profilematRadY  = squeeze(resultGUI.matRadRecalc(ixProfileY,:,ixProfileZ));
profileOrginalZ = squeeze(resultGUI.RBExDose(ixProfileY,ixProfileX,:));
profilematRadZ  = squeeze(resultGUI.matRadRecalc(ixProfileY,ixProfileX,:));
%
figure
subplot(2,2,1)
plot(profileOrginalX,'LineWidth',2),grid on,hold on,
plot(profilematRadX,'LineWidth',2),legend({'syngo','matRad'}),
xlabel('mm'),ylabel('Gy(RBE)'),title('left-right')
subplot(2,2,2)
plot(profileOrginalY,'LineWidth',2),grid on,hold on,
plot(profilematRadY,'LineWidth',2),legend({'syngo','matRad'}),
xlabel('mm'),ylabel('Gy(RBE)'),title('anterior-posterior')
subplot(2,2,3)
plot(profileOrginalZ,'LineWidth',2),grid on,hold on,
plot(profilematRadZ,'LineWidth',2),legend({'syngo','matRad'}),
xlabel('mm'),ylabel('Gy(RBE)'),title('inferior-superior')
subplot(2,2,4)


%%

matRad_gammaIndex(resultGUI.RBExDose,resultGUI.matRadRecalc, ...
[ct.resolution.x ct.resolution.y ct.resolution.z],[3 3],ixProfileZ);
%
set(gca,'XTick',[0:25:400])
set(gca,'XTickLabel',2*[0:25:400])
xlabel('[mm]')
set(gca,'YTick',[0:25:400])
set(gca,'YTickLabel',2*[0:25:400])
axis equal
axis([100 200 50 150])
ylabel('[mm]')

%%
figure
medianPTV = median(resultGUI.RBExDose(cst{14,4}{1}));
relDiff = (resultGUI.RBExDose-resultGUI.matRadRecalc)/medianPTV*100;
matRad_plotSliceWrapper(gca,ct,cst,1,relDiff,3,ixProfileZ+2);%,[],[],colorcube,[],doseWindow,[]);
title('rel. difference [%] (normalized to median PTV dose)')
axis([100 200 50 150])

%%
figure
%matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.RBExDose,3,ixProfileZ);%,[],[],colorcube,[],doseWindow,[]);
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.RBExDose,3,ixProfileZ, ...
                        [],[],[],[],[1 13],12*[10 30 50 70 90 95 107]/100,logical([0 0 1 1 1 1 1 0 0 0 0 0 0 0 0 1 1 1]),'[Gy (RBE)]',1,'LineWidth',2);%,colorBarLabel,boolPlotLegend,varargin)
plot(1/ct.resolution.y*pln.isoCenter(1,1)*[1 1],[1 ct.cubeDim(2)],'w--','LineWidth',2)
axis([100 200 50 150])
title('syngo')

%
figure
%matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.matRadRecalc,3,ixProfileZ);%,[],[],colorcube,[],doseWindow,[]);
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.matRadRecalc,3,ixProfileZ, ...
                        [],[],[],[],[1 13],12*[10 30 50 70 90 95 107]/100,logical([0 0 1 1 1 1 1 0 0 0 0 0 0 0 0 1 1 1]),'[Gy (RBE)]',1,'LineWidth',2);%,colorBarLabel,boolPlotLegend,varargin)
%plot(1/ct.resolution.y*pln.isoCenter(1,1)*[1 1],[1 ct.cubeDim(2)],'w--','LineWidth',2)
axis([100 200 50 150])
title('matRad')
%%
figure
plot(profileOrginalY,'LineWidth',2),grid on,hold on,
plot(profilematRadY,'LineWidth',2),legend({'syngo','matRad'}),
axis([100 200 0 13])
xlabel('mm'),ylabel('Gy(RBE)'),title('anterior-posterior')


