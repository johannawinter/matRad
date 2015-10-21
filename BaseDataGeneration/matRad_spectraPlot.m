%%
% This function is for plotting some curves out of spc files. In addition
% the RBE_inital data as well as the Stopping Power curves can be loaded
% and plotted. 

%% read spc file
clc
clear
close all

vEnergy = 280;
if ~ismac
    pathSpec = '\\Mac\Home\Documents\Heidelberg\TRiP98DATA\SPC\12C\RF3MM\FLUKA_NEW3_12C.H2O.MeV28000.mat';
else
    pathSpec = '\\Mac\Home\Documents\Heidelberg\TRiP98DATA\SPC\12C\RF3MM\FLUKA_NEW3_12C.H2O.MeV28000.xlsx';
end


load(pathSpec);
sParticles = MetaSPC.particles;
vDepth = [SPC.depths];

%% plot fluence 
sColor={'red','green','blue','red','green','blue','black'};
sLineSpec={'--','--','--','-','-' ,'-' ,'-'};

figure,
set(gcf,'Color',[1 1 1]);
for j = 1:SPC(1).numParticles
    vY = zeros(length(vDepth),1);
    for i = 1:length(vDepth)
        vY(i) = sum(SPC(i).(sParticles{j}).N);
        vY(i) = sum(SPC(i).(sParticles{j}).dNdE.*SPC(i).(sParticles{j}).dE);
    end
    subplot(231),plot(vDepth,vY,[sLineSpec{j} sColor{j}],'Linewidth',3),hold on
end

legend(sParticles),grid on, xlabel('depth in [cm]','FontSize',14),ylabel('rel. particle number per primary','FontSize',14),
title(['Energy = ' num2str(vEnergy) ' MeV/u'],'FontSize',14);
set(gca,'FontSize',14');
set(gca,'YScale','log');
set(gca,'YLim',[1E-5,2]);


%% plot energy spectra direct at peak

[~,idx]= min(abs([SPC.depths]-SPC(1).peakPos));

for i = 1:length(sParticles)
 subplot(232),plot(SPC(idx).(sParticles{i}).Emid,SPC(idx).(sParticles{i}).dNdE,[sLineSpec{i} sColor{i}],'Linewidth',3),hold on
end
legend(sParticles),grid on, xlabel('Energy in [MeV/u]','FontSize',14),ylabel('rel. number of particles per primary','FontSize',14);
title('energy spectra at peak','FontSize',14),
grid minor
set(gca,'FontSize',14');
set(gca,'YScale','log');
set(gca,'YLim',[.5E-5,0.1]);


%% load and display stopping powers
if ~ismac
    path = '\\Mac\Home\Documents\Heidelberg\TRiP98DATA\DEDX\dEdxFLUKAxTRiP.dedx';
else
    path = '\\psf\Home\Documents\Heidelberg\TRiP98DATA\DEDX\dEdxFLUKAxTRiP.dedx';
end

[ MetadEdx, dEdx ] = matRad_readdEdx( path );

for i = 1:length(sParticles)
    subplot(233),plot(dEdx.(sParticles{i}).Energy,dEdx.(sParticles{i}).dEdx,[sLineSpec{i} sColor{i}],'Linewidth',3),hold on
end
legend(sParticles),grid on
set(gca,'YScale','log','XScale','log'),xlabel('Energy in [MeV/u]'),ylabel('stopping power in [MeVcm^2/g]'),
title('stopping powers');

%% load depth dose distributions
cd('../')
load('carbonBaseDataHIT.mat');
cd('BioDataGeneration\')
[~,Idx]=min(abs([baseData.energy]-vEnergy));

DoseAccum = zeros(length(vDepth),1);
DoseParticle= zeros(length(vDepth),1);

for i = 1:length(sParticles)
    for x = 1:length(vDepth);  
        dEdxInterp = interp1(dEdx.(sParticles{i}).Energy,dEdx.(sParticles{i}).dEdx,SPC(x).(sParticles{i}).Emid,'linear','extrap')';
        DoseParticle(x) = SPC(x).(sParticles{i}).N*dEdxInterp; 
    end
    subplot(234),plot(vDepth,DoseParticle,[sLineSpec{i} sColor{i}],'Linewidth',3),hold on
    DoseAccum = DoseAccum+DoseParticle;
end
plot(vDepth,DoseParticle,[sLineSpec{i+1} sColor{i+1}],'LineWidth',3)
set(gca,'YScale','log')
set(gca,'YLim',[0.1 1000])
set(gca,'XLim',[0 30])
xlabel('depth in [cm]')
ylabel('MeVcm^2/g')
title('particle dose distributions')
sParticles{1,7}='total dose';
legend(sParticles);
sParticles = sParticles(1:6);
set(gca,'FontSize',14);
grid on

%% compare depth dose curves
subplot(235),plot(baseData(Idx).depths,baseData(Idx).Z,'r','LineWidth',4),hold on,grid on
subplot(235),plot(vDepth*10,DoseAccum,'k','LineWidth',4)
      
xlabel('depth in [cm]')
ylabel('[MeVcm^2/g]')       
legend({'ddd orginal','ddd calculated'})
title('comparison of calculated ddd and loaded ddd(TRiP)')
set(gca,'FontSize',14);

%% plot dose averaged LET

LETaccum = zeros(length(vDepth),1);

for i = 1:length(sParticles)
    
    LETmax = 0;
    LETnumerator = zeros(length(vDepth),1); 
    Dose          = zeros(length(vDepth),1); 
    
    for x = 1:length(vDepth); 
        
        dEdxInterp = interp1(dEdx.(sParticles{i}).Energy,dEdx.(sParticles{i}).dEdx,SPC(x).(sParticles{i}).Emid,'linear','extrap')';  
        LETnumerator(x)=(SPC(x).(sParticles{i}).N*(dEdxInterp.^2));
        Dose(x)=(SPC(x).(sParticles{i}).N*(dEdxInterp));

    end
    
    LETaccum = LETaccum + LETnumerator;
    subplot(236),plot(vDepth,LETnumerator./(Dose*10),[sLineSpec{i} sColor{i}],'Linewidth',3),hold on
    LET.(sParticles{i}).vDepth = vDepth;
    LET.(sParticles{i}).LET = LETnumerator./(Dose*10);
end
LET_tot = LETaccum./(DoseAccum*10);
subplot(236),plot(vDepth,LET_tot,[sLineSpec{i+1} sColor{i+1}],'Linewidth',3)
set(gca,'YScale','log')
set(gca,'YLim',[1 500]),set(gca,'XLim',[0 45])
xlabel('depth in cm')
ylabel('LET in [keV/µm]')
title('particle LET distributions')
sParticles{1,7}='total LET';
legend(sParticles);
sParticles = sParticles(1:6);
set(gca,'FontSize',14);
grid on


%% load RBE spc files

if ~ismac
    path = '\\Mac\Home\Documents\Heidelberg\TRiP98DATA\RBE';
else
    path = '\\psf\Home\Documents\Heidelberg\TRiP98DATA\RBE';
end

RBE = matRad_readRBE(path);

%% plot RBE spectra of specific cell type;
CellType = 1;
figure,
set(gcf,'Color',[1 1 1]);
for i = 1:length(RBE(1).particle)
    vX = [RBE(CellType).(RBE(1).particle{i}).Energy];
    vY = [RBE(CellType).(RBE(1).particle{i}).RBE];
    plot(vX,vY,'Linewidth',3),hold on;
end

str = sprintf('celltype: alpha_x: %f and beta_x: %f',RBE(CellType).alpha,RBE(CellType).beta);
legend(sParticles),xlabel('energy [MeV/u]'),ylabel('RBE'),grid on;
title(str);
set(gca,'FontSize',16)

%% asses alpha_p and beta_p
celltype = 13;
particle = 'carbon';
load('alphaEnergyInfnab2.mat')
load('alphaEnergyInfnabNB1RGB.mat')
sLineSpec2 = {':',':',':','-.','-.','-.'};
% extract meta data for current cell line
sParticles= MetaSPC.particles;
RBEcellLine = RBE(celltype);
alpha_x = RBE(celltype).alpha;
beta_x = RBE(celltype).beta;
Dcut = RBE(celltype).cut;
Smax = alpha_x+(2*beta_x)*Dcut;
Anuc = pi*(RBE(celltype).rnucleus^2); %µm^2
Anuc = Anuc/(10000^2);

% allocate variables
CntDepth   = length(vDepth);
a_num      = zeros(CntDepth,1);
a_numINFN = zeros(CntDepth,1);
a_numRapid = zeros(CntDepth,1);
b_numRapid = zeros(CntDepth,1);

figure(3),set(gcf,'Color',[1 1 1]);

for i = 1:length(MetaSPC.particles)
    
    Particle = RBEcellLine.particle{1,i};
    vRBE = RBEcellLine.(Particle).RBE;
    a_ion = (vRBE.*alpha_x)';
    subplot(121),
    plot(alphaInfn.(Particle).Energy,alphaInfn.(Particle).alpha,...
        [sLineSpec{i} sColor{i}],'LineWidth',2),hold on,grid on, grid minor,...
        title('alpha_{D} vs energy from RBE_{inital} using rapidScholz and INFN'),...
        set(gca,'XScale','log'),xlabel('Energy in MeV'),ylabel('raw alpha in Gy-1'),...
        set(gca,'FontSize',14);
    b_ion = (Smax-a_ion)./(2*Dcut);
    
    % rapid calculation according to Krämer
    LTE_RBE = interp1(dEdx.(Particle).Energy,dEdx.(Particle).dEdx,...
        RBEcellLine.(Particle).Energy,'pchip','extrap');
    d1 = ((1.602189e-10 .* LTE_RBE )/ Anuc);
    % S1 is the surviving fraction for a single particle traversal
    S1 = exp(-a_ion'.*d1);
    a_ion_rapid = (1-S1)./d1;
    % calculate scaling factor
    f = a_ion_rapid./a_ion';
    b_ion_rapid=(f.^2).*b_ion';
    subplot(121),plot(RBEcellLine.(Particle).Energy,...
        a_ion_rapid,[sLineSpec2{i} sColor{i}],'LineWidth',2),hold on,...
        grid on, grid minor, ...
        title('alpha_{D} vs energy from RBE_{inital} using rapidScholz and INFN'),...
        set(gca,'XScale','log'),xlabel('Energy in MeV'),ylabel('raw alpha in Gy-1'),...
        set(gca,'FontSize',14);
    
    % initialize some vectors
    a_Z       = zeros(CntDepth,1);
    a_Z_INFN  = zeros(CntDepth,1);
    a_Z_rapid = zeros(CntDepth,1);
    b_Z_rapid = zeros(CntDepth,1);
    dose_Z    = zeros(CntDepth,1);
    
    for depth = 1:CntDepth; 
        
        Fluence = (SPC(depth).(Particle).N);
        dEdx_interp = interp1(dEdx.(Particle).Energy,dEdx.(Particle).dEdx,...
            SPC(depth).(Particle).Emid,'pchip','extrap')';
        dose_Z(depth)= (Fluence*dEdx_interp);            
        
        a_ion_INFN_  = interp1(alphaInfn.(Particle).Energy,...
            alphaInfn.(Particle).alpha,SPC(depth).(Particle).Emid,'pchip','extrap');
        a_ion_rapid_ =interp1(RBEcellLine.(Particle).Energy,...
            a_ion_rapid,SPC(depth).(Particle).Emid,'pchip','extrap');
        b_ion_rapid_ =interp1(RBEcellLine.(Particle).Energy,b_ion_rapid,...
            SPC(depth).(Particle).Emid,'pchip','extrap');
        
        a_Z_INFN(depth) = (a_ion_INFN_.*dEdx_interp')*Fluence'; 
        a_Z_rapid(depth)= (a_ion_rapid_.*dEdx_interp')*Fluence';
        b_Z_rapid(depth)= (sqrt(b_ion_rapid_).*dEdx_interp')*Fluence';
        
    end
    
    subplot(122),plot(vDepth,a_Z_rapid./dose_Z,[sLineSpec{i} sColor{i}],'Linewidth',3),hold on
    
    a_numINFN  = a_numINFN + a_Z_INFN;
    a_numRapid = a_numRapid +a_Z_rapid;
    b_numRapid = b_numRapid+ b_Z_rapid;
end


subplot(122),plot(vDepth,(a_numRapid./DoseAccum),[sLineSpec{i+1} sColor{i+1}],'Linewidth',5),grid minor,grid on;

title('alphas contributions from a mono-energetic carbon ion beam with 350MeV/u');
sParticles{1,7}='mixed field alpha';
subplot(122),legend(sParticles)
xlabel('depth in [cm]');
ylabel('alpha in Gy^-1');
set(gca,'FontSize',14);
set(gca,'YLim',[0 2]),set(gca,'XLim',[0 30])
subplot(121),
legend({'H_{INFN}','H_{rapidScholz}','He_{INFN}','He_{rapidScholz}','Li_{INFN}','Li_{rapidScholz}','Be_{INFN}','Be_{rapidScholz}'...
    ,'B_{INFN}','B_{rapidScholz}','C_{INFN}','C_{rapidScholz}'})

load('carbonBaseData.mat');
load (['baseDataHIT' filesep 'C12_280MeVAlpha01.mat']);
load (['baseDataHIT' filesep 'refCNAOAlpha05E280.mat']);
load (['baseDataHIT' filesep 'refCNAObeta001E280.mat']);

[~,EnergyIdx] =(min(abs([baseData(:).energy]-vEnergy)));

%% alpha depth curves
figure,grid on,grid minor ,hold on,title('comparison of alpha-depth curves - 280MeV alpha_x = 0.1Gy^-1'),xlabel('depth in [cm]'),ylabel('alpha in Gy^-1')
       plot(baseData(EnergyIdx).depths./10,baseData(EnergyIdx).alpha(:,1),'k','Linewidth',3)
       plot(vDepth,(alpha_numeratorInfn./dose_accum),[sLineSpec{4} sColor{2}],'Linewidth',3),
       plot(vDepth,(alpha_numeratorRapid./dose_accum),[sLineSpec{4} sColor{3}],'Linewidth',3),
       
       vT = (0:0.005:30)*10;
       Z_interp=interp1(baseData(EnergyIdx).depths, baseData(EnergyIdx).Z,vT','pchip');
       AlphadEdx_interp=interp1(refCNAOALPHA.depth*10,refCNAOALPHA.alphadEdx,vT,'pchip');
       plot(vT/10,AlphadEdx_interp'./Z_interp,[sLineSpec{4} sColor{4}],'Linewidth',3),
       
       legend({'A.Mairani-LEM4 & CNAO data','from INFN & my SPC data','rapidScholz & my SPC data','A.Mairani-LEM1 & CNAO data'})
       set(gca,'FontSize',14),set(gca,'XLim',[0 30])
       
%% alpha*dEdx depth curves       
figure,grid on,grid minor,hold on,title('comparison of alpha-dose-depth curves - 280MeV alpha_x = 0.1Gy^-1'),xlabel('depth in [cm]'),ylabel('dEdx * alpha')
     plot(refCNAOALPHA.depth,refCNAOALPHA.alphadEdx,'r','Linewidth',3),set(gca,'FontSize',16)
     plot(vDepth,alpha_numeratorRapid,'b','Linewidth',3)
     legend({'RBE_{initial} & rapidScholz & my SPC','reference curve from A.Mairani(LEM1)'})
      set(gca,'FontSize',14),set(gca,'XLim',[0 30])
         

 %% beta depth curve     
figure,plot(vDepth,(beta_numeratorRapid./dose_accum).^2,'Linewidth',3),hold on,title('comparison of dose averaged beta depth curves')
     BetadEdx_interp=interp1(refCNAOALPHA.depth*10,refCNAOBETA.sqBetadEdx,vT,'pchip')';      
     plot(vT./10,(BetadEdx_interp./Z_interp).^2,'LineWidth',3)
      grid on, grid minor,xlabel('depth in cm'),ylabel('beta in Gy^-2'),set(gca,'Fontsize',14), legend({'beta from rapidScholz','beta from A.Mairani'})
%% beta * dEdx
figure,plot(refCNAOBETA.depth,refCNAOBETA.sqBetadEdx,'r','LineWidth',3),hold on
       plot(vDepth,beta_numeratorRapid,'b','Linewidth',3),hold on,title('comparison of dose averaged beta depth curves')
       grid on, grid minor,xlabel('depth in cm'),ylabel('sqrt(Beta)*dEdx in Gy^-2'),set(gca,'Fontsize',14),  legend({'RBE_{initial} & rapidScholz & my SPC','reference curve from A.Mairani(LEM1)'})
       set(gca,'XLim',[0 30])


vContribSPC=sum(ContribDepthSPC,1);
vContribINFN=sum(ContribDepthINFN,1);

figure,subplot(211),bar([ContribDepthSPC(:,1),ContribDepthINFN(:,1)]','stacked'),legend(sParticles(1:6)),xlabel({'(left) RBE_{ini}-rapidScholz-mySPC','(right) INFN-mySPC'}),ylabel('alpha contributions Gy^-1'),title('depth= 8.26cm, considering the whole energy range'),set(gca,'FontSize',13),grid on, grid minor
       subplot(212),bar([ContribDepthSPC(:,3),ContribDepthINFN(:,3)]','stacked'),legend(sParticles(1:6)),xlabel({'(left) RBE_{ini}-rapidScholz-mySPC','(right) INFN-mySPC'}),ylabel('alpha contributions Gy^-1'),title('depth= 8.26cm, considering only contributions from energies < 11MeV'),set(gca,'FontSize',13),grid on, grid minor     

figure, bar([ContribDepthSPC(:,2) ContribDepthINFN(:,2)]','stacked') ,legend(sParticles(1:6)),xlabel({'(left) RBE_{ini}-rapidScholz-mySPC','(right) INFN-mySPC'}),ylabel('dose in cGy'),title('dose contribution at depth= 8.26cm'),set(gca,'FontSize',13),grid on, grid minor      
       
figure,subplot(211),bar([ContribDepthSPC(:,1)/vContribSPC(2) ContribDepthINFN(:,1)/vContribINFN(2)]','stacked'),legend(sParticles(1:6)),xlabel({'(left) RBE_{ini}-rapidScholz-mySPC','(right) INFN-mySPC'}),ylabel('dose averaged alpha contributions Gy^-1'),title('depth= 8.26cm, considering the whole energy range'),set(gca,'FontSize',13),grid on, grid minor
       subplot(212),bar([ContribDepthSPC(:,3)/vContribSPC(4) ContribDepthINFN(:,3)/vContribINFN(4)]','stacked'),legend(sParticles(1:6)),xlabel({'(left) RBE_{ini}-rapidScholz-mySPC','(right) INFN-mySPC'}),ylabel('dose averaged alpha contributions Gy^-1'),title('depth= 8.26cm,  considering only contributions from energies < 11MeV'),set(gca,'FontSize',13),grid on, grid minor


       
%% plot GSI data
load('C:\Users\wieserh\Documents\matRad\BioDataGeneration\GSI_Chardoma_Carbon_BioData.mat')
figure,plot(str2num(SPC(1).peakPos)-stBioData{1,1}(3).Depths,stBioData{1,1}(3).Alpha,'Linewidth',3),grid on, grid minor, hold on
      plot(vDepth,(alpha_numeratorRapid./dose_accum),'Linewidth',3),
      legend({'from MTPS','alpha from SPC rapidScholz'}),xlabel('depth in [cm]'), ylabel('alpha in Gy^-1'), set(gca,'FontSize',14)
      set(gca,'Xlim',[-20,50])


