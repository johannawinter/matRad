%%
% This function is for plotting some curves out of spc files. In addition
% the RBE_inital data as well as the Stopping Power curves can be loaded
% and plotted. 

%% read spc file
% clc
% clear
% close all
%%

pathToTrip = 'E:\TRiP98DATA_HIT-20131120';

if ~ismac
    pathSpec = [pathToTrip filesep 'SPC' '\12C\RF3MM\FLUKA_NEW3_12C.H2O.MeV35000.mat'];
else
    pathSpec = '\\Mac\Home\Documents\Heidelberg\TRiP98DATA\SPC\12C\RF3MM\FLUKA_NEW3_12C.H2O.MeV28000.xlsx';
end

defaultFontSize = 14;
load(pathSpec);
sParticles = SPC.meta.particles;
vDepth = [SPC.data.depths];
vEnergy =  SPC.meta.energy;

%% plot fluence 
sColor   = {'red','green','blue','red','green','blue','black','red','green','blue','red','green','blue','black'};
sLineSpec= {'--' ,'--'   ,'--'  ,'-'  ,'-'    ,'-'   ,'-'    ,':'  ,':'    ,':'   ,'-.' ,'-.'   ,'-.'  ,'-.'};
sLineSpec2= {'-.' ,'-.' ,'-.'  ,':'  ,':'   ,':' };

figure,
set(gcf,'Color',[1 1 1]); 
for idxPart = 1:length(sParticles)
    vFlux = zeros(length(vDepth),1);
    for idxDepth = 1:length(vDepth)
        vFlux(idxDepth) = sum(SPC.data(idxDepth).(sParticles{idxPart}).N);
        vFlux(idxDepth) = sum(SPC.data(idxDepth).(sParticles{idxPart}).dNdE.*SPC.data(idxDepth).(sParticles{idxPart}).dE);
    end
    subplot(231),plot(vDepth,vFlux,[sLineSpec{idxPart} sColor{idxPart}],'Linewidth',3),hold on
end

legend(sParticles),grid on, xlabel('depth in [cm]','Interpreter','Latex','FontSize',defaultFontSize),ylabel('rel. particle number (p.p.)','Interpreter','Latex','FontSize',defaultFontSize),
title(['Energy = ' num2str(vEnergy) ' MeV/u'],'Interpreter','Latex','FontSize',defaultFontSize);
set(gca,'FontSize',defaultFontSize), grid minor
set(gca,'YScale','log');
set(gca,'YLim',[1E-5,2]);



%% plot energy spectra direct at peak

[~,idx]= min(abs([SPC.data.depths]-SPC.data(1).peakPos));

for i = 1:length(sParticles)
 subplot(232),plot(SPC.data(idx).(sParticles{i}).Emid,SPC.data(idx).(sParticles{i}).dNdE,[sLineSpec{i} sColor{i}],'Linewidth',3),hold on
end
legend(sParticles),grid on, xlabel('Energy in [MeV/u]','Interpreter','Latex','FontSize',defaultFontSize),ylabel('rel. number of particles per energy','Interpreter','Latex','FontSize',defaultFontSize);
title('energy spectra at peak','Interpreter','Latex','FontSize',defaultFontSize),
grid minor
set(gca,'FontSize',defaultFontSize');
set(gca,'YScale','log');

%% load and display stopping powers
[ MetadEdx, dEdx ] = matRad_readdEdx(pathToTrip);

for i = 1:length(sParticles)
    subplot(233),plot(dEdx.(sParticles{i}).Energy,dEdx.(sParticles{i}).dEdx,[sLineSpec{i} sColor{i}],'Linewidth',3),hold on
end
legend(sParticles),grid on
set(gca,'YScale','log','XScale','log'),xlabel('Energy in [MeV/u]','Interpreter','Latex','FontSize',defaultFontSize),ylabel('stopping power in $[\frac{MeVcm^2}{g}$]','Interpreter','Latex','FontSize',defaultFontSize),
title('stopping power','Interpreter','Latex','FontSize',defaultFontSize);

%% load depth dose distributions
if length(sParticles) > 1
    load('carbon_HIT.mat');
else 
    load('protons_HIT.mat');
end
[~,Idx]=min(abs([machine.data.energy]-vEnergy));

vDoseAccum    =  zeros(length(vDepth),1);
DoseParticle  =  zeros(length(vDepth),1);

for idxPart = 1:length(sParticles)
    for idxDepth = 1:length(vDepth);  
        dEdxInterp             = interp1(dEdx.(sParticles{idxPart}).Energy,dEdx.(sParticles{idxPart}).dEdx,...
                                         SPC.data(idxDepth).(sParticles{idxPart}).Emid,'linear','extrap')';
        DoseParticle(idxDepth) = SPC.data(idxDepth).(sParticles{idxPart}).N*dEdxInterp; 
    end
    subplot(234),plot(vDepth,DoseParticle,[sLineSpec{idxPart} sColor{idxPart}],'Linewidth',3),hold on
    vDoseAccum = vDoseAccum + DoseParticle;
end
plot(vDepth,vDoseAccum,[sLineSpec{i+1} sColor{idxPart+1}],'LineWidth',3)
set(gca,'YScale','log')
xlabel('depth in [cm]','Interpreter','Latex','FontSize',defaultFontSize)
ylabel('$[\frac{MeVcm^2}{g}$]','Interpreter','Latex','FontSize',defaultFontSize)
title('dose distributions per particle','Interpreter','Latex','FontSize',defaultFontSize)
sParticles{1,end+1}='total dose';
legend(sParticles);
sParticles = sParticles(1:end-1);
set(gca,'FontSize',defaultFontSize);
grid on, grid minor

%% compare depth dose curves
subplot(235),plot(vDepth,vDoseAccum,'k','LineWidth',4),hold on,grid on
subplot(235),plot(machine.data(Idx).depths/10,machine.data(Idx).Z,'r--','LineWidth',4)

      
xlabel('depth in [cm]','Interpreter','Latex','FontSize',defaultFontSize)
ylabel('$[\frac{MeVcm^2}{g}$]','Interpreter','Latex','FontSize',defaultFontSize)       
legend({'ddd orginal','ddd calculated'},'Interpreter','Latex','FontSize',defaultFontSize)
title('comparison of calculated ddd and loaded ddd(TRiP)','Interpreter','Latex','FontSize',defaultFontSize)
set(gca,'FontSize',defaultFontSize),grid on, grid minor;

%% plot dose averaged LET

vLET_Accum = zeros(length(vDepth),1);

for idxPart = 1:length(sParticles)
    
    LETmax          = 0;
    vLET_Numerator  = zeros(length(vDepth),1); 
    vDose           = zeros(length(vDepth),1); 
    
    for idxDepth = 1:length(vDepth); 
        
        dEdxInterp               = interp1(dEdx.(sParticles{idxPart}).Energy,dEdx.(sParticles{idxPart}).dEdx,SPC.data(idxDepth).(sParticles{idxPart}).Emid,'linear','extrap')';  
        vLET_Numerator(idxDepth) = (SPC.data(idxDepth).(sParticles{idxPart}).N*(dEdxInterp.^2));
        vDose(idxDepth)          = (SPC.data(idxDepth).(sParticles{idxPart}).N*(dEdxInterp));

    end
    
    vLET_Accum = vLET_Accum + vLET_Numerator;
    subplot(236),plot(vDepth,vLET_Numerator./(vDose*10),[sLineSpec{idxPart} sColor{idxPart}],'Linewidth',3),hold on
    LET.(sParticles{idxPart}).vDepth = vDepth;
    LET.(sParticles{idxPart}).LET = vLET_Numerator./(vDose*10);
end

LET_tot = vLET_Accum./(vDoseAccum*10);
subplot(236),plot(vDepth,LET_tot,[sLineSpec{idxPart+1} sColor{idxPart+1}],'Linewidth',3)
set(gca,'YScale','log')
xlabel('depth in [cm]','Interpreter','Latex','FontSize',defaultFontSize)
ylabel('LET in $[\frac{keV}{\mu}$]','Interpreter','Latex','FontSize',defaultFontSize)
title('particle LET distributions','Interpreter','Latex','FontSize',defaultFontSize)
sParticles{1,end+1}='total LET';
legend(sParticles);
sParticles = sParticles(1:end-1);
set(gca,'FontSize',defaultFontSize);
grid on,grid minor


%% load RBE spc files
RBE = matRad_readRBE(pathToTrip);

%% plot RBE spectra of specific cell type;
CellType = 13;
figure(2),
set(gcf,'Color',[1 1 1]);
for i = 1:length(RBE(1).particle)
    vX = [RBE(CellType).(RBE(1).particle{i}).Energy];
    vY = [RBE(CellType).(RBE(1).particle{i}).RBE];
    h1 = subplot(231);plot(vX,vY,[sLineSpec{i} sColor{i}],'Linewidth',2),hold on,grid on,grid minor;
end
grid on, grid minor,
xlabel('[MeV/u]','Interpreter','Latex','FontSize',defaultFontSize);
ylabel(['LEM I - $RBE_{ini}$ - celltype: ' num2str(CellType)],'Interpreter','Latex','FontSize',defaultFontSize);
title('LEM I - expected S of one cell by one random travelersal' ,'Interpreter','Latex','FontSize',defaultFontSize);
set(gca,'XScale','log');
set(gca,'XLim',[.9E-1,1000]);
legend(RBE(CellType).particle)

for i = 1:6
    vX = [RBE(CellType).(RBE(1).particle{i}).Energy];
    vY = [RBE(CellType).(RBE(1).particle{i}).RBE];
    h2 = subplot(232);plot(vX,vY*(RBE(CellType).alpha),[sLineSpec{i} sColor{i}],'Linewidth',3),hold on;
end

xlabel('[MeV/u]','Interpreter','Latex','FontSize',defaultFontSize);
ylabel(['LEM I - $\alpha_{local}$ - celltype: ' num2str(CellType)],'Interpreter','Latex','FontSize',defaultFontSize);
title('radiosensitivity $\alpha_{local}$','Interpreter','Latex','FontSize',defaultFontSize);
set(gca,'XScale','log'),grid on, grid minor;
set(gca,'XLim',[.9E-1,1000]);
legend(RBE(CellType).particle(1:6))


particle = 'carbon';
load('alphaEnergyInfnab2.mat')
% extract meta data for current cell line
sParticles= MetaSPC.particles;
RBEcellLine = RBE(CellType);
alpha_x = RBE(CellType).alpha;
beta_x = RBE(CellType).beta;
Dcut = RBE(CellType).cut;
Smax = alpha_x+(2*beta_x)*Dcut;
Anuc = pi*(RBE(CellType).rnucleus^2); %µm^2
Anuc = Anuc/(10000^2);

% allocate variables
CntDepth   = length(vDepth);
a_num      = zeros(CntDepth,1);
a_numINFN = zeros(CntDepth,1);
a_mixedField= zeros(CntDepth,1);
b_mixedField = zeros(CntDepth,1);


for i = 1:length(MetaSPC.particles)
    
    Particle = RBEcellLine.particle{1,i};
    vRBE = RBEcellLine.(Particle).RBE;
    a_local = (vRBE.*alpha_x)';
    subplot(233),
    plot(alphaInfn.(Particle).Energy,alphaInfn.(Particle).alpha,...
        [sLineSpec2{i} sColor{i}],'LineWidth',2),hold on,grid on, grid minor,...
        title('$\alpha_{lowDose}$ vs energy from $RBE_{inital}$ using rapidScholz and INFN','Interpreter','Latex'),...
        set(gca,'XScale','log'),xlabel('Energy in [MeV/u]','Interpreter','Latex'),...
        ylabel('$\alpha_{lowDose}$ in $Gy^{-1}$','Interpreter','Latex'),...
        set(gca,'FontSize',defaultFontSize),grid on,grid minor;
    b_local = (Smax-a_local)./(2*Dcut);
    
    % rapid calculation according to Krämer
    LTE_RBE = interp1(dEdx.(Particle).Energy,dEdx.(Particle).dEdx,...
        RBEcellLine.(Particle).Energy,'pchip','extrap');
    d1 = ((1.602189e-10 .* LTE_RBE )/ Anuc);
    % S1 is the surviving fraction for a single particle traversal
    S1 = exp(-a_local'.*d1);
    a_lowDose = (1-S1)./d1;
    % calculate scaling factor
    f = a_lowDose./a_local';
    b_lowDose=(f.^2).*b_local';
    h3 = subplot(233);plot(RBEcellLine.(Particle).Energy,...
        a_lowDose,[sLineSpec{i} sColor{i}],'LineWidth',2);hold on,
        grid on, grid minor

    % initialize some vectors
    a_Ztot       = zeros(CntDepth,1);
    a_Z_INFN  = zeros(CntDepth,1);
    a_Z       = zeros(CntDepth,1);
    b_Z       = zeros(CntDepth,1);
    dose_Z    = zeros(CntDepth,1);
    
    for depth = 1:CntDepth; 
        
        Fluence = (SPC(depth).(Particle).N);
        dEdx_interp = interp1(dEdx.(Particle).Energy,dEdx.(Particle).dEdx,...
            SPC(depth).(Particle).Emid,'pchip','extrap')';
        dose_Z(depth)= (Fluence*dEdx_interp);            
        
        a_ion_INFN_  = interp1(alphaInfn.(Particle).Energy,...
            alphaInfn.(Particle).alpha,SPC(depth).(Particle).Emid,'pchip','extrap');
        a_lowDoseInterp =interp1(RBEcellLine.(Particle).Energy,...
            a_lowDose,SPC(depth).(Particle).Emid,'pchip','extrap');
        b_lowDoseInterp =interp1(RBEcellLine.(Particle).Energy,b_lowDose,...
            SPC(depth).(Particle).Emid,'pchip','extrap');
        
        a_Z_INFN(depth) = (a_ion_INFN_.*dEdx_interp')*Fluence'; 
        a_Z(depth)= (a_lowDoseInterp.*dEdx_interp')*Fluence';
        b_Z(depth)= (sqrt(b_lowDoseInterp).*dEdx_interp')*Fluence';
        
    end
    
    h5 = subplot(235);plot(vDepth,a_Z./dose_Z,[sLineSpec{i} sColor{i}],'Linewidth',3),hold on
    grid on, grid minor, 
    xlabel('depth in [cm]','Interpreter','Latex','FontSize',defaultFontSize)
    ylabel('$\alpha_{E0_Z}$ in $[Gy^{-1}]$','Interpreter','Latex','FontSize',defaultFontSize)
    title(['$\alpha_{E0_Z}$ radiosensitivity per fragmentation of C-beam of ' num2str(vEnergy) ' MeV'],'Interpreter','Latex','FontSize',defaultFontSize)
    set(gca,'FontSize',defaultFontSize);
    h4 = subplot(234);plot(vDepth,dose_Z,[sLineSpec{i} sColor{i}],'Linewidth',3),hold on,set(gca,'YScale','log')
    grid on, grid minor, 
    xlabel('depth in [cm]','Interpreter','Latex','FontSize',defaultFontSize)
    ylabel('$[\frac{MeVcm^2}{g}$]','Interpreter','Latex','FontSize',defaultFontSize)
    title('dose distributions per particle','Interpreter','Latex','FontSize',defaultFontSize)
    set(gca,'FontSize',defaultFontSize);

    a_numINFN  = a_numINFN + a_Z_INFN;
    a_mixedField = a_mixedField +a_Z;
    b_mixedField = b_mixedField+ b_Z;
end

LegendString = {'H_{INFN}','H_{lowDose}','He_{INFN}','He_{lowDose}','Li_{INFN}','Li_{lowDose}','Be_{INFN}','Be_{lowDose}'...
    ,'B_{INFN}','B_{lowDose}','C_{INFN}','C_{lowDose}'};
legend(h3,LegendString),grid on, grid minor;
legend(h4,sParticles),set(h4,'YLim',[1e-2 1e4]),
legend(h5,sParticles),set(h5,'XLim',[0 max(vDepth)]),grid(h5,'on'),grid minor;
subplot(236),[AX,H1,H2] = plotyy(vDepth,vDoseAccum./max(vDoseAccum),vDepth,(a_mixedField./vDoseAccum));
grid on, grid minor
set(H1,'LineWidth',3);
set(H2,'LineWidth',3);
ylabel(AX(1),'relative dose','Interpreter','Latex');
ylabel(AX(2),'dose averaged alpha - $\alpha_{E0}$ in $[Gy^{-1}]$','Interpreter','Latex');
xlabel('depth in [mm]','Interpreter','Latex');
title(['radiosensitivy $\alpha_{E0}$ of mono energetic C beam of'  num2str(vEnergy) ' MeV'],'Interpreter','Latex');

% alpha depth curves
% 
% load('carbonmachine.data.mat');
% load (['machine.dataHIT' filesep 'C12_280MeVAlpha01.mat']);
% load (['machine.dataHIT' filesep 'refCNAOAlpha05E280.mat']);
% load (['machine.dataHIT' filesep 'refCNAObeta001E280.mat']);
% 
% [~,EnergyIdx] =(min(abs([machine.data(:).energy]-vEnergy)));
% 
% 
% figure,grid on,grid minor ,hold on,title('comparison of alpha-depth curves - 280MeV alpha_x = 0.1Gy^-1'),xlabel('depth in [cm]'),ylabel('alpha in Gy^-1')
%        plot(machine.data(EnergyIdx).depths./10,machine.data(EnergyIdx).alpha(:,1),'k','Linewidth',3)
%        plot(vDepth,(alpha_numeratorInfn./dose_accum),[sLineSpec{4} sColor{2}],'Linewidth',3),
%        plot(vDepth,(alpha_numeratorRapid./dose_accum),[sLineSpec{4} sColor{3}],'Linewidth',3),
%        
%        vT = (0:0.005:30)*10;
%        Z_interp=interp1(machine.data(EnergyIdx).depths, machine.data(EnergyIdx).Z,vT','pchip');
%        AlphadEdx_interp=interp1(refCNAOALPHA.depth*10,refCNAOALPHA.alphadEdx,vT,'pchip');
%        plot(vT/10,AlphadEdx_interp'./Z_interp,[sLineSpec{4} sColor{4}],'Linewidth',3),
%        
%        legend({'A.Mairani-LEM4 & CNAO data','from INFN & my SPC data','rapidScholz & my SPC data','A.Mairani-LEM1 & CNAO data'})
%        set(gca,'FontSize',14),set(gca,'XLim',[0 30])
%        
% %% alpha*dEdx depth curves       
% figure,grid on,grid minor,hold on,title('comparison of alpha-dose-depth curves - 280MeV alpha_x = 0.1Gy^-1'),xlabel('depth in [cm]'),ylabel('dEdx * alpha')
%      plot(refCNAOALPHA.depth,refCNAOALPHA.alphadEdx,'r','Linewidth',3),set(gca,'FontSize',16)
%      plot(vDepth,alpha_numeratorRapid,'b','Linewidth',3)
%      legend({'RBE_{initial} & rapidScholz & my SPC','reference curve from A.Mairani(LEM1)'})
%       set(gca,'FontSize',14),set(gca,'XLim',[0 30])
%          
% 
%  %% beta depth curve     
% figure,plot(vDepth,(beta_numeratorRapid./dose_accum).^2,'Linewidth',3),hold on,title('comparison of dose averaged beta depth curves')
%      BetadEdx_interp=interp1(refCNAOALPHA.depth*10,refCNAOBETA.sqBetadEdx,vT,'pchip')';      
%      plot(vT./10,(BetadEdx_interp./Z_interp).^2,'LineWidth',3)
%       grid on, grid minor,xlabel('depth in cm'),ylabel('beta in Gy^-2'),set(gca,'Fontsize',14), legend({'beta from rapidScholz','beta from A.Mairani'})
% %% beta * dEdx
% figure,plot(refCNAOBETA.depth,refCNAOBETA.sqBetadEdx,'r','LineWidth',3),hold on
%        plot(vDepth,beta_numeratorRapid,'b','Linewidth',3),hold on,title('comparison of dose averaged beta depth curves')
%        grid on, grid minor,xlabel('depth in cm'),ylabel('sqrt(Beta)*dEdx in Gy^-2'),set(gca,'Fontsize',14),  legend({'RBE_{initial} & rapidScholz & my SPC','reference curve from A.Mairani(LEM1)'})
%        set(gca,'XLim',[0 30])
% 
% 
% vContribSPC=sum(ContribDepthSPC,1);
% vContribINFN=sum(ContribDepthINFN,1);
% 
% figure,subplot(211),bar([ContribDepthSPC(:,1),ContribDepthINFN(:,1)]','stacked'),legend(sParticles(1:6)),xlabel({'(left) RBE_{ini}-rapidScholz-mySPC','(right) INFN-mySPC'}),ylabel('alpha contributions Gy^-1'),title('depth= 8.26cm, considering the whole energy range'),set(gca,'FontSize',13),grid on, grid minor
%        subplot(212),bar([ContribDepthSPC(:,3),ContribDepthINFN(:,3)]','stacked'),legend(sParticles(1:6)),xlabel({'(left) RBE_{ini}-rapidScholz-mySPC','(right) INFN-mySPC'}),ylabel('alpha contributions Gy^-1'),title('depth= 8.26cm, considering only contributions from energies < 11MeV'),set(gca,'FontSize',13),grid on, grid minor     
% 
% figure, bar([ContribDepthSPC(:,2) ContribDepthINFN(:,2)]','stacked') ,legend(sParticles(1:6)),xlabel({'(left) RBE_{ini}-rapidScholz-mySPC','(right) INFN-mySPC'}),ylabel('dose in cGy'),title('dose contribution at depth= 8.26cm'),set(gca,'FontSize',13),grid on, grid minor      
%        
% figure,subplot(211),bar([ContribDepthSPC(:,1)/vContribSPC(2) ContribDepthINFN(:,1)/vContribINFN(2)]','stacked'),legend(sParticles(1:6)),xlabel({'(left) RBE_{ini}-rapidScholz-mySPC','(right) INFN-mySPC'}),ylabel('dose averaged alpha contributions Gy^-1'),title('depth= 8.26cm, considering the whole energy range'),set(gca,'FontSize',13),grid on, grid minor
%        subplot(212),bar([ContribDepthSPC(:,3)/vContribSPC(4) ContribDepthINFN(:,3)/vContribINFN(4)]','stacked'),legend(sParticles(1:6)),xlabel({'(left) RBE_{ini}-rapidScholz-mySPC','(right) INFN-mySPC'}),ylabel('dose averaged alpha contributions Gy^-1'),title('depth= 8.26cm,  considering only contributions from energies < 11MeV'),set(gca,'FontSize',13),grid on, grid minor
% 
% 
%        
% %% plot GSI data
% load('C:\Users\wieserh\Documents\matRad\BioDataGeneration\GSI_Chardoma_Carbon_BioData.mat')
% figure,plot(str2num(SPC(1).peakPos)-stBioData{1,1}(3).Depths,stBioData{1,1}(3).Alpha,'Linewidth',3),grid on, grid minor, hold on
%       plot(vDepth,(alpha_numeratorRapid./dose_accum),'Linewidth',3),
%       legend({'from MTPS','alpha from SPC rapidScholz'}),xlabel('depth in [cm]'), ylabel('alpha in Gy^-1'), set(gca,'FontSize',14)
%       set(gca,'Xlim',[-20,50])


