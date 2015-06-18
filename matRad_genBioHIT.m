
%% read spc file
clc
clear
close all


Spectra = {'hydrogen','helium','lithium','beryllium','bor','carbon','nitrogen','oxygen','fluor','neon'};
sParticles = {'H','He','Li','Be','B','C'};
sParticleLong = {'hydrogen','helium','lithium','beryllium','bor','carbon'};
sColor={'red','green','blue','red','green','blue','black'};
sLineSpec={'--','--','--','-','-' ,'-' ,'-'};
sLineSpec2 = {':',':',':','-.','-.','-.'};
DimEnergy = 37;
CntEnergy = 1;

%% load data

load(['baseDataHIT2' filesep 'dEdx.mat']);
load(['baseDataHIT2' filesep 'initialRBE.mat']);
load(['baseDataHIT2' filesep 'carbonBaseDataHIT.mat']);

mAlpha = zeros(79,length(initialRBE),DimEnergy);
mBeta = zeros(79,length(initialRBE),DimEnergy);
mDepth = zeros(79,length(initialRBE),DimEnergy);
%% get alphas and betas
path = ['baseDataHIT2' filesep];
for i=80:10:440
   
 for j = 1:length(initialRBE)
  if i/10<=9
        filename = ['C12spc' '0' num2str(i)];
  else
        filename = ['C12spc' num2str(i)];
  end
     
  load([path filename]);
  fName=fieldnames(SPC);
  SPC = SPC.(fName{1});
  vDepth = [SPC(:).depth];
  % extract meta data for current cell line
  RBE_ini = initialRBE(j);
  alpha_x = initialRBE(j).alpha;
  beta_x = initialRBE(j).beta;
  Dcut = initialRBE(j).cut;
  Smax = alpha_x+(2*beta_x)*Dcut;
  Anuc = pi*(initialRBE(j).rnucleus^2); %µm^2
  Anuc = Anuc/(10000^2);
  % allocate variables
  alpha_numeratorRapid = zeros(length(vDepth),1);
  beta_numeratorRapid = zeros(length(vDepth),1);
  dose_accum = zeros(length(vDepth),1);
  
  for k = 1:length(sParticles)
      
      RBE_ini_z = RBE_ini.(sParticleLong{k}){1,2};
      alpha_ion = ([RBE_ini_z.RBE].*alpha_x)';
      beta_ion = (Smax-alpha_ion)./(2*Dcut);
      % rapid calculation according to Krämer and Scholz 2006
      dEdx_interp_RBE = interp1(dEdx.(sParticles{k}).energy,dEdx.(sParticles{k}).dEdx,[initialRBE(j).(sParticleLong{k}){1,2}.Energy],'pchip','extrap');
      d1 = ((1.602189e-10 .* dEdx_interp_RBE )/ Anuc);
      % S1 is the surviving fraction for a single particle traversal
      S1 = exp(-alpha_ion'.*d1);
      alpha_ion_rapid = (1-S1)./d1;
      % calculate scaling factor
      f = alpha_ion_rapid./alpha_ion';
      beta_ion_rapid=(f.^2).*beta_ion';
      %initialize some vectors
      alpha_Z_rapid = zeros(length(vDepth),1);
      beta_Z_rapid = zeros(length(vDepth),1);
      dose_Z = zeros(length(vDepth),1);
        
      
      for depth = 1:length(vDepth); 
        
        Fluence = (SPC(depth,1).(sParticles{k}).N);
        SP_interp = interp1(dEdx.(sParticles{k}).energy,dEdx.(sParticles{k}).dEdx,SPC(depth,1).(sParticles{k}).Emid,'pchip','extrap')';
        dose_Z(depth)= (Fluence*SP_interp);            
        
        alpha_ion_rapid_interp=interp1([RBE_ini_z.Energy],alpha_ion_rapid,SPC(depth,1).(sParticles{k}).Emid,'pchip','extrap');
        beta_ion_rapid_interp=interp1([RBE_ini_z.Energy],beta_ion_rapid,SPC(depth,1).(sParticles{k}).Emid,'pchip','extrap');
        
        alpha_Z_rapid(depth)= (alpha_ion_rapid_interp.*SP_interp')*Fluence';
        beta_Z_rapid(depth)= (sqrt(beta_ion_rapid_interp).*SP_interp')*Fluence';
          
      
      end
      
      alpha_numeratorRapid = alpha_numeratorRapid +alpha_Z_rapid;
      beta_numeratorRapid  = beta_numeratorRapid+ beta_Z_rapid;
      dose_accum = dose_accum + dose_Z;
      
  end
     
 
  mAlpha(:,j,CntEnergy) = alpha_numeratorRapid./dose_accum;
  mBeta(:,j,CntEnergy) = (beta_numeratorRapid./dose_accum).^2;
  mDepth(:,j,CntEnergy)= vDepth;
%   figure,plot(mDepth(:,j,CntEnergy),mAlpha(:,j,CntEnergy)),grid on
%   figure,plot(mDepth(:,j,CntEnergy),mBeta(:,j,CntEnergy)),grid on
 end
 
 CntEnergy = CntEnergy+1;
 
end
% 
% 
% 
% %% asses alpha_p and beta_p
% vDepth = [SPC.depth];
% celltype = 1;
% particle = 'carbon';
% load('baseDataHIT/alphaEnergyInfnab2.mat')
% %load('baseDataHIT/alphaEnergyInfnabNB1RGB.mat')
% 
% % extract meta data for current cell line
% sParticles=sParticles(1:6);
% RBE_ini = RBE(celltype);
% alpha_x = RBE(celltype).alpha;
% beta_x = RBE(celltype).beta;
% Dcut = RBE(celltype).cut;
% Smax = alpha_x+(2*beta_x)*Dcut;
% Anuc = pi*(RBE(celltype).rnucleus^2); %µm^2
% Anuc = Anuc/(10000^2);
% 
% % allocate variables
% alpha_numerator = zeros(length(vDepth),1);
% alpha_numeratorInfn = zeros(length(vDepth),1);
% alpha_numeratorRapid = zeros(length(vDepth),1);
% beta_numeratorRapid = zeros(length(vDepth),1);
% 
% figure(9)
% figure(10)
% for i = 1:length(sParticles)
%     
%     RBE_ini_z = RBE_ini.(sParticleLong{i}){1,2};
%     alpha_ion = ([RBE_ini_z.RBE].*alpha_x)';
%     %figure(10),plot([RBE_ini_z.Energy],alpha_ion,[sLineSpec{i} sColor{i}],'LineWidth',4),hold on,grid on, grid minor, title('alpha_{ion} vs energy'),set(gca,'XScale','log'),xlabel('Energy in MeV'),ylabel('raw alpha in Gy-1'),set(gca,'FontSize',14);
%     figure(10),plot(alphaInfn.(sParticles{i}).Energy,alphaInfn.(sParticles{i}).alpha,[sLineSpec{i} sColor{i}],'LineWidth',2),hold on,grid on, grid minor, title('alpha_{ion} vs energy from RBE_{inital} and INFN'),set(gca,'XScale','log'),xlabel('Energy in MeV'),ylabel('raw alpha in Gy-1'),set(gca,'FontSize',14);;
%     
%     beta_ion = (Smax-alpha_ion)./(2*Dcut);
%     
%         
%     % rapid calculation according to Krämer
%     SP_interp_RBE = interp1(SP.(sParticles{i}).energy,SP.(sParticles{i}).dEdx,[RBE(celltype).(sParticleLong{i}){1,2}.Energy],'pchip','extrap');
%     d1 = ((1.602189e-10 .* SP_interp_RBE )/ Anuc);
%     % S1 is the surviving fraction for a single particle traversal
%     S1 = exp(-alpha_ion'.*d1);
%     alpha_ion_rapid = (1-S1)./d1;
%     % calculate scaling factor
%     f = alpha_ion_rapid./alpha_ion';
%     beta_ion_rapid=(f.^2).*beta_ion';
%     figure(10),plot([RBE(celltype).(sParticleLong{i}){1,2}.Energy],alpha_ion_rapid,[sLineSpec2{i} sColor{i}],'LineWidth',2),hold on,grid on, grid minor, title('alpha_{ion} vs energy from RBE_{inital} and INFN'),set(gca,'XScale','log'),xlabel('Energy in MeV'),ylabel('raw alpha in Gy-1'),set(gca,'FontSize',14);;
%    
%     % initialize some vectors
%     alpha_Z = zeros(length(vDepth),1);
%     alpha_Z_Infn = zeros(length(vDepth),1);
%     alpha_Z_rapid = zeros(length(vDepth),1);
%     beta_Z_rapid = zeros(length(vDepth),1);
%     dose_Z = zeros(length(vDepth),1);
%     
%     for depth = 1:length(vDepth); 
%         
%         Fluence = (SPC(depth,1).(sParticles{i}).N);
%         SP_interp = interp1(SP.(sParticles{i}).energy,SP.(sParticles{i}).dEdx,SPC(depth,1).(sParticles{i}).Emid,'pchip','extrap')';
%         dose_Z(depth)= (Fluence*SP_interp);            
%         
%         alpha_ion_interp = interp1([RBE_ini_z.Energy],alpha_ion,SPC(depth,1).(sParticles{i}).Emid,'pchip','extrap');
%         alpha_ion_Infn_interp = interp1(alphaInfn.(sParticles{i}).Energy,alphaInfn.(sParticles{i}).alpha,SPC(depth,1).(sParticles{i}).Emid,'pchip','extrap');
%         alpha_ion_rapid_interp=interp1([RBE_ini_z.Energy],alpha_ion_rapid,SPC(depth,1).(sParticles{i}).Emid,'pchip','extrap');
%         beta_ion_rapid_interp=interp1([RBE_ini_z.Energy],beta_ion_rapid,SPC(depth,1).(sParticles{i}).Emid,'pchip','extrap');
%         
%         alpha_Z(depth)      = (alpha_ion_interp.*SP_interp')*Fluence'; 
%         alpha_Z_Infn(depth) = (alpha_ion_Infn_interp.*SP_interp')*Fluence'; 
%         alpha_Z_rapid(depth)= (alpha_ion_rapid_interp.*SP_interp')*Fluence';
%         beta_Z_rapid(depth)= (sqrt(beta_ion_rapid_interp).*SP_interp')*Fluence';
%     end
%     
%     figure(9),plot(vDepth,alpha_Z_rapid./dose_Z,[sLineSpec{i} sColor{i}],'Linewidth',3),hold on
%     
%     alpha_numerator      =  alpha_numerator+alpha_Z;
%     alpha_numeratorInfn  = alpha_numeratorInfn + alpha_Z_Infn;
%     alpha_numeratorRapid = alpha_numeratorRapid +alpha_Z_rapid;
%     beta_numeratorRapid  = beta_numeratorRapid+ beta_Z_rapid;
% end
% 
% figure(9),plot(vDepth,(alpha_numeratorRapid./dose_accum),[sLineSpec{i+1} sColor{i+1}],'Linewidth',5),grid minor,grid on;
% title('alphas contributions from a mono-energetic carbon ion beam with 350MeV/u');
% sParticles{1,7}='mixed field alpha';
% legend(sParticles)
% xlabel('depth in [cm]');
% ylabel('alpha in Gy^-1');
% set(gca,'FontSize',14);
% set(gca,'YLim',[0 3]),set(gca,'XLim',[0 30])
% figure(10),legend({'H_{Infn}','H_{rapidScholz}','He_{Infn}','He_{rapidScholz}','Li_{Infn}','Li_{rapidScholz}','Be_{Infn}','Be_{rapidScholz}'...
%     ,'B_{Infn}','B_{rapidScholz}','C_{Infn}','C_{rapidScholz}'})
% 
% load('carbonBaseData.mat');
% load (['baseDataHIT' filesep 'C12_280MeVAlpha01.mat']);
% [~,EnergyIdx] =(min(abs([baseData(:).energy]-vEnergy)));
% figure,grid on,grid minor ,hold on,title('comparison of alpha vs depth'),xlabel('depth in [cm]'),ylabel('alpha in Gy^-1')
%        plot(baseData(EnergyIdx).depths/10,baseData(EnergyIdx).alpha(:,1),[sLineSpec{1} sColor{1}],'Linewidth',3)
%        plot(vDepth,(alpha_numeratorInfn./dose_accum),[sLineSpec{1} sColor{2}],'Linewidth',3),
%        plot(vDepth,(alpha_numeratorRapid./dose_accum),[sLineSpec{3} sColor{3}],'Linewidth',3),
%        plot(C12_280MeVAlpha01(:,1),C12_280MeVAlpha01(:,2),[sLineSpec{4} sColor{4}],'Linewidth',3),
%        legend({'from Andrea Mairani-LEM4','from INFN','from spc file rapidScholz','from INFN Sarah Brueningk'})
%        set(gca,'FontSize',14),set(gca,'XLim',[0 30])
% 
% figure,plot(vDepth,(beta_numeratorRapid./dose_accum).^2,'Linewidth',3),hold on,title('comparison of dose averaged beta depth curves')
%       plot(baseData(EnergyIdx).depths/10,baseData(EnergyIdx).beta(:,1),'LineWidth',3)
%       grid on, grid minor,xlabel('depth in cm'),ylabel('beta in Gy^-2'),set(gca,'Fontsize',14), legend({'beta from rapidScholz','beta from A.Mairani'})
%     
% %% plot GSI data
% % 
% % load('C:\Users\wieserh\Documents\matRad\GSI_Chardoma_Carbon_BioData.mat')
% % figure,plot(str2num(SPC(1).peakPos)-stBioData{1,1}(3).Depths,stBioData{1,1}(3).Alpha,'Linewidth',3),grid on, grid minor, hold on
% %       plot(vDepth,(alpha_numeratorRapid./dose_accum),'Linewidth',3),
% %       legend({'from MTPS','alpha from SPC rapidScholz'}),xlabel('depth in [cm]'), ylabel('alpha in Gy^-1'), set(gca,'FontSize',14)
% %       set(gca,'Xlim',[-20,50])
% 
% 
% 
% %%
% 
% % figure,plot(LET_Carbon.Infn(:,1),LET_Carbon.Infn(:,2)),hold on
% % plot(LET_Carbon.HIT.energy,LET_Carbon.HIT.dEdx./10)
% % legend({'LET from INFN','LET from dEdx file'})
