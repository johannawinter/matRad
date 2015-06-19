
%% read spc file
clc
clear
close all
vEnergy = 280;
if ~ismac
    pathSpec = 'E:\TRiP98DATA_HIT-20131120\SPC\12C\RF3MM\FLUKA_NEW3_12C.H2O.MeV28000.xlsx';
else
    pathSpec = '\\psf\Home\Documents\Heidelberg\TRiP98DATA\SPC\12C\RF3MM\FLUKA_NEW3_12C.H2O.MeV28000.xlsx';
end
READXLS = true;

%% read spc files from xls files or from .mat

if READXLS
    [~,~,raw] = xlsread(pathSpec);
    % skip the first column and first row
    raw = raw(2:end,2:end);
    % convert the data to numbers
    rawNum = str2double(raw(:,1:10));
    Cnt = 1;
    for iDepth=1:79
        SPC(iDepth,1).depthStep = iDepth;
        SPC(iDepth,1).depth = rawNum(Cnt,2);
        SPC(iDepth,1).projectile = '12C';
        SPC(iDepth,1).target = 'H2O';
        SPC(iDepth,1).energy = raw{2,13};
        SPC(iDepth,1).peakPos = raw{2,14};
        sParticles = {'H','He','Li','Be','B','C'};
        sParticlesNo = {'1002','2004','3006','4008','5010','6012'};
        for iPart = 1:length(sParticles)
           InnerCnt = 1;

            while true

                if Cnt > length(raw)
                    break;
                end

               if strcmp(sParticlesNo(iPart),raw{Cnt,3})
                    Elow{InnerCnt} = rawNum(Cnt,4);
                    Emid{InnerCnt} = rawNum(Cnt,5);
                    Ehigh{InnerCnt} = rawNum(Cnt,6);
                    dE{InnerCnt} = rawNum(Cnt,7);
                    dNdE{InnerCnt} = rawNum(Cnt,8);
                    N{InnerCnt} = rawNum(Cnt,9);
                    InnerCnt = InnerCnt+1;
                    Cnt = Cnt +1;
               else
                   SPC(iDepth,1).(sParticles{iPart}).Elow = cell2mat(Elow);
                   Elow = [];
                   SPC(iDepth,1).(sParticles{iPart}).Emid = cell2mat(Emid);
                   Emid = [];
                   SPC(iDepth,1).(sParticles{iPart}).Ehigh = cell2mat(Ehigh);
                   Ehigh=[];
                   SPC(iDepth,1).(sParticles{iPart}).dE = cell2mat(dE);
                   dE=[];
                   SPC(iDepth,1).(sParticles{iPart}).dNdE = cell2mat(dNdE);
                   dNdE=[];
                   SPC(iDepth,1).(sParticles{iPart}).N = cell2mat(N);
                   N=[];
                   % stop while loop
                   break
               end
            end

        end 
    end

    else

        path = ['baseDataHIT' filesep];

        for i = 80:10:440
            if i/10<=9
                name = ['C12spc' '0' num2str(i)];
            else
                name = ['C12spc' num2str(i)];
            end
            SPC = load([path name]);
            fName=fieldnames(SPC);
            SPC.(fName{1})(79).C.Elow = 0;
            SPC.(fName{1})(79).C.Emid =0;
            SPC.(fName{1})(79).Ehigh =0;
            SPC.(fName{1})(79).dE =0;
            SPC.(fName{1})(79).dNdE =0;
            SPC.(fName{1})(79).C.N =0;

            if i/10<=9
                save([name '.mat'],'SPC')
            else
                save([name '.mat'],'SPC')
            end

        end

end

% set last values to zero
SPC(79).C.Elow = 0;
SPC(79).C.Emid =0;
SPC(79).C.Ehigh =0;
SPC(79).C.dE =0;
SPC(79).C.dNdE =0;
SPC(79).C.N =0;
clearvars dE dNdE Ehigh Elow Emid iDepth Cnt InnerCnt iPart raw rawNum


vDepth = [SPC.depth];
%% delete carbon ionas behind peak
% for i = 67:length(vDepth)
%     SPC(i).C.N(:) = 0;
%      SPC(i).C.dNdE(:) = 0;
% end
%% plot fluence 
sColor={'red','green','blue','red','green','blue','black'};
sLineSpec={'--','--','--','-','-' ,'-' ,'-'};
figure,
for j = 1:length(sParticles)
    vY = zeros(length(vDepth),1);
    for i = 1:length(vDepth)
        vY(i) = sum(SPC(i).(sParticles{j}).N);
        %vY(i) = sum(SPC(i).(sParticles{j}).dNdE.*SPC(i).(sParticles{j}).dE);
    end
    plot(vDepth,vY,[sLineSpec{j} sColor{j}],'Linewidth',3),hold on
end

legend(sParticles),grid on, xlabel('depth in [cm]','FontSize',14),ylabel('rel. particle number; rel. fluence','FontSize',14),
title('Energy = 350 MeV/u','FontSize',14);
set(gca,'FontSize',14');
set(gca,'YScale','log');
set(gca,'YLim',[1E-5,2]);


%% plot energy spectra at position 50--> direct before brag peak
depth = 50;
data = SPC(depth,1);
h=figure;
for i = 1:length(sParticles)
 plot(data.(sParticles{i}).Emid,data.(sParticles{i}).dNdE,[sLineSpec{i} sColor{i}],'Linewidth',3),hold on
end
legend(sParticles),grid on, xlabel('Energy in [MeV/u]','FontSize',14),ylabel('rel. number of particles per energy','FontSize',14);
title('energy spectra right before bragg peak','FontSize',14),
grid minor
set(gca,'FontSize',14');
set(gca,'YScale','log');
set(gca,'YLim',[.5E-5,0.1]);


%% load and display stopping powers
if ~ismac
    path = 'E:\TRiP98DATA_HIT-20131120\DEDX\dEdxFLUKAxTRiP.dedx';
else
    path = '\\psf\Home\Documents\Heidelberg\TRiP98DATA\DEDX\dEdxFLUKAxTRiP.dedx';
end
fileID = fopen(path);
data = textscan(fileID,'%s');
data = data{1,1};
%skip RBE info
data = data(28:end,:);
Cnt = 1;
InnerCnt=1;
CntParticle = 1;
for i = 1:length(data)
   if Cnt>length(data)
        break; 
   end
   
   if strcmp(data(Cnt,1),'!dedx')
        Cnt=Cnt+1;
        
        while ~strcmp(data(Cnt,1),'!projectile')
            SP.(sParticles{CntParticle}).energy{InnerCnt} = str2double(data(Cnt,1));
            SP.(sParticles{CntParticle}).dEdx{InnerCnt} = str2double(data(Cnt+1,1));
            SP.(sParticles{CntParticle}).range{InnerCnt} = str2double(data(Cnt+2,1));
            InnerCnt = InnerCnt+1;
            Cnt = Cnt+3;
            if Cnt>length(data)
               break; 
            end
        end
        SP.(sParticles{CntParticle}).energy =(cell2mat(SP.(sParticles{CntParticle}).energy))';
        SP.(sParticles{CntParticle}).dEdx =(cell2mat(SP.(sParticles{CntParticle}).dEdx))';
        SP.(sParticles{CntParticle}).range =(cell2mat(SP.(sParticles{CntParticle}).range))';
        CntParticle = CntParticle+1;
        InnerCnt = 1;
   else
       Cnt = Cnt+1;
   end
end

figure,
for i = 1:length(sParticles)
    plot(SP.(sParticles{i}).energy,SP.(sParticles{i}).dEdx,[sLineSpec{i} sColor{i}],'Linewidth',3),hold on
end
legend(sParticles),grid on
set(gca,'YScale','log','XScale','log'),xlabel('Energy in [MeV/u]'),ylabel('stopping power in [MeVcm^2/g]'),
title('stopping powers of different particles');

%% load depth dose distributions
load(['baseDataHIT' filesep 'carbonBaseDataHIT.mat']);
if isnumeric(SPC(1).energy)
    targetEnergy = SPC(1).energy;
else
    targetEnergy =str2num(SPC(1).energy);
end
[~,idx]=min(abs([carbonBaseDataHIT.energy]-targetEnergy));
baseData = carbonBaseDataHIT(idx);
dose_accum = zeros(length(vDepth),1);
dose_particle= zeros(length(vDepth),1);
sParticles=sParticles(1:6);
figure,
for i = 1:length(sParticles)
    for depth = 1:length(vDepth);  
        SP_interp = interp1(SP.(sParticles{i}).energy,SP.(sParticles{i}).dEdx,SPC(depth,1).(sParticles{i}).Emid,'linear','extrap')';
        dose_particle(depth) = SPC(depth,1).(sParticles{i}).N*SP_interp; 
    end
    plot(vDepth,dose_particle,[sLineSpec{i} sColor{i}],'Linewidth',3),hold on
    dose_accum = dose_accum+dose_particle;
end
plot(vDepth,dose_accum,[sLineSpec{i+1} sColor{i+1}],'Linewidth',3)
set(gca,'YScale','log')
set(gca,'YLim',[0.1 1000])
set(gca,'XLim',[0 30])
xlabel('depth in [cm]')
ylabel('dose in [cGy]')
title('particle dose distributions')
sParticles{1,7}='total dose';
legend(sParticles);
set(gca,'FontSize',14);
grid on

%% compare depth dose curves
vD = baseData.depth/100;
figure,plot(vD,baseData.Z,'r','LineWidth',4),hold on,grid on
       plot(vDepth,dose_accum,'k','LineWidth',4)
      
xlabel('depth in [cm]')
ylabel('dose in [cGy]')       
legend({'ddd orginal','ddd calculated'})
title('comparison of calculated ddd based on simulated data and "full" simulated ddd')

set(gca,'FontSize',14);

%% plot dose averaged LET

LET_accum = zeros(length(vDepth),1);
sParticles=sParticles(1:6);
figure,

for i = 1:length(sParticles)
    
    LETmax = 0;
    LET_numerator = zeros(length(vDepth),1); 
    dose = zeros(length(vDepth),1); 
    for depth = 1:length(vDepth); 
        
        SP_interp = interp1(SP.(sParticles{i}).energy,SP.(sParticles{i}).dEdx,SPC(depth,1).(sParticles{i}).Emid,'linear','extrap')';  
        LET_numerator(depth)=(SPC(depth,1).(sParticles{i}).N*(SP_interp.^2));
        dose(depth)=(SPC(depth,1).(sParticles{i}).N*(SP_interp));
         % this can be excluded
%         if strcmp(sParticles{i},'C') &&  LET_numerator(depth)>LETmax
%            LETmax =  LET_numerator(depth);
%         elseif strcmp(sParticles{i},'C') && depth>10 && LET_numerator(depth)<LETmax
%            LET_numerator(depth)=0;
%         end

    end
    
    LET_accum = LET_accum + LET_numerator;
    plot(vDepth,LET_numerator./(dose*10),[sLineSpec{i} sColor{i}],'Linewidth',3),hold on
    LET.(sParticles{i}).vDepth = vDepth;
    LET.(sParticles{i}).LET = LET_numerator./(dose*10);
end
LET_tot = LET_accum./(dose_accum*10);
plot(vDepth,LET_tot,[sLineSpec{i+1} sColor{i+1}],'Linewidth',3)
set(gca,'YScale','log')
set(gca,'YLim',[1 500]),set(gca,'XLim',[0 45])
xlabel('depth in cm')
ylabel('LET in [keV/�m]')
title('particle LET distributions')
sParticles{1,7}='total LET';
legend(sParticles);
set(gca,'FontSize',14);
grid on



%% double lateral gaussian
% vX = -10:0.1:10;
% sigma1 = 1;
% sigma2 = 10;
% w = 0.15;
% vY_narr = 1/(sqrt(2*pi*(sigma1^2)))*exp(-(vX).^2/(2*(sigma1^2)));
% vY_bro = 1/(sqrt(2*pi*(sigma2^2)))*exp(-(vX).^2/(2*(sigma2^2)));
% vY = (1-w)*vY_narr+w*vY_bro;
% figure,plot(vX,vY)



%% load RBE spc files
Spectra = {'hydrogen','helium','lithium','beryllium','bor','carbon','nitrogen','oxygen','fluor','neon'};
if ~ismac
    path = 'E:\TRiP98DATA_HIT-20131120\RBE';
else
    path = '\\psf\Home\Documents\Heidelberg\TRiP98DATA\RBE';
end
folderInfo = dir(path);
CntFiles = 1;

for i = 1:length(folderInfo)
   
    if folderInfo(i).isdir ~=1 && folderInfo(i).bytes>0
       
        fileID = fopen([path filesep folderInfo(i).name]);
        data = textscan(fileID,'%s');
        data = data{1,1};
      
        FlagParse = true;
        Cnt = 1;
        % parse header
        while FlagParse
           
           if ~isempty(strfind(data{Cnt,1},'projectile'))
               FlagParse = false;
               break;
           elseif ~isempty(strfind(data{Cnt,1},'mean'))
               RBE(CntFiles).(strrep(data{Cnt,1},'!','')) = '';
                 Cnt = Cnt +1;
           else
               % parse input as number or as string
               if isempty(str2num(data{Cnt+1,1}))
                   val = data{Cnt+1,1};
               else
                   val = str2num(data{Cnt+1,1});
               end
               RBE(CntFiles).(strrep(data{Cnt,1},'!','')) = val;
               Cnt = Cnt +2;
           end

        end
      
        % parse projectile headers
        for SpecCnt = 1:length(Spectra)
            %RBE(CntFiles).(strrep(data{Cnt,1},'!','')){1,1} = data{Cnt+1,1
            RBE(CntFiles).(Spectra{SpecCnt}){1,1}.projectile = data{Cnt+1,1};
            RBE(CntFiles).(Spectra{SpecCnt}){1,1}.A = SpecCnt*2;
            RBE(CntFiles).(Spectra{SpecCnt}){1,1}.Z = SpecCnt;
            isNumber = true;
            CntDat = Cnt + 6;
            CntStruct = 1;
            while isNumber
                 if CntDat>length(data)
                    break;
                 end
                 if isempty(str2num(data{CntDat,1})) || isempty(str2num(data{CntDat+1,1})) 
                     isNumber = false;
                     break;
                 else
                     rbeDat(CntStruct).Energy = str2num(data{CntDat,1});
                     rbeDat(CntStruct).RBE = str2num(data{CntDat+1,1});
                     CntStruct = CntStruct+1;
                     CntDat = CntDat+2;
                 end
            end

            RBE(CntFiles).(Spectra{SpecCnt}){1,2} = rbeDat;
            Cnt = CntDat;
        end
        
     CntFiles = CntFiles+1;
    end
end


%% plot RBE spectra of specific cell type;

    CellType = 1;
    figure,

    for i = 1:length(Spectra)
        vX = [RBE(CellType).(Spectra{i}){1,2}(:).Energy];
        vY = [RBE(CellType).(Spectra{i}){1,2}(:).RBE];
        plot(vX,vY,'Linewidth',3),hold on;
    end

    str = sprintf('celltype: alpha_x: %f and beta_x: %f',RBE(CellType).alpha,RBE(CellType).beta);
    legend(Spectra),xlabel('energy [MeV/u]'),ylabel('RBE'),grid on;
    title(str);
    set(gca,'FontSize',16)



%% asses alpha_p and beta_p
celltype = 13;
particle = 'carbon';
load('baseDataHIT/alphaEnergyInfnab2.mat')
%load('baseDataHIT/alphaEnergyInfnabNB1RGB.mat')
sLineSpec2 = {':',':',':','-.','-.','-.'};
% extract meta data for current cell line
sParticles=sParticles(1:6);
RBE_ini = RBE(celltype);
alpha_x = RBE(celltype).alpha;
beta_x = RBE(celltype).beta;
Dcut = RBE(celltype).cut;
Smax = alpha_x+(2*beta_x)*Dcut;
Anuc = pi*(RBE(celltype).rnucleus^2); %�m^2
Anuc = Anuc/(10000^2);
sParticleLong = {'hydrogen','helium','lithium','beryllium','bor','carbon'};
% allocate variables
alpha_numerator = zeros(length(vDepth),1);
alpha_numeratorInfn = zeros(length(vDepth),1);
alpha_numeratorRapid = zeros(length(vDepth),1);
beta_numeratorRapid = zeros(length(vDepth),1);
% verify depth step
ContribDepthINFN = zeros(length(sParticles),4);
ContribDepthSPC = zeros(length(sParticles),4);

figure(9)
figure(10)
for i = 1:length(sParticles)
    
    RBE_ini_z = RBE_ini.(sParticleLong{i}){1,2};
    alpha_ion = ([RBE_ini_z.RBE].*alpha_x)';
    figure(10),plot(alphaInfn.(sParticles{i}).Energy,alphaInfn.(sParticles{i}).alpha,[sLineSpec{i} sColor{i}],'LineWidth',2),hold on,grid on, grid minor, title('alpha_{D} vs energy from RBE_{inital} using rapidScholz and INFN'),set(gca,'XScale','log'),xlabel('Energy in MeV'),ylabel('raw alpha in Gy-1'),set(gca,'FontSize',14);;
    
    beta_ion = (Smax-alpha_ion)./(2*Dcut);
    
        
    % rapid calculation according to Kr�mer
    LTE_RBE = interp1(SP.(sParticles{i}).energy,SP.(sParticles{i}).dEdx,[RBE(celltype).(sParticleLong{i}){1,2}.Energy],'pchip','extrap');
    d1 = ((1.602189e-10 .* LTE_RBE )/ Anuc);
    % S1 is the surviving fraction for a single particle traversal
    S1 = exp(-alpha_ion'.*d1);
    alpha_ion_rapid = (1-S1)./d1;
    % calculate scaling factor
    f = alpha_ion_rapid./alpha_ion';
    beta_ion_rapid=(f.^2).*beta_ion';
    figure(10),plot([RBE(celltype).(sParticleLong{i}){1,2}.Energy],alpha_ion_rapid,[sLineSpec2{i} sColor{i}],'LineWidth',2),hold on,grid on, grid minor, title('alpha_{D} vs energy from RBE_{inital} using rapidScholz and INFN'),set(gca,'XScale','log'),xlabel('Energy in MeV'),ylabel('raw alpha in Gy-1'),set(gca,'FontSize',14);;
    % initialize some vectors
    alpha_Z = zeros(length(vDepth),1);
    alpha_Z_Infn = zeros(length(vDepth),1);
    alpha_Z_rapid = zeros(length(vDepth),1);
    beta_Z_rapid = zeros(length(vDepth),1);
    dose_Z = zeros(length(vDepth),1);
    
    for depth = 1:length(vDepth); 
        
        Fluence = (SPC(depth,1).(sParticles{i}).N);
        SP_interp = interp1(SP.(sParticles{i}).energy,SP.(sParticles{i}).dEdx,SPC(depth,1).(sParticles{i}).Emid,'pchip','extrap')';
        dose_Z(depth)= (Fluence*SP_interp);            
        
        alpha_ion_Infn_interp = interp1(alphaInfn.(sParticles{i}).Energy,alphaInfn.(sParticles{i}).alpha,SPC(depth,1).(sParticles{i}).Emid,'pchip','extrap');
        alpha_ion_rapid_interp=interp1([RBE_ini_z.Energy],alpha_ion_rapid,SPC(depth,1).(sParticles{i}).Emid,'pchip','extrap');
        beta_ion_rapid_interp=interp1([RBE_ini_z.Energy],beta_ion_rapid,SPC(depth,1).(sParticles{i}).Emid,'pchip','extrap');
        

        alpha_Z_Infn(depth) = (alpha_ion_Infn_interp.*SP_interp')*Fluence'; 
        alpha_Z_rapid(depth)= (alpha_ion_rapid_interp.*SP_interp')*Fluence';
        beta_Z_rapid(depth)= (sqrt(beta_ion_rapid_interp).*SP_interp')*Fluence';
        
        if depth == 12
            Idx = SPC(depth,1).(sParticles{i}).Emid<11;
            ContribDepthSPC(i,1)= alpha_Z_rapid(depth);
            ContribDepthSPC(i,2)= dose_Z(depth);
            % cut energy
            SP_interpCut = interp1(SP.(sParticles{i}).energy,SP.(sParticles{i}).dEdx,SPC(depth,1).(sParticles{i}).Emid(Idx),'pchip','extrap')';
            alpha_Cut=interp1([RBE_ini_z.Energy],alpha_ion_rapid,SPC(depth,1).(sParticles{i}).Emid(Idx),'pchip','extrap');
            ContribDepthSPC(i,3)=((alpha_Cut.*SP_interpCut')*Fluence(Idx)'/dose_Z(depth));
            ContribDepthSPC(i,4)=((alpha_Cut.*SP_interpCut')*Fluence(Idx)'/(Fluence(Idx)*SP_interp(Idx)));
            
            ContribDepthINFN(i,1)= alpha_Z_Infn(depth);
            ContribDepthINFN(i,2)= dose_Z(depth);
            % cut energy
            SP_interpCut = interp1(SP.(sParticles{i}).energy,SP.(sParticles{i}).dEdx,SPC(depth,1).(sParticles{i}).Emid(Idx),'pchip','extrap')';
            alpha_Cut=interp1(alphaInfn.(sParticles{i}).Energy,alphaInfn.(sParticles{i}).alpha,SPC(depth,1).(sParticles{i}).Emid(Idx),'pchip','extrap');
            ContribDepthINFN(i,3)=((alpha_Cut.*SP_interpCut')*Fluence(Idx)'/dose_Z(depth));
            ContribDepthINFN(i,4)=((alpha_Cut.*SP_interpCut')*Fluence(Idx)'/(Fluence(Idx)*SP_interp(Idx)));
            
        end
        
    end
    
    figure(9),plot(vDepth,alpha_Z_rapid./dose_Z,[sLineSpec{i} sColor{i}],'Linewidth',3),hold on
    
    alpha_numeratorInfn  = alpha_numeratorInfn + alpha_Z_Infn;
    alpha_numeratorRapid = alpha_numeratorRapid +alpha_Z_rapid;
    beta_numeratorRapid  = beta_numeratorRapid+ beta_Z_rapid;
end


figure(9),plot(vDepth,(alpha_numeratorRapid./dose_accum),[sLineSpec{i+1} sColor{i+1}],'Linewidth',5),grid minor,grid on;

title('alphas contributions from a mono-energetic carbon ion beam with 350MeV/u');
sParticles{1,7}='mixed field alpha';
legend(sParticles)
xlabel('depth in [cm]');
ylabel('alpha in Gy^-1');
set(gca,'FontSize',14);
set(gca,'YLim',[0 3]),set(gca,'XLim',[0 30])
figure(10),legend({'H_{Infn}','H_{rapidScholz}','He_{Infn}','He_{rapidScholz}','Li_{Infn}','Li_{rapidScholz}','Be_{Infn}','Be_{rapidScholz}'...
    ,'B_{Infn}','B_{rapidScholz}','C_{Infn}','C_{rapidScholz}'})

load('carbonBaseData.mat');
load (['baseDataHIT' filesep 'C12_280MeVAlpha01.mat']);
[~,EnergyIdx] =(min(abs([baseData(:).energy]-vEnergy)));
figure,grid on,grid minor ,hold on,title('comparison of alpha-depth curves - 280MeV alpha_x = 0.1Gy^-1'),xlabel('depth in [cm]'),ylabel('alpha in Gy^-1')
       plot(baseData(EnergyIdx).depths/10,baseData(EnergyIdx).alpha(:,1),'k','Linewidth',3)
       plot(vDepth,(alpha_numeratorInfn./dose_accum),[sLineSpec{4} sColor{2}],'Linewidth',3),
       plot(vDepth,(alpha_numeratorRapid./dose_accum),[sLineSpec{4} sColor{3}],'Linewidth',3),
       %plot(C12_280MeVAlpha01(:,1),C12_280MeVAlpha01(:,2),[sLineSpec{4} sColor{4}],'Linewidth',3),
       legend({'A.Mairani-LEM4 & CNAO data','from INFN & my SPC data','rapidScholz & my SPC data','from INFN & Sarah Bruenings SPC data'})
       set(gca,'FontSize',14),set(gca,'XLim',[0 30])

figure,plot(vDepth,(beta_numeratorRapid./dose_accum).^2,'Linewidth',3),hold on,title('comparison of dose averaged beta depth curves')
      plot(baseData(EnergyIdx).depths/10,baseData(EnergyIdx).beta(:,1),'LineWidth',3)
      grid on, grid minor,xlabel('depth in cm'),ylabel('beta in Gy^-2'),set(gca,'Fontsize',14), legend({'beta from rapidScholz','beta from A.Mairani'})




vContribSPC=sum(ContribDepthSPC,1);
vContribINFN=sum(ContribDepthINFN,1);

figure,subplot(211),bar([ContribDepthSPC(:,1),ContribDepthINFN(:,1)]','stacked'),legend(sParticles(1:6)),xlabel({'(left) RBE_{ini}-rapidScholz-mySPC','(right) INFN-mySPC'}),ylabel('alpha contributions Gy^-1'),title('depth= 8.26cm, considering the whole energy range'),set(gca,'FontSize',13),grid on, grid minor
       subplot(212),bar([ContribDepthSPC(:,3),ContribDepthINFN(:,3)]','stacked'),legend(sParticles(1:6)),xlabel({'(left) RBE_{ini}-rapidScholz-mySPC','(right) INFN-mySPC'}),ylabel('alpha contributions Gy^-1'),title('depth= 8.26cm, considering only contributions from energies < 11MeV'),set(gca,'FontSize',13),grid on, grid minor     

figure, bar([ContribDepthSPC(:,2) ContribDepthINFN(:,2)]','stacked') ,legend(sParticles(1:6)),xlabel({'(left) RBE_{ini}-rapidScholz-mySPC','(right) INFN-mySPC'}),ylabel('dose in cGy'),title('dose contribution at depth= 8.26cm'),set(gca,'FontSize',13),grid on, grid minor      
       
figure,subplot(211),bar([ContribDepthSPC(:,1)/vContribSPC(2) ContribDepthINFN(:,1)/vContribINFN(2)]','stacked'),legend(sParticles(1:6)),xlabel({'(left) RBE_{ini}-rapidScholz-mySPC','(right) INFN-mySPC'}),ylabel('dose averaged alpha contributions Gy^-1'),title('depth= 8.26cm, considering the whole energy range'),set(gca,'FontSize',13),grid on, grid minor
       subplot(212),bar([ContribDepthSPC(:,3)/vContribSPC(4) ContribDepthINFN(:,3)/vContribINFN(4)]','stacked'),legend(sParticles(1:6)),xlabel({'(left) RBE_{ini}-rapidScholz-mySPC','(right) INFN-mySPC'}),ylabel('dose averaged alpha contributions Gy^-1'),title('depth= 8.26cm,  considering only contributions from energies < 11MeV'),set(gca,'FontSize',13),grid on, grid minor


       
%% plot GSI data

load('C:\Users\wieserh\Documents\matRad\GSI_Chardoma_Carbon_BioData.mat')
figure,plot(str2num(SPC(1).peakPos)-stBioData{1,1}(3).Depths,stBioData{1,1}(3).Alpha,'Linewidth',3),grid on, grid minor, hold on
      plot(vDepth,(alpha_numeratorRapid./dose_accum),'Linewidth',3),
      legend({'from MTPS','alpha from SPC rapidScholz'}),xlabel('depth in [cm]'), ylabel('alpha in Gy^-1'), set(gca,'FontSize',14)
      set(gca,'Xlim',[-20,50])


