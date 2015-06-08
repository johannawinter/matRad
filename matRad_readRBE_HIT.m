
%% read spc file
clc
clear
close all
load('vDepthRel.mat');
pathSpec = 'E:\TRiP98DATA_HIT-20131120\SPC\12C\RF3MM\FLUKA_NEW3_12C.H2O.MeV35000.xlsx';
%pathSpec = '\\psf\Home\Documents\Heidelberg\TRiP98DATA\SPC\12C\RF3MM\FLUKA_NEW3_12C.H2O.MeV35000.xlsx';
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
% set last values to zero
SPC(79).C.Elow = 0;
SPC(79).C.Emid =0;
SPC(79).C.Ehigh =0;
SPC(79).C.dE =0;
SPC(79).C.dNdE =0;
SPC(79).C.N =0;

clearvars dE dNdE Ehigh Elow Emid iDepth Cnt InnerCnt iPart raw rawNum


%% plot fluence 
sColor={'red','green','blue','red','green','blue','black'};
sLineSpec={'--','--','--','-','-' ,'-' ,'-'};
vDepth = [SPC.depth];
figure,
for j = 1:length(sParticles)
    vY = zeros(length(vDepth),1);
    for i = 1:length(vDepth)
        vY(i) = sum(SPC(i).(sParticles{j}).N);
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
set(gca,'FontSize',14');
set(gca,'YScale','log');
set(gca,'YLim',[.5E-5,0.1]);


%% load and display stopping powers
path = 'E:\TRiP98DATA_HIT-20131120\DEDX\dEdxFLUKAxTRiP.dedx';
%path = '\\psf\Home\Documents\Heidelberg\TRiP98DATA\DEDX\dEdxFLUKAxTRiP.dedx';
fileID = fopen(path);
data = textscan(fileID,'%s');
data = data{1,1};
%skip meta info
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
sParticles=sParticles(1:6);
figure,
for i = 1:length(sParticles)
    for depth = 1:length(vDepth);  
        SP_interp = interp1(SP.(sParticles{i}).energy,SP.(sParticles{i}).dEdx,SPC(depth,1).(sParticles{i}).Emid,'linear','extrap')';
        dose{depth} = SPC(depth,1).(sParticles{i}).N*SP_interp; 
    end
    plot(vDepth,cell2mat(dose),[sLineSpec{i} sColor{i}],'Linewidth',3),hold on
    dose_accum = dose_accum+cell2mat(dose)';
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
       plot(vDepth,dose_accum,'LineWidth',4)
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
   
end
plot(vDepth,LET_accum./(dose_accum*10),[sLineSpec{i+1} sColor{i+1}],'Linewidth',3)
set(gca,'YScale','log')
set(gca,'YLim',[1 2000]),set(gca,'XLim',[0 45])
xlabel('depth in cm')
ylabel('LET in [keV/µm]')
title('particle LET distributions')
sParticles{1,7}='total LET';
legend(sParticles);
set(gca,'FontSize',14);
grid on



%% double lateral gaussian
vX = -10:0.1:10;
sigma1 = 1;
sigma2 = 10;
w = 0.15;
vY_narr = 1/(sqrt(2*pi*(sigma1^2)))*exp(-(vX).^2/(2*(sigma1^2)));
vY_bro = 1/(sqrt(2*pi*(sigma2^2)))*exp(-(vX).^2/(2*(sigma2^2)));
vY = (1-w)*vY_narr+w*vY_bro;
figure,plot(vX,vY)



%% load RBE spc files
Spectra = {'hydrogen','helium','lithium','beryllium','bor','carbon','nitrogen','oxygen','fluor','neon'};

path = 'E:\TRiP98DATA_HIT-20131120\RBE';
%path = '\\psf\Home\Documents\Heidelberg\TRiP98DATA\RBE';
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
               meta(CntFiles).(strrep(data{Cnt,1},'!','')) = '';
                 Cnt = Cnt +1;
           else
               % parse input as number or as string
               if isempty(str2num(data{Cnt+1,1}))
                   val = data{Cnt+1,1};
               else
                   val = str2num(data{Cnt+1,1});
               end
               meta(CntFiles).(strrep(data{Cnt,1},'!','')) = val;
               Cnt = Cnt +2;
           end

        end
      
        % parse projectile headers
        for SpecCnt = 1:length(Spectra)
            %meta(CntFiles).(strrep(data{Cnt,1},'!','')){1,1} = data{Cnt+1,1
            meta(CntFiles).(Spectra{SpecCnt}){1,1}.projectile = data{Cnt+1,1};
            meta(CntFiles).(Spectra{SpecCnt}){1,1}.A = SpecCnt*2;
            meta(CntFiles).(Spectra{SpecCnt}){1,1}.Z = SpecCnt;
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

            meta(CntFiles).(Spectra{SpecCnt}){1,2} = rbeDat;
            Cnt = CntDat;
        end
        
     CntFiles = CntFiles+1;
    end
end


%% plot RBE spectra of specific cell type;

    CellType = 1;
    figure,

    for i = 1:length(Spectra)
        vX = [meta(CellType).(Spectra{i}){1,2}(:).Energy];
        vY = [meta(CellType).(Spectra{i}){1,2}(:).RBE];
        plot(vX,vY,'Linewidth',3),hold on;
    end

    str = sprintf('celltype: alpha_x: %f and beta_x: %f',meta(CellType).alpha,meta(CellType).beta);
    legend(Spectra),xlabel('energy [MeV/u]'),ylabel('RBE'),grid on;
    title(str);
    set(gca,'FontSize',16)



%% asses alpha_p and beta_p
celltype = 1;
particle = 'carbon';

sParticles=sParticles(1:6);
RBE_ini = meta(celltype);
alpha_x = meta(celltype).alpha;
beta_x = meta(celltype).beta;
Dcut = meta(celltype).cut;
Smax = alpha_x+(2*beta_x)*Dcut;
sParticleLong = {'hydrogen','helium','lithium','beryllium','bor','carbon'};
alpha_numerator = zeros(length(vDepth),1);
beta_numerator = zeros(length(vDepth),1);
dose_denominator = zeros(length(vDepth),1);
alpha_rapid = zeros(length(vDepth),1);

figure,

for i = 1:length(sParticles)
    
    RBE_ini_z = RBE_ini.(sParticleLong{i})(2);
    alpha_ion = ([RBE_ini_z{1,1}.RBE].*alpha_x)';
    beta_ion = (Smax-alpha_ion)./(2*Dcut);
   
    for depth = 1:length(vDepth); 
        
        Fluence = (SPC(depth,1).(sParticles{i}).N);
        SP_interp = interp1(SP.(sParticles{i}).energy,SP.(sParticles{i}).dEdx,SPC(depth,1).(sParticles{i}).Emid,'linear','extrap')';
        dose{depth}= (Fluence*SP_interp);            
        alpha_ion_interp = interp1([RBE_ini_z{1,1}.Energy],alpha_ion,SPC(depth,1).(sParticles{i}).Emid,'linear','extrap');
        numeratorA{depth} = (alpha_ion_interp.*SP_interp')*Fluence'; 
        
        beta_ion_interp = interp1([RBE_ini_z{1,1}.Energy],beta_ion,SPC(depth,1).(sParticles{i}).Emid,'linear','extrap');
        numeratorB{depth} = (Fluence.*sqrt(beta_ion_interp)*SP_interp);
        
        alpha_ion_avg{depth} = numeratorA{depth}/dose{depth};
        
        Fluence2 = (SPC(depth,1).(sParticles{i}).dNdE).*(SPC(depth,1).(sParticles{i}).dE);
        dose{depth}= (Fluence2*SP_interp);  
        alpha_rapidA{depth} = (alpha_ion_interp.*SP_interp')*Fluence2';         
        
    end
    plot(vDepth,cell2mat(alpha_ion_avg),[sLineSpec{i} sColor{i}],'Linewidth',3),hold on
    alpha_numerator =  alpha_numerator+cell2mat(numeratorA)';
    beta_numerator =  beta_numerator+cell2mat(numeratorB)';
    dose_denominator = dose_denominator+cell2mat(dose)';  
    alpha_rapid = alpha_rapid+cell2mat(alpha_rapidA)';
end

plot(vDepth,alpha_numerator./dose_denominator,[sLineSpec{i+1} sColor{i+1}],'Linewidth',5),grid minor,grid on;
title('alphas contributions from a mono-energetic carbon ion beam with 350MeV/u');
sParticles{1,7}='mixed field alpha';
legend(sParticles)
xlabel('depth in [cm]');
ylabel('alpha in Gy^-1');
set(gca,'FontSize',14);
set(gca,'YLim',[0 3]),set(gca,'XLim',[0 30])

load('carbonBaseData.mat');
figure,plot(vDepth,(alpha_numerator./dose_denominator),'k','Linewidth',3),grid on,hold on,title('comparison of alphas obtained from A.Maraini and spc-files')
       plot(baseData(96).depths/10,baseData(61).alpha(:,1),'Linewidth',3), legend({'from spc file','from Mairani'}),xlabel('depth in [cm]'),ylabel('alpha in Gy^-1')
set(gca,'FontSize',14);
       
figure,plot(vDepth,(sqrt(beta_numerator./dose_denominator)),'k','Linewidth',3),title('comparison of betas obtained from A.Maraini and spc-files'),grid on, hold on,
       plot(baseData(96).depths/10,baseData(61).beta(:,1),'Linewidth',3), legend({'from spc file','from Mairani'}),xlabel('depth in [cm]'),ylabel('beta in Gy^-2')
set(gca,'FontSize',14);





