
%%
clc
clear
close all
load('vDepthRel.mat');
%pathSpec = 'E:\TRiP98DATA_HIT-20131120\SPC\12C\RF3MM\FLUKA_NEW3_12C.H2O.MeV35000.xlsx';
pathSpec = '\\psf\Home\Documents\Heidelberg\TRiP98DATA\SPC\12C\RF3MM\FLUKA_NEW3_12C.H2O.MeV35000.xlsx';
[~,~,raw] = xlsread(pathSpec);
raw = raw(2:end,2:end);
rawNum = str2double(raw(:,1:10));
Cnt = 1;
for iDepth=1:79
    s(iDepth,1).depthStep = iDepth;
    s(iDepth,1).depth = rawNum(Cnt,2);
    s(iDepth,1).projectile = '12C';
    s(iDepth,1).target = 'H2O';
    s(iDepth,1).energy = raw{2,13};
    s(iDepth,1).peakPos = raw{2,14};
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
               s(iDepth,1).(sParticles{iPart}).Elow = cell2mat(Elow);
               Elow = [];
               s(iDepth,1).(sParticles{iPart}).Emid = cell2mat(Emid);
               Emid = [];
               s(iDepth,1).(sParticles{iPart}).Ehigh = cell2mat(Ehigh);
               Ehigh=[];
               s(iDepth,1).(sParticles{iPart}).dE = cell2mat(dE);
               dE=[];
               s(iDepth,1).(sParticles{iPart}).dNdE = cell2mat(dNdE);
               dNdE=[];
               s(iDepth,1).(sParticles{iPart}).N = cell2mat(N);
               N=[];
               % stop while loop
               break
           end
        end
        
    end 
end




%% 
sColor={'cyan','magenta','green','black','blue','red'};
vDepth = [s.depth];
figure,
for j = 1:length(sParticles)
    vY = zeros(78,1);
    for i = 1:78
        vY(i) = sum(s(i).(sParticles{j}).N);
    end
    plot(vDepth(1:78),vY,sColor{j},'Linewidth',3),hold on
end

legend(sParticles),grid on, xlabel('depth in cm','FontSize',14),ylabel('rel. particle fraction - fluence in cm^2','FontSize',14),
title('Energy = 350 MeV/u','FontSize',14);
set(gca,'FontSize',14');
set(gca,'YScale','log');
set(gca,'YLim',[1E-3,10]);


%% position 50 is right before brag peak
depth = 50;
data = s(depth,1);
h=figure,
for i = 1:length(sParticles)
 plot(data.(sParticles{i}).Emid,data.(sParticles{i}).dNdE,sColor{i},'Linewidth',3),hold on
end
legend(sParticles),grid on, xlabel('Energy','FontSize',14),ylabel('rel. number of particles per energy','FontSize',14);
title('depth = 20,9cm','FontSize',14),
set(gca,'FontSize',14');
set(gca,'YScale','log');
set(gca,'YLim',[.5E-5,0.1]);

%% load stopping powers
%path = 'E:\TRiP98DATA_HIT-20131120\DEDX\dEdxFLUKAxTRiP.dedx';
path = '\\psf\Home\Documents\Heidelberg\TRiP98DATA\DEDX\dEdxFLUKAxTRiP.dedx';
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
    plot(SP.(sParticles{i}).energy,SP.(sParticles{i}).dEdx,sColor{i},'Linewidth',3),hold on
end
legend(sParticles),grid on
set(gca,'YScale','log','XScale','log'),xlabel('Energy in MeV/u'),ylabel('stopping power in MeVcm^2/g'),
title('stopping powers of different particles');

%% load depth dose distributions
load(['baseDataHIT' filesep 'carbonBaseDataHIT.mat']);

[~,idx]=min(abs([carbonBaseDataHIT.energy]-s(1).energy));
baseData = carbonBaseDataHIT(idx);
D_accum = zeros(78,1);
sParticles=sParticles(1:6);
figure,
for i = 1:length(sParticles)
    for depth = 1:78;  
        D{depth} = s(depth,1).(sParticles{i}).N* ...
            interp1(SP.(sParticles{i}).energy,SP.(sParticles{i}).dEdx,s(depth,1).(sParticles{i}).Emid)'; 
    end
    plot(vDepth(1:78),cell2mat(D),sColor{i},'Linewidth',3),hold on
    D_accum = D_accum+cell2mat(D)';
end
plot(vDepth(1:78),D_accum,'Linewidth',3)
set(gca,'YScale','log')
set(gca,'YLim',[0.1 1000])
xlabel('depth in cm')
ylabel('dose in Gy')
title('particle dose distributions')
sParticles{1,7}='total dose';
legend(sParticles);
set(gca,'FontSize',14);
grid on

%% plot LET
LET = zeros(78,1);
sParticles=sParticles(1:6);
figure,
for i = 1:length(sParticles)
    for depth = 1:78; 
        
        SP_p = interp1(SP.(sParticles{i}).energy,SP.(sParticles{i}).dEdx,s(depth,1).(sParticles{i}).Emid)'; 
        dose{depth} = (s(depth,1).(sParticles{i}).N* SP_p);
            
        
        let_p{depth}=(s(depth,1).(sParticles{i}).N*(SP_p.^2))./dose{depth};
           
        
    end
    plot(vDepth(1:78),cell2mat(let_p),sColor{i},'Linewidth',3),hold on
    LET = LET + (cell2mat(let_p))';
end
plot(vDepth(1:78),LET,'Linewidth',3)
set(gca,'YScale','log')
set(gca,'YLim',[1 10000])
xlabel('depth in cm')
ylabel('LET')
title('particle LET distributions')
sParticles{1,7}='total let';
legend(sParticles);
set(gca,'FontSize',14);
grid on

%% compare depth dose curves
vD = baseData.depth/100;
figure,plot(vD,baseData.Z,'r','LineWidth',3),hold on,grid on
       plot(vDepth(1:78),D_accum,'LineWidth',3)
legend({'ddd orginal','ddd calculated'})

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

%path = 'E:\TRiP98DATA_HIT-20131120\RBE';
path = '\\psf\Home\Documents\Heidelberg\TRiP98DATA\RBE';
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
for j = 1:23
    CellType = j;
    figure,

    for i = 1:length(Spectra)
        vX = [meta(CellType).(Spectra{i}){1,2}(:).Energy];
        vY = [meta(CellType).(Spectra{i}){1,2}(:).RBE];
        plot(vX,vY,'Linewidth',3),hold on
    end

    str = sprintf('celltype: alpha_x: %f and beta_x: %f',meta(CellType).alpha,meta(CellType).beta);
    legend(Spectra),xlabel('energy [MeV/u]'),ylabel('RBE'),grid on
    title(str);
    set(gca,'FontSize',16)
end


%% asses alpha_p and beta_p
celltype = 5;
particle = 'carbon';

RBE_ini = meta(celltype).(particle)(2);
RBE_ini = RBE_ini{1,1};
alpha_x = meta(celltype).alpha;
alpha_p_ini = [RBE_ini.RBE]./alpha_x;
figure,plot([RBE_ini.Energy],alpha_p_ini)
% make it again from scretch
% convert RBE_initals to alpha_initals and apply fluence and SP according
% to formel in script to get depth dose alphas






