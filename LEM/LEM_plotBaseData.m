clc,clear,close all,
%load('carbon_Generic.mat');
load('carbon_HIT.mat');

extraAxisAptions = ['title style={font=\normalsize},'...
                    'xlabel style={font=\normalsize},'...
                    'ylabel style={font=\normalsize},',...
                    'legend style={font=\normalsize},',...
                    'ticklabel style={font=\small},'];

if strcmp(machine.meta.machine,'Generic') 
    ixAlphas = 1:1:2;
else
    ixAlphas = 1:1:10;
end

FLAG_SAVE       = true;    
defaultFontSize = 16;
foldername      = 'plotBaseData';
exportPath      = [pwd filesep 'LEM' filesep 'expFigures' filesep foldername];   
color           = colorspecs();   
addpath(exportPath) ; showInfoFlag = false; addpath('LEM'); addpath(['APM' filesep 'expFigures']); 

%% plot various alpha_p curves
filename  = ['LEM_alphaParticlesProfiles_' machine.meta.machine '.tex'];
latexPath = [exportPath filesep filename];
ix = 95;

figure(),set(gcf,'Color',[1 1 1],'Units', 'normalized', 'Position', [0.1,0.1,0.8,0.8] ),
plot(machine.data(ix).depths,machine.data(ix).alpha(:,ixAlphas),'LineWidth',2),grid on

xlabel('z[mm]','Interpreter','Latex');
ylabel('$Gy^{-1}$','Interpreter','Latex');
title(['depth-dependend dose-averaged $\alpha$-profiles;',' C12, $E_0$=' num2str(machine.data(ix).energy)],'Interpreter','Latex')

string = {strsplit(num2str(machine.data(ix).alphaBetaRatio(ixAlphas)),' ')};
for i = 1:numel(string{1})
   string{1}{1,i} = ['$\frac{\alpha_x}{\beta_x}$ ratio = '  string{1}{1,i}] ;
end

legend(string{1},'Interpreter','latex','FontSize',defaultFontSize)
set(gca,'Fontsize',defaultFontSize)

if FLAG_SAVE && exist('matlab2tikz','file') == 2
    cleanfigure;
    matlab2tikz([latexPath],'height', '12cm', 'width', '21cm','showInfo',showInfoFlag,'standalone', true,...
        'extraaxisoptions',extraAxisAptions);         

    currPath = pwd; cd(exportPath);
    if ispc  
         command = sprintf('pdflatex %s',filename);
    elseif ismac
         command = sprintf('/Library/Tex/texbin/pdflatex %s',[filename]);
    end
    [status,cmdout] = system(command);delete('*.aux');delete('*.log'); cd(currPath);
    if status > 0
        warning(['couldnt compile pdf: ' cmdout]);
    end
end


%% plot various alpha_p curves divided by the next lower alpha_p curve
filename   = ['LEM_alphaParticlesProfiles2_' machine.meta.machine '.tex'];
latexPath  = [exportPath filesep filename];

figure(),set(gcf,'Color',[1 1 1],'Units', 'normalized', 'Position', [0.1,0.1,0.8,0.8] ),
plot(machine.data(ix).depths,machine.data(ix).alpha(:,ixAlphas(2:end))./machine.data(ix).alpha(:,ixAlphas(1:end-1)),'LineWidth',2),grid on

xlabel('z[mm]','Interpreter','Latex');
ylabel('$dimensionless$','Interpreter','Latex');
title(['relation between depth-dependend dose-averaged $\alpha$-profiles;',' C12, $E_0$=' num2str(machine.data(ix).energy)],'Interpreter','Latex')

string1 = {strsplit(num2str(machine.data(ix).alphaBetaRatio(ixAlphas(2:end))),' ')};
string2 = {strsplit(num2str(machine.data(ix).alphaBetaRatio(ixAlphas(1:end-1))),' ')};

for i = 1:numel(string1{1})
   string3{1}{1,i} = ['$\frac{\alpha_p' string1{1}{1,i} '}{' '\alpha_p' string2{1}{1,i} '}$'] ;
end

legend(string3{1},'Interpreter','latex','FontSize',defaultFontSize)
set(gca,'Fontsize',defaultFontSize)

if FLAG_SAVE && exist('matlab2tikz','file') == 2
    cleanfigure;
    matlab2tikz([latexPath],'height', '12cm', 'width', '21cm','showInfo',showInfoFlag,'standalone', true,...
        'extraaxisoptions',extraAxisAptions);         

    currPath = pwd; cd(exportPath);
    if ispc  
         command = sprintf('pdflatex %s',filename);
    elseif ismac
         command = sprintf('/Library/Tex/texbin/pdflatex %s',[filename]);
    end
    [status,cmdout] = system(command);  delete('*.aux');delete('*.log'); cd(currPath);
    if status > 0
        warning(['couldnt compile pdf: ' cmdout]);
    end
end




%% plot alpha beta ratio relations
filename   = ['LEM_alphaBetaRelation_' machine.meta.machine '.tex'];
latexPath  = [exportPath filesep filename];

figure(),set(gcf,'Color',[1 1 1],'Units', 'normalized', 'Position', [0.1,0.1,0.8,0.8] ),
subplot(221),plot(machine.data(ix).alphaBetaRatio(ixAlphas),machine.data(ix).alphaX(ixAlphas),'LineWidth',2),grid on
xlabel(['$\frac{\alpha_x}{\beta_x}$'],'Interpreter','Latex'),ylabel(['$\alpha_x$'],'Interpreter','Latex')
set(gca,'FontSize',defaultFontSize);

subplot(222),plot(machine.data(ix).alphaBetaRatio(ixAlphas),machine.data(ix).betaX(ixAlphas),'LineWidth',2),grid on
xlabel(['$\frac{\alpha_x}{\beta_x}$'],'Interpreter','Latex'),ylabel(['$\beta_x$'],'Interpreter','Latex')
set(gca,'FontSize',defaultFontSize);

subplot(223),plot(machine.data(ix).depths,machine.data(ix).alpha(:,ixAlphas)./machine.data(ix).alphaX(ixAlphas),'LineWidth',2),grid on
xlabel(['depth'],'Interpreter','Latex'),ylabel(['$\frac{\alpha_p}{\alpha_x}$'],'Interpreter','Latex')
set(gca,'FontSize',defaultFontSize);
legend(string{1},'Interpreter','latex','FontSize',defaultFontSize)
set(gca,'Fontsize',defaultFontSize)

subplot(224),plot(machine.data(ix).depths,machine.data(ix).alpha(:,ixAlphas)./machine.data(ix).betaX(ixAlphas),'LineWidth',2),grid on
xlabel(['depth'],'Interpreter','Latex'),ylabel(['$\frac{\alpha_p}{\beta_x}$'],'Interpreter','Latex')
set(gca,'FontSize',defaultFontSize);
legend(string{1},'Interpreter','latex','FontSize',defaultFontSize)
set(gca,'Fontsize',defaultFontSize)


if FLAG_SAVE && exist('matlab2tikz','file') == 2
    cleanfigure;
    matlab2tikz([latexPath],'height', '16cm', 'width', '21cm','showInfo',showInfoFlag,'standalone', true,...
        'extraaxisoptions',extraAxisAptions);         

    currPath = pwd; cd(exportPath);
    if ispc  
         command = sprintf('pdflatex %s',filename);
    elseif ismac
         command = sprintf('/Library/Tex/texbin/pdflatex %s',[filename]);
    end
    [status,cmdout] = system(command);delete('*.aux');delete('*.log'); cd(currPath);
    if status > 0
        warning(['couldnt compile pdf: ' cmdout]);
    end
end

%%
load('RBE.mat');

ixTargetTissue = 1;
ixNormalTissue = 13;

figure(),set(gcf,'Color',[1 1 1],'Units', 'normalized', 'Position', [0.1,0.1,0.8,0.8] ),
plot(RBE(ixTargetTissue).C.Energy,RBE(ixTargetTissue).C.RBE * RBE(ixTargetTissue).alpha,'LineWidth',2),set(gca,'xScale','log'),grid on, hold 
plot(RBE(ixNormalTissue).C.Energy,RBE(ixNormalTissue).C.RBE * RBE(ixNormalTissue).alpha,':','LineWidth',2)


%%
machineOrg = load('carbon_HIT2.mat');
machineOrg = machineOrg.machine;
load('carbon_HIT_MYLEM.mat');

ix = 198;

ixAlphas = 1:1:10;
numAlpha = numel(ixAlphas);

figure(),set(gcf,'Color',[1 1 1],'Units', 'normalized', 'Position', [0.1,0.1,0.8,0.8] ),
for i = 1:numAlpha
    plot(machineOrg.data(ix).depths,machineOrg.data(ix).alpha(:,ixAlphas(i)),'LineWidth',2,'color',color.cube(i,:));grid on; hold on
    plot(machine.data(ix).depths(1:2:end),machine.data(ix).alpha(1:2:end,ixAlphas(i)),'x','LineWidth',2,'color',color.cube(i,:));
end

xlabel('z[mm]','Interpreter','Latex');
ylabel('$Gy^{-1}$','Interpreter','Latex');
title(['comparison $\alpha$-profiles; HIT vs. myLEM1',' C12, $E_0$=' num2str(machine.data(ix).energy)],'Interpreter','Latex')
set(gca,'FontSize',defaultFontSize)

filename   = ['LEM_alphaProfilesComparison_' machine.meta.machine '.tex'];
latexPath  = [exportPath filesep filename];

if FLAG_SAVE && exist('matlab2tikz','file') == 2
    cleanfigure;
    matlab2tikz([latexPath],'height', '16cm', 'width', '21cm','showInfo',showInfoFlag,'standalone', true,...
        'extraaxisoptions',extraAxisAptions);         

    currPath = pwd; cd(exportPath);
    if ispc  
         command = sprintf('pdflatex %s',filename);
    elseif ismac
         command = sprintf('/Library/Tex/texbin/pdflatex %s',[filename]);
    end
    [status,cmdout] = system(command);delete('*.aux');delete('*.log'); cd(currPath);
    if status > 0
        warning(['couldnt compile pdf: ' cmdout]);
    end
end



%% 
filename  = ['LEM_comparisonRBEini_' machine.meta.machine '.tex'];
latexPath = [exportPath filesep filename];

RBE = load('RBE');
RBE_HIT = RBE.RBE;
load('RBE_LEM1.mat');

part     = 'C';

figure(),set(gcf,'Color',[1 1 1],'Units', 'normalized', 'Position', [0.1,0.1,0.8,0.8] )
ixTissue = 1;
h1 = plot(RBE(ixTissue).(part).Energy,RBE(ixTissue).(part).RBE,'LineWidth',2,'color',color.dkfzdB);hold on
h2 = plot(RBE_HIT(ixTissue).(part).Energy,RBE_HIT(ixTissue).(part).RBE,'x','LineWidth',2,'color',color.dkfzdB); set(gca,'xScale','log'),grid on,grid minor,

ixTissue = 3;
h3 = plot(RBE(ixTissue).(part).Energy,RBE(ixTissue).(part).RBE,'LineWidth',2,'color',color.dre);hold on
h4 = plot(RBE_HIT(ixTissue).(part).Energy,RBE_HIT(ixTissue).(part).RBE,'x','LineWidth',2,'color',color.dre); set(gca,'xScale','log'),grid on,grid minor,

ixTissue = 5;
h5 = plot(RBE(ixTissue).(part).Energy,RBE(ixTissue).(part).RBE,'LineWidth',2,'color',color.gre);,hold on
h6 = plot(RBE_HIT(ixTissue).(part).Energy,RBE_HIT(ixTissue).(part).RBE,'x','LineWidth',2,'color',color.gre); set(gca,'xScale','log'),grid on,grid minor,

legend([h1 h2 h3 h4 h5 h6],{'$\frac{\alpha_x}{\beta_x}=1$','$\frac{\alpha_x}{\beta_x}=1$; HIT',...
                            '$\frac{\alpha_x}{\beta_x}=3$','$\frac{\alpha_x}{\beta_x}=3$; HIT',...
                            '$\frac{\alpha_x}{\beta_x}=5$','$\frac{\alpha_x}{\beta_x}=5$; HIT'},'Interpreter','Latex')

xlabel('energy [MeV]','Interpreter','Latex');
ylabel('$RBE_{ini}$','Interpreter','Latex');
title('$RBE_{ini}$ of carbon particles for different tissu types','Interpreter','Latex')
set(gca,'FontSize',20)

if FLAG_SAVE && exist('matlab2tikz','file') == 2
    cleanfigure;
    matlab2tikz([latexPath],'height', '16cm', 'width', '21cm','showInfo',showInfoFlag,'standalone', true,...
        'extraaxisoptions',extraAxisAptions);         

    currPath = pwd; cd(exportPath);
    if ispc  
         command = sprintf('pdflatex %s',filename);
    elseif ismac
         command = sprintf('/Library/Tex/texbin/pdflatex %s',[filename]);
    end
    [status,cmdout] = system(command);delete('*.aux');delete('*.log'); cd(currPath);
    if status > 0
        warning(['couldnt compile pdf: ' cmdout]);
    end
end



%%





