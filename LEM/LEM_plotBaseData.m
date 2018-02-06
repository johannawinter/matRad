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
clc,clear,close all
machineOrg = load('carbon_HIT.mat');
machineOrg = machineOrg.machine;
load('carbon_HIT_MYLEM.mat');

ix = 198;

ixAlphas = 1:1:10;
numAlpha = numel(ixAlphas);

figure(),set(gcf,'Color',[1 1 1],'Units', 'normalized', 'Position', [0.1,0.1,0.8,0.8] ),
for i = 1:numAlpha
    plot(machineOrg.data(ix).depths,machineOrg.data(ix).alpha(:,ixAlphas(i)),'LineWidth',2,'color',color.cube(i,:));grid on; hold on
    plot(machine.data(ix).depths,machine.data(ix).alpha(:,ixAlphas(i)),'x','LineWidth',1,'color',color.cube(i,:));
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









