% figure for MA
% Bragg peak heterogeneity
clear, close all
addpath(genpath('submodules'))
addpath(genpath('phantomAnalysis'))

load('RicData')

energyStep = 112;
pEnergyStep0exp = eval(['pE' num2str(energyStep) '0exp']);
pEnergyStepAexp = eval(['pE' num2str(energyStep) 'Aexp']);
WET_A = 8.9; %[mm]

BPheteroFig = figure;
hold on
plot(pEnergyStep0exp(:,1),       pEnergyStep0exp(:,2),'-', 'DisplayName','pristine')
plot(pEnergyStepAexp(:,1)+WET_A, pEnergyStepAexp(:,2),'--', 'DisplayName','degraded')
% xlim([55 100])
xlabel('depth in water [mm]')
ylabel('dose [a.u.]')
legend('show','location','northwest')

savefig(BPheteroFig,'X:\Masterarbeit\figures\BraggPeakHetero.fig')
matlab2tikz('X:\Masterarbeit\figures\BraggPeakHetero.tex','width','\fwidth')

