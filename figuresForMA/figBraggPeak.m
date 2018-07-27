% figure for MA
% Bragg peak
clear, close all
addpath(genpath('submodules'))

load('protons_HIT_APMgantry')

energyIx1 = 1;
energyIx2 = 128;
energyIx3 = 255;

depths1 = machine.data(energyIx1).depths;
dose1 = machine.data(energyIx1).Z.doseORG;
depths2 = machine.data(energyIx2).depths;
dose2 = machine.data(energyIx2).Z.doseORG;
depths3 = machine.data(energyIx3).depths;
dose3 = machine.data(energyIx3).Z.doseORG;

BPfig = figure;
hold on
plot(depths1,dose1,'DisplayName',['E_0 = ' num2str(machine.data(energyIx1).energy) ' MeV'])
plot(depths2,dose2,'DisplayName',['E_0 = ' num2str(machine.data(energyIx2).energy) ' MeV'])
plot(depths3,dose3,'DisplayName',['E_0 = ' num2str(machine.data(energyIx3).energy) ' MeV'])
xlabel('depth in water [mm]')
ylabel('dose per 10^6 primaries [Gy]')
legend('show')

savefig(BPfig,'X:\Masterarbeit\figures\BraggPeak.fig')
matlab2tikz('X:\Masterarbeit\figures\BraggPeak.tex','width','\fwidth')
