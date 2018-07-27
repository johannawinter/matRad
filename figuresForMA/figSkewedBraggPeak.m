% figure for MA
% skewed Bragg peak
% data from U:\E0408\ImPlan\Porcine lung project\2016-11-27_MainRun\PeakFinderData\LP_2314.mcc
clear, close all
addpath(genpath('submodules'))

depths = [22:.5:58];
dose = [1.2801, 1.2901, 1.3009, 1.3127, 1.3250, 1.3364, 1.3509, 1.3630, 1.3773,...
    1.3916, 1.4066, 1.4228, 1.4400, 1.4559, 1.4748, 1.4950, 1.5151, 1.5356, ...
    1.5590, 1.5823, 1.6088, 1.6378, 1.6670, 1.6997, 1.7333, 1.7707, 1.8113, ...
    1.8562, 1.9018, 1.9536, 2.0101, 2.0647, 2.1204, 2.1674, 2.2090, 2.2387,...
    2.2567, 2.2616, 2.2603, 2.2411, 2.2127, 2.1752, 2.1262, 2.0708, 2.0102,...
    1.9474, 1.8800, 1.8140, 1.7368, 1.6353, 1.5218, 1.3861, 1.2417, 1.0792,...
    0.91504, 0.75681, 0.61483, 0.48828, 0.39463, 0.32032, 0.26655, 0.23248,...
    0.20868, 0.19336, 0.18323, 0.17630, 0.17070, 0.16720, 0.16368, 0.16072,...
    0.15827, 0.15561, 0.15244];
doseShifted = dose - 0.15244;
    
myFig = figure;
hold on
plot(depths,doseShifted)
xlim([22 58])
xlabel('depth in water [mm]')
ylabel('dose [a.u.]')

savefig(myFig,'X:\Masterarbeit\figures\skewedBraggPeak.fig')
matlab2tikz('X:\Masterarbeit\figures\skewedBraggPeak.tex','width','\fwidth')