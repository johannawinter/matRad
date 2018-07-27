% figures for MA
% basic code for different export options
clear, close all
addpath(genpath('submodules'))

myFig = figure;

% different saving options
savefig(myFig,'X:\Masterarbeit\figures\basic.fig')
saveas(myFig,'X:\Masterarbeit\figures\basic','eps')
matlab2tikz('X:\Masterarbeit\figures\basic.tex','width','\fwidth')
matlab2tikz('X:\Masterarbeit\figures\basic.tex','width','\fwidth','height','\fheight')
matlab2tikz('X:\Masterarbeit\figures\basic.tex','standalone',true)
