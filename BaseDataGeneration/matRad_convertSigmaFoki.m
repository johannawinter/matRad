clc
clear
close all
%load('protons_Generic.mat');
load('carbon_Generic.mat');

for i = 1:length(machine.data)
    
    SigmaIni = machine.data(i).sigma(1);
    
    for j = 1:length(machine.data(i).sigma)
        machine.data(i).sigma(j) = sqrt(machine.data(i).sigma(j)^2 - SigmaIni^2);
    end
    
    machine.data(i).initFocus(1).dist  = [0 machine.meta.SAD];
    machine.data(i).initFocus(1).sigma = [SigmaIni SigmaIni];
    
end

