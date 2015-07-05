function [ baseData ] = extrapDeeper(baseData,visBool)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if isempty(visBool)
    visBool=0;
end

MethodPhysicalData = 'linear';  % 'const'
MethodBioData = 'const';      % 'const'

% add deltaX milimeter
deltaX = 15; 
% number of values which should be added
NumValues2Add = 3;

for i = 1:length(baseData)
    %% extrapolate depth
    vXold = baseData(i).depths;
    OldLength = length(vXold);
    vXdelta = linspace(baseData(i).depths(end),baseData(i).depths(end) + deltaX,NumValues2Add+1);
    vXnew = zeros(length( baseData(i).depths)+length(vXdelta)-1,1);
    vXnew(1:length(vXold))=vXold;
    vXnew(length(vXold)+1:length(vXold)+NumValues2Add)=vXdelta(2:end);
    baseData(i).depths=vXnew;
    %% extrapolate depth depended values per energy
    switch MethodPhysicalData
        case 'linear' 
            baseData(i).Z  = interp1(vXold,baseData(i).Z,vXnew,'linear','extrap');
            baseData(i).sigma1 = interp1(vXold,baseData(i).sigma1,vXnew,'linear','extrap');
            baseData(i).sigma2 = interp1(vXold,baseData(i).sigma2,vXnew,'linear','extrap');
            baseData(i).weight = interp1(vXold,baseData(i).weight,vXnew,'linear','extrap');
        case'const'
            baseData(i).Z(OldLength+1:OldLength+NumValues2Add)=baseData(i).Z(end);
            baseData(i).sigma1(OldLength+1:OldLength+NumValues2Add)=baseData(i).sigma1(end);
            baseData(i).sigma2(OldLength+1:OldLength+NumValues2Add)=baseData(i).sigma2(end);
            baseData(i).weigth(OldLength+1:OldLength+NumValues2Add)=baseData(i).weight(end);
    end
    %% extrapolate biological data
    if isfield(baseData(i),'alpha')
         CntCellLines = size(baseData(i).alpha,2);
        switch MethodBioData
            case 'linear'
                for j=1:CntCellLines
                    baseData(i).alpha(:,j)=interp1(vXold,baseData(i).alpha(:,j),vXnew,'linear','extrap');
                    baseData(i).beta(:,j)= interp1(vXold,baseData(i).beta(:,j),vXnew,'linear','extrap');
                end
            case 'const'
               for j=1:CntCellLines
                    baseData(i).alpha(OldLength+1:OldLength+NumValues2Add,j)= baseData(i).alpha(OldLength,j);
                    baseData(i).beta(OldLength+1:OldLength+NumValues2Add,j)= baseData(i).beta(OldLength,j);
               end
        end
    end
    
    if visBool
        figure,
        subplot(121),plot(baseData(i).depths,baseData(i).Z,'b')
        subplot(122),plot(baseData(i).depths,baseData(i).alpha,'r')
    end
end




end

