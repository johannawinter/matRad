function plotIndividualGaussComponents(baseData,ellSq,heteroCorrType)
% function to compare convolution (complete and depth based) of APM 
%   Gauss curves and of Bragg peak (sum of them) by plotting.
%   Function should be called from within matRad_calcParticleDoseBixel to
%   get base data, ellSq and heteroCorrType.

% plot pristine Gauss curves
if strcmp(heteroCorrType,'complete')
    GaussFig = figure;
end
hold on
title(['Inidivdual Gauss curves contribution to BP, convolution type: ' ...
    heteroCorrType ', energy ' num2str(baseData.energy) ' MeV'])
xlabel('rad. depth [mm]')
for i = 1:10
    GaussPlot(i,:) = baseData.Z.weight(i)/(baseData.Z.width(i)*sqrt(2*pi)).* ...
        exp(-(baseData.depths-baseData.Z.mean(i)).^2 / (2*baseData.Z.width(i).^2));
    GaussLine(i) = plot(baseData.depths, GaussPlot(i,:), 'b');
end
% add pristine Bragg peak
sumLine = plot(baseData.depths, sum(GaussPlot,1),'g');

% plot convoluted Gauss curves
for i = 1:10
    GaussConvPlot(i,:) = baseData.Z.weight(i)/(sqrt(ellSq(1,i))*sqrt(2*pi)).* ...
        exp(-(baseData.depths-baseData.Z.mean(i)).^2 / (2*sqrt(ellSq(1,i)).^2));
    if strcmp(heteroCorrType,'complete')
        GaussConvLine(i) = plot(baseData.depths, GaussConvPlot(i,:), 'r');
    elseif strcmp(heteroCorrType,'depthBased')
        GaussConvLine(i) = plot(baseData.depths, GaussConvPlot(i,:), 'c');
    end
end
% add convoluted Bragg peak
if strcmp(heteroCorrType,'complete')
    sumConvLine = plot(baseData.depths, sum(GaussConvPlot,1), '--r');
elseif strcmp(heteroCorrType,'depthBased')
    sumConvLine = plot(baseData.depths, sum(GaussConvPlot,1), '--c');
end

% add lung boundaries
lungBoundary1 = 30 - baseData.offset;
lungBoundary2 = lungBoundary1 + 50*.297;
lungLine(1) = plot([lungBoundary1 lungBoundary1],[0 55],'-.k');
lungLine(2) = plot([lungBoundary2 lungBoundary2],[0 55],'-.k');

% add legend
if strcmp(heteroCorrType,'complete')
    legend([GaussLine(1),sumLine,GaussConvLine(1),sumConvLine,lungLine(1)], ...
        'pristine Gauss curves','pristine BP',...
        ['convoluted Gauss curves (' heteroCorrType ')'],['convoluted BP (' heteroCorrType ')'],...
        'lung boundaries','location','northwest')
end


% savefig(GaussFig, ...
%     ['C:\Matlab\Analysis phantom degradation\implementationComparison\Gauss_' heteroCorrType])



