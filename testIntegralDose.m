function [differenceDose2Dose1,differenceDose3Dose1,differenceDose4Dose1] = ...
    testIntegralDose(resultGUI,doseCube)

% Test conservation of integral dose of doseCubes in resultGUI.
%   doseCubes must be given in cell arrays, e.g., doseCube{1} = resultGUI.RBExDose


% define dose cubes for comparison if not specified
if ~exist('doseCube','var') || isempty(doseCube)
    doseCube{1} = resultGUI.RBExDose;
    doseCube{2} = resultGUI.matRadRecalc_RBExDose;
    doseCube{3} = resultGUI.matRadHetero_RBExDose;
    
    doseCubeName{1} = 'original RBExDose';
    doseCubeName{2} = 'matRad-recalculated RBExDose';
    doseCubeName{3} = 'heterogeneity RBExDose';
else
    doseCubeName = ' ';
end

% integrate dose over complete dose cube
integralDose1 = sum(sum(sum(doseCube{1})));
integralDose2 = sum(sum(sum(doseCube{2})));
try    
    integralDose3 = sum(sum(sum(doseCube{3})));
catch
end
try    
    integralDose4 = sum(sum(sum(doseCube{4})));
catch
end

% calculate differences
differenceDose2Dose1 = (integralDose2-integralDose1)/integralDose1;
try
    differenceDose3Dose1 = (integralDose3-integralDose1)/integralDose1;
catch
end
try
    differenceDose4Dose1 = (integralDose4-integralDose1)/integralDose1;
catch
end

% print differences
fprintf(['Integral difference ' doseCubeName{2} ' - ' doseCubeName{1} ': ' ...
    num2str(differenceDose2Dose1*100,2) '%%. \n'])
try
    fprintf(['Integral difference ' doseCubeName{3} ' - ' doseCubeName{1} ': ' ...
        num2str(differenceDose3Dose1*100,2) '%%. \n'])
catch
end
try
    fprintf(['Integral difference ' doseCubeName{4} ' - ' doseCubeName{1} ': ' ...
        num2str(differenceDose4Dose1*100,2) '%%. \n'])
catch
end

