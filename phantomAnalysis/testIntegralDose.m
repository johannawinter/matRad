function [differenceDose2Dose1,differenceDose3Dose1,differenceDose3Dose2] = ...
    testIntegralDose(resultGUI,dose1,dose2,dose3,dose4)
% function to test conservation of integral dose in resultGUI of 
%   dose1 and dose2, and optionally dose3.

% integrate dose over complete dose cube
integralDose1 = sum(sum(sum(resultGUI.(dose1))));
integralDose2 = sum(sum(sum(resultGUI.(dose2))));
try    
    integralDose3 = sum(sum(sum(resultGUI.(dose3))));
catch
end
try    
    integralDose4 = sum(sum(sum(resultGUI.(dose4))));
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
fprintf(['Integral difference ' dose2 ' - ' dose1 ': ' num2str(differenceDose2Dose1*100,2) '%%. \n'])
try
    fprintf(['Integral difference ' dose3 ' - ' dose1 ': ' num2str(differenceDose3Dose1*100,2) '%%. \n'])
catch
end
try
    fprintf(['Integral difference ' dose4 ' - ' dose1 ': ' num2str(differenceDose4Dose1*100,2) '%%. \n'])
catch
end

