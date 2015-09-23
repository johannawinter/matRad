function [Dose_Gy] = LET2Dose(LET_MeVcm2_g, CntHits, Area_cm2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

C = 1.602176487e-13;
Fluence_cm2 = CntHits/Area_cm2;
Dose_Gy = C * LET_MeVcm2_g * Fluence_cm2*1000;

end

