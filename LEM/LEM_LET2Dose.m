function [Dose_Gy] = LEM_LET2Dose(LET_MeVcm2_g, CntHits, Area_cm2)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LEM_LET2Dose
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2016, Hans-Peter Wieser
%
% h.wieser@dkfz.de
%
% This file is NOT part of the official matRad release. 
% This file has to be used only for internal purposes! 
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function determines the dose in Gy given the LET, the number of hits
% and the area

C = 1.6021766208e-10;
Fluence_cm2 = CntHits/Area_cm2;
Dose_Gy = C * LET_MeVcm2_g * Fluence_cm2;

end

