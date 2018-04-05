function [Sigma_SIS,vEnergySIS] = matRad_getSigmaSIS(Identifier,basePath,FocusIdx)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad_getSigmaSIS script 
% This script can be used to parse all foci or one specific foci of a sis file. 
% Please note that it returns SIGMA and not FWHM. Values in the SIS file
% determine the beam width at iso-center position in FWHM.
%
% call
%   [Sigma_SIS,vEnergySIS] = matRad_getSigmaSIS(Identifier,basePath,FocusIdx)
%
% input
%   Identifier:      machine base data file
%   basePath:        path to TRiP folder
%   FocusIdx:        focus index  (optional)
%
% output
%   machine:         machine base data file including an additional field representing the lateral
%                    cutoff
%   
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
%     ASCI image to explain relevant distances for beam widing                                                              .:irB2LuuLvi, 
%                                                                                                                     .,bsghrt<c      hcghhfgbf
%                                                                                                                   b v                       grst
%                                                                                                                 ds                              ds  
%                                                                                                               iG@                                u@   
% :                                                                                                           .Su.                                  22  
% OB1                                                                                                        .@                                      @  
%   F@S                                                                                                      LO                                      @  
%     UBG                                                                                                     @            PATIENT                   u0  
%       L@8.                                                                                                  jB                                    Bv  
%         7@B.                                                                                                 B                                   Y@   
%  SOURCE   :@BN0OZ8G8G8GOG8Z8G8GOGOGOZ8G8Z8ZOG8ZOG8GOZ8GOG8GOGOG8G8ZOG8Z8GOENBk                               @:                                 v@    
%            OU                                                               ZM                             7M2                                  @     
%    LSu.    BU                                                               08                           uBu                .GB@j               B:    
%  .@B@B@L   B5                                                               GZ XX_______________________@5_________________d@@B@@d              uO    
%  FB@B@B@   @1                      BEAM APPLICATION SYSTEM                  0G XX                      @:                  rB@B@B@               @.   
%   @B@@@,   B1                                                               EZ                        LO                    uB@B@                :@F  
%     8:     BL                                                               SG                        i@                      k   ISO CENTER       OP 
%     F      B@i777777v7v777777777777777777777777777v777v7777777777777777777ri@N                         Bu                     P                    @i 
%     0    i@Oi                                                               ZL                        .1iXP.                  G                    @  
%     N  7@M,                                                                 ;7                        ,k  :B                  O                   5M  
%     5FB8.                                                                   7L                        ,N   SM                 G                   B   
%    u@1                                                                      Lv                        ,N    @                 O                  Bi   
% ,XB7X                                                                       vL                        .0    i@                G                iB7    
% @7  U                       <-   SSD = 5383 mm + X mm  ->                   vv                         U     7MUULi.          O             .7MM.     
%     M:                                                                      rr                        v@       ..:72qFi.      P          iu0Nv        
%   ,B@BuvYLYLYLYLJLYLJLJYJYJLJLJLJLJLJYYLYYJLJLYLYLJLYLYLJLJLYLJLJLYLJLYLYLLrOG7vJLJYJLJYYLJLJLJLYLJLL1@B@7j          ,71qS2L;.M,...:iv2X0j:           
%    .Bv                                                                      7v                        :i                  .:irB2LuuLvi,               
%     F                                                                       v7                                                P                       
%     E.                                                                      vv                                                E                       
%     8.                                                                      7r                                                P                       
%     G.                                                                      ir             <-  BAMStoISO 1126 mm ->           E                       
%     8.                                                                     Z@B2:rrrrr;r;r;r;r;r;r;rrrrr;rrr;rrrrrrr;rrrrrrrii@@@:ir               
%     G.                                                                     PB@,.                                              PBB          
%     O.                                                                                                                        5                       
%     Z.                                                                                                                        0                       
%     G                                                                                                                         Z                       
%     2                                                <- SAD 6509 mm ->                                                        5                       
%    ;@S                                                                                                                       ,@L                      
%   ;@B@vvvLvLvLvLvLvLvLLLvLvLLLvLLLvLvLvLvLLLvLvLLLvLvLvLvLvLvLvLvLvLvLvLLLvLLLvLvLvLLLLLvLvLLLLLvLvLvLvLvLvLLLvLvLvLvLvLLLvv;@B@LL.                   
%     Z.                                                                                                                        @,                      
%     1                                                                                                                                                 
%     

%#ok<*AGROW>

if nargin < 3
    FocusIdx = 1:7;
end

switch Identifier
    case {'c','C'}
        FullFilePath = [basePath filesep 'SIS' filesep '12C_1.6.2008.sis'];
    case {'p','P','h''H'}
        FullFilePath = [basePath filesep 'SIS' filesep '1H_1.1.2009.sis'];
    case {'o','O'}
        FullFilePath = [basePath filesep 'SIS' filesep '16O_20.12.2010'];
    otherwise
        error('unkown particle type')
end

%open file handle
fHandle = fopen(FullFilePath,'r');

if fHandle < 0
    display('Could not open SIS file');
end
    
Cnt=1;
currentline = fgetl(fHandle);

while ischar(currentline)
    % each data row starts with 'energy'
    if(strfind(currentline,'energy'))
        vEnergySIS(Cnt,1)  = str2double(currentline(8:13));     
        vFoci = cell2mat(textscan(currentline,'energy %*f focus %f %f %f %f %f %f %f')); 
        FWHM(Cnt,:) = vFoci(FocusIdx);
        Cnt = Cnt+1;
    end
    currentline = fgetl(fHandle);       
end

Sigma_SIS = (FWHM./(2*sqrt(2*log(2))));
%close file handle
fclose(fHandle);

end




