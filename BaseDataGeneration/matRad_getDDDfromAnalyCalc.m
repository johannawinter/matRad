function [ddd] = matRad_getDDDfromAnalyCalc(Identifier,R0, ZSec, ASec, vDepth)

switch Identifier
    case {'p','P','h''H'}
           getDose = @(R0,depth) getProtonDose(R0,vDepth);
    case {'c','C'}
           getDose = @(R0,depth) getCarbonDose(R0,vDepth);
    otherwise
        error('unkown particle type')
end


ddd = getDose(R0,vDepth);

end


function [DDD] = getProtonDose(R0,depth)
% Caclulate the depth dose curve for a proton beam with given R0 and depth
% values.
% Based on T. Bortfeld code
% Version 0.1
p = 1.77;        % no units
%alpha = 0.0022* 3^(1-p)/2^2;  % cm/MeV^(-p)
alpha = 0.0022;
beta  = 0.012;    % 1/cm
rho   = 1.0;       % density of water g/cm^3
%dosefrac = 1.0;  % ???
RT       = 0.0002;     % ???
epsilon  = 0.2;   % no tail in the initial energy spectrum
vargamma = 0.6;      % Necessary?

% Specifiy special parameters for this case
%R0 = 14.64;      % cm (0.8 * max dose)
%peak = 17.47499; % cm (max dose)
sigmaE = 0.8;    % energy invariance of machine

% the original file specifies a filename at this point - not required here

%% Calculation of dose, LET_d and LET_t

% Calculate sigma
sigmaMono = 0.012 * R0^(0.935);  % ???
E0 = (R0/alpha)^(1/p);          % MeV (initial energy)
sigma = sqrt( sigmaMono^2+(sigmaE*alpha*p*E0^(p-1))^2); % ???

% Some more precalculations
z = depth;
zeta = (z-R0)./sigma;   % benannt wie im Paper
xi = (z-RT-R0)./sigma;  % ???

% Calculate physical DOSE
dose = sigma.^(1/p).*gamma(1/p).*(cyl_gaussj(-1./p,zeta)./sigma ...
    + cyl_gaussj(-(1/p)-1,zeta).*(beta/p+vargamma*beta+epsilon/R0) )/ ...
	(sqrt(2*pi)*rho*p*alpha^(1/p)*(1+beta*R0));

	
%	calculate LETd
DDD.LETd_RT = ( 0.1 .* (sigma^(2./p).*gamma(2./p) ...
    .* ( cyl_gaussj(-2./p,xi) - cyl_gaussj(-2./p,zeta) ) ...  
    - 2 .* (RT/2.) .^ (2./p) .* exp(-(zeta+xi).^2./8)) ...
	./ (p^2.*alpha^(1./p).*(2./p-1.).*(sigma^(1./p+1.).*gamma(1./p+1) ...
    * ( cyl_gaussj(-1./p-1.,xi) -  cyl_gaussj(-1./p-1.,zeta) )  ...
    - 2.* (RT/2.).^(1./p+1.).*exp(-(zeta+xi).^2./8))));


%DDD.dose = dose./max(dose);
DDD.dose = dose;
 
end

function [DDD] = getCarbonDose(R0,depth)
        error('not implemented');
end
