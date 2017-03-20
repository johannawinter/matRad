function result = matRad_cylGauss(a,zin)
% Implementation according to the one found in T Bortfelds (adapted by JJ
% Wilkens) code for proton dose and LET calculation

%%  Original comments
%	This is a modified version of cyl_gauss. Modified by
%	Jan Wilkens, DKFZ, February 2, 2002
%	It avoids some numerical instabilities and is valid for 
%	0 > a > -3.
%
% Date written:	August 1, 1996
% Revisions:
% Author:       Thomas Bortfeld, DKFZ Heidelberg
% Keywords:     special functions, parabolic cylinder function,
%               degenerate hypergeometric function
% Description:	This is a straight-forward implementation of the product 
%               of the parabolic cylinder function D_a(z) with 
%               exp(-z**2/4), which
%               is valid for real arguments and negative a, 0 > a > -2.
%               It uses a Taylor series for z close to 0, (cf. [1], 
%               sections 9.240 and 9.210),
%               and an asymptotic expansion for larger z
%               (cf. [1], section 9.246). 
% References:	[1] Gradshteyn & Ryzhik; Table of integrals, series and 
%               products; Academic Press; London; 1980.

%% Code

% Set result vector
result = zeros(size(zin));

% Parameters (used by JJ Wilkens):
branch = 3.5; % For a<-2 it might be better to use branch 3.5 to avoid numerical instability
nterm = 10;

% Alternative parameters:
% branch = 4.0;
% nterm = 6.0;

% Check if a is in range
if (-3.0 > a > 0.0)
    disp('a is out of range');
    return;
end

for j=1:numel(zin)
    
    z = zin(j);
    
    if (abs(z) > branch)
        z_square = z.^2;
        nterm_2 = fix(nterm.*branch./abs(z)+0.5);
    
        asympt = 1;
        term_i = 1;
    
        u = a;
    
        for i=1:1:nterm_2
            term_i = term_i .* u .* (1-u) / (2.*i.*z_square);
            asympt = asympt + term_i;
            u = u - 2;
        end
    
        if (z > 0)
            result(j) = exp(-z_square./2).*z.^a.*asympt;
        else
            asympt_2 = 1;
            term_i = 1;
            u = a + 1;
            for i = 1:1:nterm_2
                term_i = term_i * u * (u+1)/(2*i*z_square);
                asympt_2 = asympt_2 + term_i;
                u = u+2;
            end
            
            result(j) =  cos(pi.*a).*exp(-z_square./2) .* (-z).^a .* asympt + sqrt(2.*pi) ./ gamma(-a) .* (-z).^(-a-1.) * asympt_2;
            
        end
        
    else
        c1 = 0.5;
        c2 = 1.5;
        t1 = sqrt(pi)/gamma((1-a)/2)*hypergd(c1+a./2,c1,-z.^2./2);
        t2 = sqrt(2.*pi).*z./gamma(-a/2.)* hypergd(c2+(a-1)./2,c2,-z.^2./2);
        diff = t1-t2;
        result(j) = 2^(a/2).*diff;
    end
end
end


function result = hypergd(a,b,z)
% Implementation according to the one found in T Bortfelds (adapted by JJ
% Wilkens) code for proton dose and LET calculation

%%  Original comments
%
%	This is a modified version of hyperg, working with
%	double precision, as needed in cyl_gaussj.
%	Modified by Jan Wilkens, DKFZ, Feb 8, 2002
%
%
% Date written:	August 1, 1996
% Revisions:
% Author:       Thomas Bortfeld, DKFZ Heidelberg
% Keywords:     special functions, degenerate hypergeometric function
% Description:	Taylor series expansion of the degenerate hypergeometric
%               function Phi(alpha,gamma,z) (cf. [1], section 
%       		9.210), useful for small values of |z|. The value of eps
%               can be adjusted to calculate more or less significant 
%       		digits.
% References:	[1] Gradshteyn & Ryzhik; Table of integrals, series and 
%               products; Academic Press; London; 1980.	 

%% Code

% Parameter
eps =  0.00000000001;

result = 1;
term_i = 1;

for i = 1:1:1000
	  term_i = term_i * (a+i-1)*z/(b+i-1)/i;
	  result = result + term_i;
      
	  if (abs(term_i) < eps*result)
         return
      end
end

disp('No convergence in HYPERG !!!');
return 
end