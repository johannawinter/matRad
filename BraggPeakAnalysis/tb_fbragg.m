function y = tb_fbragg(x,maxdepth,unc)

% data from Bortfeld's range to energy fit R = alpha*E^p
p     = 1.77;
alpha = .0022;

E_0    = (maxdepth/alpha).^(1/p);

% call fbragg(x,maxdepth)
% x = depth in cm - or specify E0 = initial Energy of beam in MeV
% maxdepth = range of proton beam, maxdepth=alpha*E0^p
% epsilon = fraction of primary fluence contribution to the "tail" of the energy spectrum
% sigma_E0 in MeV , energy spectrum of beam
% good beams: gaussian E-spectrum -> epsilon = 0

% ---------- Data used in Bortfeld, med. phys., 1997 ----------------------
% sigma includes both, sigma_E0 and sigma_mono due to straggling. see below (line 57)

% HCL measured data
epsilon= 0.20;
sigma= 0.27;

if nargin > 2
    sigma = sqrt(sigma^2 + unc^2);
end

% TRIUMF_70
%maxdepth= 3.17;
%epsilon= 0.04;
%sigma= 0.06;

% TRIUMF_70_NOZZLE
%maxdepth= 3.17;
%epsilon= 0.09;
%sigma= 0.06;

% TRIUMF_116
%maxdepth= 8.8;
%epsilon= 0.0;
%sigma= 0.13;

% TRIUMF_116_NOZZLE
%maxdepth= 8.8;
%epsilon= 0.13;
%sigma= 0.13;

% ---------- end of example data ------------------------------------------


q          = 1.0;
beta       = .012;
gammaValue = .6;

  
    % sigma^2 = sigma_mono^2 + sigma_E0^2*alpha^2*p^2*E0^(2p-2) (see bortfeld, 1997)
    %si_r=sqrt(si.^2*alpha^2*p^2*E0.^(2*p-2));
    %sigma=sqrt((.012*maxdepth.^.935).^2 + si_r^2);


    %------------------ Main program -----------------------------------
    % original file:bragg_analyt.f
    % Author:	Thomas Bortfeld, DKFZ Heidelberg, August 1, 1996
    % Revisions: 	converted from fortran to matlab, modified to
    % functional form (Dezember 2004, Daniel Pflugfelder)

    % format bp:
    % bp(:,1) -> x-value
    % bp(:,2) -> depth-dose
    % bp(:,3) -> depth-dose due to nuclear reactions

bp(:,1) = x';

v = maxdepth - bp(:,1);

bp(:,2) = sigma^(1/p) * gamma(1/p) / ((2*pi)^.5 * p * alpha^(1/p) * (1+beta*maxdepth)) * ...
           ( cyl_gauss(-1/p,-v./sigma) / sigma + ...
             cyl_gauss(-1/p-q,-v./sigma) * (beta/p + gammaValue*beta + epsilon/maxdepth) );

y = bp(:,2);

% plot bragg peak
%figure
%hold on
%plot(bp(:,1),bp(:,2),'r')
%plot(bp(:,1),bp(:,3),'r--')
%title(['Pristine Bragg peak with range R = ' num2str(maxdepth) ' cm (' num2str(E_0) ' MeV)'])
   


%------------------ functions required to calculate dose ----------

    function y = cyl_gauss(a,z)
    %C Date written:	August 1, 1996
    %C Revisions:	converted to matlab, 2.12.2004 (Daniel Pflugfelder)
    %C Author:	Thomas Bortfeld, DKFZ Heidelberg
    %C Keywords:	special functions, parabolic cylinder function,
    %C		degenerate hypergeometric function
    %C Description:	This is a straight-forward implementation of the product of
    %C		the parabolic cylinder function D_a(z) with exp(-z**2/4),
    %C		which is valid for real arguments and negative a, 0 > a > -2.
    %C		It uses a Taylor series for z close to 0, (cf. [1], sections
    %C		9.240 and 9.210), and an asymptotic expansion for larger z
    %C		(cf. [1], section 9.246). 
    %C References:	[1] Gradshteyn & Ryzhik; Table of integrals, series and 
    %C		products; Academic Press; London; 1980.	 

        branch = 4.0;
        nterm = 6 ;
        
        if ((a>=0)||(a<-2))
            'A is out of range'
        else

            y=zeros(size(z));
            zgb=find(abs(z)>branch);  %indexes for which abs(z(i))>branch

            if length(zgb)
                % use asymptotic expansion
                z_square(zgb) = z(zgb).^2;
                z_square=z_square';
                nterm_2_max = round(max(nterm*branch./abs(z(zgb)))+0.5);

                in_zero=find(z==0);
                z(in_zero)=.000000001; % avoid division by zero.
                nterm_2=round(nterm*branch./abs(z)+0.5);
                z(in_zero)=0;

                asympt = ones(size(z));
                term_i = ones(size(z));
                u = a;

                for i = 1:nterm_2_max
                    term_i(zgb) = term_i(zgb) * u * (1-u) ./ (2*i*z_square(zgb));
                    asympt(zgb,i+1) = asympt(zgb,i) + term_i(zgb);
                    u = u - 2;
                end


                zgbgn=find((z>0) & (abs(z)>branch));	%indexes for which abs(z(i))>branch and z(i)>0 
                v_asym=asympt(:,nterm_2_max+1);
                for k=1:length(zgbgn)
                    v_asym(zgbgn(k))=asympt(zgbgn(k),nterm_2(zgbgn(k))+1);
                end

                y(zgbgn) = exp(-z_square(zgbgn)/2) .* z(zgbgn).^a .* v_asym(zgbgn);%asympt(nterm_2(zgbgn),zgbgn)

                zgbsn=find((z<=0) & (abs(z)>branch));	% indexes for which abs(z(i))>branch and z(i)<=0
                %                                  asymptotic expansion for z < 0
                if ~isempty(zgbsn)
                    clear asympt;
                    asympt=v_asym;

                    asympt_2 = ones(size(z));
                    term_i = ones(size(z));
                    u = a + 1;
                    for i = 1:nterm_2_max
                        term_i(zgbsn) = term_i(zgbsn) * u * (u+1) ./ (2*i*z_square(zgbsn));
                        asympt_2(zgbsn,i+1) = asympt_2(zgbsn,i) + term_i(zgbsn); 
                        u = u + 2;
                    end

                    for k=1:length(zgbsn)
                        v_asym(zgbsn(k))=asympt_2(zgbsn(k),nterm_2(zgbsn(k))+1);
                    end

                    y(zgbsn) = cos(pi*a) .* exp(-z_square(zgbsn)/2) .* (-z(zgbsn)).^a .* asympt(zgbsn) + sqrt(2*pi) / exp(gammaln(-a)) .* (-z(zgbsn)).^(-a-1) .* v_asym(zgbsn);
                end
            end


            zsb=find(abs(z)<=branch);  %indexes for which abs(z(i))<=branch
            %                                   use Taylor series
            if (zsb)
                c1 = 0.5;
                c2 = 1.5;
                term1=ones(size(z));
                term2=ones(size(z));
                term1(zsb) = sqrt(pi)/exp(gammaln((1-a)/2)) .* hypergeom(c1+a/2,c1,-z(zsb).^2/2);
                term2(zsb) = sqrt(2*pi).*z(zsb)/exp(gammaln(-a/2)) .* hypergeom(c2+(a-1)/2,c2,-z(zsb).^2/2);
                %term1, term2
                y(zsb) = 2^(a/2) * (term1(zsb)-term2(zsb));
            end
        end
    end

end



