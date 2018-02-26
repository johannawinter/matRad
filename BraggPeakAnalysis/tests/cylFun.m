%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function u=cylFun(a,x)
%       ===================================================================
%       Purpose: Compute parabolic cylinder function U(a,x)
%       Input  : a --- Parameter (|a| < 5)
%                x --- Argument (|x| < 5)
%       Output : u ------ U(a,x)                
%       ===================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  E. Cojocaru, January 2009
%  Observations, suggestions and recommendations are welcome at e-mail:
%   ecojocaru@theory.nipne.ro
% Downloaded from http://de.mathworks.com/matlabcentral/fileexchange/63405-bragg-peak-analysis?focused=8040338&tab=function
%   on 01.12.2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c=zeros(1,100);
d=zeros(1,100);
eps=1.0e-15;

c0=1;
c1=a;
c(1)=a;
for k1=4:2:200
    m=k1/2;
    c(m)=a*c1+((k1-2)*(k1-3)*c0)/4;
    c0=c1;
    c1=c(m);
end
y1=1;
r=1;
for k=1:100
    r=0.5*r.*x.^2/(k*(2*k-1));
    r1=c(k)*r;
    y1=y1+r1;
    if(all(abs(r1./y1))<=eps && k>30)
        break
    end
end

d1=1;
d2=a;
d(1)=1.0;
d(2)=a;
for  k2=5:2:160
    m=floor((k2+1)./2);
    d(m)=a*d2+(k2-2)*(k2-3)*d1/4;
    d1=d2;
    d2=d(m);
end
y2=1;
r=1;
for  k=1:100
    r=0.5*r.*x.^2/(k*(2*k+1));
    r1=d(k+1).*r;
    y2=y2+r1;
    if(all(abs(r1./y2))<= eps && k>30)
        break
    end
end
y2=x.*y2;

if a<0 && (a+1/2)==fix(a+1/2)
    ar=pi*(1/4+a/2);
    f1=gamma(1/4-a/2)/(sqrt(pi)*2^(a/2+1/4));
    f2=gamma(3/4-a/2)/(sqrt(pi)*2^(a/2-1/4));
    u=cos(ar)*f1*y1-sin(ar)*f2*y2;
else
    p0=sqrt(pi)/(2^(a/2+1/4));
    g1=gamma(1/4+a/2);
    g3=gamma(3/4+a/2);
    u=p0*(y1/g3-sqrt(2).*y2/g1);
end