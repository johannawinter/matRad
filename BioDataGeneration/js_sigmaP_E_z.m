function [data, sigmaSISsq] = js_sigmaP_E_z

data=js_getProton_Data; 
sigmaSISsq = tk_P_sigma_SISsq;          

global perc1      % percentage of integral dose to be cut! 0.3
global perc2      % percentage at which sampling starts! 5
waitb=waitbar(0,'Plase wait calculation in progress');

for i = 1:1:255
    clear wN
    clear s1
    clear s2
    clear depthq
    clear ion
    clear latcut
    
voxelsize=1;
waitbar(i/255,waitb,['Plase wait calculation in progress: step ' num2str(i) ' of 255'] );
energy=data(i).energy;
depth=data(i).depth/10;
FWHM1=data(i).FWHM1;
weight=data(i).weight;
FWHM2=data(i).FWHM2;
ion=data(i).DDD;

latcut=zeros(size(depth));




h = sigmaSISsq( i,1);

wN = abs(weight);
s11 = 1/(2*sqrt(2*log(2)))*abs(FWHM1);
s22 = 1/(2*sqrt(2*log(2)))*abs(FWHM2);
s1(:,1) = sqrt( h + s11(:,1).^2);
s2(:,1) = sqrt( h + s22(:,1).^2);

a=round(5*max(s2));
x=-a:voxelsize:a;
z=-a:voxelsize:a;
N=numel(x);

% endcal=max(find(34<depth ,1,'first'),find(depth<=37,1,'last'));
% if size(endcal,1)<1
    endcal=size(depth,1);
% end
endcal=min(endcal,size(depth,1));
endp=min(depth(endcal),38);

[~,o]=max(ion .* ( (1-wN(:,1))./s1(:,1).^2 + wN(:,1)./s2(:,1).^2 ) );

Dmax = (1-wN(o,1))  ./ (s1(o,1).^2 ) + wN(o,1)  ./ (s2(o,1).^2 ) ;
Dmaxi = @(io) (Dmax * io - 65500).^2;
io=fminsearch(Dmaxi,1000000);
CalFac = max(ion)*data(i).absCalFac;

ion(:)=round(io.*ion(:)/max(ion));
data(i).DDD = ion;
data(i).absCalFac = CalFac ./ io;




depthq=0:voxelsize/10:endp;
L2 = zeros(endcal,N,N);
D = zeros(endcal,N,N);

     
% Double Gaussian

for k=1:endcal
for j=1:N
    L2(k,:,j) = (1-wN(k,1)) .* exp(-.5 .* ((x )./s1(k,1)) .^ 2) .* exp(-.5 .*((z(j) )./s1(k,1)) .^ 2) ./ (s1(k,1).^2 ) ...
                  + wN(k,1) .* exp(-.5 .* ((x )./s2(k,1)) .^ 2) .* exp(-.5 * ((z(j) )./s2(k,1)) .^ 2) ./ (s2(k,1).^2 ) ;
end
    D(k,:,:)=L2(k,:,:)*ion(k,1);
end   



D2compl=D(:,:,round(max(size(z))/2));
D2compl(D2compl<1)=0;
D = interp1(depth(1:endcal),D,depthq,'linear');
D(D<1)=0;


              Dtemp = squeeze(D(:));
              F1 = @(Thres1) (sum(Dtemp(Dtemp<=Thres1))/sum(Dtemp) - perc1/100).^2;
              F2 = @(Thres) (sum(Dtemp(Dtemp<=Thres))/sum(Dtemp) - perc2/100).^2;
              Thresh = fminsearch(F1,20);
              Thresh2 = fminsearch(F2,200);
   data(i).DoseCutValue=Thresh;
   data(i).DoseSampleValue=Thresh2;
   Dcut=D;
   Dcut(Dcut<Thresh)=0;
   Dcut=Dcut(:,:,round(max(size(z))/2));
   [~,temp2]=find(Dcut>0,1,'last');
%    D=D(:,:,round(max(size(z))/2));
%    [~,temp2]=find(D>0,1,'last');

D2compl(D2compl<Thresh)=0;
for k=1:endcal
    D2=D2compl(k,:);
    [~,temp1]=find(D2>0,1,'last');
    if sum(temp1)==0
        temp1=round(size(x,2)/2);
    end
    latcut(k,1)=abs(x(max(temp1)))+3;
end
   
%    data(i).LateralCutValue=abs(x(max(temp2))+3);
   data(i).LateralCutValue=abs(x(max(temp2)))+3;
   data(i).LatCutVal=latcut;
   

end
delete(waitb)
