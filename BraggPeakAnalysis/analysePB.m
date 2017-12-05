function results = analysePB(z, dose, varargin)
% ANALYSEBP the function performs analysis of Bragg Peak Curve giving the
% parameters of the depth dose distribution. Three methods of analysis are
% possible using spline fit, 3rd polynomial fit to the peak [1] and bortfeld
% fit [2]. All parameters of Bragg curve are calculated using normalisation 
% to maximum of spline, poly3 fit or bortfeld fit for spline, poly3 and
% bortfeld analysis method respectively. 
%
%   results = ANALYSEBP(z,dose) analyse a single bragg curve using all
%   possible methods (spline, poly3 and bortfeld). Inputs are:
%       z - vector of depth positions
%       dose - value of dose at depth (it is required that
%       lenght(dose)=length(z))
%   As a result the function returns a structure with fields: spline, poly3
%   and bortfeld containing following parameters of analysis. Common
%   information is:
%       D100    - dose at maximum
%       R100    - position of maximum
%       R90P    - position of 90% proximal build-up
%       R90D    - position of 90% distal fall-off
%       R80P    - position of 80% proximal build-up
%       R80D    - position of 80% distal fall-off
%       R20     - position of 20% distal fall-off
%       R10     - position of 10% distal fall-off
%       R50P    - position of 50% proximal build-up
%       R50D    - position of 50% distal fall-off
%       FWHM    - full width at half maximum
%       DFO2080 - distal fall-off 20%-80%
%       DFO1090 - distal fall-off 10%-90%
%   Additional information are:
%       spline:
%           zSpline - vector of depth positions after spline
%           doseSpline - vector of dose at depth after spline
%           doseSplineNorm - vector of dose at depth after spline normalised
%                            to maximum
%       poly3:
%           fitResults - fitobject of 3rd polynomial fit
%           GOF - goodness-of-fit
%           z -  vector of depth positions (copy of input)
%           dose - vector of dose at depth (copy of input)
%           doseNorm - vector of dose at depth normalised to maximum of fit
%       bortfeld:
%           z -  vector of depth positions (copy of input) [mm]
%           dose - vector of dose at depth (copy of input)
%           doseNorm - vector of dose at depth normalised to maximum of fit
%           fitResultBraggPeak - fitobject of bortfeld fit to peak part
%           GOFBraggPeak - goodness-of-fit to peak part
%           fitResultPlateau - fitobject of bortfeld fit to plateau part
%           GOFPlateau - goodness-of-fit to plateau part
%           E0 - initial energy calculated based on bortfeld fit [MeV]
%           sigmaMono - range uncertainty resulting from natural spread
%                       calculated based on bortfeld fit [mm]
%           sigmaE0 - energy spread calculated based on bortfeld fit [MeV]
%
%  results = ANALYSEBP(z,dose,PARAM1,VAL1,...) additional parameters can be
%  as follows: 
%       'method' - method of analysis, it takes values (def. 'full'):
%             'full' - full analysis using all possible methods
%             'spline' - use spline method of analysis
%             'poly3' - perform 3rd polynomial fit to peak [1]
%             'bortfeld' - perform bortfeld fit [2]
%       'AccuracySpline' - step for spline analisis (def. 10% of minimum
%                          step in z) 
%       'AccuracyBortfeld' - accuracy of bortfeld analisis (def. 0.001 [mm]) 
%       'levelPoly' - level above which the 3rd polynomial will be fitted.
%                     (def. 0.8 meaning 80%) 
%   
%   NOTE:
%   For borthelfd analysis the user must assure that the z parameter is in
%   [mm] unit. 
 
%   [1] Parodi, K., et al. (2012). "Monte Carlo simulations to support start-up 
%   and treatment planning of scanned proton and carbon ion therapy at a 
%   synchrotron-based facility", Physics in Medicine and Biology, 57(12), 3759–84
%
%   [2] Bortfeld, T. (1997). "An analytical approximation of the Bragg curve
%   for therapeutic proton beams", Medical Physics, 24(12), 2024. 
%
%   See also fit, 
%
% Created by Jan Gajewski (jan.gajewski@ifj.edu.pl) on 14.06.2017
% Modified by Jan Gajewski (jan.gajewski@ifj.edu.pl) on 26.07.2017
% Prepared in Matlab R2016a
% Downloaded from http://de.mathworks.com/matlabcentral/fileexchange/63405-bragg-peak-analysis?focused=8040338&tab=function
%   on 29.11.2017


%% parser for inputs
% initialize parser
p=inputParser;
p.CaseSensitive=false;
p.FunctionName='analysePB';
% required inputs
addRequired(p,'z', @(x) validateattributes(x,{'numeric'},{'vector'})); % is numeric and is 1D array
addRequired(p,'dose', @(x) validateattributes(x,{'numeric'},{'vector'})); % is numeric and is 1D array
% optional inputs
% parameter inputs
addParameter(p,'method','full',@(x) any(validatestring(x,{'full','poly3','bortfeld','spline'}))) % choose method of analysis
addParameter(p,'AccuracySpline',0.1*min(diff(z)),@(x) validateattributes(x,{'numeric'},{'scalar'})); % accuracy of spline analysis
addParameter(p,'AccuracyBortfeld',0.001,@(x) validateattributes(x,{'numeric'},{'scalar'})); % accuracy of bortfeld parameters calculation (assuming the input is in mm)
addParameter(p,'levelPoly',0.8,@(x) validateattributes(x,{'numeric'},{'scalar'})); % level above which the polynominals 3rd is fitted
% parse inputs
parse(p, z, dose, varargin{:});
% get parameters values
AccuracySpline=p.Results.AccuracySpline;
AccuracyBortfeld=p.Results.AccuracyBortfeld;
levelPoly=p.Results.levelPoly;
method=validatestring(p.Results.method,{'full','poly3','bortfeld','spline'});

%% perform spline analysis 
if any(strcmpi(method,{'spline','full'}))
    % Perform spline
    results.spline.zSpline=min(z):AccuracySpline:max(z);
    results.spline.doseSpline=spline(z,dose, results.spline.zSpline);
    results.spline.doseSplineNorm=results.spline.doseSpline/max(results.spline.doseSpline);
    % Calculate parameters with normalisation to maximum of spline
    results.spline.D100=max(results.spline.doseSpline);
    results.spline.R100=results.spline.zSpline(results.spline.doseSpline==results.spline.D100);
    results.spline.R90P=getLineIntValue(results.spline.zSpline,results.spline.doseSplineNorm,0.9,'prox');
    results.spline.R90D=getLineIntValue(results.spline.zSpline,results.spline.doseSplineNorm,0.9,'dist');
    results.spline.R80P=getLineIntValue(results.spline.zSpline,results.spline.doseSplineNorm,0.8,'prox');
    results.spline.R80D=getLineIntValue(results.spline.zSpline,results.spline.doseSplineNorm,0.8,'dist');
    results.spline.R20=getLineIntValue(results.spline.zSpline,results.spline.doseSplineNorm,0.2,'dist');
    results.spline.R10=getLineIntValue(results.spline.zSpline,results.spline.doseSplineNorm,0.1,'dist');
    results.spline.R50P=getLineIntValue(results.spline.zSpline,results.spline.doseSplineNorm,0.5,'prox');
    results.spline.R50D=getLineIntValue(results.spline.zSpline,results.spline.doseSplineNorm,0.5,'dist');
    results.spline.FWHM=results.spline.R50D-results.spline.R50P;
    results.spline.DFO2080=results.spline.R20-results.spline.R80D;
    results.spline.DFO1090=results.spline.R10-results.spline.R90D;
end
%% perform fit of 3rd polynominal
if any(strcmpi(method,{'poly3','full'}))
    % Fit options
    [results.poly3.z,results.poly3.dose]=prepareCurveData(z,dose);
    poly3FitType=fittype('poly3');
    poly3FitOpt=fitoptions(poly3FitType);
    poly3FitOpt.Robust = 'Bisquare';
    poly3FitOpt.Lower = [-Inf 0 -Inf 0];
    poly3FitOpt.Upper = [0 Inf 0 Inf];
    poly3FitOpt.Normalize = 'off';
    poly3FitOpt.Exclude = (dose <=(max(dose)*levelPoly));
    % 3rd polynominal fit
    [results.poly3.fitResult,results.poly3.GOF]=fit(results.poly3.z,results.poly3.dose, poly3FitType, poly3FitOpt);
    % Calculate parameters with normalisation to maximum of fit
    coefficients=coeffvalues(results.poly3.fitResult)';
    results.poly3.R100=(-2*coefficients(2)-sqrt(4*coefficients(2)^2-(12*coefficients(1)*coefficients(3))))/(6*coefficients(1));
    results.poly3.D100=results.poly3.fitResult(results.poly3.R100);
    results.poly3.doseNorm=results.poly3.dose/results.poly3.D100;
    results.poly3.R90P=getLineIntValue(results.poly3.z,results.poly3.doseNorm,0.9,'prox');
    results.poly3.R90D=getLineIntValue(results.poly3.z,results.poly3.doseNorm,0.9,'dist');
    results.poly3.R80P=getLineIntValue(results.poly3.z,results.poly3.doseNorm,0.8,'prox');
    results.poly3.R80D=getLineIntValue(results.poly3.z,results.poly3.doseNorm,0.8,'dist');
    results.poly3.R20=getLineIntValue(results.poly3.z,results.poly3.doseNorm,0.2,'dist');
    results.poly3.R10=getLineIntValue(results.poly3.z,results.poly3.doseNorm,0.1,'dist');
    results.poly3.R50P=getLineIntValue(results.poly3.z,results.poly3.doseNorm,0.5,'prox');
    results.poly3.R50D=getLineIntValue(results.poly3.z,results.poly3.doseNorm,0.5,'dist');
    results.poly3.FWHM=results.poly3.R50D-results.poly3.R50P;
    results.poly3.DFO2080=results.poly3.R20-results.poly3.R80D;
    results.poly3.DFO1090=results.poly3.R10-results.poly3.R90D;
    results.poly3.z=results.poly3.z(results.poly3.dose>=(max(dose)*levelPoly));
    results.poly3.dose=results.poly3.fitResult(results.poly3.z);
    results.poly3.doseNorm=results.poly3.dose/results.poly3.D100;

end

%% perform bortfeld fit
if any(strcmpi(method,{'bortfeld','full'}))
    % Find starting parameters and bounds
    bortParam.p=1.77; % [-]
    bortParam.alpha=0.022; % [mm*MeV^(-p)]
    bortParam.zSpline=min(z):AccuracySpline:max(z);
    bortParam.doseSpline=spline(z,dose, bortParam.zSpline);
    bortParam.doseSplineNorm=bortParam.doseSpline/max(bortParam.doseSpline);
    bortParam.D100=max(bortParam.doseSpline);
    bortParam.R0=getLineIntValue(bortParam.zSpline,bortParam.doseSplineNorm,0.8,'dist'); % [mm] range of protons (approx. 80% of distal fall-off)
    bortParam.E0=(bortParam.R0/bortParam.alpha)^(1/bortParam.p); % [MeV] initial energy
    bortParam.sigmaMono=(0.012*bortParam.R0^0.935)/10; % [mm] width of Gaussian range straggling
    bortParam.sigmaE0=0.01*bortParam.E0; % [MeV] width of Gaussian energy spectrum
    bortParam.sigma=sqrt(bortParam.sigmaMono^2+bortParam.sigmaE0^2*bortParam.alpha^2*bortParam.p^2*bortParam.E0^(2*bortParam.p-2)); % [mm] proton range disparsion
    bortParam.epsilon=0.1; % [-] Fraction of primary fluence contribut- ing to the "tail" of the energy spectrum
    bortParam.phi=0.1; % 1/mm^2 initial fluence (normalisation factor)
    
    bortParam.Exclude=bortParam.R0+[-3 5]*bortParam.sigma;
    bortParam.dR0=bortParam.R0+[-0.3 0.3]*bortParam.sigma;
    bortParam.dsigma=bortParam.sigma*[0.5 3];
    bortParam.depsilon=[-10 10];
    bortParam.dphi=[0 Inf];
    % Fit options BraggPeak
    [results.bortfeld.z,results.bortfeld.dose] = prepareCurveData(z,dose);
    bortfeldFitTypeBraggPeak=fittype('phi*((exp(-((R0-z).^2)/(4*sigma^2))*sigma^0.565)/(1+0.012*R0)).*((11.26/sigma).*cylFun(0.565-0.5,(-(R0-z)./sigma))+(0.157+11.26*(epsilon/R0)).*cylFun(1.565-0.5,(-(R0-z)./sigma)))',...
                                     'independent','z','coefficients',{'phi','R0','epsilon','sigma'});
    bortfeldFitOptBraggPeak=fitoptions(bortfeldFitTypeBraggPeak);
    bortfeldFitOptBraggPeak.Lower=[bortParam.dphi(1) bortParam.dR0(1) bortParam.depsilon(1) bortParam.dsigma(1)];
    bortfeldFitOptBraggPeak.StartPoint=[bortParam.phi bortParam.R0 bortParam.epsilon bortParam.sigma];
    bortfeldFitOptBraggPeak.Upper=[bortParam.dphi(2) bortParam.dR0(2) bortParam.depsilon(2) bortParam.dsigma(2)];
    bortfeldFitOptBraggPeak.Exclude=(results.bortfeld.z<bortParam.Exclude(1)) | (results.bortfeld.z>bortParam.Exclude(2));
    bortfeldFitTypePlateau=fittype('(phi/(1+0.012*R0))*(17.93*(R0-z).^-0.435+(0.444+31.7*(epsilon/R0))*(R0-z).^0.565)',...
                                     'independent','z','coefficients',{'phi','R0','epsilon'});
    bortfeldFitOptPlateau=fitoptions(bortfeldFitTypePlateau);
    bortfeldFitOptPlateau.Lower=[bortParam.dphi(1) bortParam.dR0(1) bortParam.depsilon(1)];
    bortfeldFitOptPlateau.StartPoint=[bortParam.phi bortParam.R0 bortParam.epsilon];
    bortfeldFitOptPlateau.Upper=[bortParam.dphi(2) bortParam.dR0(2) bortParam.depsilon(2)];
    bortfeldFitOptPlateau.Exclude=(results.bortfeld.z<bortParam.dR0(1)*0.7 | results.bortfeld.z>=bortParam.Exclude(1));
    
    % Bortfeld fit
    [results.bortfeld.fitResultBraggPeak,results.bortfeld.GOFBraggPeak]= fit(results.bortfeld.z,results.bortfeld.dose, bortfeldFitTypeBraggPeak, bortfeldFitOptBraggPeak);
    [results.bortfeld.fitResultPlateau,results.bortfeld.GOFPlateau]= fit(results.bortfeld.z,results.bortfeld.dose, bortfeldFitTypePlateau, bortfeldFitOptPlateau);
    
    % Calculate parameters with normalisaation to maximum of bortfeld fit
    results.bortfeld.z=(bortParam.dR0(1)*0.7):AccuracyBortfeld:bortParam.Exclude(2);
    results.bortfeld.dose=[results.bortfeld.fitResultPlateau(results.bortfeld.z(results.bortfeld.z<bortParam.Exclude(1)))*(results.bortfeld.fitResultBraggPeak(bortParam.Exclude(1))/results.bortfeld.fitResultPlateau(bortParam.Exclude(1)));...
                           results.bortfeld.fitResultBraggPeak(results.bortfeld.z(results.bortfeld.z>=bortParam.Exclude(1)))]';
    results.bortfeld.doseNorm=results.bortfeld.dose/max(results.bortfeld.dose);
    results.bortfeld.D100=max(results.bortfeld.dose);
    results.bortfeld.R100=results.bortfeld.z(results.bortfeld.dose==results.bortfeld.D100);
    results.bortfeld.R90P=getLineIntValue(results.bortfeld.z,results.bortfeld.doseNorm,0.9,'prox');
    results.bortfeld.R90D=getLineIntValue(results.bortfeld.z,results.bortfeld.doseNorm,0.9,'dist');
    results.bortfeld.R80P=getLineIntValue(results.bortfeld.z,results.bortfeld.doseNorm,0.8,'prox');
    results.bortfeld.R80D=getLineIntValue(results.bortfeld.z,results.bortfeld.doseNorm,0.8,'dist');
    results.bortfeld.R20=getLineIntValue(results.bortfeld.z,results.bortfeld.doseNorm,0.2,'dist');
    results.bortfeld.R10=getLineIntValue(results.bortfeld.z,results.bortfeld.doseNorm,0.1,'dist');
    results.bortfeld.R50P=getLineIntValue(results.bortfeld.z,results.bortfeld.doseNorm,0.5,'prox');
    results.bortfeld.R50D=getLineIntValue(results.bortfeld.z,results.bortfeld.doseNorm,0.5,'dist');
    results.bortfeld.FWHM=results.bortfeld.R50D-results.bortfeld.R50P;
    results.bortfeld.DFO2080=results.bortfeld.R20-results.bortfeld.R80D;
    results.bortfeld.DFO1090=results.bortfeld.R10-results.bortfeld.R90D;  
    results.bortfeld.E0=(results.bortfeld.fitResultBraggPeak.R0/bortParam.alpha)^(1/bortParam.p);
    results.bortfeld.sigmaMono=(0.012*results.bortfeld.fitResultBraggPeak.R0^0.935);
    results.bortfeld.sigmaE0=sqrt((results.bortfeld.fitResultBraggPeak.sigma^2-results.bortfeld.sigmaMono^2)/(bortParam.alpha^2*bortParam.p^2*results.bortfeld.E0^(2*bortParam.p-2)));
end

function zVal=getLineIntValue(z,doseNorm,level,site)
    switch site
        case 'prox'
            z=z(find(doseNorm>=level,1,'first')+[-1 0]);
            doseNorm=doseNorm(find(doseNorm>=level,1,'first')+[-1 0]);
            lineInt=polyfit(z, doseNorm, 1);
            zVal=(level-lineInt(2))/lineInt(1);
        case 'dist'
            z=z(find(doseNorm>=level,1,'last')+[0 1]);
            doseNorm=doseNorm(find(doseNorm>=level,1,'last')+[0 1]);
            lineInt=polyfit(z, doseNorm, 1);
            zVal=(level-lineInt(2))/lineInt(1);
    end
end

%% DEBUG
% bortfeldFitOpt
% results.bortfeld.fitResultBraggPeak
% results.bortfeld.GOFBraggPeak.rsquare
% results.bortfeld.fitResultPlateau
% results.bortfeld.GOFPlateau.rsquare
% results.bortfeld

% plot(z,dose/results.bortfeld.D100,'.-b')
% hold on
% plot([1 1]*results.bortfeld.R50P,[0 1],'-g')
% plot([1 1]*bortParam.Exclude(1),[0 0.7],'-g')
% plot(results.bortfeld.z,results.bortfeld.doseNorm,'.-r');
% hold off
% xlim([0 bortParam.Exclude(2)])
% ylim([0 1.1])
% grid on

% plot(z,dose/results.poly3.D100,'.-b')
% hold on
% plot([1 1]*results.poly3.R50P,[0 1],'-g')
% hold off
% xlim([0 300])
% ylim([0 1.1])
% grid on

end

