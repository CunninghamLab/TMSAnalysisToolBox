function [finIndex,finVal,time,PL,Icoarse,PLcoarse] = prolific(data,SamplingRate,fCutLow,fCutHigh,distr,gridcoarse,varargin)
% PROLIFIC - PROfile LIkelihood based on FIbonaCci search
% Find the onset of an EMG signal by Fibonacci maximization of the Profile Likelihood
%
% INPUT:
% data         Vector containing the EMG signal
% SamplingRate Sampling rate (observations per second).
%              Mandatory argument, no default
% fCutLow      Low-pass frequency of low-pass filter (fCutHigh < 0 or undefined) or
%              lower band-pass frequency of band-bass filter (fCutHigh > 0)
%              Unit: Hz, 0 <= fCutLow < SamplingRate/2
%              Default: 0
% fCutHigh     High-pass frequency of high-pass filter (fCutLow < 0 or undefined) or
%              high band-pass frequency in band-bass filter (fCutLow >= 0)
%              Unit: Hz, 0 < fCutHigh <= SamplingRate/2
%              Default: 0
%              If both fcutLow and fCutHigh are missing no filter is applied
% distr        Probability distribution of the Profile Likelihood
%              Possible options: 'Weibull' (default), 'Gauss', 'Laplace',
%              'Lognormal', 'Exponential', 'Cauchy', 'Gamma', 'Logistic',
%              'Extreme value', 'Birnbaum-Sanders', 'Burr'
%              distr is a cell array of two elements (distribution at left
%              and right of the onset). If only one distribution is provided,
%              it is intended to be the same at left and right of the
%              onset.
% gridcoarse   Width of the uniform initial coarse grid
%              Unit: ms
%              Default: 150 ms
% varargin     Optional arguments in the form of the couple ('option', value),
%              where possible options are:
%               'verbose' Amount of printouts:
%                         0: No printout (default)
%                         1: Moderate printouts (initial and final solution)
%                         2: Extensive printouts (include iterations)
%              'plotflag' Plot produced:
%                         0: no plot produced (default)
%                         >= 1: Mixed plot of PL and TKEO-conditioned signal
%                         >= 2: Also plots of PL, signal and TKEO-conditioned signal
%                'unit_t' Time unit:
%                         'ms'  millisecond (default)
%                         's'  second
%                         'index'  index of the array
%                'unit_s' Unit of the input EMG signal:
%                         'mV'  milliVolt (default)
%                         'V'   Volt
%             'maintitle' Main title of the plots (string)
%              'saveplot' Whether to save the plots on files eventually
%                         0: Do not save plots (figures not closed)
%                         1: Save plots (figures closed, default)
%            'folderplot' Subfolder where to save plots eventually
%                         (default current directory)
%            'formatplot' Format of the saved plots eventually
%                         'pdf': PDF (default)
%                         'eps': EPS
%                         'jpg': JPEG
%                         'png': PNG
%                'fullPL' Compute PL on the full grid:
%                         0: PL is not computed (default)
%                         1: PL is computed
%
% OUTPUT:
% finIndex: Index of the Onset vector
%   finVal: Value of the Profile Likelihood at the estimated Onset
%     time: Elapsed time (s)
%       PL: Profile likelihood on the whole grid (if plotflag>0 or fullPL=1)
%
% EXAMPLES
% All default values
% x = randn(5000,1) + [repmat(0,1500,1); repmat(3,2500,1); repmat(0,1000,1)];
% [Onset, Value] = prolific(x,2048);
%
% Specify distribution and width of the coarse grid
% [Onset, Value] = prolific(x,2048,[],[],'Gauss',150);
%
% Extensive output
% [Onset, Value] = prolific(x,2048,[],[],'Gauss',150,'verbose',2,'plotflag',2,'saveplot',0)
%
% Save output plot
% [Onset, Value] = prolific(x,2048,[],[],'Gauss',150,'verbose',2,'plotflag',1,'saveplot',1,'formatplot','pdf')
% 
% References:
% S.E. Selvan, D. Allexandre, U. Amato, B. Della Vecchia, G.H. Yue:
% A Fast and Robust Profile-Likelihood-Based Muscle Onset Detection in EMG
% Using Discrete Fibonacci Search. To be published in IEEE Access
%
% S.E. Selvan, D. Allexandre, U. Amato, G.H. Yue (2018): Unsupervised Stochastic
% Strategies for Robust Detection of Muscle Activation Onsets in Surface
% Electromyogram. IEEE Trans. Neur. Syst. Rehabil. Engin. 26(6), 1279-1291
%
% Requires MATLAB Statistics and Machine Learning Toolbox
% Uses parts of code eLattice by John L. Weatherwax
% http://waxworksmath.com/Authors/N_Z/Wilde/Code/eLattice.m
%
% Version: 0.99.1 - June 2020
%
% Copyright: S.E. Selvan, D. Allexandre, U. Amato, B. Della Vecchia, G.H. Yue
%

% Default values of the input arguments
fCutLow_def = 0;
fCutHigh_def = 0;
gridcoarse_def = 150;
distr_def = 'Weibull';

verbose_def = 0;
plotflag_def = 0;
unit_t_def = 'ms';
unit_s_def = 'mV';
maintitle_def = '';
saveplot_def = 1;
folderplot_def = 'plot';
formatplot_def = 'pdf';
fullPL_def = 0;

fontsize = 14;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse variable input arguments of the function %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Count input arguments
nArgs = length(varargin);
if round(nArgs/2) ~= nArgs/2
    error('PROLIFIC requires an even number of parameters')
end

% Check and process fixed arguments
if ~exist('data','var') || isempty(data)
    error('PROLIFIC: data not provided and no default. Exit program')
end
if ~isnumeric(data)
    error('PROLIFIC: data not numeric. Exit program')
end
if length(data) <= 1
    error('PROLIFIC: data not array. Exit progran')
end
if ~exist('SamplingRate','var') || isempty(SamplingRate)
    error('PROLIFIC: SamplingRate non provided and no default. Exit program')
end
if ~exist('fCutLow','var') || isempty(fCutLow), fCutLow = fCutLow_def; end
if ~isnumeric(fCutLow)
    warning('PROLIFIC: fCutLow not numeric. Setting default')
    fCutLow = fCutLow_def;
end
if length(fCutLow) > 1
    warning('PROLIFIC: fCutLow not a scalar. Setting default')
    fCutLow = fCutLow_def;
end
if fCutLow < 0 || fCutLow > SamplingRate/2
    warning('PROLIFIC: fCutLow out of the range [0 - SamplingRate/2]. Setting default')
end
if ~exist('fCutHigh','var') || isempty(fCutHigh), fCutHigh = fCutHigh_def; end
if ~isnumeric(fCutHigh)
    warning('PROLIFIC: fCutHigh not numeric. Setting default')
    fCutHigh = fCutHigh_def;
end
if length(fCutHigh) > 1
    warning('PROLIFIC: fCutHigh not a scalar. Setting default')
    fCutHigh = fCutHigh_def;
end
if fCutHigh < 0 || fCutHigh > SamplingRate/2
    warning('PROLIFIC: fCutHigh out of the range [0 - SamplingRate/2]. Setting default')
end
if ~exist('distr','var') || isempty(distr), distr = distr_def; end
if ~exist('gridcoarse','var') || isempty(gridcoarse), gridcoarse = gridcoarse_def; end
if ~isnumeric(gridcoarse)
    warning('PROLIFIC: gridcoarse not numeric. Setting default')
    gridcoarse = gridcoarse_def;
end
if length(gridcoarse) > 1
    warning('PROLIFIC: gridcoarse not a scalar. Setting default')
    gridcoarse = gridcoarse_def;
end
if length(data) <= 1
    error('PROLIFIC: data not array. Exit progran')
end
% Make sure distr is cell array with 2 distributions (for Left and Right)
if iscell(distr)
    if length(distr) == 1
        if ~ischar(distr{1})
            error('PROLIFIC: distr not a string. Exit program')
        end
        distr = [distr distr];
    end
    if length(distr) >= 2
        if ~ischar(distr{1}) || ~ischar(distr{2})
            error('PROLIFIC: distr not a string. Exit program')
        end
        distr = distr(1:2);
    end
else
    if ~ischar(distr)
        error('PROLIFIC: distr not a string. Exit program')
    end
    distr = {distr distr};
end
if ~any(strcmpi({'Weibull','Gauss','Laplace','Lognormal','Exponential',...
        'Cauchy', 'Gamma', 'Logistic','Extreme value','Birnbaum-Sanders',...
        'Burr'},distr{1}))
    warning('PROLIFIC: distr not an allowed distribution - setting default')
    distr{1} = distr_def;
end
if ~any(strcmpi({'Weibull','Gauss','Laplace','Lognormal','Exponential',...
        'Cauchy', 'Gamma', 'Logistic','Extreme value','Birnbaum-Sanders',...
        'Burr'},distr{2}))
    warning('PROLIFIC: distr not an allowed distribution - setting default')
    distr{2} = distr_def;
end

if ~isnumeric(SamplingRate)
    error('PROLIFIC: SamplingRate not numeric. Exit program')
end
if SamplingRate <= 0
    error('PROLIFIC: SamplingRate negative. Exit program')
end
if ~isnumeric(fCutHigh)
    error('PROLIFIC: fCutLow not numeric. Exit program')
end
if ~isnumeric(fCutHigh)
    error('PROLIFIC: fCutHigh not numeric. Exit program')
end
if fCutHigh < 0 || fCutHigh > SamplingRate/2
    error('PROLIFIC: fCutHigh out of the range [0 - SamplingRate/2]. Exit program')
end

% Parse variable input arguments
for n = 1:2:nArgs
    arg = lower(varargin{n});
    switch arg
        case 'verbose', verbose = varargin{n+1};
        case 'plotflag', plotflag = varargin{n+1};
        case 'unit_t', unit_t = varargin{n+1};
        case 'unit_s', unit_s = varargin{n+1};
        case 'maintitle', maintitle = varargin{n+1};
        case 'saveplot', saveplot = varargin{n+1};
        case 'folderplot', folderplot = varargin{n+1};
        case 'formatplot', formatplot = varargin{n+1};
        case 'fullpl', fullPL = varargin{n+1};
        case 'gridcoarse', gridcoarse = varargin{n+1};
        otherwise, warning('%s is not a recognized argument name',varargin{n})
    end
end

% Check if variable input arguments exist
if ~exist('verbose','var'), verbose = verbose_def; end
if ~exist('plotflag','var'), plotflag = plotflag_def; end
if ~exist('unit_t','var'), unit_t = unit_t_def; end
if ~exist('unit_s','var'), unit_s = unit_s_def; end
if ~exist('maintitle','var'), maintitle = maintitle_def; end
if ~exist('saveplot','var'), saveplot = saveplot_def; end
if ~exist('folderplot','var'), folderplot = folderplot_def; end
if ~exist('formatplot','var'), formatplot = formatplot_def; end
if ~exist('fullPL','var'), fullPL = fullPL_def; end
if ~exist('gridcoarse','var'), gridcoarse = gridcoarse_def; end

% Check and process variable input arguments

% Check and process verbose
if ~isnumeric(verbose)
    warning(['PROLIFIC: value for verbose (' verbose ') not recognized - setting default'])
    verbose = verbose_def;
end
verbose = round(verbose);
if verbose < 0, verbose = verbose_def; end

% Check and process plotflag
if ~isnumeric(plotflag)
    warning(['PROLIFIC: value for plotflag (' plotflag ') not recognized - setting default'])
    plotflag = plotflag_def;
end
plotflag = round(plotflag);
if plotflag < 0, plotflag = plotflag_def; end
if plotflag > 0
    % Check and process unit_t
    % c: label of the time unit for the plots
    if ~ischar(unit_t)
        warning(['PROLIFIC: value for unit_t (' unit_t ') not recognized - setting default'])
        unit_t = unit_t_def;
    end
    if ~any(strcmpi({'i','s','m'},unit_t(1:1)))
        warning(['PROLIFIC: value for unit_t (' unit_t ') not recognized - setting default'])
        unit_t = unit_t_def;
    end
    switch lower(unit_t(1:1))
        case 'i', c = 1; xlab = '(Index)';
        case 's', c = 1/SamplingRate; xlab = '(sec)';
        case 'm', c = 1000/SamplingRate; xlab = '(ms)';
        otherwise
            warning(['PROLIFIC: value for unit_t (' unit_t ') not recognized - setting default'])
            xlab = '(ms)';
    end
    xlab = ['EMG signal duration ' xlab];
    
    % Check and process unit_s
    % Label of the signal for the plots
    if ~ischar(unit_s)
        warning(['PROLIFIC: value for unit_s (' unit_s ') not recognized - setting default'])
        unit_s = unit_s_def;
    end
    if ~any(strcmpi({'v','m'},unit_s(1:1)))
        warning(['PROLIFIC: value for unit_s (' unit_t ') not recognized - setting default'])
        unit_s = unit_s_def;
    end
    switch lower(unit_s(1:1))
        case 'v', ylab = '(V)';
        case 'm', ylab = '(mV)';
        otherwise
            warning(['PROLIFIC: value for unit_s (' unit_s ') not recognized - setting default'])
            ylab = '(mV)';
    end
    
    % Check and process maintitle
    if ~ischar(maintitle)
        warning(['PROLIFIC: value for maintitle (' num2str(maintitle) ') not recognized - setting default'])
        maintitle = maintitle_def;
    end
    form3 = maintitle(~isspace(maintitle));
    
    % Check and process saveplot
    if ~isnumeric(saveplot)
        warning(['PROLIFIC: value for saveplot (' saveplot ') not recognized - setting default'])
        saveplot = saveplot_def;
    end
    saveplot = round(saveplot);
    if saveplot < 0 || saveplot > 1, saveplot = saveplot_def; end
    
    if saveplot
        % Check and process folderplot
        if ~ischar(folderplot)
            warning(['PROLIFIC: value for folderplot (' num2str(folderplot) ') not recognized - setting default'])
            folderplot = folderplot_def;
        end
        if ~exist(folderplot,'dir')
            warning(['PROLIFIC: folder (' folderplot ') not existing - creating'])
            mkdir(folderplot);
        end
        % Check and process formatplot
        if ~ischar(formatplot)
            warning(['PROLIFIC: value for formatplot (' num2str(formatplot) ') not recognized - setting default'])
            formatplot = formatplot_def;
            
        end
        if ~any(strcmpi({'pdf','jpeg','jpg','eps','png'},formatplot))
            warning(['PROLIFIC: value for formatplot (' formatplot ') not admittable - setting default'])
            formatplot = formatplot_def;
        end
        switch lower(formatplot)
            case {'jpg','jpeg'}
                form1 = '-djpeg100'; form2 = '.jpg';
            case 'pdf'
                form1 = '-dpdf'; form2 = '.pdf';
            case 'eps'
                form1 = '-depsc'; form2 = '.eps';
            case 'png'
                form1 = '-dpng'; form2 = '.png';
        end
    end
    
    % col: Colors for the plots
    col=[0         0.4470    0.7410
        0.8500    0.3250    0.0980
        1.0000    0.0000    1.0000
        0.9290    0.6940    0.1250
        0.4660    0.6740    0.1880
        0.3010    0.7450    0.9330
        0.6350    0.0780    0.1840];
end

tstart = tic;

%%%%%%%%%%%%%%%%%%%%
% FIBONACCI search %
%%%%%%%%%%%%%%%%%%%%

% TKEO-preprocess the signal
EMG = EnergyTKEO(data,SamplingRate,fCutLow,fCutHigh);
N = length(data);

% The FIBONACCI NUMBERS:
F_NUMBERS = zeros(100,1);
F_NUMBERS(1) = 1; F_NUMBERS(2) = 1;
i = 2;
while F_NUMBERS(i) <= N
    i = i+1;
    F_NUMBERS(i) = F_NUMBERS(i-1) + F_NUMBERS(i-2);
end
F_NUMBERS(i+1:100) = [];

% Set the coarse uniform grid
Ncoarse = round(N/SamplingRate/gridcoarse*1000);
Icoarse = round(linspace(1,N,Ncoarse))';
% Icoarse=1 failes because it is the first data; set to half grid
Icoarse(1) = round(Icoarse(2)/2);
% Analogously in the case of right boundary
Icoarse(Ncoarse) = round(0.5*(Icoarse(Ncoarse)+Icoarse(Ncoarse-1)));

% STEP 1: Estimate PL on the coarse uniform grid
PLcoarse = zeros(Ncoarse,1);
if verbose >= 2
    fprintf('*** FIBONACCI ONSET SEARCH - distribution')
    fprintf(' %s-%s ',distr{1},distr{2})
    fprintf('\nStep 1 - Search on a coarse unform grid of length %i:\n',Ncoarse)
end
% Compute PL on the uniform grid
for ii = 1:Ncoarse
    PLcoarse(ii) = ProfLik(Icoarse(ii),distr,0);
end
% Find all maxima of PL
indmax = find(PLcoarse(2:Ncoarse-1) > PLcoarse(1:Ncoarse-2) & ...
    PLcoarse(2:Ncoarse-1) > PLcoarse(3:Ncoarse)) + 1;
% Eventually add the first and/or last point if larger than the nearest ones
if PLcoarse(Ncoarse) > PLcoarse(Ncoarse-1), indmax = [indmax; Ncoarse]; end
if PLcoarse(1) > PLcoarse(2), indmax = [1; indmax]; end
% Estimate points to the left and right of the first maximum (if within the sample)
PLindex1 = Icoarse(max(indmax(1)-1,1));
PLindex2 = Icoarse(min(indmax(1)+1,Ncoarse));
if verbose > 2, fprintf('PLindex1 - PLindex2 %i %i\n',PLindex1,PLindex2), end

time = toc(tstart);

if verbose >= 2
    fprintf('Found %i maxima: ',length(indmax))
    for ii = 1:length(indmax)
        fprintf('%i (%i) ',indmax(ii),Icoarse(indmax(ii)))
    end
    fprintf('\n')
    fprintf('Selected subregion: %i-%i\n',PLindex1,PLindex2)
end

if fullPL > 0 || plotflag > 0
    % Compute PL on the full grid for plotting
    % Allow at least 5 points for estimating parameters
    PL = NaN(N,1);
    for ii = 5:N-5
        PL(ii)  =ProfLik(ii,distr);
    end
else
    PL = [];
end

if plotflag > 1
    fig1 = figure('NumberTitle', 'off', 'Name', 'EMG signal'); clf
    plot((1:N)*c,data)
    axis tight
    xlabel(xlab)
    ylabel(['EMG Signal ' ylab])
    title(maintitle,'Interpreter', 'none')
    set(gca,'fontsize', fontsize)
    if saveplot > 0
        print(form1,fullfile(folderplot,['plot_emg_' form3 form2]))
        close(fig1.Number)
    end
    
    fig2 = figure('NumberTitle', 'off', 'Name', 'TKEO-conditioned EMG amplitude'); clf
    plot((1:N)*c,EMG)
    axis tight
    xlabel(xlab)
    ylabel('TKEO-conditioned EMG amplitude')
    title(maintitle,'Interpreter', 'none')
    set(gca,'fontsize', fontsize)
    if saveplot > 0
        print(form1,fullfile(folderplot,['plot_TKEO_' form3 form2]))
        close(fig2.Number)
    end
    
    fig3 = figure('NumberTitle', 'off', 'Name', 'PL'); clf
    h1 = zeros(5,1);
    h1(1) = plot((1:N)*c,PL,'Color',[0.50 0.50 0.50],'LineWidth',2);
    hold on
    h1(2) = plot(Icoarse*c,PLcoarse,'o','MarkerFaceColor',col(1,:),'MarkerEdgeColor',col(1,:));
    h1(3) = plot(Icoarse(max(indmax(1)-1,1):min(indmax(1)+1,Ncoarse))*c, ...
        PLcoarse(max(indmax(1)-1,1):min(indmax(1)+1,Ncoarse)), ...
        '-','LineWidth',2,'Color',col(2,:));
end

if plotflag > 0
    fig4 = figure('NumberTitle', 'off', 'Name', 'PL and Signal'); clf
    h = zeros(5+plotflag>1,1);
    set(fig4,'defaultAxesColorOrder',[0 0 0]);
    % Double y axes handled manually instead of yyaxis because yyaxis
    % mandatorily puts the left axis in back. Here define ax2 before ax1 so
    % that ax2 is put in back
    ax2 = axes; ax1 = axes;
    line(ax2,(1:N)*c,EMG,'Color',col(5,:));
    % Plot again the first point just to make a ticker line in the legend
    h(6) = line(ax2,c,EMG(1),'Color',col(5,:),'LineWidth',2);
    cc=3; % Reduction factor for the right plot
    tk1=min(EMG); tk2=max(EMG);
    dtk = tk2 - tk1;
    y2min = min(tk1,-dtk*cc*0.1); % 0.075
    y2max = y2min + dtk*cc;
    ax2.YLim = [y2min,y2max];
    ax2.XLim = [1 N]*c;
    ax2.FontSize = fontsize;
    ylabel(ax2,'TKEO-conditioned EMG amplitude')
    h(1) = plot(ax1,(1:N)*c,PL,'Color',[0.50 0.50 0.50],'LineWidth',2);
    hold on
    h(2) = plot(ax1,Icoarse*c,PLcoarse,'o','MarkerFaceColor',col(1,:),'MarkerEdgeColor',col(1,:));
    h(3) = plot(ax1,Icoarse(max(indmax(1)-1,1):min(indmax(1)+1,Ncoarse))*c, ...
        PLcoarse(max(indmax(1)-1,1):min(indmax(1)+1,Ncoarse)), ...
        '-','LineWidth',2,'Color',col(2,:));
    plot(ax1,[c,N*c],[ax1.YLim(2),ax1.YLim(2)],'k')
end

tstart = tic;

% Find the largest Fibonacci number than the size of the subregion
fIndex = find(F_NUMBERS >= PLindex2-PLindex1+1,1,'first');
% xArray starts from the first index of the subregion
xArray = PLindex1+(0:F_NUMBERS(fIndex))-1;
PLindex2 = xArray(end);
% Fix initial right boundary of the Fibonacci search
lBIndex = 0;
hBIndex = F_NUMBERS(fIndex);

%------------------------------------------------------------------------
% EVALUATE STARTING TWO POINTS OF SEARCH ( "X1" and "X2" ):
%------------------------------------------------------------------------

testLBIndex = F_NUMBERS(fIndex-2);
testLB      = xArray(testLBIndex);
lBVal       = ProfLik(testLB, distr);

time = time + toc(tstart);

if plotflag > 1, figure(fig3.Number), h1(4) = plot(testLB*c,lBVal,'*','Color',col(3,:)); end
if plotflag > 0, figure(fig4.Number), h(4) = plot(ax1,testLB*c,lBVal,'*','Color',col(3,:)); end

tstart = tic;

testHBIndex = F_NUMBERS( fIndex-1 );
testHB      = xArray(testHBIndex);
hBVal       = ProfLik(testHB, distr);

time = time + toc(tstart);

if plotflag > 1, figure(fig3.Number), plot(testHB*c,hBVal,'*','Color',col(3,:)), end
if plotflag > 0, figure(fig4.Number), plot(ax1,testHB*c,hBVal,'*','Color',col(3,:)), end

%------------------------------------------------------------------------
% START LATTICE FIBONACCI SEARCH
%------------------------------------------------------------------------
if verbose >= 2, fprintf('Step 2 - Fine search:\nit:  0 - a: %4i - S1: %4i (%f) - S2: %4i (%f) - b: %4i\n', ...
        lBIndex,testLBIndex,lBVal,testHBIndex,hBVal,hBIndex), end

tstart = tic;

for ii = 1:fIndex-4
    if verbose>2, fprintf('%i - lBIndex %i, testLBIndex %i,lBVal %f ,testHBIndex %i, hBVal %f,hBIndex %i', ...
            ii,lBIndex,testLBIndex,lBVal,testHBIndex,hBVal,hBIndex), end
    if( lBVal < hBVal )
        %'caso lBVal < hBval'
        lBIndex = testLBIndex;
        testIndex = lBIndex + (hBIndex - testHBIndex);
        testPt    = xArray(testIndex);
        if verbose > 2, fprintf('testIndex %i - testPt %i\n',testIndex,testPt), end
        if( testIndex < testHBIndex )
            %'caso testIndex < testHBIndex'
            testLBIndex = testIndex;
            testLB      = testPt;
            if verbose > 2, fprintf('chiama Profilik per testLB %i\n',teslLB), end
            lBVal       = ProfLik(testLB, distr);
            
            time = time + toc(tstart);
            
            if plotflag > 1, figure(fig3.Number), plot(testLB*c,lBVal,'*','Color',col(3,:)), end
            if plotflag > 0, figure(fig4.Number), plot(ax1,testLB*c,lBVal,'*','Color',col(3,:)), end
            
            tstart = tic;
        else
            % 'caso testIndex !< testHBIndex'
            testLBIndex = testHBIndex;
            testLB      = testHB;
            lBVal       = hBVal;
            testHBIndex = testIndex;
            testHB      = testPt;
            if verbose > 2, fprintf('chiama Profilik per testHB %i\n',testHB), end
            hBVal       = ProfLik(testHB, distr);
            
            time = time + toc(tstart);
            
            if plotflag > 1, figure(fig3.Number), plot(testHB*c,hBVal,'*','Color',col(3,:)), end
            if plotflag > 0, figure(fig4.Number), plot(ax1,testHB*c,hBVal,'*','Color',col(3,:)), end
            
            tstart = tic;
        end
    else
        %'caso lBVal !< hBval'
        hBIndex = testHBIndex;
        testIndex = hBIndex - (testLBIndex - lBIndex);
        testPt = xArray( testIndex );
        if( testIndex > testLBIndex )
            %'caso testIndex > testLBIndex'
            testHBIndex = testIndex;
            testHB      = testPt;
            hBVal       = ProfLik(testHB, distr);
        else
            %'caso testIndex !> testLBIndex'
            testHBIndex = testLBIndex;
            testHB      = testLB;
            hBVal       = lBVal;
            testLBIndex = testIndex;
            testLB      = testPt;
            lBVal       = ProfLik(testLB, distr);
            
            time = time + toc(tstart);
            
            if plotflag > 1, figure(fig3.Number), plot(testLB*c,lBVal,'*','Color',col(3,:)), end
            if plotflag > 0, figure(fig4.Number), plot(ax1,testLB*c,lBVal,'*','Color',col(3,:)), end
            
            tstart = tic;
            
        end
    end
    
    if verbose >= 2, fprintf('it: %2i - a: %4i - S1: %4i (%f) - S2: %4i (%f) - b: %4i\n', ...
            ii,lBIndex,testLBIndex,lBVal,testHBIndex,hBVal,hBIndex), end
end

if( lBVal < hBVal )
    finIndex = testHB;
    finVal = hBVal;
else
    finIndex = testLB;
    finVal = lBVal;
end

time = time + toc(tstart);

if verbose >= 1, fprintf('PROLIFIC maximum: %i (%f)\n',...
        finIndex,finVal),
end

if plotflag > 1
    figure(fig3.Number)
    h1(5) = plot(finIndex*c,finVal,'o','MarkerFaceColor',col(4,:),'MarkerEdgeColor','k','MarkerSize',10);
    xlabel(xlab,'FontSize',fontsize)
    ylabel('Profile Likelihood','FontSize',fontsize)
    axis tight
    legend(h1,{'Profile Likelihood','Coarse grid PL','Selected range','Trial points','Solution'}, ...
        'Location','Best','FontSize',fontsize)
    title(maintitle,'Interpreter', 'none')
    set(gca,'fontsize', fontsize)
    if saveplot
        print(form1,fullfile(folderplot,['plot_PL_' form3 form2]))
        close(fig3.Number)
    end
end
if plotflag > 0
    figure(fig4.Number)
    h(5) = plot(ax1,finIndex*c,finVal,'o','MarkerFaceColor',col(4,:),'MarkerEdgeColor','k','MarkerSize',10);
    ylabel(ax1,'Profile likelihood','FontSize',fontsize)
    xlabel(ax1,xlab,'FontSize',fontsize)
    ax1.XLim = [1 N]*c;
    ax1.FontSize = fontsize;
    % legend(ax2,h,{'Exhaustive search','Coarse grid search','Selected range','DFS trial points','Solution','TKEO output'},'Location','Best')
    legend(ax2,h,{'Profile likelihood','Coarse grid PL','Selected range','Trial points','Solution','TKEO amplitude'}, ...
        'Location','Best','FontSize',fontsize)
    
    if length(indmax) > 1
        lastindmax = Icoarse(indmax(end));
    else
        lastindmax = N;
    end
    arrowY = -0.02*dtk*cc; % Y coordinate in the ax2 axis
    arrowYlab = -0.040*dtk*cc; % Y coordinate in the ax2 axis
    % Plot labels in the axis ax1 instead of ax2 to make them appear front
    % Therefore change Y coordinate
    tmp = ax1.YLim; y1min = tmp(1); y1max = tmp(2);
    % Compute arrowYlab in ax1 coordinates
    arrowYlab = y1min+(arrowYlab-y2min)/(y2max-y2min)*(y1max-y1min);
    h=annotation('arrow');
    set(h,'parent', ax2, 'position', [finIndex*c arrowY (lastindmax-finIndex)*c 0])
    h=annotation('arrow');
    set(h,'parent', ax2, 'position', [lastindmax*c arrowY -(lastindmax-finIndex)*c 0])
    h=annotation('arrow');
    set(h,'parent', ax2, 'position', [0 arrowY finIndex*c 0])
    h=annotation('arrow');
    set(h,'parent', ax2, 'position', [finIndex*c arrowY -finIndex*c 0])
    text(ax1,finIndex/2*c,arrowYlab,'Baseline','HorizontalAlignment','center', ...
        'FontSize',fontsize)
    text(ax1,(finIndex+(lastindmax-finIndex)/2)*c,arrowYlab,'Muscular activity', ...
        'HorizontalAlignment','center','FontSize',fontsize)
    
    title(maintitle,'Interpreter', 'none','FontSize',fontsize)
    
    ax2.XTick = []; % Plot ticks only for axis ax1
    ax2.XLabel = []; % Plot ticks only for axis ax1
    ax1.Box = 'off'; % Avoid mixing ax1 ticks on the axis ax2
    ax2.Box = 'off'; % Force filling the whole box in the figure
    ax2.YAxisLocation = 'right'; % move ax2 at the right
    set(gca, 'Color', 'None') % make ax1 front and ax2 back
    
    % linkaxes([ax1 ax2]) % Link axes in case of zooming
    
    if saveplot
        print(form1,fullfile(folderplot,['plot_PL_TKEO_' form3 form2]))
        close(fig4.Number)
    end
end % of if plotflag > 0

%------------------------------------------------------------------------
% END LATTICE FIBONACCI SEARCH
%------------------------------------------------------------------------

    function elbowPt = ProfLik(elbowInd, distr, flag_short)
        if ~exist('flag_short','var'), flag_short = 0; end
        if ~isnumeric('flag_short'), flag_short = 0; end
        if flag_short < 0 || flag_short > 1, flag_short = 0; end
        % Function to compute Profile Likelihood with respect to index elbowInd
        % Distributions1 = {'Gauss','Laplace','Cauchy','Logistic','Extreme value'}; % Distributions on -inf,+inf
        Distributions2 = {'Lognormal','Weibull','Gamma','BirnbaumSaunders','Burr','Exponential'}; % Distributions on 0,+inf
        nbDataPoints = length(EMG);
        if elbowInd > nbDataPoints
            elbowPt = -1e30*(elbowInd-nbDataPoints);
        else
            % Compute parameters for both sides
            elbowPt = 0;
            for ipiece = 1:2 % Left and right of the candidate onset
                % Put in smpl the sample at the left (ipiece=1) or right
                % (ipiece=2) of the candidate onset
                switch ipiece
                    case 1
                        if flag_short == 0
                            if verbose>2, fprintf('ProfLik - ipiece %i,flag_short %i\n',ipiece,flag_short), end
                            if verbose>2, fprintf('ProfLik - max(1,PLindex1-10) elbowInd %i %i \n',max(1,PLindex1-10),elbowInd), end
                            smpl = EMG(1:elbowInd);
                        else
                            smpl = EMG(max(1,PLindex1-10):elbowInd);
                        end
                    case 2
                        if flag_short == 0
                            smpl = EMG(elbowInd+1:nbDataPoints);
                        else
                            if verbose>2, fprintf('ProfLik - ipiece %i,flag_short %i\n',ipiece,flag_short), end
                            if verbose>2, fprintf('ProfLik - elbowInd+1 PLindex2+10 %i %i \n',elbowInd+1,PLindex2+10), end
                            smpl = EMG(elbowInd+1:PLindex2+10);
                        end
                end
                % Rectify the sample if the distribution is defined on [0,+inf[
                if ismember(distr{ipiece},Distributions2), smpl = abs(smpl); end
                Ndata = length(smpl);
                switch lower(distr{ipiece}(1:5))
                    case 'gauss'
                        Mean = mean(smpl);
                        Var = var(smpl);
                        Llik = -Ndata*log(2*pi*Var)/2 - sum((smpl-Mean).^2)/2/Var;
                    case 'lapla'
                        smpl = sort(smpl);
                        LocationLaplace = smpl(round(Ndata/2));
                        ScaleLaplace = mean(abs(smpl-LocationLaplace));
                        Llik = -Ndata*log(2*ScaleLaplace) - sum(abs(smpl-LocationLaplace))/ScaleLaplace;
                    case 'cauch'
                        smpl = sort(smpl);
                        locCauchy = mean(smpl(round(0.38*Ndata):round(0.62*Ndata)));
                        scaleCauchy = 0.5*iqr(smpl);
                        Llik = Ndata*log(scaleCauchy/pi) - sum(log((scaleCauchy^2 + (smpl-locCauchy).^2)));
                        % parCau = cauchyfit(datasort(datasort ~= 0));
                        % Llik = Ndata*log/(parCau(2)/pi) - sum(log((parCau(2)^2 + (smpl-parCau(1)).^2)));
                    case 'logno'
                        smpl(smpl == 0) = []; % Remove 0 to avoid log 0
                        smpl = log(smpl);
                        Mean = mean(smpl);
                        Var = var(smpl);
                        Llik = -Ndata*log(2*pi*Var)/2 - sum(smpl) - sum((smpl-Mean).^2)/2/Var;
                    case 'weibu'
                        smpl(smpl == 0) = [];
                        parWeibull = fitdist(smpl,'Weibull');
                        Llik = length(smpl)*log(parWeibull.B/parWeibull.A^parWeibull.B) + ...
                            (parWeibull.B-1)*sum(log(smpl)) - ...
                            sum(smpl.^parWeibull.B)/parWeibull.A^parWeibull.B;
                        % refpdf = pdf(parWeibull,smpl); % CHECK OK!
                    case 'gamma'
                        smpl(smpl == 0) = [];
                        parGamma = fitdist(smpl,'Gamma');
                        Llik = length(smpl)*(-log(gamma(parGamma.a))-parGamma.a*log(parGamma.b)) + ...
                            (parGamma.a-1)*sum(log(smpl))-sum(smpl)/parGamma.b;
                        % refpdf = pdf(parGamma,smpl); % CHECK OK!
                    case 'expon' % Exponential
                        parExp = 1/mean(smpl);
                        % refcdf = 1-exp(parExp*smpl);
                        Llik = Ndata*log(parExp) - parExp*sum(smpl);
                    case 'extre' % Extreme value
                        % Extreme value of maximum instead of minimum (use -smpl)
                        smpl = -smpl;
                        parExtr = fitdist(smpl,'ExtremeValue');
                        dum = (smpl-parExtr.mu)/parExtr.sigma;
                        Llik = -Ndata*log(parExtr.sigma) + sum(dum) - sum(exp(dum)) ;
                        % refpdf = pdf(parExtr,smpl); % CHECK OK!
                    case 'logis'
                        parLogis = fitdist(smpl,'Logistic');
                        dum = (smpl-parLogis.mu)/parLogis.sigma;
                        Llik = -sum(dum) - Ndata*log(parLogis.sigma) - 2*sum(log(1+exp(-dum)));
                        % refpdf = pdf(parLogis,smpl); % CHECK OK!
                    case 'birnb'
                        smpl(smpl==0) = [];
                        Ndatadum = length(smpl);
                        parBirnbaum = fitdist(smpl,'BirnbaumSaunders');
                        dum1 = sqrt(smpl/parBirnbaum.beta); dum2 = sqrt(parBirnbaum.beta./smpl);
                        Llik = -Ndatadum*log(2*pi)/2 -sum((dum1-dum2).^2)/2/parBirnbaum.gamma^2 ...
                            + sum(log(dum1+dum2)) - Ndatadum*log(2*parBirnbaum.gamma) ...
                            - sum(log(smpl));
                        % refpdf = pdf(parBirnbaum,smpl); % CHECK OK!
                    case 'burr ' % Burr removed because of error that Weibull is better
                        % smpl = random('Burr',1,2,3,1000,1); % TRIAL DATA FOR CHECK
                        smpl(smpl==0) = [];
                        parBurr = fitdist(smpl,'Burr');
                        dum = smpl/parBurr.alpha;
                        Llik = length(smpl)*log(parBurr.k*parBurr.c/parBurr.alpha) ...
                            + (parBurr.c-1)*sum(log(dum)) ...
                            - (parBurr.k+1)*sum(log(1+dum.^parBurr.c));
                        % refpdf = pdf(parBurr,smpl); % CHECK OK!
                end
                elbowPt = elbowPt + Llik;
            end % of for ipiece = 1:2
        end % of if elbowInd > nbDataPoints
        
    end % of function elbowPt

end % of function prolific