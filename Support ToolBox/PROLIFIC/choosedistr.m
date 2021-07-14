function distr = choosedistr(data,test,SamplingRate,fCutLow,fCutHigh,varargin)
% function choosedistr
% Search the best distribution at the left and right of an Onset (estimated
% by PROLIFIC using distr_fg distribution) according to the
% Kolmogorov-Smirnov (KS), Lilliefors (L) and Anderson-Darling (AD) test
%
% Distribution      | Allowed test |      Range
% ------------------+--------------+------------------
% Gauss             |    KS L AD   |  [-infty,+infty]
% Extreme Value     |    KS L AD   |  [-infty,+infty]
% Laplace           |    KS        |  [-infty,+infty]
% Cauchy            |    KS        |  [-infty,+infty]
% Logistic          |    KS        |  [-infty,+infty]
% Lognormal         |    KS L AD   |  [0,+infty]
% Weibull           |    KS L AD   |  ]0,+infty]
% Gamma             |    KS        |  ]0,+infty]
% Birnbaum-Saunders |    KS        |  [0,+infty]
% Exponential       |    KS L AD   |  [0,+infty]
% Burr              |    KS        |  ]0,+infty]
%
% INPUT:
% data         Vector containing the EMG signal
%              Mandatory argument, no default
% test         Test for finding the best distributions:
%              KS: Kolmogorov-Smirnov test
%              L : Lilliefors test (default)
%              AD: Anderson-Darling test
%              Mandatory argument, no default
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
% varargin     Optional arguments in the form of the couple ('option', value),
%              where possible options are:
%               'dist_fg' Distribution to estimate initial Onset (string)
%                         One of the allowed distributions (default: 'Weibull')
%               'verbose' Amount of printout:
%                         0: No printout (default)
%                         1: Minimal printout
%                         >1: Extensive printout
%              'plotflag' Plot produced:
%                         0:  no plot produced (default)
%                         >0: Plots produced
%             'maintitle' Main title of the plots (string)
%              'saveplot' Whether to save the plots on files eventually
%                         0:  Do not save plots (figures not closed)
%                         >0: Save plots (figures closed, default)
%            'folderplot' Subfolder where to save plots eventually (string)
%                         (default current directory)
%            'formatplot' Format of the saved plots eventually (string)
%                         'pdf': PDF (default)
%                         'eps': EPS
%                         'jpg': JPEG
%             'removedis' Distributions to remove from the search
%                         String (if only 1 distribution to remove) or
%                         cell of strings for 2 or more distributions,
%                         e.g., 'Burr' or {'Burr','Logistic'}
%                         Default: '' (no distribution to remove)
%
% OUTPUT:
% distr Cell array of 2 elements containing distributions chosen by the test at
%       the left and right of the Onset, respectively
%
% Requires Matlab Toolbox Statistics and Machine Learning
%
% EXAMPLES:
% Default values
% x = randn(5000,1) + [repmat(0,1500,1); repmat(3,3500,1)];
% distr = choosedistr(x,'L',1500);
%
% Example with extensive output
% distr = choosedistr(x,'AD',1500,0,0,'verbose',1,'plotflag',1,...
%                             'saveplot',1,'formatplot','pdf')
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
%
% Known issues:
% A 1-second pause has been inserted between printouts of plots because
%   Matlab does not correctly flushes the previous plot
% lillietest and adtest p-values are not accurate for small p (they are
%   indeed tabulated so that ties often occur if no distribution well fits
%   data). In case of ties then the distribution with the smallest
%   statistic is chosen
%
% Version: 0.99 - June 2020
%
% Copyright: S.E. Selvan, D. Allexandre, U. Amato, B. Della Vecchia, G.H. Yue
%

fCutLow_def = 0;
fCutHigh_def = 0;

verbose_def = 0;
plotflag_def = 0;
maintitle_def = '';
saveplot_def = 0;
folderplot_def = 'plot';
formatplot_def = 'pdf';
test_def = 'L';
distr_fg_def = 'Weibull';
removedis_def = '';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse variable input arguments of the function %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Count input arguments
nArgs = length(varargin);
if round(nArgs/2) ~= nArgs/2
    error('choosedistr: requires an even number of parameters')
end

if ~exist('data','var') || isempty(data)
    warning('choosedistr: data non provided and no default. Exit program')
    return
end
if ~isnumeric(data)
    error('PROLIFIC: data not numeric. Exit program')
end
if length(data) <= 1
    error('PROLIFIC: data not array. Exit program')
end
if ~exist('test','var') || isempty(test)
    warning('choosedistr: test not assigned - setting default')
    test = test_def;
end
if ~ischar(test)
    warning('choosedistr: test not a string - setting default')
    test = test_def;
end
if ~any(contains({'KS','L','AD'},upper(test)))
    warning(['choosedistr: value for test (' test ') not recognized - setting default'])
    test = test_def;
end
% Check and process fixed arguments
if ~exist('SamplingRate','var') || isempty(SamplingRate)
    error('choosedistr: SamplingRate non provided and no default. Exit program')
end
if ~isnumeric(SamplingRate)
    error('choosedistr: SamplingRate not numeric. Exit program')
end
if SamplingRate <= 0
    error('choosedistr: SamplingRate negative. Exit program')
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

switch upper(test)
    case 'KS'
        Distributions1 = {'Gauss','Laplace','Cauchy','Logistic','Extreme Value'}; % Distributions on -inf,+inf
        Distributions2 = {'Lognormal','Weibull','Gamma','Exponential','Birnbaum-Saunders','Burr'}; % Distributions on 0,+inf
    case 'L'
        Distributions1 = {'Gauss','Extreme Value'}; % Distributions on -inf,+inf
        Distributions2 = {'Lognormal','Weibull','Exponential'}; % Distributions on 0,+inf
    case 'AD'
        Distributions1 = {'Gauss','Extreme Value'}; % Distributions on -inf,+inf
        Distributions2 = {'Lognormal','Weibull','Exponential'}; % Distributions on 0,+inf
end
Distributions = [Distributions1, Distributions2];

% Parse variable input arguments
for n = 1:2:nArgs
    arg = lower(varargin{n});
    switch arg
        case 'verbose', verbose = varargin{n+1};
        case 'plotflag', plotflag = varargin{n+1};
        case 'maintitle', maintitle = varargin{n+1};
        case 'saveplot', saveplot = varargin{n+1};
        case 'folderplot', folderplot = varargin{n+1};
        case 'formatplot', formatplot = varargin{n+1};
        case 'distr_fg', distr_fg = varargin{n+1};
        case 'removedis', removedis = varargin{n+1};
        otherwise, warning('%s is not a recognized argument name',varargin{n})
    end
end

% Check if variable input arguments exist
if ~exist('verbose','var'), verbose = verbose_def; end
if ~exist('plotflag','var'), plotflag = plotflag_def; end
if ~exist('maintitle','var'), maintitle = maintitle_def; end
if ~exist('saveplot','var'), saveplot = saveplot_def; end
if ~exist('folderplot','var'), folderplot = folderplot_def; end
if ~exist('formatplot','var'), formatplot = formatplot_def; end
if ~exist('distr_fg','var'), distr_fg = distr_fg_def; end
if ~exist('removedis','var'), removedis = removedis_def; end

% Check and process variable input arguments

% Check and process verbose
if ~isnumeric(verbose)
    warning(['choosedistr: value for verbose (' verbose ') not recognized - setting default'])
    verbose = verbose_def;
end
verbose = round(verbose);
if verbose < 0, verbose = verbose_def; end

% Check and process removedis
if ~isempty(removedis)
    if ~ischar(removedis) && ~iscell(removedis)
        warning('choosedistr: value for removedis (%f) not recognized - setting default',removedis)
        removedis = [];
    elseif ischar(removedis)
        if ~any(strcmpi(Distributions,removedis))
            warning(['choosedistr: value for removedis (' removedis ') not recognized - setting default'])
            removedis = [];
        else
            removedis = {removedis}; % Make a cell array with one element
        end
    else % removedis is a cell array
        Nrem = length(removedis);
        indtodrop = zeros(Nrem,1);
        for nrem = 1:Nrem
            if ~ischar(removedis{nrem})
                warning(['choosedistr: value for removedis (' removedis{nrem} ') not recognized - dropped'])
                indtodrop(nrem) = 1;
            else
                if ~any(strcmpi(Distributions,removedis{nrem}))
                    warning(['choosedistr: value for removedis (' removedis{nrem} ') not recognized - dropped'])
                    indtodrop(nrem) = 1;
                end
            end
        end
        removedis(indtodrop==1) = [];
    end
end
if ~isempty(removedis)
    for distr = removedis
        Distributions1(strcmpi(Distributions1,char(distr))) = [];
        Distributions2(strcmpi(Distributions2,char(distr))) = [];
    end
end


% Check and process plotflag
if ~isnumeric(plotflag)
    warning(['choosedistr: value for plotflag (' plotflag ') not recognized - setting default'])
    plotflag = plotflag_def;
end
plotflag = round(plotflag);
if plotflag < 0, plotflag = plotflag_def; end
if plotflag > 0
    % Check and process maintitle
    if ~ischar(maintitle)
        warning(['choosedistr: value for maintitle (' num2str(maintitle) ') not recognized - setting default'])
        maintitle = maintitle_def;
    end
    
    % Check and process saveplot
    if ~isnumeric(saveplot)
        warning(['choosedistr: value for saveplot (' saveplot ') not recognized - setting default'])
        saveplot = saveplot_def;
    end
    saveplot = round(saveplot);
    if saveplot < 0, saveplot = saveplot_def; end
    
    % Check and process folderplot
    if ~ischar(folderplot)
        warning(['choosedistr: value for folderplot (' num2str(folderplot) ') not recognized - setting default'])
        folderplot = folderplot_def;
    end
    if ~exist(folderplot,'dir')
        warning(['choosedistr: folder (' folderplot ') not existing - creating'])
        mkdir(folderplot);
    end
    % Check and process formatplot
    if ~ischar(formatplot)
        warning(['choosedistr: value for formatplot (' num2str(formatplot) ') not recognized - setting default'])
        formatplot = formatplot_def;
        
    end
    if ~any(strcmpi({'pdf','jpeg','jpg','eps'},formatplot))
        warning(['choosedistr: value for formatplot (' formatplot ') not admittable - setting default'])
        formatplot = formatplot_def;
    end
    switch lower(formatplot)
        case {'jpg','jpeg'}
            form1 = '-djpeg100'; form2 = '.jpg';
        case 'pdf'
            form1 = '-dpdf'; form2 = '.pdf';
        case 'eps'
            form1 = '-depsc'; form2 = '.eps';
    end
    
    hist_gray = [0.75 0.75 0.75];
    cmap = lines(100); % Default parula color map for lines (only 7 different)
    fontsize = 13;
    labDir = {'Left','Right'}; % Labels for plot titles
end

Ndistr1 = length(Distributions1);
Ndistr2 = length(Distributions2);
Distributions = [Distributions1, Distributions2];
Ndistr = Ndistr1 + Ndistr2;
distr_type = [ones(1,Ndistr1) ones(1,Ndistr2)+1];
ind1 = find(distr_type == 1);
ind2 = find(distr_type == 2);

if verbose > 0
    fprintf('Search of optimal dsitribution by %s test\n',test)
    fprintf('Candidate distributions: ')
    for i=1:length(Distributions1), fprintf('%s ',Distributions1{i}), end
    for i=1:length(Distributions2), fprintf('%s ',Distributions2{i}), end
    fprintf('\n')
end

N = length(data);
[finIndex,~,~,~,Icoarse,PLcoarse] = prolific(data,SamplingRate,fCutLow,fCutHigh,distr_fg);

% Find the last maximum lastindmax (possibly offset)
% Find all maxima of PL on the coarse grid
Ncoarse = length(Icoarse);
indmax = find(PLcoarse(2:Ncoarse-1) > PLcoarse(1:Ncoarse-2) & ...
    PLcoarse(2:Ncoarse-1) > PLcoarse(3:Ncoarse)) + 1;
% Eventually add the first and/or last point if larger than the nearest ones
if PLcoarse(Ncoarse) > PLcoarse(Ncoarse-1), indmax = [indmax; Ncoarse]; end
if PLcoarse(1) > PLcoarse(2), indmax = [1; indmax]; end
if length(indmax) > 1
    lastindmax = Icoarse(indmax(end));
else
    lastindmax = N;
end

% Remove first and last values from EMG because 0 due the TKEO algorithm
ind = [2,max(2,finIndex-1); finIndex,min(N-1,max(lastindmax,finIndex+10))];

% TKEO-preprocess the signal
EMG = EnergyTKEO(data,SamplingRate,fCutLow,fCutHigh);

distr = cell(2,1);
if plotflag > 0
    fig1 = figure('NumberTitle', 'off', 'Name', [test ' - CDF - All Real']); clf
    fig2 = figure('NumberTitle', 'off', 'Name', [test ' - CDF - Positive']); clf
    fig3 = figure('NumberTitle', 'off', 'Name', [test ' - PDF - All Real']); clf
    fig4 = figure('NumberTitle', 'off', 'Name', [test ' - PDF - Positive']); clf
end
for ipiece=1:2
    datasort1 = sort(EMG(ind(ipiece,1):ind(ipiece,2)));
    datasort2 = sort(abs(EMG(ind(ipiece,1):ind(ipiece,2))));
    Ni = ind(ipiece,2)-ind(ipiece,1)+1;
    stat = zeros(Ndistr,1);
    pv = zeros(Ndistr,1);
    labdistr = cell(Ndistr+1,1);
    labdistr{1} = 'Data';
    if plotflag > 0
        figure(fig1)
        subplot(2,1,ipiece)
        h1 = zeros(Ndistr1+1,1);
        h1(1) = plot(datasort1,(1:Ni)/Ni,'--k','LineWidth',2);
        hold on
        
        figure(fig2)
        subplot(2,1,ipiece)
        h2 = zeros(Ndistr2+1,1);
        h2(1) = plot(datasort2,(1:Ni)/Ni,'--k','LineWidth',2);
        hold on
        
        hist1_Nbin = round(length(datasort1)/50);
        hist2_Nbin = round(length(datasort2)/50);
        histmax = zeros(2,1);
        figure(fig3)
        subplot(2,1,ipiece)
        h3 = zeros(Ndistr1+1,1);
        hist1 = histogram(datasort1,hist1_Nbin,'FaceColor',hist_gray,'EdgeColor',hist_gray,'FaceAlpha',1);
        histmax(1) = max(hist1.BinCounts);
        hold on
        h3(1) = plot(NaN,NaN,'LineWidth',5,'Color',hist_gray);
        
        figure(fig4)
        h4 = zeros(Ndistr2+1,1);
        subplot(2,1,ipiece)
        hist2 = histogram(datasort2,hist2_Nbin,'FaceColor',hist_gray,'EdgeColor',hist_gray,'FaceAlpha',1);
        histmax(2) = max(hist2.BinCounts);
        hold on
        h4(1) = plot(NaN,NaN,'LineWidth',5,'Color',hist_gray);
    end
    pdfmax = [-1,-1];
    ndistr1 = 0; ndistr2 = 0;
    for ndistr = 1:Ndistr
        if distr_type(ndistr) == 1, ndistr1 = ndistr1 + 1; datasort = datasort1; end
        if distr_type(ndistr) == 2, ndistr2 = ndistr2 + 1; datasort = datasort2; end
        switch lower(Distributions{ndistr}(1:4))
            case 'gaus'
                switch upper(test)
                    case 'KS'
                        Mean = mean(datasort);
                        Std = std(datasort);
                        refcdf = normcdf(datasort,Mean,Std);
                        refpdf = normpdf(datasort,Mean,Std);
                        [~,pv(ndistr),stat(ndistr)] = kstest(datasort,'CDF',[datasort refcdf]);
                    case 'L'
                        [~, pv(ndistr),stat(ndistr)] = lillietest(datasort,'Distribution','normal');
                        params = [mean(datasort), std(datasort)];
                        refcdf = normcdf(datasort, params(1), params(2));
                        refpdf = normpdf(datasort, params(1), params(2));
                    case 'AD'
                        [~, pv(ndistr),stat(ndistr)] = adtest(datasort,'Distribution','normal');
                        params = [mean(datasort), std(datasort)];
                        refcdf = normcdf(datasort, params(1), params(2));
                        refpdf = normpdf(datasort, params(1), params(2));
                end
            case 'logn'
                switch upper(test)
                    case 'KS'
                        Mean = mean(log(datasort(datasort>0)));
                        Std = std(log(datasort(datasort>0)));
                        refcdf = normcdf(log(datasort),Mean,Std);
                        refpdf = normpdf(log(datasort),Mean,Std)./datasort;
                        [~,pv(ndistr),stat(ndistr)] = kstest(datasort,'CDF',[datasort refcdf]);
                    case 'L' % log(Lognormal) ~ Guassian
                        [~, pv(ndistr),stat(ndistr)] = lillietest(log(datasort(datasort>0)),'Distribution','normal');
                        params = [mean(log(datasort(datasort>0))), std(log(datasort(datasort>0)))];
                        refcdf = normcdf(log(datasort),params(1), params(2));
                        % Fy(y)=Fx(g^-1(y)) = Fx(log y)
                        refpdf = normpdf(log(datasort),params(1), params(2));
                        % fy(y)=fx(g^-1(y))d/dy g^-1(y)
                        % with y=g(x)=exp(x); d/dy g^-1(y) = d/dy log y = 1/y
                        refpdf = refpdf./datasort;
                    case 'AD' % log(Lognormal) ~ Guassian
                        [~, pv(ndistr),stat(ndistr)] = adtest(log(datasort(datasort>0)),'Distribution','normal');
                        params = [mean(log(datasort(datasort>0))), std(log(datasort(datasort>0)))];
                        refcdf = normcdf(log(datasort),params(1), params(2));
                        % Fy(y)=Fx(g^-1(y)) = Fx(log y)
                        refpdf = normpdf(log(datasort),params(1), params(2));
                        % fy(y)=fx(g^-1(y))d/dy g^-1(y)
                        % with y=g(x)=exp(x); d/dy g^-1(y) = d/dy log y = 1/y
                        refpdf = refpdf./datasort;
                end
            case 'weib'
                switch upper(test)
                    case 'KS'
                        params = fitdist(datasort(datasort>0),'Weibull');
                        refcdf = cdf(params,datasort);
                        refpdf = pdf(params,datasort);
                        [~,pv(ndistr),stat(ndistr)] = kstest(datasort,'CDF',[datasort refcdf]);
                    case 'L' % log(Weibull) ~ Extreme Value
                        [~, pv(ndistr),stat(ndistr)] = lillietest(log(datasort(datasort>0)),'Distribution','ev');
                        params = evfit(log(datasort(datasort>0)));
                        refcdf = evcdf(log(datasort),params(1),params(2));
                        refpdf = evpdf(log(datasort),params(1),params(2)); % because parameters estimated on log(datasort)
                        refpdf = refpdf./datasort; % because parameters estimated on log(datasort)
                    case 'AD' % log(Weibull) ~ Extreme Value
                        [~, pv(ndistr),stat(ndistr)] = adtest(log(datasort(datasort>0)),'Distribution','ev');
                        params = evfit(log(datasort(datasort>0)));
                        refcdf = evcdf(log(datasort),params(1),params(2));
                        refpdf = evpdf(log(datasort),params(1),params(2)); % because parameters estimated on log(datasort)
                        refpdf = refpdf./datasort; % because parameters estimated on log(datasort)
                end
            case 'extr'
                datasort = -datasort(Ni:-1:1); % EV of maximum instead of minimum
                switch upper(test)
                    case 'KS'
                        params = fitdist(datasort,'ev');
                        refcdf = cdf(params,datasort);
                        refpdf = pdf(params,datasort);
                        [~,pv(ndistr),stat(ndistr)] = kstest(datasort,'CDF',[datasort refcdf]);
                    case 'L'
                        [~, pv(ndistr),stat(ndistr)] = lillietest(datasort,'Distribution','ev');
                        params = evfit(datasort);
                        refcdf = evcdf(datasort,params(1),params(2));
                        refpdf = evpdf(datasort,params(1),params(2));
                    case 'AD'
                        [~, pv(ndistr),stat(ndistr)] = adtest(datasort,'Distribution','ev');
                        params = evfit(datasort);
                        refcdf = evcdf(datasort,params(1),params(2));
                        refpdf = evpdf(datasort,params(1),params(2));
                end
                datasort = -datasort(Ni:-1:1); refpdf = refpdf(Ni:-1:1); refcdf = 1 - refcdf(Ni:-1:1); % EV of maximum instead of minimum
            case 'expo'
                switch upper(test)
                    case 'KS'
                        params = fitdist(datasort,'exponential');
                        refcdf = cdf(params,datasort);
                        refpdf = pdf(params,datasort);
                        [~,pv(ndistr),stat(ndistr)] = kstest(datasort,'CDF',[datasort refcdf]);
                    case 'L'
                        [~, pv(ndistr),stat(ndistr)] = lillietest(datasort(datasort>0),'Distribution','exponential');
                        params = mean(datasort(datasort>0));
                        refcdf = expcdf(datasort, params(1));
                        refpdf = exppdf(datasort, params(1));
                    case 'AD'
                        [~, pv(ndistr),stat(ndistr)] = adtest(datasort(datasort>0),'Distribution','exponential');
                        params = mean(datasort(datasort>0));
                        refcdf = expcdf(datasort, params(1));
                        refpdf = exppdf(datasort, params(1));
                end
            case 'gamm'
                switch upper(test)
                    case 'KS'
                        params = fitdist(datasort(datasort>0),'Gamma');
                        refcdf = cdf(params,datasort);
                        refpdf = pdf(params,datasort);
                        [~,pv(ndistr),stat(ndistr)] = kstest(datasort,'CDF',[datasort refcdf]);
                end
            case 'lapl'
                switch upper(test)
                    case 'KS'
                        Median = datasort(round(Ni/2));
                        Location = mean(abs(datasort-Median));
                        refcdf = 0.5*exp(-abs(datasort - Median)/Location);
                        indtmp = find(datasort - Median > 0);
                        refcdf(indtmp) = 1 - refcdf(indtmp);
                        refpdf = 0.5*exp(-abs(datasort - Median)/Location)/Location;
                        [~,pv(ndistr),stat(ndistr)] = kstest(datasort,'CDF',[datasort refcdf]);
                end
            case 'cauc'
                switch upper(test)
                    case 'KS'
                        locCauchy = mean(datasort(round(0.38*Ni):round(0.62*Ni)));
                        scaleCauchy = 0.5*iqr(datasort);
                        refcdf = 0.5 + atan((datasort-locCauchy)./scaleCauchy)/pi;
                        refpdf = 1./(pi*scaleCauchy*(1+((datasort-locCauchy)./scaleCauchy).^2));
                        % parCau = cauchyfit(datasort(datasort ~= 0)); % ,'info2');
                        % refcdf = 0.5 + atan((datasort-parCau(1))./parCau(2))/pi;
                        [~,pv(ndistr),stat(ndistr)] = kstest(datasort,'CDF',[datasort refcdf]);
                end
            case 'logi'
                switch upper(test)
                    case 'KS'
                        params = fitdist(datasort,'Logistic');
                        refcdf = cdf(params,datasort);
                        refpdf = pdf(params,datasort);
                        [~,pv(ndistr),stat(ndistr)] = kstest(datasort,'CDF',[datasort refcdf]);
                end
            case 'birn'
                switch upper(test)
                    case 'KS'
                        params = fitdist(datasort(datasort>0),'BirnbaumSaunders');
                        refcdf = cdf(params,datasort);
                        refpdf = pdf(params,datasort);
                        [~,pv(ndistr),stat(ndistr)] = kstest(datasort,'CDF',[datasort refcdf]);
                end
            case 'burr'
                switch upper(test)
                    case 'KS'
                        params = fitdist(datasort(datasort>0),'Burr');
                        refcdf = cdf(params,datasort);
                        refpdf = pdf(params,datasort);
                        [~,pv(ndistr),stat(ndistr)] = kstest(datasort,'CDF',[datasort refcdf]);
                end
        end
        if plotflag > 0
            Apdf = 1;
            % Value of Apdf computed from the data is approximately 1 and fails if
            % distribution tends to infinity
            % Apdf = 0.5*(sum(refpdf(1:Ni-1).*datasort(2:Ni))-sum(refpdf(2:Ni).*datasort(1:Ni-1)) + ...
            %        refpdf(Ni)*datasort(Ni) - refpdf(1)*datasort(1));
            pdfmax(distr_type(ndistr)) = max(pdfmax(distr_type(ndistr)), max(refpdf)/Apdf*length(datasort)*hist1.BinWidth);
            switch distr_type(ndistr)
                case 1
                    figure(fig1)
                    h1(ndistr1+1) = plot(datasort,refcdf,'LineWidth',1.5,'Color',cmap(ndistr1,:));
                    figure(fig3)
                    h3(ndistr1+1) = plot(datasort,refpdf/Apdf*length(datasort)*hist1.BinWidth, ...
                        'LineWidth',1.5,'Color',cmap(ndistr1,:));
                case 2
                    figure(fig2)
                    h2(ndistr2+1) = plot(datasort,refcdf,'LineWidth',1.5,'Color',cmap(ndistr2,:));
                    figure(fig4)
                    h4(ndistr2+1) = plot(datasort,refpdf/Apdf*length(datasort)*hist2.BinWidth, ...
                        'LineWidth',1.5,'Color',cmap(ndistr2,:));
            end
            labdistr{ndistr+1} = [Distributions{ndistr} sprintf(' %5.3f',stat(ndistr))];
        end
    end
    if plotflag > 0
        figure(fig1), legend(h1,labdistr{[1 ind1+1]},'Location','Best','FontSize',fontsize)
        xlabel('x','interpreter','none'), ylabel('CDF(x)','interpreter','none')
        set(gca,'FontSize',fontsize)
        title([maintitle ' ' labDir{ipiece} ' of the onset'])
        figure(fig2), legend(h2,labdistr{[1 ind2+1]},'Location','Best','FontSize',fontsize)
        xlabel('|x|','interpreter','none'), ylabel('CDF(|x|)','interpreter','none')
        set(gca,'FontSize',fontsize)
        title([maintitle ' ' labDir{ipiece} ' of the onset'])
        figure(fig3)
        if pdfmax(1) > 1.5*histmax(1), ylim([0,histmax(1)*1.3]), end
        legend(h3,labdistr{[1 ind1+1]},'Location','Best','FontSize',fontsize)
        xlabel('x','interpreter','none'), ylabel('PDF(x)','interpreter','none')
        set(gca,'FontSize',fontsize)
        title([maintitle ' ' labDir{ipiece} ' of the onset'])
        figure(fig4)
        if pdfmax(2) > 1.5*histmax(2), ylim([0,histmax(2)*1.3]), end
        legend(h4,labdistr{[1 ind2+1]},'Location','Best','FontSize',fontsize)
        xlabel('|x|','interpreter','none'), ylabel('PDF(|x|)','interpreter','none')
        set(gca,'FontSize',fontsize)
        title([maintitle ' ' labDir{ipiece} ' of the onset'])
    end
    switch upper(test)
        case 'KS'
            % Sort by ascending stat and by descending pv in case of ties
            [~,ibest] = sortrows([pv stat],[2,-1]); ibest = ibest(1);
            % Sort by descending pv and by ascending stat in case of ties
        case {'L','AD'}
            [~,ibest] = sortrows([pv stat],[-1,2]); ibest = ibest(1);
    end
    if verbose > 1
        for ndistr = 1:Ndistr
            signp = '';
            if pv(ndistr) <= 0.1, signp = '.'; end
            if pv(ndistr) <= 0.05, signp = '*'; end
            if pv(ndistr) <= 0.01, signp = '**'; end
            if pv(ndistr) <= 0.001, signp = '***'; end
            fprintf('Piece %i: Distribution %18s (Stat: %5.3f; p-v: %f %3s)\n',ipiece,Distributions{ndistr},stat(ndistr),pv(ndistr),signp)
        end
    end
    if verbose == 1
        signp = '';
        if pv(ibest) <= 0.1, signp = '.'; end
        if pv(ibest) <= 0.05, signp = '*'; end
        if pv(ibest) <= 0.01, signp = '**'; end
        if pv(ibest) <= 0.001, signp = '***'; end
        fprintf('Piece %i: Distribution %s (Stat: %5.3f; p-v: %f %3s) - ',ipiece,Distributions{ibest},stat(ibest),pv(ibest),signp)
    end
    distr{ipiece} = Distributions{ibest};
end
if verbose > 0, fprintf('\n'), end

if plotflag > 0 && saveplot > 0
    figure(fig1)
    print(form1,fullfile(folderplot,['plot_' test '_TKEOcdf' form2]))
    close(fig1)
    pause(1)
    figure(fig2)
    print(form1,fullfile(folderplot,['plot_' test '_absTKEOcdf' form2]))
    close(fig2)
    pause(2)
    figure(fig3)
    print(form1,fullfile(folderplot,['plot_' test '_TKEOpdf' form2]))
    close(fig3)
    pause(1)
    figure(fig4)
    print(form1,fullfile(folderplot,['plot_' test '_absTKEOpdf' form2]))
    close(fig4)
    pause(1)
end

end