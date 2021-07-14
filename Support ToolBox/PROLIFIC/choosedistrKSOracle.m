function distr = choosedistrKSOracle(data,OnsetGold,SamplingRate,fCutLow,fCutHigh,varargin)
% function choosedistrKSOracle
% Search the best distribution at the left and right of an input Onset
%
% INPUT:
% data         Vector containing the EMG signal
%              Mandatory argument, no default
% OnsetGold    True Onset (visually detected by experts)
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
%            'folderplot' Subfolder where to save plots eventually
%                         (default current directory)
%            'formatplot' Format of the saved plots eventually
%                         'pdf': PDF (default)
%                         'eps': EPS
%                         'jpg': JPEG
%                         'png': PNG
%
% OUTPUT:
% distr: Cell array of 2 elements containing distributions chosen by KS at
%        the left and right of the Onset, respectively
%
% EXAMPLES:
% Default values
% x = randn(5000,1) + [repmat(0,1500,1); repmat(3,3500,1)];
% distr = choosedistrKSOracle(x,1500,1000);
%
% Extensive output
% distr = choosedistrKSOracle(x,1500,2048,0,0,'verbose',1,'plotflag',1,...
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
% A pause has been inserted between printouts of plots because Matlab does
%   not correctly flushes the previous plot
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse variable input arguments of the function %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Count input arguments
nArgs = length(varargin);
if round(nArgs/2) ~= nArgs/2
    error('choosedistrKSOracle: requires an even number of parameters')
end

if ~exist('data','var') || isempty(data)
    error('choosedistrKSOracle: data non provided and no default. Exit program')
end
if ~isnumeric(data)
    error('PROLIFIC: data not numeric. Exit program')
end
if length(data) <= 1
    error('PROLIFIC: data not array. Exit program')
end
if ~exist('OnsetGold','var') || isempty(OnsetGold)
    error('choosedistrKSOracle: OnsetGold non provided and no default. Exit program')
end
% Check and process fixed arguments
if ~exist('SamplingRate','var') || isempty(SamplingRate)
    error('choosedistrKSOracle: SamplingRate non provided and no default. Exit program')
end
if ~isnumeric(SamplingRate)
    error('choosedistrKSOracle: SamplingRate not numeric. Exit program')
end
if SamplingRate <= 0
    error('choosedistrKSOracle: SamplingRate negative. Exit program')
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

% Check and process variable input arguments

% Check and process verbose
if ~isnumeric(verbose)
    warning(['choosedistrKSOracle: value for verbose (' verbose ') not recognized - setting default'])
    verbose = verbose_def;
end
verbose = round(verbose);
if verbose < 0, verbose = verbose_def; end

% Check and process plotflag
if ~isnumeric(plotflag)
    warning(['choosedistrKSOracle: value for plotflag (' plotflag ') not recognized - setting default'])
    plotflag = plotflag_def;
end
plotflag = round(plotflag);
if plotflag < 0, plotflag = plotflag_def; end
if plotflag > 0
    % Check and process maintitle
    if ~ischar(maintitle)
        warning(['choosedistrKSOracle: value for maintitle (' num2str(maintitle) ') not recognized - setting default'])
        maintitle = maintitle_def;
    end
    
    % Check and process saveplot
    if ~isnumeric(saveplot)
        warning(['choosedistrKSOracle: value for saveplot (' saveplot ') not recognized - setting default'])
        saveplot = saveplot_def;
    end
    saveplot = round(saveplot);
    if saveplot < 0, saveplot = saveplot_def; end
    
    if saveplot > 0
        % Check and process folderplot
        if ~ischar(folderplot)
            warning(['choosedistrKSOracle: value for folderplot (' num2str(folderplot) ') not recognized - setting default'])
            folderplot = folderplot_def;
        end
        if ~exist(folderplot,'dir')
            warning(['choosedistrKSOracle: folder (' folderplot ') not existing - creating'])
            mkdir(folderplot);
        end
        % Check and process formatplot
        if ~ischar(formatplot)
            warning(['choosedistrKSOracle: value for formatplot (' num2str(formatplot) ') not recognized - setting default'])
            formatplot = formatplot_def;
            
        end
        if ~any(strcmpi({'pdf','jpeg','jpg','eps','png'},formatplot))
            warning(['choosedistrKSOracle: value for formatplot (' formatplot ') not admittable - setting default'])
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
    
    hist_gray = [0.75 0.75 0.75];
    cmap = lines(100); % Default parula color map for lines (only 7 different)
    fontsize = 13;
    labDir = {'Left','Right'}; % Labels for plot titles
end

Distributions1 = {'Gauss','Laplace','Cauchy','Extreme Value'}; % ,'Logistic'}; % Distributions on -inf,+inf
Distributions2 = {'Lognormal','Weibull','Gamma'}; % ,'BirnbaumSaunders'}; % Distributions on 0,+inf ,'Burr '
Ndistr1 = length(Distributions1);
Ndistr2 = length(Distributions2);
Distributions = [Distributions1, Distributions2];
Ndistr = Ndistr1 + Ndistr2;
distr_type = [ones(1,Ndistr1) repmat(2,1,Ndistr2)];
ind1 = find(distr_type == 1);
ind2 = find(distr_type == 2);

if verbose > 0
    fprintf('ORACLE search of optimal dsitribution for Gold onset %i\n',OnsetGold)
    fprintf('Candidate distributions: ')
    for i=1:length(Distributions1), fprintf('%s ',Distributions1{i}), end
    for i=1:length(Distributions2), fprintf('%s ',Distributions2{i}), end
    fprintf('\n')
end

N = length(data);
ind = [1,OnsetGold-1; OnsetGold,N];

% TKEO-preprocess the signal
EMG = EnergyTKEO(data,SamplingRate,fCutLow,fCutHigh);

if plotflag > 0
    fig1 = figure('NumberTitle', 'off', 'Name', 'OracleKS - CDF - All Real'); clf
    fig2 = figure('NumberTitle', 'off', 'Name', 'OracleKS - CDF - Positive'); clf
    fig3 = figure('NumberTitle', 'off', 'Name', 'OracleKS - PDF - All Real'); clf
    fig4 = figure('NumberTitle', 'off', 'Name', 'OracleKS - PDF - Positive'); clf
end
distr = cell(2,1);

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
        % Ndata = length(datasort);
        % hgrid = linspace(datasort(1),datasort(end),Ngrid);
        switch lower(Distributions{ndistr}(1:5))
            case 'gauss'
                Mean = mean(datasort);
                Std = std(datasort);
                refcdf = normcdf(datasort,Mean,Std);
                refpdf = normpdf(datasort,Mean,Std);
            case 'lapla'
                Median = datasort(round(Ni/2));
                Location = mean(abs(datasort-Median));
                refcdf = 0.5*exp(-abs(datasort - Median)/Location);
                indtmp = find(datasort - Median > 0);
                refcdf(indtmp) = 1 - refcdf(indtmp);
                refpdf = 0.5*exp(-abs(datasort - Median)/Location)/Location;
            case 'cauch'
                locCauchy = mean(datasort(round(0.38*Ni):round(0.62*Ni)));
                scaleCauchy = 0.5*iqr(datasort);
                refcdf = 0.5 + atan((datasort-locCauchy)./scaleCauchy)/pi;
                refpdf = 1./(pi*scaleCauchy*(1+((datasort-locCauchy)./scaleCauchy).^2));
                % parCau = cauchyfit(datasort(datasort ~= 0)); % ,'info2');
                % refcdf = 0.5 + atan((datasort-parCau(1))./parCau(2))/pi;
            case 'logno'
                Mean = mean(log(datasort(datasort>0)));
                Std = std(log(datasort(datasort>0)));
                refcdf = normcdf(log(datasort),Mean,Std);
                refpdf = normpdf(log(datasort),Mean,Std)./datasort;
            case 'weibu'
                parWeibull = fitdist(datasort(datasort>0),'Weibull');
                refcdf = cdf(parWeibull,datasort);
                refpdf = pdf(parWeibull,datasort);
            case 'gamma'
                parGamma = fitdist(datasort(datasort>0),'Gamma');
                refcdf = cdf(parGamma,datasort);
                refpdf = pdf(parGamma,datasort);
            case 'logis'
                parLogis = fitdist(datasort,'Logistic');
                refcdf = cdf(parLogis,datasort);
                refpdf = pdf(parLogis,datasort);
            case 'birnb'
                parBirnbaum = fitdist(datasort(datasort>0),'BirnbaumSaunders');
                refcdf = cdf(parBirnbaum,datasort);
                refpdf = pdf(parBirnbaum,datasort);
            case 'burr ' % Burr removed because of error that Weibull is better
                parBurr = fitdist(datasort(datasort>0),'Burr');
                refcdf = cdf(parBurr,datasort);
                refpdf = pdf(parBurr,datasort);
            case 'extre' % Extreme Value
                datasort = -datasort; % Change sign for EV of maximum and not minimum
                parExtr = fitdist(datasort,'ExtremeValue');
                refcdf = 1-cdf(parExtr,datasort);
                refpdf = pdf(parExtr,datasort);
                datasort = -datasort; % Revert to original datasort
        end
        [~,pv(ndistr),stat(ndistr)] = kstest(datasort,'CDF',[datasort refcdf]);
        if plotflag > 0
            Apdf = 1;
            % Value computed from the data is approximately 1 and fails if
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
        xlabel('$x$','interpreter','latex'), ylabel('$CDF(x)$','interpreter','latex')
        set(gca,'FontSize',fontsize)
        title([maintitle ' ' labDir{ipiece} ' of the onset'])
        figure(fig2), legend(h2,labdistr{[1 ind2+1]},'Location','Best','FontSize',fontsize)
        xlabel('$|x|$','interpreter','latex'), ylabel('$CDF(|x|)$','interpreter','latex')
        set(gca,'FontSize',fontsize)
        title([maintitle ' ' labDir{ipiece} ' of the onset'])
        figure(fig3)
        if pdfmax(1) > 1.5*histmax(1), ylim([0,histmax(1)*1.3]), end
        legend(h3,labdistr{[1 ind1+1]},'Location','Best','FontSize',fontsize)
        xlabel('$x$','interpreter','latex'), ylabel('$PDF(x)$','interpreter','latex')
        set(gca,'FontSize',fontsize)
        title([maintitle ' ' labDir{ipiece} ' of the onset'])
        figure(fig4)
        if pdfmax(2) > 1.5*histmax(2), ylim([0,histmax(2)*1.3]), end
        legend(h4,labdistr{[1 ind2+1]},'Location','Best','FontSize',fontsize)
        xlabel('$|x|$','interpreter','latex'), ylabel('$PDF(|x|)$','interpreter','latex')
        set(gca,'FontSize',fontsize)
        title([maintitle ' ' labDir{ipiece} ' of the onset'])
    end
    
    [~,ibest] = min(stat);
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
    print(form1,fullfile(folderplot,['plot_OracleKS_TKEOcdf' form2]))
    close(fig1)
    pause(1)
    figure(fig2)
    print(form1,fullfile(folderplot,['plot_OracleKS_absTKEOcdf' form2]))
    close(fig2)
    pause(1)
    figure(fig3)
    print(form1,fullfile(folderplot,['plot_OracleKS_TKEOpdf' form2]))
    close(fig3)
    pause(1)
    figure(fig4)
    print(form1,fullfile(folderplot,['plot_OracleKS_absTKEOpdf' form2]))
    close(fig4)
    pause(1)
end

end