%{
ZScoreOnsetOffsetDetect - onset/offset detection by z-scored envelope.

HOW TO USE THIS TEMPLATE
  Edit only the three ZONE blocks below:
    ZONE 1  set numVar and your parameter labels
    ZONE 2  unpack your parameters from UserVar
    ZONE 3  write your detection method and fill the outputs
  Everything else is handled for you by createPluginFigure.

Adapted from the FieldTrip example ft_trialfun_emgdetect (FieldTrip-independent):
  https://www.fieldtriptoolbox.org/example/preproc/trialdef_emg/
For each trial it: (1) high-pass filters, (2) takes the Hilbert envelope,
(3) boxcar-smooths, (4) z-standardizes, (5) thresholds and keeps the first
crossing that PERSISTS (>50% of a window stays past the threshold) as the onset,
and the next qualifying down-crossing as the offset.

INPUTS (provided by the app - do not change)
  app                - handle to the main app (e.g. app.Time, in seconds)
  existFig           - whether the parameter pop-up is already open
  PluginsFolderName  - folder for this plugin's settings file
  AnalyzeSampleRate  - sample rate after processing, Hz
  PreStimData        - baseline window, numSamples x numTrials
  SelectedTrialsData - trial signals in VOLTS, numTrials x 1 cell

OUTPUTS (shapes and units)
  AllOnOffsetTime - numTrials x 2: onset (col 1) & offset (col 2), SECONDS
  OnsetLimit      - 1 x numTrials, onset threshold,  VOLTS   (or [])
  OffsetLimit     - 1 x numTrials, offset threshold, VOLTS   (or [])
  meanPreStimData - 1 x numTrials baseline mean,     VOLTS   (or [])
  CustomAnalysisOpts - the pop-up object, returned untouched
%}
function [AllOnOffsetTime, OnsetLimit, OffsetLimit, meanPreStimData, CustomAnalysisOpts]=ZScoreOnsetOffsetDetect(app, existFig, PluginsFolderName, AnalyzeSampleRate, PreStimData, SelectedTrialsData)

% ======================= ZONE 1: your parameters =======================
numVar = 5;                                   % how many parameters you need
ListofVariableLabels = {'High-pass cutoff (Hz)', 'Smoothing window (ms)', 'Threshold (z-score)', 'Persistence window (ms)', 'Plot trial (0=off)'};
DefaultValues        = [10, 5, 1, 20, 0];     % first-run defaults, same order/units as the labels
% =======================================================================
assert(numel(ListofVariableLabels)==numVar, 'numVar must equal the number of labels.');

% --------------------------- DO NOT EDIT -------------------------------
% Builds the parameter pop-up and loads/saves this plugin's settings.
% mfilename tells the helper which settings file belongs to this plugin.
[CustomAnalysisOpts, UserVar] = createPluginFigure(app, existFig, ...
    app.CustomAnalysisOpts, PluginsFolderName, numVar, ListofVariableLabels, mfilename, DefaultValues);
% -----------------------------------------------------------------------

% ======================= ZONE 2: unpack parameters =====================
% UserVar is (numVar+1) x 2; column 2 holds the values.
% Index 1 is ALWAYS the auto-added Start Time. YOUR parameters start at 2.
StartTime  = UserVar{1,2}*0.001;   % ms -> s, to match app.Time
HPcutoff   = UserVar{2,2};         % Hz  (high-pass cutoff)
SmoothWin  = UserVar{3,2}*0.001;   % ms -> s (boxcar smoothing window)
ZThreshold = UserVar{4,2};         % z-score units (activity threshold)
PersistWin = UserVar{5,2}*0.001;   % ms -> s (a crossing must persist this long)
PlotTrial  = round(UserVar{6,2});  % trial number to diagnose (0 = no plot)
% =======================================================================

% ======================= ZONE 3: detection method ======================
Fs        = AnalyzeSampleRate;   % Hz
TrialTime = app.Time;            % seconds, one time value per sample

% Snap the requested Start Time onto the nearest sample.
[~, StartIndx] = min(abs(TrialTime - StartTime));
StartTime      = TrialTime(StartIndx);

% Echo the snapped Start Time back into the pop-up (field found by grid row 1).
ef = CustomAnalysisOpts.Children.Children;
ef = ef(arrayfun(@(c) string(class(c))=="matlab.ui.control.NumericEditField", ef));
ef(arrayfun(@(c) c.Layout.Row==1, ef)).Value = StartTime*1000;

% Boxcar length in samples, and a Butterworth high-pass designed once.
boxN      = max(1, round(SmoothWin*Fs));
persistN  = max(1, round(PersistWin*Fs));   % persistence window in samples
minFrac   = 0.5;                            % >this fraction of the window must agree
Nyq       = Fs/2;
filtOrder = 4;
useFilter = HPcutoff > 0 && HPcutoff < Nyq;
if useFilter
    [bHP, aHP] = butter(filtOrder, HPcutoff/Nyq, 'high');
end

nTrials         = numel(SelectedTrialsData);
AllOnOffsetTime = nan(nTrials, 2);
OnsetLimit      = nan(1, nTrials);
OffsetLimit     = nan(1, nTrials);
meanPreStimData = [];   % not used by this method

% Diagnostic-plot setup: only active when PlotTrial points at a real trial.
doPlot   = PlotTrial >= 1 && PlotTrial <= nTrials;
plotData = [];

for i = 1:nTrials %for each trial
    x = SelectedTrialsData{i,1}(:).';   % row vector, volts

    %1) High-pass filter (zero-phase). filtfilt needs length > 3*order;
    %   fall back to mean-removal for very short trials.
    if useFilter && numel(x) > 3*filtOrder
        xf = filtfilt(bHP, aHP, x);
    else
        xf = x - mean(x);
    end

    %2) Hilbert envelope (rectified outline of the burst)
    env = abs(hilbert(xf));

    %3) Boxcar smoothing (moving average)
    sm = conv(env, ones(1,boxN)/boxN, 'same');

    %4) z-standardize: mean 0, std 1
    smStd = std(sm);
    if smStd == 0
        smStd = eps;
    end
    z = (sm - mean(sm)) ./ smStd;

    %5) Threshold, then require each crossing to PERSIST: more than minFrac of
    %   the persistence window must stay on the correct side of the threshold.
    above = z > ZThreshold;
    above(1:StartIndx-1) = false;       % ignore anything before Start Time

    % Onset: first up-crossing whose window is ABOVE the threshold.
    on = [];
    for c = find(diff([false, above]) == 1)   % first sample of each run above
        w = c:min(numel(z), c+persistN-1);
        if mean(z(w) > ZThreshold) > minFrac
            on = c; break;
        end
    end

    if ~isempty(on)
        AllOnOffsetTime(i,1) = TrialTime(on);

        % Offset: first down-crossing after onset whose window is BELOW.
        off   = [];
        downs = find(diff([false, above]) == -1);   % first sample back below
        for c = downs(downs > on)
            w = c:min(numel(z), c+persistN-1);
            if mean(z(w) <= ZThreshold) > minFrac
                off = c; break;
            end
        end
        if ~isempty(off)
            AllOnOffsetTime(i,2) = TrialTime(off);
        end

        %Convert the z-score threshold back to volts for the plotted line.
        OnsetLimit(i)  = mean(sm) + ZThreshold*smStd;
        OffsetLimit(i) = OnsetLimit(i);
    end

    % Stash this trial's data for the optional diagnostic plot.
    if doPlot && i == PlotTrial
        plotData = struct('t',TrialTime, 'raw',x, 'sm',sm, 'z',z, ...
                          'onT',AllOnOffsetTime(i,1), 'offT',AllOnOffsetTime(i,2));
    end
end %end for each trial

% Optional diagnostic plot: shows the envelope, threshold, and detected
% onset/offset for one trial, so you can see how the inputs change detection.
if doPlot
    plotZScoreDiagnostic(plotData, ZThreshold, StartTime, PlotTrial);
elseif PlotTrial >= 1
    warning('Plot trial %d exceeds the number of trials (%d).', PlotTrial, nTrials);
end
% =======================================================================

end %end function

% ===================== DO NOT EDIT BELOW =====================
function plotZScoreDiagnostic(D, zThresh, startT, trialNum)
% One-trial diagnostic. Reuses a tagged figure so re-runs update in place.
figTag = 'ZScoreDetectDiagnostic';
f = findobj('Type','figure','Tag',figTag);
if isempty(f)
    f = figure('Tag',figTag,'Color','w','Name','Z-score onset/offset diagnostic');
else
    figure(f); clf(f);
end

t = D.t(:).'; raw = D.raw(:).'; sm = D.sm(:).'; z = D.z(:).';
onT = D.onT; offT = D.offT;

smStd = std(sm); if smStd==0, smStd = eps; end
vThresh = mean(sm) + zThresh*smStd;   % z threshold expressed in volts

% Panel 1: raw signal + smoothed envelope (volts)
ax1 = subplot(2,1,1); hold(ax1,'on'); box(ax1,'on');
plot(ax1, t, raw, 'Color',[0.75 0.75 0.75]);
plot(ax1, t, sm,  'Color',[0 0.2 0.8], 'LineWidth',1.5);
yline(ax1, vThresh, '--', 'threshold', 'Color',[0.85 0.2 0.2]);
markOnOff(ax1, startT, onT, offT);
ylabel(ax1,'volts');
title(ax1, 'raw signal and smoothed envelope');
legend(ax1, {'raw','envelope'}, 'Location','best'); hold(ax1,'off');

% Panel 2: z-scored envelope (the space detection actually works in)
ax2 = subplot(2,1,2); hold(ax2,'on'); box(ax2,'on');
zpad = 0.05*(max(z)-min(z)) + eps;
yl   = [min(z)-zpad, max(z)+zpad];
active = z > zThresh; active(t < startT) = false;   % what detection uses
shadeRuns(ax2, t, active, yl);
plot(ax2, t, z, 'Color',[0 0.2 0.8], 'LineWidth',1.5);
yline(ax2, zThresh, '--', sprintf('z = %.2f', zThresh), 'Color',[0.85 0.2 0.2]);
markOnOff(ax2, startT, onT, offT);
ylim(ax2, yl); xlabel(ax2,'time (s)'); ylabel(ax2,'z-score');
title(ax2,'z-scored envelope - detection space'); hold(ax2,'off');
end

function markOnOff(ax, startT, onT, offT)
xline(ax, startT, ':', 'Color',[0 0.5 0]);
if ~isnan(onT),  xline(ax, onT,  '-', 'Color',[0 0 0.8]);  end
if ~isnan(offT), xline(ax, offT, '-', 'Color',[0.6 0 0.6]); end
end

function shadeRuns(ax, t, mask, yl)
% Shade each contiguous run where mask is true.
d = diff([false, logical(mask), false]);
s = find(d==1); e = find(d==-1)-1;
for k = 1:numel(s)
    x0 = t(s(k)); x1 = t(e(k));
    patch(ax, [x0 x1 x1 x0], [yl(1) yl(1) yl(2) yl(2)], [1 0.85 0.85], ...
          'EdgeColor','none', 'FaceAlpha',0.5);
end
end
% =======================================================================
