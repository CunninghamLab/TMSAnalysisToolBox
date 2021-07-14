function out = EnergyTKEO(inp, SamplingRate, fCutLow, fCutHigh, order)
%global time
% Performs Teager-Kaiser energy operator (TKEO) pre-conditioning to improve
% contrast between the baseline and the signal and applies a low-, high- or
% band-pass filter to the signal
%
% INPUT
% inp          Input signal
% SamplingRate Sampling rate (observations per second)
%              Values greater than 0; mandatory
% fCutLow      Low-pass frequency of low-pass filter (fCutHigh < 0 or undefined) or
%              lower band-pass frequency of band-bass filter (fCutHigh > 0)
%              Unit: Hz, 0 <= fCutLow < SamplingRate/2
% fCutHigh     High-pass frequency of high-pass filter (fCutLow < 0 or undefined) or
%              high band-pass frequency in band-bass filter (fCutLow >= 0)
%              Unit: Hz, 0 < fCutHigh <= SamplingRate/2
% order        Order of the filter
%              order > 0 (default: 2)
%
% OUTPUT
% out          Filtered and TKEO-preconditioned signal
%
% References:
% S.E. Selvan, D. Allexandre, U. Amato, B. Della Vecchia, G.H. Yue:
% A Fast and Robust Profile-Likelihood-Based Muscle Onset Detection in EMG
% Using Discrete Fibonacci Search. To be published in IEEE Access
%
% S.E. Selvan, D. Allexandre, U. Amato, G.H. Yue: Unsupervised Stochastic
% Strategies for Robust Detection of Muscle Activation Onsets in Surface
% Electromyogram. IEEE Trans. Neur. Syst. Rehabil. Engin. 26(6), 1279-1291
%
% Version: 0.99 - June 2020
%
% Copyright: S.E. Selvan, D. Allexandre, U. Amato, G.H. Yue

% Matlab function butter requires normalized frequency strictly in the
% range (0,1), therefore a small delta is introduced to keep within limits
delta = 1e-5;

%%%%%%%%%%%%%%%%%%%%%%
% Checking arguments %
%%%%%%%%%%%%%%%%%%%%%%

% Sampling rate must be positive otherwise skipping filtering
if ~exist('SamplingRate','var')
    error('EnergyTKEO: Sampling rate not provided')
end
if ~isnumeric(SamplingRate)
    error('EnergyTKEO: Sampling rate not numeric')
end
if SamplingRate <= 0
    error('EnergyTKEO: Not positive sampling rate')
end

if ~exist('order','var'), order = 2; end
if ~isnumeric(order)
    warning('EnergyTKEO: order not numeric - set to default')
end
if order <= 0
    warning('EnergyTKEO: order not positive - skip filtering')
    out = inp;
else
    
    %%%%%%%%%%%%%
    % Filtering %
    %%%%%%%%%%%%%
    
    if ~exist('fCutLow','var'), fCutLow = -1; end
    if isempty(fCutLow), fCutLow = -1; end
    if ~isnumeric(fCutLow), fCutLow = -1; end
    fCutLow = min(fCutLow,SamplingRate/2);
    
    if ~exist('fCutHigh','var'), fCutHigh = -1; end
    if isempty(fCutHigh), fCutHigh = -1; end
    if ~isnumeric(fCutHigh), fCutHigh = -1; end
    fCutHigh = min(fCutHigh,SamplingRate/2);
    
    if (fCutLow < 0 && fCutHigh < 0) || fCutLow == fCutHigh
        % No fitering; return the same signal
        a=NaN; b=NaN;
    elseif fCutLow < 0
        % High-pass filter
        if fCutHigh == (SamplingRate/2)
            % No filtering if high-pass frequency maximum
            a = NaN; b = NaN;
        else
            [b,a] = butter(order,fCutHigh/(SamplingRate/2)-delta,'high');
        end
    elseif fCutHigh < 0
        % Low-pass filter
        if fCutLow == 0
            % No filtering if low-pass frequency is 0
            a = NaN; b = NaN;
        else
            [b,a] = butter(order,fCutLow/(SamplingRate/2)-delta,'low');
        end
    elseif fCutHigh > fCutLow
        % Band-pass filter
        [b,a] = butter(order,[fCutLow/(SamplingRate/2)+delta,fCutHigh/(SamplingRate/2)-delta],'bandpass');
    else
        % Band-stop filter
        [b,a] = butter(order,[fCutHigh/(SamplingRate/2)+delta,fCutLow/(SamplingRate/2)-delta],'stop');
    end
    
    % Filtering the input signal
    if isnan(a) | isnan(b), out = inp; else, out = filter(b,a,inp); end
end
%plot(time,out, 'k','LineWidth',1.25);
% TKEO pre-processing
N = length(inp);
out = [0; out(2:N-1).^2-(out(1:N-2).*out(3:N)); 0];


