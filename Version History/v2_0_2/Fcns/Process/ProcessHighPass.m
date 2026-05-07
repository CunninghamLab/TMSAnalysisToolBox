%Apply high pass filter
%Inputs: Data, sample rate, cutoff frequency
%Outpus: filtered data

function filteredData=ProcessHighPass(app,Data,SampleRate,Cutoff)

%check for invalid values
if Cutoff> SampleRate/2
    msgbox('Invalid cutoff frequency');
    Cutoff=SampleRate/2-1;
    app.HighPassEditField.Value=Cutoff;
elseif Cutoff == 0
    Cutoff=1;
    app.HighPassEditField.Value=Cutoff;
end

newData=nan(size(Data))';

%Filter
if strcmp(app.HighPassFilterType.Value, "Butterworth") % butterworth
    [b,a] = butter(app.HighPassFilterOrder.Value,Cutoff/(SampleRate/2),'high');
elseif strcmp(app.HighPassFilterType.Value, "Chebyshev I") % chebyshev type I
    [b,a] = cheby1(app.HighPassFilterOrder.Value, app.HighPassFilterRp.Value, Cutoff/(SampleRate/2),'high');
elseif strcmp(app.HighPassFilterType.Value, "Chebyshev II") % chebyshev type Ii
    [b,a] = cheby2(app.HighPassFilterOrder.Value, app.HighPassFilterRs.Value, Cutoff/(SampleRate/2),'high');
elseif strcmp(app.HighPassFilterType.Value, "Elliptic") % elliptic
    [b,a] = ellip(app.HighPassFilterOrder.Value, app.HighPassFilterRp.Value, app.HighPassFilterRs.Value, Cutoff/(SampleRate/2),'high');
end

for i=1:length(Data(:,1)) %for each trial
    TrialData=Data(i,:);
    if any(isnan(TrialData)) %if there are any nans
        NumLoc=find(~isnan(TrialData)); %locations of real numbers
        TrialData=TrialData(NumLoc);
    else
        NumLoc=1:length(TrialData);
    end

    newData(NumLoc,i)=filtfilt(b,a,TrialData);
end

filteredData=newData';

end