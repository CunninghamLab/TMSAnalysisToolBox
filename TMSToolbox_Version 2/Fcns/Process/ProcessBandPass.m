%Apply band pass filter
%Inputs: Data, sample rate, lower and upper bound of filter
%Outputs: filtered data

function filteredData=ProcessBandPass(app,Data,SampleRate,StartHz,EndHz)

%Cutoff
a=StartHz;
b=EndHz;

%check for invalid values
if a > SampleRate/2
    msgbox('Invalid start frequency');
    a=SampleRate/2-2; %default to second highest acceptable value
    app.BandPassStartEditField.Value=a;
elseif a == 0
    a=1;
    app.BandPassStartEditField.Value=a;
end
if b > SampleRate/2
    msgbox('Invalid end frequency');
    b=SampleRate/2-1; %default to highest acceptable value
    app.BandPassEndEditField.Value=b;
elseif b == 0
    b=1;
    app.BandPassEndEditField.Value=b;
end

newData=nan(size(Data))';

%Filter
[bpb,bpa] = butter(4,[a b]./(SampleRate/2),'bandpass');
for i=1:length(Data(:,1)) %for each trial
    TrialData=Data(i,:);
    if any(isnan(TrialData)) %if there are any nans
        NumLoc=find(~isnan(TrialData)); %locations of real numbers
        TrialData=TrialData(NumLoc);
    else 
        NumLoc=1:length(TrialData);
    end

    newData(NumLoc,i)=filtfilt(bpb,bpa,TrialData);
end

filteredData=newData';

end