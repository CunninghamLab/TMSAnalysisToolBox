%Apply low pass filter
%Inputs: Data, sample rate, cutoff frequency 
%Outpus: filtered data

function filteredData=ProcessLowPass(app,Data,SampleRate,Cutoff)

%check for invalid values
if Cutoff> SampleRate/2
    msgbox('Invalid cutoff frequency');
    Cutoff=SampleRate/2-1; 
    app.LowPassEditField.Value=Cutoff;
elseif Cutoff == 0
    Cutoff=1;
    app.LowPassEditField.Value=Cutoff;
end

newData=nan(size(Data))';

%Filter
[b,a] = butter(6,Cutoff/(SampleRate/2),'low');
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