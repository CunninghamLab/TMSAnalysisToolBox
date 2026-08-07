%Apply Notch filter 
%Inputs: Data, sample rate, frequency to filter out
%Outputs: filtered data

function filteredData=ProcessNotch(app,Data, SampleRate,Freq)

%Check for invalid values
if (Freq-1) > (SampleRate/2) || (Freq+1) > (SampleRate/2)
    msgbox('Invalid frequency');
    Freq=(SampleRate/2)-2;
    app.NotchEditField.Value=Freq;

end

newData=nan(size(Data))';

%Cuttoff
a=Freq-1;
b=Freq+1;

%Filter
[nb,na] = butter(2, [a b]./(SampleRate/2),'stop');
for i=1:length(Data(:,1)) %for each trial
    TrialData=Data(i,:);
    if any(isnan(TrialData)) %if there are any nans
        NumLoc=find(~isnan(TrialData)); %locations of real numbers
        TrialData=TrialData(NumLoc);
    else 
        NumLoc=1:length(TrialData);
    end

    newData(NumLoc,i)=filtfilt(nb,na,TrialData);
end

filteredData=newData';


end