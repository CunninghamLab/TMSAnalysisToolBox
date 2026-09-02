%check override in custom analysis function
%Inputs:current trial, app variable, current MissingAnalyze value, current Start and End values
%Outputs: Start and End values, new MissingAnalyze value, Analyze value


function [Start,End,MissingAnalyze,Analyze]=checkOverride(i,app,MissingAnalyze,Start,End)
if app.OverrideUsed == 0 %Override is not used
    %Determine the index for the onset and offset
    OnsetTime=app.AllOnOffsetTime(i,1);
    OffsetTime=app.AllOnOffsetTime(i,2);
    Tol=eps("double");
    Start=find(abs(app.Time - OnsetTime) < Tol);
    End=find(abs(app.Time - OffsetTime) < Tol);
    %if either are not found, skip this trial
    if isempty(Start) || isempty(End)
        MissingAnalyze=MissingAnalyze+1; %add to the Missing Analyze counter, MissingAnalyze is initialized at 0 in the main app
        Analyze=0; %don't analyze this trial, fill with nans
    else
        Analyze=1; %analyze this trial
    end
else
    Analyze=1;

end


end