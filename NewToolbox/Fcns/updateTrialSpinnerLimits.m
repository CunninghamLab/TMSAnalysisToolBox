%Update Trial spinner limits
%Inputs: Public app properties
%Outputs: Updates Trial Spinner limits

function updateTrialSpinnerLimits(app)

    app.TrialSpinner.Limits = [1 app.numTrials];  %update trial spinner limits
    if app.Trial > app.numTrials %if current trial number is larger than new numTrials set trial to last trial
        app.TrialSpinner.Value=app.numTrials;
        app.Trial=app.numTrials;
    end

end






