%Update Trial Select List Box
%Inputs: Public app properties, number of trials
%Outputs: Updates trials listed in TrialSelectListBox

function updateMultiTrialList(app,numTrials)

%update number of trials in list
app.TrialSelectListBox.Items=string(1:numTrials);




end


