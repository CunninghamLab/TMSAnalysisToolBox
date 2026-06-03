%{
Custom Onset/Offset Detection Method 
The user should not edit any inputs or outputs of the function
Inputs: These are all from the main app data
app - allows access to main app figure
existFig - is the pop-up fig open 
PluginsFolderName - name of location of the analysis plugins folder
AnalyzeSampleRate - the sample rate of the data after processing
PreStimData - data within the prestim range set in the main app figure, numSamples x numTrials matrix, each column is the preStimData for a trial
SelectedTrialsData - data in the selected trials, unit=volts, numTrials x 1 cell, column 1 contains the data for each trial, if using the average each condition is a trial  

Outputs: These are all expected outputs, ensure that the formats match what is expected
AllOnOffsetTime - the onset and offset times found for each trial, format numTrials x 2 double matrix, column 1 are the onset times, column 2 are the offset times, time is in seconds
CustomAnalysisOpts - the pop-up figure object, user shouldn't edit this
meanPreStimData - mean of the rectified prestim data for each trial, 1 x numTrials double vector
OnsetLimit - the Onset threshold limit for each trial, 1 x numTrials double vector
OffsetLimit - the offset threshold limit for each trial, 1 x numTrials double vector

%Note: if any of the outputs are values that are not needed for your method,
please create empty variables for them. (i.e. meanPrestimData=[])


%}
function [AllOnOffsetTime, OnsetLimit, OffsetLimit, meanPreStimData, CustomAnalysisOpts]=AnalysisPluginTemplate(app, existFig, PluginsFolderName, AnalyzeSampleRate, PreStimData, SelectedTrialsData)

%% Create figure
%User Inputs============================================================================================================

%variables
numVar=3; %how many variables

%Give each variable a label
ListofVariableLabels={'Variable 1 Label', 'Variable 2 Label', 'Variable 3 Label'};
%============================================================================================================

%Do not edit this line
[CustomAnalysisOpts,UserVar]=createFigure(app,existFig,app.CustomAnalysisOpts,PluginsFolderName,numVar,ListofVariableLabels); %Create pop-up figure, do not edit

%% Detection Method
%Define your user variables
%UserVar is a numVarx2 cell that contains the labels and values in the pop-up figure
%Column 1: labels  Column 2 values
%Then write your detection method

StartTime=UserVar{1,2}*0.001; %Convert from ms to seconds, app.Time unit is seconds %Time in data to start the search for an onset, this is automatically the first value
Variable1=UserVar{2,2};
Variable2=UserVar{3,2};
Variable3=UserVar{4,2};

% VV Write your detection method here VV 



%If these variables are not created in the custom detection method, create empty vectors for them
% AllOnOffsetTime=[]; %return in unit seconds
% OnsetLimit=[];      %return in unit volts
% OffsetLimit=[];     %return in unit volts
% meanPreStimData=[]; %return in unit volts

end %end function

function [CustomAnalysisOpts,UserVar]=createFigure(app,existFig,appCustomAnalysisOpts,PluginsFolderName,numVar,ListofVariableLabels)


existSettings=isfile([PluginsFolderName '\' mfilename '_Settings.mat']);
ListofVariableLabels(2:end+1)=ListofVariableLabels;
ListofVariableLabels{1}='Start Time (s)';
CustomAnalysisOpts=appCustomAnalysisOpts;
numRows=numVar+1;
if existFig == 0 %if the figure doesn't already exist
    numColumns=2; %number of columns, one for the label one for the variable

    CustomAnalysisOpts.Position(3:4)=[130*numColumns 40*(numRows+1)];
    Grid=uigridlayout(CustomAnalysisOpts,[numRows+1 numColumns]);
    Grid.RowHeight(:,1:numRows)={'fit'};
    Grid.ColumnWidth(:,1:numColumns)={'fit'};

    for i=1:numRows %Create the edit field for each variable
        Var.(['Var' num2str(i) 'Label'])=uilabel(Grid,'text',ListofVariableLabels{i});
        Var.(['Var' num2str(i) 'Label']).Layout.Row=i; Var.(['Var' num2str(i) 'Label']).Layout.Column=1;
        Var.(['Var' num2str(i) 'EditField'])=uieditfield(Grid, 'numeric');
        Var.(['Var' num2str(i) 'EditField']).Layout.Row=i; Var.(['Var' num2str(i) 'EditField']).Layout.Column=2;
    end

    UpdateButton=uibutton(Grid,'text','Update','ButtonPushedfcn', @(src,event) UpdateButtonPushed(CustomAnalysisOpts,Var));
    UpdateButton.Enable=1;
    UpdateButton.Layout.Row=numRows+1; %last row
    UpdateButton.Layout.Column=1;

    if existSettings == 1 %fill in with settings
        load(string([PluginsFolderName '\' mfilename '_Settings.mat']),"UserVar");
        for i=1:numRows %Create the edit field for each variable
            Var.(['Var' num2str(i) 'EditField']).Value=UserVar{i,2};
        end

    end

    uiwait(CustomAnalysisOpts);

else %if the pop-up exists
    UserVar=cell(numRows,2);
    UserVar(:,1)=ListofVariableLabels';
    e=1; Values=nan(1,numRows);
    for i=1:length(appCustomAnalysisOpts.Children.Children)
        ValuesLoc = string(class(appCustomAnalysisOpts.Children.Children(i))) == "matlab.ui.control.NumericEditField";
        if ValuesLoc == 1
            Values(e)=appCustomAnalysisOpts.Children.Children(i).Value;
            e=e+1;
        end
    end
    Values(isnan(Values))=[]; 
    if length(Values) ~= numRows %if there are different amount of values for how many rows there are throw an error
        error('Number of figure values does not match number of rows');
    end
    UserVar(:,2)=num2cell(Values');


end

    %Function for figure
    function UpdateButtonPushed(CustomAnalysisOpts,Var)

        %When the update button is clicked update the settings file
        uiresume(CustomAnalysisOpts);

        %Define User Variables
        UserVar=cell(numRows,2);
        for v=1:numRows
            UserVar{v,2}=Var.(['Var' num2str(v) 'EditField']).Value;
            UserVar{v,1}=Var.(['Var' num2str(v) 'Label']).Text;
        end

        %Update/create settings file
        save(string([PluginsFolderName '\' mfilename '_Settings.mat']), "UserVar");

    end

app.StartTimemsEditField.Value=UserVar{1,2}; %for plotting in main graph
end %end createFigure()