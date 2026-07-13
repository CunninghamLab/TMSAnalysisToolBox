%{
Note: this is an example meant to show the user how to use the processing plugin, it is not meant to be used as a proper processing function for data
Custom Processing function
The user should not edit any inputs or outputs of the function
Inputs: These are all from the main app data
app - allows access to main app figure
existFig - is the pop-up fig open 
CustomFunctionName - name of custom function
PluginsFolderName - name of plugins folder
OrigData - data from app, numTrials x numSamples double matrix, processes for all the trials of all the conditions

Outputs:
Data - processed data, numTrials x numSamples double matrix
CustomProcessingOpts - the pop-up figure object, user shouldn't edit this

%}

function [Data, CustomProcessingOpts]=ProcessingPluginExample2(app, existFig,CustomFunctionName, PluginsFolderName, OrigData)

%% Create figure
%User Inputs============================================================================================================

%variables
numVar=2; %how many variables

%Give each variable a label
ListofVariableLabels={'Divider 1', 'Divider 2'};
%============================================================================================================

%Do not edit this line
[CustomProcessingOpts,UserVar]=createProcessingPluginFigure(existFig,app.CustomProcessingOpts.(CustomFunctionName),PluginsFolderName,numVar,ListofVariableLabels,mfilename); %Create pop-up figure, do not edit

%% Processing Code
%Define your user variables
%UserVar is a numVarx2 cell that contains the labels and values in the pop-up figure
%Column 1: labels  Column 2 values
%Then write your detection method

Divider1=UserVar{1,2};
Divider2=UserVar{2,2};

% VV Write your processing code here VV 

Data=OrigData/Divider1/Divider2;

end



