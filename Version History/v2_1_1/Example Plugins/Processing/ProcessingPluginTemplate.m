%{
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

function [Data, CustomProcessingOpts]=ProcessingPluginTemplate(app, existFig,CustomFunctionName, PluginsFolderName, OrigData)

%% Create figure
%User Inputs============================================================================================================

%variables
numVar=3; %how many variables

%Give each variable a label
ListofVariableLabels={'Variable 1 Label', 'Variable 2 Label', 'Variable 3 Label'};
%============================================================================================================

%Do not edit this line
[CustomProcessingOpts,UserVar]=createProcessingPluginFigure(existFig,app.CustomProcessingOpts.(CustomFunctionName),PluginsFolderName,numVar,ListofVariableLabels,mfilename); %Create pop-up figure, do not edit

%% Processing Code
%Define your user variables
%UserVar is a numVarx2 cell that contains the labels and values in the pop-up figure
%Column 1: labels  Column 2 values
%Then write your detection method

Variable1=UserVar{1,2};
Variable2=UserVar{2,2};
Variable3=UserVar{3,2}; 

% VV Write your processing code here VV 


end



