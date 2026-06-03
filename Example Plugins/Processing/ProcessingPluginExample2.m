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

function [Data, CustomProcessingOpts]=ProcessingPluginExample2(app, existFig,CustomFunctionName, PluginsFolderName, OrigData)

%% Create figure
%User Inputs============================================================================================================

%variables
numVar=3; %how many variables

%Give each variable a label
ListofVariableLabels={'Multiple Factor', 'Variable 2 Label', 'Variable 3 Label'};
%============================================================================================================

%Do not edit this line
[CustomProcessingOpts,UserVar]=createFigure(existFig,app.CustomProcessingOpts.(CustomFunctionName),PluginsFolderName,numVar,ListofVariableLabels); %Create pop-up figure, do not edit

%% Detection Method
%Define your user variables
%UserVar is a numVarx2 cell that contains the labels and values in the pop-up figure
%Column 1: labels  Column 2 values
%Then write your detection method

Variable1=UserVar{1,2};
Variable2=UserVar{2,2};
Variable3=UserVar{3,2}; 

% VV Write your processing code here VV 
Data=OrigData*Variable1;

end


function [CustomProcessingOpts,UserVar]=createFigure(existFig,appCustomProcessingOpts,PluginsFolderName,numVar,ListofVariableLabels)


existSettings=isfile([PluginsFolderName '\' mfilename '_Settings.mat']);
CustomProcessingOpts=appCustomProcessingOpts;
numRows=numVar;
if existFig == 0 %if the figure doesn't already exist
    numColumns=2; %number of columns, one for the label one for the variable

    CustomProcessingOpts.Position(3:4)=[140*numColumns 40*(numRows+1)];
    Grid=uigridlayout(CustomProcessingOpts,[numRows+1 numColumns]);
    Grid.RowHeight(:,1:numRows)={'fit'};
    Grid.ColumnWidth(:,1:numColumns)={'fit'};

    for i=1:numRows %Create the edit field for each variable
        Var.(['Var' num2str(i) 'Label'])=uilabel(Grid,'text',ListofVariableLabels{i});
        Var.(['Var' num2str(i) 'Label']).Layout.Row=i; Var.(['Var' num2str(i) 'Label']).Layout.Column=1;
        Var.(['Var' num2str(i) 'EditField'])=uieditfield(Grid, 'numeric');
        Var.(['Var' num2str(i) 'EditField']).Layout.Row=i; Var.(['Var' num2str(i) 'EditField']).Layout.Column=2;
    end

    UpdateButton=uibutton(Grid,'text','Update','ButtonPushedfcn', @(src,event) UpdateButtonPushed(CustomProcessingOpts,Var));
    UpdateButton.Enable=1;
    UpdateButton.Layout.Row=numRows+1; %last row
    UpdateButton.Layout.Column=1;

    if existSettings == 1 %fill in with settings
        load(string([PluginsFolderName '\' mfilename '_Settings.mat']),"UserVar");
        for i=1:numRows %Create the edit field for each variable
            Var.(['Var' num2str(i) 'EditField']).Value=UserVar{i,2};
        end

    end

    uiwait(CustomProcessingOpts);

else %if the pop-up exists
    UserVar=cell(numRows,2);
    UserVar(:,1)=ListofVariableLabels';
    e=1; Values=nan(1,numRows);
    for i=1:length(appCustomProcessingOpts.Children.Children)
        ValuesLoc = string(class(appCustomProcessingOpts.Children.Children(i))) == "matlab.ui.control.NumericEditField";
        if ValuesLoc == 1
            Values(e)=appCustomProcessingOpts.Children.Children(i).Value;
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
    function UpdateButtonPushed(CustomProcessingOpts,Var)

        %When the update button is clicked update the settings file
        uiresume(CustomProcessingOpts);

        %Define User Variables
        UserVar=cell(numRows,2);
        for v=1:numRows
            UserVar{v,2}=Var.(['Var' num2str(v) 'EditField']).Value;
            UserVar{v,1}=Var.(['Var' num2str(v) 'Label']).Text;
        end

        %Update/create settings file
        save(string([PluginsFolderName '\' mfilename '_Settings.mat']), "UserVar");

    end

end %end createFigure()

