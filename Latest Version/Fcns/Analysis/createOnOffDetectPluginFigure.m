%{
createPluginFigure - shared parameter pop-up builder for onset/offset detection plugins.

DO NOT EDIT. Every onset/offset detection plugin calls this. It builds the numeric
parameter pop-up, auto-prepends a "Start Time (ms)" field, and loads/saves
the calling plugin's settings file.

Inputs:
  app, existFig, appCustomOnOffDetectOpts, PluginsFolderName, numVar,
  ListofVariableLabels - passed straight through from the plugin.
  PluginName - pass mfilename from the plugin so each plugin gets its own
               <PluginName>_Settings.mat (do NOT hardcode a name here).

Outputs:
  CustomOnOffDetectOpts - the pop-up figure object.
  UserVar - (numVar+1) x 2 cell: column 1 = labels, column 2 = values.
            Row 1 is always the Start Time field.
%}
function [CustomOnOffDetectOpts,UserVar]=createOnOffDetectPluginFigure(app,existFig,appCustomOnOffDetectOpts,PluginsFolderName,numVar,ListofVariableLabels,PluginName,DefaultValues)

settingsFile  = fullfile(PluginsFolderName, [PluginName '_Settings.mat']);
existSettings = isfile(settingsFile);

ListofVariableLabels(2:end+1)=ListofVariableLabels;
ListofVariableLabels{1}='Start Time (ms)';
CustomOnOffDetectOpts=appCustomOnOffDetectOpts;
numRows=numVar+1;

% Optional first-run default field values (used only before a settings file
% exists). Accept numVar entries (your parameters, Start Time defaults to 0)
% or numRows entries (Start Time included).
if nargin < 8 || isempty(DefaultValues)
    DefaultValues = [];
elseif numel(DefaultValues)==numVar
    DefaultValues = [0, DefaultValues(:).'];      % prepend Start Time default
elseif numel(DefaultValues)==numRows
    DefaultValues = DefaultValues(:).';
else
    error('DefaultValues must have %d or %d elements.', numVar, numRows);
end

if existFig == 0 %pop-up doesn't exist yet - build it
    numColumns=2; %one column for the label, one for the value

    CustomOnOffDetectOpts.Position(3:4)=[130*numColumns 40*(numRows+1)];
    Grid=uigridlayout(CustomOnOffDetectOpts,[numRows+1 numColumns]);
    Grid.RowHeight(:,1:numRows)={'fit'};
    Grid.ColumnWidth(:,1:numColumns)={'fit'};

    for i=1:numRows %create a label + numeric field for each row
        Var.(['Var' num2str(i) 'Label'])=uilabel(Grid,'text',ListofVariableLabels{i});
        Var.(['Var' num2str(i) 'Label']).Layout.Row=i; Var.(['Var' num2str(i) 'Label']).Layout.Column=1;
        Var.(['Var' num2str(i) 'EditField'])=uieditfield(Grid, 'numeric');
        Var.(['Var' num2str(i) 'EditField']).Layout.Row=i; Var.(['Var' num2str(i) 'EditField']).Layout.Column=2;
    end

    UpdateButton=uibutton(Grid,'text','Update','ButtonPushedfcn', @(src,event) UpdateButtonPushed(CustomOnOffDetectOpts,Var));
    UpdateButton.Enable=1;
    UpdateButton.Layout.Row=numRows+1; %last row
    UpdateButton.Layout.Column=1;

    if existSettings %pre-fill fields from the saved settings
        load(settingsFile,"UserVar");
        if numRows ~= length(UserVar(:,1)) % the number of variables has changed
            %enter default values
            for i=1:numRows
                Var.(['Var' num2str(i) 'EditField']).Value=DefaultValues(i);
            end
            UpdateButtonPushed(CustomOnOffDetectOpts,Var);
        else
            for i=1:numRows
                Var.(['Var' num2str(i) 'EditField']).Value=UserVar{i,2};
            end
        end
    elseif ~isempty(DefaultValues) %first run: show the defaults
        for i=1:numRows
            Var.(['Var' num2str(i) 'EditField']).Value=DefaultValues(i);
        end
    end

    uiwait(CustomOnOffDetectOpts); %wait for the user to click Update

else %pop-up already exists - read the current field values
    UserVar=cell(numRows,2);
    UserVar(:,1)=ListofVariableLabels';
    e=1; Values=nan(1,numRows);
    for i=1:length(appCustomOnOffDetectOpts.Children.Children)
        if string(class(appCustomOnOffDetectOpts.Children.Children(i)))=="matlab.ui.control.NumericEditField"
            Values(e)=appCustomOnOffDetectOpts.Children.Children(i).Value;
            e=e+1;
        end
    end
    Values(isnan(Values))=[];
    if length(Values) ~= numRows
        error('Number of figure values does not match number of rows');
    end
    UserVar(:,2)=num2cell(Values');
end

    %Nested callback: save settings when Update is clicked
    function UpdateButtonPushed(CustomOnOffDetectOpts,Var)
        uiresume(CustomOnOffDetectOpts);
        UserVar=cell(numRows,2);
        for v=1:numRows
            UserVar{v,2}=Var.(['Var' num2str(v) 'EditField']).Value;
            UserVar{v,1}=Var.(['Var' num2str(v) 'Label']).Text;
        end
        save(settingsFile, "UserVar");
    end

app.StartTimemsEditField.Value=UserVar{1,2}; %for plotting in main graph
end %end createPluginFigure()
