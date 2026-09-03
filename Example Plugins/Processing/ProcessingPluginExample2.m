%{
ProcessingPluginExample2 - Simple example file for processing plugin functions

HOW TO USE THIS TEMPLATE
  Edit only the three ZONE blocks below:
    ZONE 1  set numVar, your parameter labels, and optional default values
    ZONE 2  unpack your parameters from UserVar
    ZONE 3  write your processing code 
  Everything else is handled for you by createPluginFigure.

INPUTS (provided by the app - do not change)
  app                  - handle to the main app (e.g. app.Time, in seconds)
  existFig             - whether the parameter pop-up is already open
  CustomFunctionName   - name of custom function
  PluginsFolderName    - folder for this plugin's settings file 
  OrigData             - data from app, numTrials x numSamples double matrix, processes for all the trials of all the conditions

OUTPUTS (you must return these shapes and units)
  Data                  - processed data, numTrials x numSamples double matrix
  CustomOnOffDetectOpts - the pop-up object, returned untouched
  Return [] for any output your method does not compute.
%}
function [Data, CustomProcessingOpts]=ProcessingPluginExample2(app, existFig,CustomFunctionName, PluginsFolderName, OrigData)

% ======================= ZONE 1: your parameters =======================
numVar=3;       % how many parameters you need
ListofVariableLabels={'Divider 1', 'Divider 2', 'Divider 3'};
% =======================================================================
assert(numel(ListofVariableLabels)==numVar, 'numVar must equal the number of labels.');

% --------------------------- DO NOT EDIT -------------------------------
% Builds the parameter pop-up and loads/saves this plugin's settings.
% mfilename tells the helper which settings file belongs to this plugin.
% DefaultValues pre-fills the pop-up the first time (before settings exist).
[CustomProcessingOpts,UserVar]=createProcessingPluginFigure(existFig,app.CustomProcessingOpts.(CustomFunctionName) ... 
    ,PluginsFolderName,numVar,ListofVariableLabels,mfilename); %Create pop-up figure, do not edit
% -----------------------------------------------------------------------

% ======================= ZONE 2: unpack parameters =====================
% UserVar is (numVar) x 2; column 2 holds the values.
Divider1=UserVar{1,2};
Divider2=UserVar{2,2};
Divider3=UserVar{3,2}; 
% =======================================================================

% ======================= ZONE 3: Processing Code =======================

Data=OrigData./Divider1./Divider2./Divider3;

% =======================================================================

end



