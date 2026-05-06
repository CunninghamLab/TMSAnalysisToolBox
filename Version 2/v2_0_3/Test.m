

clear; clc;

%User inputs
numRows=1;
numColumns=2;

%Create figure
CustomOpts=uifigure();
CustomOpts.Position=[680 458 130*numColumns 75*numRows];
Grid=uigridlayout(CustomOpts,[numRows+1 numColumns]);
Grid.RowHeight(:,1:numRows)={'fit'};
Grid.ColumnWidth(:,1:numColumns)={'fit'};
mcdEditFieldLabel=uilabel(Grid,'text','MCD Constant');
mcdEditFieldLabel.Layout.Row=1; mcdEditFieldLabel.Layout.Column=1;
mcdEditField=uieditfield(Grid, 'numeric'); %,'ValueChangedFcn',@(btn,event) MCDConstantChanged());
mcdEditField.Layout.Row=1; mcdEditField.Layout.Column=2;
ContinueButton=uibutton(Grid,'text','Continue');
ContinueButton.ButtonPushedFcn='uiresume(CustomOpts)';
ContinueButton.Layout.Row=numRows+1; %last row
ContinueButton.Layout.Column=1;

uiwait(CustomOpts);



disp('Done');