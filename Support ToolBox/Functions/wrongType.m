function [A, B, C] = wrongType()
% WRONGTYPE Summary of this function goes here
%     [A,B, C] = wrongType()
% [A,B,C] refers the [app.EditField.Visible,app.EventTypePanel.Visible,
% app.SetEventWindowPanel.Visible]
%
% if the data type does not matck, reimport the data.Display the
% instruction.
global A B C

A = 'on';
B = 'off';
C = 'off';

end

