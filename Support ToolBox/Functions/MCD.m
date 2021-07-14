function [lowEndLimit,highEndLimit] = MCD(time,data,preStimDurInd,MCDconstant)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


 % Formats according to paper, first column is number of data point, second
            % column is time, third column is emg value
            dataMCD = (preStimDurInd)';
            dataMCD = [dataMCD time(preStimDurInd)'];
            dataMCD = [dataMCD data(preStimDurInd)'];
            
            % Fourth column is dn-(dn-1)
            differences = dataMCD(2:end,3) - dataMCD(1:end-1,3);
            differences = [0;differences];
            dataMCD = [dataMCD differences];
            
            % Fifth column is absolute values of differences
            differencesAbs = abs(differences);
            dataMCD = [dataMCD differencesAbs];
            
            % Finds mean emg and MCD (mean of absolute differences)
            mPreStim = mean(dataMCD(:,3));
            mcd = mean(dataMCD(:,5));
            
            % Finds low end and high end limit by formula of mean emg +- mcd*2.66
            
            % Mcd analysis is based on the number 2.66, given in paper. Otherwise, base
            % the low end limit on num of standard deviations below prestim rms
            lowEndLimit = mPreStim - mcd*MCDconstant;
            highEndLimit = mPreStim + mcd*MCDconstant;
            
  
end

