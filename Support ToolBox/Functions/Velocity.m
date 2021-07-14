function [VelocityData] = Velocity(data,time)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    VelocityData = zeros(size(data,1), size(time,2));
    for iT = 1:(size(data,1))
        for i = 1:length(time)-1
            VelocityData(iT,i) = (data(iT,i+1)-data(iT,i))/(time(i+1)-time(i)) ;
        end
    end
end

