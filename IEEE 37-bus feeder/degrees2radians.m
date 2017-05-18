function [ phaseOut ] = degrees2radians( phaseIn )
%DEGREES2RADIANS converts degrees to radians
%   degrees2radians( phaseIn ) converts the angle phaseIn in degrees to its
%   equivalent in radians.
% Created by Hafez Bazrafshan 8/26/2016
phaseOut=-pi+ (pi/180)*(phaseIn+180);
end
