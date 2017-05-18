function [ phaseOut ] = radian2degrees( phaseIn )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

phaseOut=-180+((phaseIn+pi)./(2*pi))*360;
end

