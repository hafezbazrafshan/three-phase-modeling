function [ y] = naNForEmptyCell( x )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

if ~isempty(x)
    y=x;
else
    y=NaN;
end
    
end

