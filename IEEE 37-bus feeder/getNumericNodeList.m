function [ sNodes ] = getNumericNodeList( lineBuses, busNames )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
sNodes=zeros(length(lineBuses),1); 
for jj=1:length(lineBuses) 
    sNodes(jj)=find(strcmp(busNames,lineBuses(jj))); 
    
end



end
