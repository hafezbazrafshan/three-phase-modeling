function [ yLoadImpedance ] =getYLoadImpedance(cMat, ePage, yL_load )
%GETYLOADIMPEDANCE obtains the load impedance matrix to be incorporated
%into the Z-BUS method.

J=size(yL_load,1); 

for j=1:J
    yLoadImpedance(j,:)= cMat(j,1).*yL_load(j,1).*ePage(j,:,1) + cMat(j,2).*yL_load(j,2).*ePage(j,:,2);
    
end
end

