function [ network ] = obtainVoltages( network )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

v2struct(network); 
v2struct(noLoadQuantities);
v2struct(loadQuantities);
v2struct(ZBusResults);

availableBusIndices=network.availableBusIndices;

v3phase=reshape(v,3,N+1).';



% adding n' of the regulators to the labels:
v3phaseRegs=[];


resultsVMag=[abs(v3phase);abs(v3phaseRegs)];
resultsVPhase=[radian2degrees(angle(v3phase)); radian2degrees(angle(v3phaseRegs))];

%% 

solution=v2struct(resultsVMag,resultsVPhase, v3phase,v3phaseRegs, err,success);
network.solution=solution;
end

