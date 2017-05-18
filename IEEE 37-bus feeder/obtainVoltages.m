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
v3phaseRegs=NaN(1,3);

if strcmp(regulatorTypes,'ideal')
 v3phaseRegs(1,:)=(inv(regulatorVoltageGains{1})*v3phase(find(strcmp(busNames,'799')),:).').'; % three phases
else
v3phaseRegs(1,:)=(inv(regulatorVoltageGains{1})*(v3phase(find(strcmp(busNames,'799')),:).'-...
    regulatorImpedances{1}*inv(regulatorCurrentGains{1})*(regulatorYNMn{1}*v3phase(find(strcmp(busNames,'799')),:).'-regulatorYNMm{1}*v3phase(find(strcmp(busNames,'701')),:).'))).';
end

resultsVMag=[abs(v3phase);abs(v3phaseRegs)];
resultsVPhase=[radian2degrees(angle(v3phase)); radian2degrees(angle(v3phaseRegs))];

%% 

solution=v2struct(resultsVMag,resultsVPhase, v3phase,v3phaseRegs, err,success);
network.solution=solution;
end

