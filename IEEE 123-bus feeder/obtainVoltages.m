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
v3phaseRegs=NaN(4,3);


if strcmp(regulatorTypes,'ideal')
v3phaseRegs(1,:)=(inv(regulatorVoltageGains{1})*v3phase(find(strcmp(busNames,'150')),:).').'; % three phases
v3phaseRegs(2,1)=   inv(regulatorVoltageGains{2}(1,1))*v3phase(find(strcmp(busNames,'9')),1);
v3phaseRegs(3,[1,3])=(inv(regulatorVoltageGains{3}([1,3],[1,3]))*v3phase(find(strcmp(busNames,'25')),[1,3]).').';
v3phaseRegs(4,:)=(inv(regulatorVoltageGains{4})*v3phase(find(strcmp(busNames,'160')),:).').';
else
v3phaseRegs(1,:)=(inv(regulatorVoltageGains{1})*(v3phase(find(strcmp(busNames,'150')),:).'-...
    regulatorImpedances{1}*inv(regulatorCurrentGains{1})*(regulatorYNMn{1}*v3phase(find(strcmp(busNames,'150')),:).'-regulatorYNMm{1}*v3phase(find(strcmp(busNames,'149')),:).'))).';
v3phaseRegs(2,1)=(inv(regulatorVoltageGains{2}(1,1))*(v3phase(find(strcmp(busNames,'9')),1).'-...
    regulatorImpedances{2}(1,1)*inv(regulatorCurrentGains{2}(1,1))*(regulatorYNMn{2}(1,1)*v3phase(find(strcmp(busNames,'9')),1).'-regulatorYNMm{2}(1,1)*v3phase(find(strcmp(busNames,'14')),1).'))).';
v3phaseRegs(3,[1,3])=(inv(regulatorVoltageGains{3}([1,3],[1,3]))*(v3phase(find(strcmp(busNames,'25')),[1,3]).'-...
    regulatorImpedances{3}([1,3],[1,3])*inv(regulatorCurrentGains{3}([1,3],[1,3]))*(regulatorYNMn{3}([1,3],[1,3])*v3phase(find(strcmp(busNames,'25')),[1,3]).'-regulatorYNMm{3}([1,3],[1,3])*v3phase(find(strcmp(busNames,'26')),[1,3]).'))).';
v3phaseRegs(4,:)=(inv(regulatorVoltageGains{4})*(v3phase(find(strcmp(busNames,'160')),:).'-...
    regulatorImpedances{4}*inv(regulatorCurrentGains{4})*(regulatorYNMn{4}*v3phase(find(strcmp(busNames,'160')),:).'-regulatorYNMm{4}*v3phase(find(strcmp(busNames,'67')),:).'))).';

end

resultsVMag=[abs(v3phase);abs(v3phaseRegs)];
resultsVPhase=[radian2degrees(angle(v3phase)); radian2degrees(angle(v3phaseRegs))];

%% 

solution=v2struct(resultsVMag,resultsVPhase, v3phase,v3phaseRegs, err,success);
network.solution=solution;
end

