function [network]=computeNoLoadVoltage(network, v0mags, v0phases)

if nargin <3 
    v0mags=[1;1;1]; 
    v0phases=degrees2radians([0; -120; 120]); 
end


v2struct(network); 




%% 3. Setting up slack fixed voltage and finding the no-load solution
v0phases=degrees2radians(v0phases);
v0=v0mags.*exp(sqrt(-1)*v0phases); 
w=-inv(Y)*Y_NS*v0;
wCheck=zeros(3*(N+1),1);
wCheck(availableBusIndices,1)=[w;v0];
wCheck(missingBusIndices,1)=NaN;
w3phase=reshape(wCheck,3,N+1).';


noLoadQuantities=v2struct(v0mags,v0phases,v0,w,wCheck,w3phase); 
network.noLoadQuantities=noLoadQuantities;