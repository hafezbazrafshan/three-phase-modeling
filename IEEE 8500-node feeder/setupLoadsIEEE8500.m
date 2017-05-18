function [ network ] = setupLoadsIEEE8500( network )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

v2struct(network); 
v2struct(noLoadQuantities);
availableBusIndices=network.availableBusIndices;








%% 4. Importing loads:
 sL=sparse(zeros(N*3,1));
 % import loads : 
[trName,trBusName] =importLoadXfmrs([datapath,'/LoadXfmrs.csv']);
 load([datapath,'/PQloads.mat']); 
 for ii=1:length(trBusName)
     nIdx=find(strcmp(busNames,trBusName(ii)));
     trNameForPhase=char(trName(ii));
     loadIdx= find( and(strcmp(strtrim(lower(Element)),['transformer.',lower(char(trName(ii)))]), Terminal==1));
    if ~isempty(loadIdx)
     sLoad=1000*(PkW(loadIdx(1))+sqrt(-1)*Qkvar(loadIdx(1)))/Sbase;
         phase=trNameForPhase(end); 
         switch phase
             case 'A'
                 sL((nIdx-1)*3+1)=sLoad;
             case 'B'
                  sL((nIdx-1)*3+2)=sLoad;
             case 'C'
                 sL((nIdx-1)*3+3)=sLoad;
         end
    end
         
         
     end
     
 
 
sL_load=sL(availableBusIndices(1:end-3));
sL_loadPerPhase=reshape(sL,3,N).';


loadQuantities=v2struct(sL_load, sL_loadPerPhase); 
network.loadQuantities=loadQuantities;
end

