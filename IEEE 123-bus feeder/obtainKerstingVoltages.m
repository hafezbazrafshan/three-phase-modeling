function[network]=obtainKerstingVoltages(network)

v2struct(network); 
v2struct(solution); 

kerstingpath=[pwd,'/IEEE-123 feeder data'];

% Kersting Solution
[BusKersting,pu1Kersting,Angle1Kersting,pu2Kersting, Angle2Kersting, pu3Kersting,Angle3Kersting] =...
    importKerstingSolutions([kerstingpath,'/powerflow.txt']);


KerstingVMag=NaN(size(resultsVMag)); 
KerstingVPhase=NaN(size(resultsVPhase));



for ii=1:length(busNamesWithRegs)
    
    busName=busNamesWithRegs{ii};
   
    
    busIdx=find(strcmp(lower(BusKersting),busName)); 
    if ~isempty(busIdx)
    
    vaMag=NaN;
    vbMag=NaN;
    vcMag=NaN;
    vaPhase=NaN;    
    vbPhase=NaN;
    vcPhase=NaN;
    
    
    
         vaMag=pu1Kersting{busIdx}; 
    vaPhase=Angle1Kersting{busIdx}; 
    if vaMag<0.1
    vaMag=NaN;
    vaPhase=NaN;
    end
   
    
        
    
      

        vbMag=pu2Kersting{busIdx}; 
    vbPhase=Angle2Kersting{busIdx}; 
    if vbMag<0.1
    vbMag=NaN;
    vbPhase=NaN;
    end

    
        
        
         
    vcMag=pu3Kersting{busIdx}; 
    vcPhase=Angle3Kersting{busIdx}; 
    if vcMag<0.1
    vcMag=NaN;
    vcPhase=NaN;
    end
   
    
    KerstingVMag(ii,:)=[vaMag, vbMag, vcMag]; 
    KerstingVPhase(ii,:)=[vaPhase, vbPhase,vcPhase]; 

    else 
        pause;
    end
    
end



Kerstingsolution.KerstingVMag=KerstingVMag;
Kerstingsolution.KerstingVPhase=KerstingVPhase;
network.Kerstingsolution=Kerstingsolution;
end

