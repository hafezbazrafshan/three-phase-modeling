function [network]=obtainDSSVoltages(network)

v2struct(network); 
v2struct(solution); 

dsspath=[pwd,'/IEEELVopenDSSdata'];
%% 1. DSS solution
 [Bus,Node1,Angle1,pu1,Node2,Angle2,pu2,Node3,Angle3,pu3] = importDSSvoltage([dsspath,'/LVTest_EXP_VOLTAGES.csv']);

DSSVMag=NaN(size(resultsVMag)); 
DSSVPhase=NaN(size(resultsVPhase));



for ii=1:length(busNamesWithRegs)
    
    busName=busNamesWithRegs{ii};
   
    
    busIdx=find(strcmp(lower(Bus),lower(busName))); 
    if ~isempty(busIdx)
    
    vaMag=NaN;
    vbMag=NaN;  
    vcMag=NaN;
    vaPhase=NaN;
    vbPhase=NaN;
    vcPhase=NaN;
    
    
    if Node1(busIdx)==1
         vaMag=pu1(busIdx); 
    vaPhase=Angle1(busIdx); 
    if vaMag<0.1
    vaMag=NaN;
    vaPhase=NaN;
    end
    elseif Node1(busIdx)==2
        vbMag=pu1(busIdx); 
    vbPhase=Angle1(busIdx); 
    if vbMag<0.1
    vbMag=NaN;
    vbPhase=NaN;
    end
    elseif Node1(busIdx)==3
    vcMag=pu1(busIdx); 
    vcPhase=Angle1(busIdx); 
    if vcMag<0.1
    vcMag=NaN;
    vcPhase=NaN;
    end
       
    end
    
        
    
        if Node2(busIdx)==1
         vaMag=pu2(busIdx); 
    vaPhase=Angle2(busIdx); 
    if vaMag<0.1
    vaMag=NaN;
    vaPhase=NaN;
    end
    elseif Node2(busIdx)==2
        vbMag=pu2(busIdx); 
    vbPhase=Angle2(busIdx); 
    if vbMag<0.1
    vbMag=NaN;
    vbPhase=NaN;
    end
    elseif Node2(busIdx)==3
    vcMag=pu2(busIdx); 
    vcPhase=Angle2(busIdx); 
    if vcMag<0.1
    vcMag=NaN;
    vcPhase=NaN;
    end
      
        end
    
        
        
               if Node3(busIdx)==1
         vaMag=pu3(busIdx); 
    vaPhase=Angle3(busIdx); 
    if vaMag<0.1
    vaMag=NaN;
    vaPhase=NaN;
    end
    elseif Node3(busIdx)==2
        vbMag=pu3(busIdx); 
    vbPhase=Angle3(busIdx); 
    if vbMag<0.1
    vbMag=NaN;
    vbPhase=NaN;
    end
               elseif Node3(busIdx)==3
    vcMag=pu3(busIdx); 
    vcPhase=Angle3(busIdx); 
    if vcMag<0.1
    vcMag=NaN;
    vcPhase=NaN;
    end
    
               end
    
    DSSVMag(ii,:)=[vaMag, vbMag, vcMag]; 
    DSSVPhase(ii,:)=[vaPhase, vbPhase,vcPhase]; 
    else 
        pause;
    end
    
end

DSSsolution.DSSVMag=DSSVMag;
DSSsolution.DSSVPhase=DSSVPhase;
network.DSSsolution=DSSsolution;



