function [loadBuses,loadTypes,Ph1,Ph2,Ph3,Ph4,Ph5,Ph6]  = importloads(filename1,filename2)
%only for ieee-lv

[~,Bus,phases] = importLoadBusPhase(filename1);
[~,kw,kvar]=importLoadValue(filename2);
loadTypes={};
Ph1=[];
Ph2=[];
Ph3=[];
Ph4=[];
Ph5=[];
Ph6=[];

loadBuses=Bus;
loadTypes=repmat({'Y-PQ'},size(Bus,1),1); 
Ph1=zeros(size(Bus));
Ph2=zeros(size(Bus)); 
Ph3=zeros(size(Bus));
Ph4=zeros(size(Bus));
Ph5=zeros(size(Bus));
Ph6=zeros(size(Bus));

Ph1Idx=find(strcmp(phases,'A')); 
Ph3Idx=find(strcmp(phases,'B')); 
Ph5Idx=find(strcmp(phases,'C')); 

Ph1(Ph1Idx)=kw(Ph1Idx); 
Ph2(Ph1Idx)=kvar(Ph1Idx); 

Ph3(Ph3Idx)=kw(Ph3Idx); 
Ph4(Ph3Idx)=kvar(Ph3Idx);

Ph5(Ph5Idx)=kw(Ph5Idx); 
Ph6(Ph5Idx)=kvar(Ph5Idx); 




    
    
end











