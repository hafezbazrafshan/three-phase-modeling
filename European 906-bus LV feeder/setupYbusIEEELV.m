function [ network ] = setupYbusIEEELV( regulatorTypes, epsilon )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if nargin <2 
    epsilon=1e-5;
end
if nargin<1
    regulatorTypes='non-ideal';
end

currentFolder=pwd;
datapath=[pwd,'/LV feeder data'];
[lineBuses1,lineBuses2,lineLengths,lineCodes] = importlines([datapath,'/Lines.csv']);

 %% Setting up the list of nodes, giving them a unique order from 1:N, setting the substation to node N+1

lineBuses1=strtrim(cellstr(num2str(lineBuses1)));
lineBuses2=strtrim(cellstr(num2str(lineBuses2)));
lineStatuses=repmat({''},length(lineBuses1),1);


 %% Substation transformer configuration
lineBuses1=[lineBuses1;'Sourcebus'];
lineBuses2=[lineBuses2; '1']; 
lineCodes=[lineCodes;0]; 
lineStatuses=[lineStatuses; {'subXFM'}]; 
lineLengths=[lineLengths;0];

%% Base voltages:
Sbase=800000; % from the substation transformer
Vbase=416/sqrt(3); % line to neutral conversion (secondary of the substation transformer)
Zbase=(Vbase^2)/(Sbase);
Ybase=1./Zbase;


%% Replacing lineCodes with their integers alternative (manual)
lineCodes(strcmp(lineCodes,'2c_.007'))={1};
lineCodes(strcmp(lineCodes,'2c_.0225'))={2};
lineCodes(strcmp(lineCodes,'2c_16'))={3};
lineCodes(strcmp(lineCodes,'35_SAC_XSC'))={4};
lineCodes(strcmp(lineCodes,'4c_.06'))={5};
lineCodes(strcmp(lineCodes,'4c_.1'))={6};
lineCodes(strcmp(lineCodes,'4c_.35'))={7};
lineCodes(strcmp(lineCodes,'4c_185'))={8};
lineCodes(strcmp(lineCodes,'4c_70'))={9};
lineCodes(strcmp(lineCodes,'4c_95_SAC_XC'))={10};


lineCodes=cell2mat(lineCodes);

%% Impedances in Ohm per mile for all the line configurations [NEED TO CHANGE]
a=exp(1i*2*pi/3); 
T=[1 1 1; 1 a.^2 a; 1 a a.^2];


Z1=[3.97 0 0; 0 3.97 0; 0 0 3.97]+sqrt(-1)*[0.099 0 0; 0 0.099 0; 0 0 0.099]; 
Z2=[1.257 0 0; 0 1.257 0; 0  0 1.257]+ sqrt(-1)*[0.085 0 0; 0 0.085 0; 0 0 0.085]; 
Z3=[1.2 0 0;  0 1.15 0; 0 0 1.15]+sqrt(-1)*[0.088 0 0; 0 0.088 0; 0 0 0.088]; 
Z4=diag([0.76, 0.868, 0.868])+sqrt(-1)*diag([0.092, 0.092, 0.092]); 
Z5=diag([1.581, 0.469, 0.469])+sqrt(-1)* diag([0.091 , 0.075, 0.075]);
Z6=diag([0.959,0.274,0.274])+sqrt(-1)*diag([0.079,0.073,0.073]);
Z7=diag([0.319,0.089,0.089])+sqrt(-1)*diag([0.076,0.0675,0.0675]);
Z8=diag([0.58,0.166,0.166])+sqrt(-1)*diag([0.078,0.068,0.068]);
Z9=diag([1.505,0.446,0.446])+sqrt(-1)*diag([0.083 ,0.071,0.071]);
Z10=diag([0.804, 0.322,0.322])+sqrt(-1)*diag([0.093,0.074 ,0.074 ]);

Zcfg(:,:,1)=T*  Z1*inv(T);
Zcfg(:,:,2)= T*Z2*inv(T);
Zcfg(:,:,3)=T*Z3*inv(T); 
Zcfg(:,:,4)=T*Z4*inv(T);
Zcfg(:,:,5)=T*Z5*inv(T); 
Zcfg(:,:,6)=T*Z6*inv(T); 
Zcfg(:,:,7)=T*Z7*inv(T); 
Zcfg(:,:,8)=T*Z8*inv(T); 
Zcfg(:,:,9)=T*Z9*inv(T); 
Zcfg(:,:,10)=T*Z10*inv(T); 

Zcfg=Zcfg./1000; % convert km to meters.

Ycfg=zeros(size(Zcfg));




%% 3.  Organizing bus names
% collecting bus names:
busNames=unique([lineBuses1;lineBuses2]);


%% Putting substation index at the very end:
% substation='799';
substation='SourceBus';
% substation='1';
substationIndex=find(strcmp(busNames,substation));

if substationIndex< length(busNames) % if substation is not the end bus
    busNames=[busNames(1:substationIndex-1); busNames(substationIndex+1:end); busNames(substationIndex)];
end


busNamesWithRegs=[busNames];
regulatorVoltageGains=[]; 
regulatorCurrentGains=[];
regulatorImpedances=[]; 
regulatorYNMn=[];
regulatorYNMm=[];


%% 5. Create Ytilde and Ybranch
% Ybranch is simplified to be without Yshunt (for testing purposes)
% Ytilde considers Yshunt and Ycaps
nBuses=length(busNames);
busSet=1:nBuses;
JbusSet=1:3*(nBuses);
nBr=length(lineBuses1);
JBranchSet=1:3*nBr;
phaseNodesLinIndices=reshape(JbusSet,3,nBuses).';
phaseNodes=repmat([1,2,3],nBuses,1);
N=nBuses-1;

CNodes=sparse(zeros(nBr,nBuses));
CIndices=sparse(zeros(nBr*3, nBuses*3));
Ytilde=sparse(zeros(3*nBuses));
Ybranch=sparse(zeros(3*nBr));

sNodes=getNumericNodeList(lineBuses1,busNames); 
rNodes=getNumericNodeList(lineBuses2,busNames); 






for ii=1:length(lineBuses1)
    jBranchIdx=(ii-1)*3+1:ii*3;
    bus1=lineBuses1(ii);
    bus2=lineBuses2(ii);
    nIdx=find(strcmp(busNames,bus1));
    mIdx=find(strcmp(busNames,bus2));
    jmIdx=(mIdx-1)*3+1:mIdx*3;
     jnIdx=(nIdx-1)*3+1:nIdx*3;
   
    
   
   
    
    CNodes(ii,nIdx)=1;
    CNodes(ii,mIdx)=-1;
    
    
   
    
    Cmat=zeros(3);
   
        YtildeNMn=zeros(3,3);
    YtildeNMm=zeros(3,3);
    YtildeMNn=zeros(3,3);
    YtildeMNm=zeros(3,3);

    
    
    
    switch lineStatuses{ii}
        
            
        case 'subXFM'
  zt=0.01*(4i)*3;
          yt=1./zt;
          Y1=diag([yt; yt; yt]); 
          Y2=(1/3) *[2* yt, -yt, -yt;
                            -yt, 2*yt, -yt; 
                            -yt, -yt, 2*yt]; 
                        
           Y3=(1/sqrt(3))*[-yt, yt, 0; 
                                    0, -yt, yt; 
                                    yt, 0, -yt]; 
                                
                
                                
             Y2hat=Y2-1i*epsilon*abs(yt)*eye(3);                    
          
                 
    
   
    
    
    % Ytilde

    YtildeNMn(availablePhases,availablePhases)=Y2hat; 
    YtildeNMm(availablePhases,availablePhases)=-Y3;
    YtildeMNn(availablePhases,availablePhases)=-Y3.';
    YtildeMNm(availablePhases,availablePhases)=Y1;
    
    
   
    
    Ytilde(jnIdx, jnIdx)= Ytilde(jnIdx, jnIdx)+YtildeNMn;
    Ytilde(jnIdx, jmIdx)=-YtildeNMm;
    Ytilde(jmIdx,jmIdx)=Ytilde(jmIdx,jmIdx)+YtildeMNm;
    Ytilde(jmIdx,jnIdx)=-YtildeMNn;
    
    
    % Ybranch
    YbranchII=zeros(3,3);
    YbranchII(availablePhases,availablePhases)=Y2hat;
    Ybranch(jBranchIdx,jBranchIdx)=YbranchII;
%     
%     
%     
        otherwise
            
               
    Zseries=Zcfg(:,:,lineCodes(ii)); 
    Yshunt=Ycfg(:,:,lineCodes(ii));
     availablePhases=find(any(Zseries));
  
 
  
     Zseries=Zseries(availablePhases,availablePhases); 
     Yshunt=Yshunt(availablePhases,availablePhases); 
    
    % length of the line:
    lineLength=lineLengths(ii);


    
    Yseries=inv(Zseries);
    % Ytilde
    
 

    YtildeNMn(availablePhases,availablePhases)=(Yseries/lineLength/Ybase+0*Yshunt*lineLength/Ybase);
    YtildeNMm(availablePhases,availablePhases)=(Yseries/lineLength/Ybase);
    YtildeMNn(availablePhases,availablePhases)=(Yseries/lineLength/Ybase);
    YtildeMNm(availablePhases,availablePhases)=(Yseries/lineLength/Ybase+0*Yshunt*lineLength/Ybase);
    
    Ytilde(jnIdx, jnIdx)= Ytilde(jnIdx, jnIdx)+YtildeNMn;
    Ytilde(jnIdx, jmIdx)=-YtildeNMm;
    Ytilde(jmIdx,jmIdx)=Ytilde(jmIdx,jmIdx)+YtildeMNm;
    Ytilde(jmIdx,jnIdx)=-YtildeMNn;
    
    
    % Ybranch
    YbranchII=zeros(3,3);
    YbranchII(availablePhases,availablePhases)=(Yseries/lineLength/Ybase);
    Ybranch(jBranchIdx,jBranchIdx)=YbranchII;
               Cmat=zeros(3);
    Cmat(availablePhases,availablePhases)=eye(length(availablePhases));
    CIndices( jBranchIdx, jnIdx) = Cmat;
    CIndices( jBranchIdx, jmIdx) = -Cmat; 
            
    end
end
    
    %% create a connected List of Nodes
connectedPath = gen_path(nBuses, CNodes.', sNodes, rNodes); 
connectedPath=[connectedPath,setdiff([1:nBuses],connectedPath)];
    
%% 6. Adding Ycap [no capacitors]
Ycap=sparse(zeros(size(Ytilde)));
    
    
    %% Finding available phases
 availableBranchIndices=find(any(CIndices,2)); 
 missingBranchIndices=find(~any(CIndices,2)); 
 
  availableBusAppIndices=find(any(CIndices,1)).'; 
 missingBusAppIndices=find(~any(CIndices,1)).'; 
 
 availableBusIndices=find(any(Ytilde,2)); 
 missingBusIndices=find(~any(Ytilde,2)); 
 
 % availableBusAppIndices and availableBusIndices should match
 
 YbusS=Ytilde( availableBusIndices, availableBusIndices);
 YcapS=Ycap(availableBusIndices, availableBusIndices); 
Ynet=YbusS+YcapS;
Y=Ynet(1:end-3,1:end-3);
Y_NS=Ynet(1:end-3,end-2:end);
Y_SS=Ynet(end-2:end,end-2:end);
Ybus=Y;  % the rank deficient admittance matrix




name='IEEE-LV';
network=v2struct(name,busNames, nBuses,busSet, JbusSet,nBr,JBranchSet, busNamesWithRegs,connectedPath,...
     lineBuses1, lineBuses2,lineLengths,lineCodes,lineStatuses,...
    CNodes,CIndices,Ytilde, Ybranch, Ycap, regulatorVoltageGains, regulatorCurrentGains,...
    regulatorImpedances,regulatorYNMn, regulatorYNMm,regulatorTypes,epsilon,...
    phaseNodesLinIndices, phaseNodes, N, datapath,...
    Sbase,Vbase, Zbase, Ybase, ...
    availableBusIndices, missingBusIndices, ...
    YbusS, YcapS, Ynet, Y, Y_NS, Y_SS, Ybus);
    
    
    
    
  
    
    
