function [ network ] = setupYBusIEEE37( regulatorTypes, epsilon )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if nargin <2 
    epsilon=1e-5;
end
if nargin<1
    regulatorTypes='non-ideal';
end

currentFolder=pwd;
datapath=[pwd,'/IEEE-37 feeder data'];
[lineBuses1,lineBuses2,lineLengths,lineCodes] = importlines([datapath,'/Line Data.xls']);

 %% Setting up the list of nodes, giving them a unique order from 1:N, setting the substation to node N+1

lineBuses1=cellstr(num2str(lineBuses1));
lineBuses2=cellstr(num2str(lineBuses2));
lineStatuses=repmat({''},length(lineBuses1),1);




%% Xfm-1 transformer configuration
XfmIdx=find(isnan(lineCodes)); 
lineCodes(XfmIdx,1)=2; % transformer is type 2
lineStatuses{XfmIdx,1}='Xfm-1'; 


%% Regulator configuration
regIdx=find(strcmp(lineBuses1,'799')); 
lineStatuses{regIdx,1}='reg1'; 

% % %% Substation transformer configuration
lineBuses1=[lineBuses1;'sourcebus'];
lineBuses2=[lineBuses2; '799']; 
lineCodes=[lineCodes;2]; 
lineStatuses=[lineStatuses; {'subXFM'}]; 
lineLengths=[lineLengths;0];



%% Base voltages:
Sbase=2500000; % from the substation transformer
Vbase=4800/sqrt(3); % line to neutral conversion (secondary of the substation transformer)
Zbase=(Vbase^2)/Sbase;
Ybase=1./Zbase;


%% Impedances in Ohm per mile for all the line configurations
Zcfg(:,:,1)=[0.2926+0.1973i, 0.0673-0.0368i, 0.0337-0.0417i; 
                0.0673-0.0368i, 0.2646+0.1900i, 0.0673-0.0368i;
                0.0337-0.0417i,  0.0673-0.0368i, 0.2926+0.1973i]/5280;
Ycfg(:,:,1)=1j*(10^(-6))*[159.8   0   0;
    0        159.8   0;
   0        0     159.8]/5280;  % for configuration 721



Zcfg(:,:,2)= [0.4751+0.2973i, 0.1629-0.0326i, 0.1234-0.0607i; 
                    0.1629-0.0326i, 0.4488+0.2678i, 0.1629-0.0326i;
                  0.1234-0.0607i, 0.1629-0.0326i, 0.4751+0.2973i ]/5280;

Ycfg(:,:,2)=1j*(10^(-6))* [127.8306  0   0;
    0       127.8306   0;
   0 0  127.8306]/5280;  % for configuration 722

Zcfg(:,:,3)=[1.2936+0.6713i, 0.4871+0.2111i, 0.4585+0.1521i; 
    0.4871+0.2111i, 1.3022+0.6326i, 0.4871+0.2111i;
    0.4585+0.1521i, 0.4871+0.2111i, 1.2936+0.6713i]/5280;
Ycfg(:,:,3)=1j*(10^(-6))*[74.8405   0   0;
    0   74.8405   0;
    0   0    74.8405]/5280;


Zcfg(:,:,4)=[2.0952+0.7758i, 0.5204+0.2738i, 0.4926+0.2123i;
    0.5204+0.2738i, 2.1068+0.7398i, 0.5204+0.2738i;
    0.4926+0.2123i,  0.5204+0.2738i ,   2.0952+0.7758i]/5280;
Ycfg(:,:,4)=1j*(10^(-6))*[60.2483  0    0;
    0    60.2483   0;
   0  0        60.2483]/5280;


%% 3.  Organizing bus names
% collecting bus names:
busNames=unique([lineBuses1;lineBuses2]);

%% Putting substation index at the very end:
% substation='799';
substation='sourcebus';
substationIndex=find(strcmp(busNames,substation));

if substationIndex< length(busNames) % if substation is not the end bus
    busNames=[busNames(1:substationIndex-1); busNames(substationIndex+1:end); busNames(substationIndex)];
end




busNamesWithRegs=[busNames;'799r'];  
regulatorVoltageGains=cell(1,1); 
regulatorCurrentGains=cell(1,1); 
regulatorImpedances=cell(1,1); 
regulatorYNMn=cell(1,1);
regulatorYNMm=cell(1,1); 



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
        case 'reg1' % reg 1
          tapAB=7;
            tapBC=4;

 
            arAB= 1-0.00625*tapAB;
            arBC=1-0.00625*tapBC;
            
             Av=[arAB, 1-arAB,0; 
                    0, 1, 0 ; 
                    0, 1-arBC,arBC];
             Ai=[1/arAB, 0,0; 
                 1-1/arAB,1,1-1/arBC;
                 0, 0, 1/arBC];
             
             
%              

    ztReg=(2.5/2)*3*0.01i;
    Zregulator=[ztReg, 0, 0; 0 0 0; 0 0 ztReg]; 
    


             
             

    
    Zseries=Zcfg(:,:,lineCodes(ii)-720); 
    Yshunt=Ycfg(:,:,lineCodes(ii)-720);
     availablePhases=find(any(Zseries));
  
     Zseries=Zseries(availablePhases,availablePhases); 
     Yshunt=Yshunt(availablePhases,availablePhases); 
    
    % length of the line:
    lineLength=lineLengths(ii);
    
    
        YtildeNMn(availablePhases,availablePhases)=( inv(Zseries*lineLength)/Ybase+0*Yshunt*lineLength/Ybase);
    YtildeNMm(availablePhases,availablePhases)=(inv(Zseries*lineLength)/Ybase);
    YtildeMNn(availablePhases,availablePhases)=(inv(Zseries*lineLength)/Ybase);
    YtildeMNm(availablePhases,availablePhases)=(inv(Zseries*lineLength)/Ybase+0*Yshunt*lineLength/Ybase);
    
    
    FSVR=eye(3)+ YtildeNMn*inv(Av)*Zregulator*Ai; 
    % Ytilde
    if strcmp(regulatorTypes,'ideal')

    Ytilde(jnIdx, jnIdx)= Ytilde(jnIdx, jnIdx)+Ai*(YtildeNMn)*inv(Av);
    Ytilde(jnIdx, jmIdx)=-Ai*YtildeNMm;
    Ytilde(jmIdx,jmIdx)=Ytilde(jmIdx,jmIdx)+YtildeMNm;
    Ytilde(jmIdx,jnIdx)=-YtildeMNn*inv(Av);
 
    else
        
    Ytilde(jnIdx, jnIdx)= Ytilde(jnIdx, jnIdx)+Ai*inv(FSVR)*(YtildeNMn)*inv(Av);
    Ytilde(jnIdx, jmIdx)=-Ai*inv(FSVR)*YtildeNMm;
    Ytilde(jmIdx,jmIdx)=Ytilde(jmIdx,jmIdx)+YtildeMNm - YtildeMNm*inv(Av)*Zregulator*Ai*inv(FSVR)*YtildeNMm;
    Ytilde(jmIdx,jnIdx)=-(YtildeMNn*inv(Av)-YtildeMNn*inv(Av)*Zregulator*Ai*inv(FSVR)*(YtildeNMn)*inv(Av));
    end
regulatorVoltageGains{1}=Av;
regulatorCurrentGains{1}=Ai; 
regulatorImpedances{1}=Zregulator; 
regulatorYNMn{1}=Ai*inv(FSVR)*(YtildeNMn)*inv(Av);
regulatorYNMm{1}=Ai*inv(FSVR)*YtildeNMm;
    
    % Ybranch
    YbranchII=zeros(3,3);
    YbranchII(availablePhases,availablePhases)=(inv(Zseries*lineLength)/Ybase);
    Ybranch(jBranchIdx,jBranchIdx)=YbranchII;
            
            
        
                Cmat=zeros(3);
    Cmat(availablePhases,availablePhases)=eye(length(availablePhases));
    CIndices( jBranchIdx, jnIdx) = Cmat;
    CIndices( jBranchIdx, jmIdx) = -Cmat;
            
        case 'subXFM'

                       zt  =0.01*[2+8i]*3;
            yt=1./zt;
            Y2= (1/3) *[ 2*yt, -yt, -yt; -yt,2*yt,-yt; -yt,-yt,2*yt];
            
            
              Y3=(1/sqrt(3))*[-yt, yt, 0; 
                                    0, -yt, yt; 
                                    yt, 0, -yt]; 
                                
            
            Y1=diag([yt;yt;yt]); 
            
            Y2hat1=Y2+abs(yt)*epsilon*eye(3);
                        Y2hat2=Y2+abs(yt)*(epsilon/2)*eye(3);

            
            
            
            availablePhases=[1;2;3];
    

    
        YtildeNMn(availablePhases,availablePhases)=Y2hat1; 
    YtildeNMm(availablePhases,availablePhases)=Y2hat2;
    YtildeMNn(availablePhases,availablePhases)=Y2hat2;
    YtildeMNm(availablePhases,availablePhases)=Y2hat1;
    

    
    
    
    Ytilde(jnIdx, jnIdx)= Ytilde(jnIdx, jnIdx)+YtildeNMn;
    Ytilde(jnIdx, jmIdx)=-YtildeNMm;
    Ytilde(jmIdx,jmIdx)=Ytilde(jmIdx,jmIdx)+YtildeMNm;
    Ytilde(jmIdx,jnIdx)=-YtildeMNn;
    
    
    % Ybranch
    YbranchII=zeros(3,3);
    YbranchII(availablePhases,availablePhases)=Y2;
    Ybranch(jBranchIdx,jBranchIdx)=YbranchII;
            
            
            
                Cmat=zeros(3);
    Cmat(availablePhases,availablePhases)=eye(length(availablePhases));
    CIndices( jBranchIdx, jnIdx) = Cmat;
    CIndices( jBranchIdx, jmIdx) = -Cmat;
%     
%     
%     
        case 'Xfm-1'
           zt=0.01*[0.09+1.81i]*((480/Vbase).^2)*Sbase./(500000); 
           yt=1./zt;
        Y2= (1/3) *[ 2*yt, -yt, -yt; -yt,2*yt,-yt; -yt,-yt,2*yt];
   
%               Y2hat1=Y2+abs(yt)*epsilon*(100/3)*(1/5)*eye(3);
%                         Y2hat2=Y2+abs(yt)*(epsilon/2)*(100/3)*eye(3);

                          Y2hat1=Y2+abs(yt)*epsilon*eye(3);
                        Y2hat2=Y2+abs(yt)*(epsilon/2);

            
 
         
          
                 
    
   availablePhases=[1;2;3];
    
    
    % Ytilde

    YtildeNMn(availablePhases,availablePhases)=Y2hat1; 
    YtildeNMm(availablePhases,availablePhases)=Y2hat2;
    YtildeMNn(availablePhases,availablePhases)=Y2hat2;
    YtildeMNm(availablePhases,availablePhases)=Y2hat1;
    
    
   
    
    Ytilde(jnIdx, jnIdx)= Ytilde(jnIdx, jnIdx)+YtildeNMn;
    Ytilde(jnIdx, jmIdx)=-YtildeNMm;
    Ytilde(jmIdx,jmIdx)=Ytilde(jmIdx,jmIdx)+YtildeMNm;
    Ytilde(jmIdx,jnIdx)=-YtildeMNn;
    
    
    % Ybranch
    YbranchII=zeros(3,3);
    YbranchII(availablePhases,availablePhases)=Y2;
    Ybranch(jBranchIdx,jBranchIdx)=YbranchII;
            
            
            
                Cmat=zeros(3);
    Cmat(availablePhases,availablePhases)=eye(length(availablePhases));
    CIndices( jBranchIdx, jnIdx) = Cmat;
    CIndices( jBranchIdx, jmIdx) = -Cmat;
    
        otherwise
            
               
    Zseries=Zcfg(:,:,lineCodes(ii)-720); 
    Yshunt=Ycfg(:,:,lineCodes(ii)-720);
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
connectedPath = gen_path(38, CNodes.', sNodes, rNodes); 
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




name='IEEE-37';
network=v2struct(name,busNames, nBuses,busSet, JbusSet,nBr,JBranchSet, busNamesWithRegs,connectedPath,...
     lineBuses1, lineBuses2,lineLengths,lineCodes,lineStatuses,...
    CNodes,CIndices,Ytilde, Ybranch, Ycap, regulatorVoltageGains, regulatorCurrentGains,...
    regulatorImpedances,regulatorYNMn, regulatorYNMm,regulatorTypes,epsilon,...
    phaseNodesLinIndices, phaseNodes, N, datapath,...
    Sbase,Vbase, Zbase, Ybase, ...
    availableBusIndices, missingBusIndices, ...
    YbusS, YcapS, Ynet, Y, Y_NS, Y_SS, Ybus);


end

