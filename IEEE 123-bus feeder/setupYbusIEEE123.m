function [ network ] = setupYbusIEEE123( regulatorTypes, epsilon )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if nargin <2 
    epsilon=1e-5;
end
if nargin<1
    regulatorTypes='non-ideal';
end

datapath=[pwd,'/IEEE-123 feeder data'];
[lineBuses1,lineBuses2,lineLengths,lineCodes] = importlines([datapath,'/line data.xls']);
 [capacitorBuses,capA,capB,capC] = importcapacitors([datapath,'/cap data']);

 %% Setting up the list of nodes, giving them a unique order from 1:N, setting the substation to node N+1

lineBuses1=strtrim(cellstr(num2str(lineBuses1)));
lineBuses2=strtrim(cellstr(num2str(lineBuses2)));
lineStatuses=repmat({''},length(lineBuses1),1);
capacitorBuses=strtrim(cellstr(num2str(capacitorBuses)));



%% Switch configuration
lineBuses1=[lineBuses1;'13';'18';'60';'61';'97';'150'];
lineBuses2=[lineBuses2;'152';'135';'160';'61s';'197';'149'];
lineLengths=[lineLengths; 1e-3;1e-3;1e-3;1e-3;1e-3;1e-3];
lineCodes=[lineCodes;1;1;1;1;1;1];% switch assumed of configuration one 
lineStatuses=[lineStatuses; {''}; {''}; {''}; {''}; {''}; {''}];

%% Xfm-1 configuration
lineBuses1=[lineBuses1;'61s']; 
lineBuses2=[lineBuses2;'610']; 
lineLengths=[lineLengths;0]; 
lineCodes=[lineCodes;0];
lineStatuses=[lineStatuses;'Xfm-1']; 

%% Substation transformer:
% lineBuses1=[lineBuses1; 'sourcebus']; 
% lineBuses2=[lineBuses2;'150']; 
% lineCodes=[lineCodes;2]; 
% lineLengths=[lineLengths;0];
% lineStatuses=[lineStatuses; {'subXFM'}]; 



%% Regulator
regidx1=find(and(strcmp( lineBuses1,'150'), strcmp(lineBuses2,'149'))); 
lineStatuses{regidx1}='reg1'; 

regidx2=find(and(strcmp(lineBuses1,'9'), strcmp(lineBuses2,'14'))); 
lineStatuses{regidx2}='reg2'; 

regidx3=find(and(strcmp(lineBuses1,'25'), strcmp(lineBuses2,'26'))); 
lineStatuses{regidx3}='reg3'; 

regidx4=find(and(strcmp(lineBuses1,'160'), strcmp(lineBuses2,'67'))); 
lineStatuses{regidx4}='reg4'; 




%% Base voltages:
Sbase=5000*1000; % from the substation transformer
Vbase=4160/sqrt(3); % line to neutral conversion (secondary of the substation transformer)
Zbase=Vbase^2/Sbase;
Ybase=1./Zbase;


%% Line Codes (must match the ones from openDSS)
% (12 configurations)
Zcfg(:,:,1)= [0.4576+1.0780i   0.1560+0.5017i   0.1535+0.3849i;
    0.1560+0.5017i     0.4666+1.0482i   0.1580+0.4236i;
    0.1535+0.3849i       0.1580+0.4236i    0.4615+1.0651i]/5280;
Ycfg(:,:,1)=1j*(10^(-6))*[5.6765   -1.8319   -0.6982;
    -1.8319        5.9809   -1.1645;
    -0.6982         -1.1645      5.3971]/5280;


Zcfg(:,:,2)= [0.4666+1.0482i   0.1580+0.4236i   0.1560+0.5017i;
    0.1580+0.4236i 0.4615+1.0651i   0.1535+0.3849i;
    0.1560+0.5017i      0.1535+0.3849i        0.4576+1.0780i]/5280;

Ycfg(:,:,2)=1j*(10^(-6))* [5.9809   -1.1645   -1.8319;
    -1.1645       5.3971   -0.6982;
    -1.8319 -0.6982 5.6765]/5280;

Zcfg(:,:,3)=[0.4615+1.0651i   0.1535+0.3849i   0.1580+0.4236i;
    0.1535+0.3849i  0.4576+1.0780i   0.1560+0.5017i;
    0.1580+0.4236i    0.1560+0.5017i  0.4666+1.0482i]/5280;
Ycfg(:,:,3)=1j*(10^(-6))*[5.3971   -0.6982   -1.1645;
    -0.6982   5.6765   -1.8319;
    -1.1645   -1.8319    5.9809]/5280;


Zcfg(:,:,4)=[0.4615+1.0651i   0.1580+0.4236i   0.1535+0.3849i;
    0.1580+0.4236i 0.4666+1.0482i   0.1560+0.5017i;
    0.1535+0.3849i        0.1560+0.5017i      0.4576+1.0780i]/5280;
Ycfg(:,:,4)=1j*(10^(-6))*[5.3971   -1.1645   -0.6982;
    -1.1645    5.9809   -1.8319;
    -0.6982  -1.8319         5.6765]/5280;


Zcfg(:,:,5)=[0.4666+1.0482i   0.1560+0.5017i   0.1580+0.4236i;
    0.1560+0.5017i  0.4576+1.0780i   0.1535+0.3849i;
    0.1580+0.4236i         0.1535+0.3849i         0.4615+1.0651i]/5280;
Ycfg(:,:,5)=1j*(10^(-6))*[5.9809   -1.8319   -1.1645
    -1.8319      5.6765   -0.6982
    -1.1645      -0.6982          5.3971]/5280;


Zcfg(:,:,6)=[ 0.4576+1.0780i   0.1535+0.3849i   0.1560+0.5017i
    0.1535+0.3849i 0.4615+1.0651i   0.1580+0.4236i
    0.1560+0.5017i       0.1580+0.4236i       0.4666+1.0482i]/5280;
Ycfg(:,:,6)=1j*(10^(-6))*[5.6765   -0.6982   -1.8319
    -0.6982    5.3971   -1.1645
    -1.8319     -1.1645       5.9809]/5280;



Zcfg(:,:,7)=[0.4576+1.0780i      0                   0.1535+0.3849i;
    0                             0                   0;
    0.1535+0.3849i       0         0.4615+1.0651i]/5280;

Ycfg(:,:,7)=1j*(10^(-6))*[5.1154    0   -1.0549;
    0    0 0;
    -1.0549 0         5.1704]/5280;

Zcfg(:,:,8)=[ 0.4576+1.0780i   0.1535+0.3849i   0.0000;
    0.1535+0.3849i 0.4615+1.0651i     0.0000 ;
    0.0000             0.0000  0.0000]/5280;
Ycfg(:,:,8)=1j*(10^(-6))*[5.1154   -1.0549    0.0000
    -1.0549   5.1704    0.0000
    0.0000          0.0000  0.0000]/5280;

Zcfg(:,:,9)=[1.3292+1.3475i   0   0;
    0 0 0;
    0 0 0] /5280;
Ycfg(:,:,9)=1j*(10^(-6))*[4.5193    0.0000    0.0000
    0  0.0000    0.0000
    0      0      0.0000]/5280;


Zcfg(:,:,10)=[0 0 0;
    0 1.3292+1.3475i    0;
    0  0 0]/5280;
Ycfg(:,:,10)=1j*(10^(-6))*[0 0 0;
    0 4.5193 0;
    0 0 0]/5280;



Zcfg(:,:,11)=[ 0   0   0;
    0  0   0;
    0 0   1.3292+1.3475i]/5280;
Ycfg(:,:,11)=1j*(10^(-6))*[ 0  0   0;
    0   0   0;
    0    0  4.5193]/5280;

Zcfg(:,:,12)=[ 1.5209+0.7521i   0.5198+0.2775i   0.4924+0.2157i;
    0.5198+0.2775i 1.5329+0.7162i   0.5198+0.2775i;
    0.4924+0.2157i        0.5198+0.2775i         1.5209+0.7521i]/5280;
Ycfg(:,:,12)=1j*(10^(-6))*[67.2242  0    0;
    0 67.2242    0;
    0       0    67.2242]/5280;



%% 3.  Organizing bus names
% collecting bus names:
busNames=unique([lineBuses1;lineBuses2]);

%% Putting substation index at the very end:
substation='150';
% substation='sourcebus';
substationIndex=find(strcmp(busNames,substation));

if substationIndex< length(busNames) % if substation is not the end bus
    busNames=[busNames(1:substationIndex-1); busNames(substationIndex+1:end); busNames(substationIndex)];
end
busNamesWithRegs=[busNames;'150r';'9r';'25r';'160r'];  
regulatorVoltageGains=cell(4,1); 
regulatorCurrentGains=cell(4,1); 
regulatorImpedances=cell(4,1); 
regulatorYNMn=cell(4,1);
regulatorYNMm=cell(4,1); 



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
    AR=zeros(3);
    
   
   
    
    CNodes(ii,nIdx)=1;
    CNodes(ii,mIdx)=-1;
    
    
   
    
    Cmat=zeros(3);
   
        YtildeNMn=zeros(3,3);
    YtildeNMm=zeros(3,3);
    YtildeMNn=zeros(3,3);
    YtildeMNm=zeros(3,3);

    
    
    
    switch lineStatuses{ii}
        case 'reg1' % reg 1
            tapA=7;
            tapB=7;
            tapC=7;
            arA=1-0.00625*tapA;
            arB=1-0.00625*tapB;
            arC=1-0.00625*tapC;
            
            
              
             Av=[arA,0,0;
                    0, arB, 0 ; 
                    0, 0,arC];
             Ai=inv(Av);
%              

    ztReg=(Sbase/5000000)*3*0.00001i;
    Zregulator=[ztReg./arA, 0, 0; 0 ztReg./arB 0; 0 0 ztReg./arC]; 
    

            
            
            
           
            

    
    Zseries=Zcfg(:,:,lineCodes(ii)); 
    Yshunt=Ycfg(:,:,lineCodes(ii));
     availablePhases=find(any(Zseries));
  
     Zseries=Zseries(availablePhases,availablePhases); 
     Yshunt=Yshunt(availablePhases,availablePhases); 
    
    % length of the line:
    lineLength=lineLengths(ii);
    
    
    YtildeNMn(availablePhases,availablePhases)=( inv(Zseries*lineLength)/Ybase+0.5*Yshunt*lineLength/Ybase);
    YtildeNMm(availablePhases,availablePhases)=(inv(Zseries*lineLength)/Ybase);
    YtildeMNn(availablePhases,availablePhases)=(inv(Zseries*lineLength)/Ybase);
    YtildeMNm(availablePhases,availablePhases)=(inv(Zseries*lineLength)/Ybase+0.5*Yshunt*lineLength/Ybase);
    
    % Ytilde
    if strcmp(regulatorTypes,'ideal')

    Ytilde(jnIdx, jnIdx)= Ytilde(jnIdx, jnIdx)+Ai*(YtildeNMn)*inv(Av);
    Ytilde(jnIdx, jmIdx)=-Ai*YtildeNMm;
    Ytilde(jmIdx,jmIdx)=Ytilde(jmIdx,jmIdx)+YtildeMNm;
    Ytilde(jmIdx,jnIdx)=-YtildeMNn*inv(Av);
 
    else
    Ytilde(jnIdx, jnIdx)= Ytilde(jnIdx, jnIdx)+Ai*inv(eye(3)+ YtildeNMn *inv(Av)*Zregulator)*(YtildeNMn)*inv(Av);
    Ytilde(jnIdx, jmIdx)=-Ai*inv(eye(3)+ YtildeNMn *inv(Av)*Zregulator)*YtildeNMm;
    Ytilde(jmIdx,jmIdx)=Ytilde(jmIdx,jmIdx)+YtildeMNm - YtildeMNm*Zregulator*Ai*inv(eye(3)+ YtildeNMn *inv(Av)*Zregulator)*YtildeNMm;
    Ytilde(jmIdx,jnIdx)=-(YtildeMNn*inv(Av)-YtildeMNn*Zregulator*Ai*inv(eye(3)+ YtildeNMn *inv(Av)*Zregulator)*(YtildeNMn)*inv(Av));
    end
regulatorVoltageGains{1}=Av;
regulatorCurrentGains{1}=Ai; 
regulatorImpedances{1}=Zregulator; 
regulatorYNMn{1}=Ai*inv(eye(3)+ YtildeNMn *inv(Av)*Zregulator)*(YtildeNMn)*inv(Av);
regulatorYNMm{1}=Ai*inv(eye(3)+ YtildeNMn *inv(Av)*Zregulator)*YtildeNMm;
    
    
    % Ybranch
    YbranchII=zeros(3,3);
    YbranchII(availablePhases,availablePhases)=(inv(Zseries*lineLength)/Ybase);
    Ybranch(jBranchIdx,jBranchIdx)=YbranchII;
            
            
        
                Cmat=zeros(3);
    Cmat(availablePhases,availablePhases)=eye(length(availablePhases));
    CIndices( jBranchIdx, jnIdx) = Cmat;
    CIndices( jBranchIdx, jmIdx) = -Cmat;
    
        case 'reg2'
              tapA=-1;
               arA=1-0.00625*tapA;
            
            
            
              
             Av=[arA, 0, 0; 0 0 0; 0 0 0]; 
            Ai = [1./arA, 0, 0; 0 0 0; 0 0 0];
        invAv=Ai;

    ztReg=(Sbase/2000000)*1*0.0001i;
    Zregulator=[ztReg./arA, 0, 0; 0 0 0; 0 0 0]; 
    

            
            
            
           
            

    
    Zseries=Zcfg(:,:,lineCodes(ii)); 
    Yshunt=Ycfg(:,:,lineCodes(ii));
     availablePhases=find(any(Zseries));
  
     Zseries=Zseries(availablePhases,availablePhases); 
     Yshunt=Yshunt(availablePhases,availablePhases); 
    
    % length of the line:
    lineLength=lineLengths(ii);
    
    
    YtildeNMn(availablePhases,availablePhases)=( inv(Zseries*lineLength)/Ybase+0.5*Yshunt*lineLength/Ybase);
    YtildeNMm(availablePhases,availablePhases)=(inv(Zseries*lineLength)/Ybase);
    YtildeMNn(availablePhases,availablePhases)=(inv(Zseries*lineLength)/Ybase);
    YtildeMNm(availablePhases,availablePhases)=(inv(Zseries*lineLength)/Ybase+0.5*Yshunt*lineLength/Ybase);
    
    % Ytilde
    if strcmp(regulatorTypes,'ideal')

    Ytilde(jnIdx, jnIdx)= Ytilde(jnIdx, jnIdx)+Ai*(YtildeNMn)*invAv;
    Ytilde(jnIdx, jmIdx)=-Ai*YtildeNMm;
    Ytilde(jmIdx,jmIdx)=Ytilde(jmIdx,jmIdx)+YtildeMNm;
    Ytilde(jmIdx,jnIdx)=-YtildeMNn*invAv;
 
    else
    Ytilde(jnIdx, jnIdx)= Ytilde(jnIdx, jnIdx)+Ai*inv(eye(3)+ YtildeNMn *invAv*Zregulator)*(YtildeNMn)*invAv;
    Ytilde(jnIdx, jmIdx)=-Ai*inv(eye(3)+ YtildeNMn *invAv*Zregulator)*YtildeNMm;
    Ytilde(jmIdx,jmIdx)=Ytilde(jmIdx,jmIdx)+YtildeMNm - YtildeMNm*Zregulator*Ai*inv(eye(3)+ YtildeNMn *invAv*Zregulator)*YtildeNMm;
    Ytilde(jmIdx,jnIdx)=-(YtildeMNn*invAv-YtildeMNn*Zregulator*Ai*inv(eye(3)+ YtildeNMn *invAv*Zregulator)*(YtildeNMn)*invAv);
    end
regulatorVoltageGains{2}=Av;
regulatorCurrentGains{2}=Ai; 
regulatorImpedances{2}=Zregulator; 
regulatorYNMn{2}=Ai*inv(eye(3)+ YtildeNMn *invAv*Zregulator)*(YtildeNMn)*invAv;
regulatorYNMm{2}=Ai*inv(eye(3)+ YtildeNMn *invAv*Zregulator)*YtildeNMm;
    
    % Ybranch
    YbranchII=zeros(3,3);
    YbranchII(availablePhases,availablePhases)=(inv(Zseries*lineLength)/Ybase);
    Ybranch(jBranchIdx,jBranchIdx)=YbranchII;
            
            
        
                Cmat=zeros(3);
    Cmat(availablePhases,availablePhases)=eye(length(availablePhases));
    CIndices( jBranchIdx, jnIdx) = Cmat;
    CIndices( jBranchIdx, jmIdx) = -Cmat;
        case 'reg3'
          tapA=0;
         tapC=-1;
             arA=1-0.00625*tapA;
             arC=1-0.00625*tapC;
            
            
            
              
             Av=[arA, 0, 0; 0 0 0; 0 0 arC]; 
            Ai = [1./arA, 0, 0; 0 0 0; 0 0 1./arC];
        invAv=Ai;

    ztReg=(Sbase/2000000)*1*0.0001i;
    Zregulator=[ztReg./arA, 0, 0; 0 0 0; 0 0 ztReg./arC]; 
    

            
            
            
           
            

    
    Zseries=Zcfg(:,:,lineCodes(ii)); 
    Yshunt=Ycfg(:,:,lineCodes(ii));
     availablePhases=find(any(Zseries));
  
     Zseries=Zseries(availablePhases,availablePhases); 
     Yshunt=Yshunt(availablePhases,availablePhases); 
    
    % length of the line:
    lineLength=lineLengths(ii);
    
    
    YtildeNMn(availablePhases,availablePhases)=( inv(Zseries*lineLength)/Ybase+0.5*Yshunt*lineLength/Ybase);
    YtildeNMm(availablePhases,availablePhases)=(inv(Zseries*lineLength)/Ybase);
    YtildeMNn(availablePhases,availablePhases)=(inv(Zseries*lineLength)/Ybase);
    YtildeMNm(availablePhases,availablePhases)=(inv(Zseries*lineLength)/Ybase+0.5*Yshunt*lineLength/Ybase);
    
    % Ytilde
    if strcmp(regulatorTypes,'ideal')

    Ytilde(jnIdx, jnIdx)= Ytilde(jnIdx, jnIdx)+Ai*(YtildeNMn)*invAv;
    Ytilde(jnIdx, jmIdx)=-Ai*YtildeNMm;
    Ytilde(jmIdx,jmIdx)=Ytilde(jmIdx,jmIdx)+YtildeMNm;
    Ytilde(jmIdx,jnIdx)=-YtildeMNn*invAv;
 
    else
    Ytilde(jnIdx, jnIdx)= Ytilde(jnIdx, jnIdx)+Ai*inv(eye(3)+ YtildeNMn *invAv*Zregulator)*(YtildeNMn)*invAv;
    Ytilde(jnIdx, jmIdx)=-Ai*inv(eye(3)+ YtildeNMn *invAv*Zregulator)*YtildeNMm;
    Ytilde(jmIdx,jmIdx)=Ytilde(jmIdx,jmIdx)+YtildeMNm - YtildeMNm*Zregulator*Ai*inv(eye(3)+ YtildeNMn *invAv*Zregulator)*YtildeNMm;
    Ytilde(jmIdx,jnIdx)=-(YtildeMNn*invAv-YtildeMNn*Zregulator*Ai*inv(eye(3)+ YtildeNMn *invAv*Zregulator)*(YtildeNMn)*invAv);
    end
regulatorVoltageGains{3}=Av;
regulatorCurrentGains{3}=Ai; 
regulatorImpedances{3}=Zregulator; 
regulatorYNMn{3}=Ai*inv(eye(3)+ YtildeNMn *invAv*Zregulator)*(YtildeNMn)*invAv;
regulatorYNMm{3}=Ai*inv(eye(3)+ YtildeNMn *invAv*Zregulator)*YtildeNMm;
    
    
    % Ybranch
    YbranchII=zeros(3,3);
    YbranchII(availablePhases,availablePhases)=(inv(Zseries*lineLength)/Ybase);
    Ybranch(jBranchIdx,jBranchIdx)=YbranchII;
            
            
        
                Cmat=zeros(3);
    Cmat(availablePhases,availablePhases)=eye(length(availablePhases));
    CIndices( jBranchIdx, jnIdx) = Cmat;
    CIndices( jBranchIdx, jmIdx) = -Cmat;
            
        case 'reg4'
                    tapA=8;
            tapB=1;
            tapC=5;
            arA=1-0.00625*tapA;
            arB=1-0.00625*tapB;
            arC=1-0.00625*tapC;
            
            
              
             Av=[arA,0,0;
                    0, arB, 0 ; 
                    0, 0,arC];
             Ai=inv(Av);
%              

    ztReg=(Sbase/2000000)*1*0.0001i;
    Zregulator=[ztReg./arA, 0, 0; 0 ztReg./arB 0; 0 0 ztReg./arC]; 
    

            
            
            
           
            

    
    Zseries=Zcfg(:,:,lineCodes(ii)); 
    Yshunt=Ycfg(:,:,lineCodes(ii));
     availablePhases=find(any(Zseries));
  
     Zseries=Zseries(availablePhases,availablePhases); 
     Yshunt=Yshunt(availablePhases,availablePhases); 
    
    % length of the line:
    lineLength=lineLengths(ii);
    
    
    YtildeNMn(availablePhases,availablePhases)=( inv(Zseries*lineLength)/Ybase+0.5*Yshunt*lineLength/Ybase);
    YtildeNMm(availablePhases,availablePhases)=(inv(Zseries*lineLength)/Ybase);
    YtildeMNn(availablePhases,availablePhases)=(inv(Zseries*lineLength)/Ybase);
    YtildeMNm(availablePhases,availablePhases)=(inv(Zseries*lineLength)/Ybase+0.5*Yshunt*lineLength/Ybase);
    
    % Ytilde
    if strcmp(regulatorTypes,'ideal')

    Ytilde(jnIdx, jnIdx)= Ytilde(jnIdx, jnIdx)+Ai*(YtildeNMn)*inv(Av);
    Ytilde(jnIdx, jmIdx)=-Ai*YtildeNMm;
    Ytilde(jmIdx,jmIdx)=Ytilde(jmIdx,jmIdx)+YtildeMNm;
    Ytilde(jmIdx,jnIdx)=-YtildeMNn*inv(Av);
 
    else
    Ytilde(jnIdx, jnIdx)= Ytilde(jnIdx, jnIdx)+Ai*inv(eye(3)+ YtildeNMn *inv(Av)*Zregulator)*(YtildeNMn)*inv(Av);
    Ytilde(jnIdx, jmIdx)=-Ai*inv(eye(3)+ YtildeNMn *inv(Av)*Zregulator)*YtildeNMm;
    Ytilde(jmIdx,jmIdx)=Ytilde(jmIdx,jmIdx)+YtildeMNm - YtildeMNm*Zregulator*Ai*inv(eye(3)+ YtildeNMn *inv(Av)*Zregulator)*YtildeNMm;
    Ytilde(jmIdx,jnIdx)=-(YtildeMNn*inv(Av)-YtildeMNn*Zregulator*Ai*inv(eye(3)+ YtildeNMn *inv(Av)*Zregulator)*(YtildeNMn)*inv(Av));
    end
regulatorVoltageGains{4}=Av;
regulatorCurrentGains{4}=Ai; 
regulatorImpedances{4}=Zregulator; 
regulatorYNMn{4}=Ai*inv(eye(3)+ YtildeNMn *inv(Av)*Zregulator)*(YtildeNMn)*inv(Av);
regulatorYNMm{4}=Ai*inv(eye(3)+ YtildeNMn *inv(Av)*Zregulator)*YtildeNMm;
    
    
    
    % Ybranch
    YbranchII=zeros(3,3);
    YbranchII(availablePhases,availablePhases)=(inv(Zseries*lineLength)/Ybase);
    Ybranch(jBranchIdx,jBranchIdx)=YbranchII;
            
            
        
                Cmat=zeros(3);
    Cmat(availablePhases,availablePhases)=eye(length(availablePhases));
    CIndices( jBranchIdx, jnIdx) = Cmat;
    CIndices( jBranchIdx, jmIdx) = -Cmat;
            
    
    
%         case 'subXFM'
%            zt=0.0001*0.01*(1+sqrt(-1)*8)*3;
%           yt=1./zt;
%           Y1=diag([yt; yt; yt]); 
%           Y2=(1/3) *[2* yt, -yt, -yt;
%                             -yt, 2*yt, -yt; 
%                             -yt, -yt, 2*yt]; 
%                         
%            Y3=(1/sqrt(3))*[-yt, yt, 0; 
%                                     0, -yt, yt; 
%                                     yt, 0, -yt]; 
%                                 
%                 
%                           
%              Y2hat=Y2+epsilon*eye(3);                    
%           
%                  
%     
%       availablePhases=[1;2;3];
%     
%     
%     % Ytilde
% 
%     YtildeNMn(availablePhases,availablePhases)=Y2hat; 
%     YtildeNMm(availablePhases,availablePhases)=-Y3;
%     YtildeMNn(availablePhases,availablePhases)=-Y3.';
%     YtildeMNm(availablePhases,availablePhases)=Y1;
%     
%     
%    
%     
%     Ytilde(jnIdx, jnIdx)= Ytilde(jnIdx, jnIdx)+YtildeNMn;
%     Ytilde(jnIdx, jmIdx)=-YtildeNMm;
%     Ytilde(jmIdx,jmIdx)=Ytilde(jmIdx,jmIdx)+YtildeMNm;
%     Ytilde(jmIdx,jnIdx)=-YtildeMNn;
%     
%     
%     % Ybranch
%     YbranchII=zeros(3,3);
%     YbranchII(availablePhases,availablePhases)=Y2hat;
%     Ybranch(jBranchIdx,jBranchIdx)=YbranchII;
%     
    
            
        case 'Xfm-1'
           zt=0.01*(1.27+2.72i)*(Sbase/150000)*3;
          yt=1./zt;
          Y1=diag([yt;yt;yt]); 
          Y2=(1/3) *[ 2*yt, -yt, -yt; -yt,2*yt,-yt; -yt,-yt,2*yt];
                        
%           
           Y2hat1=Y2+abs(yt)*epsilon*eye(3);
                        Y2hat2=Y2+abs(yt)*(epsilon/2)*eye(3);

    


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
            
               
    Zseries=Zcfg(:,:,lineCodes(ii)); 
    Yshunt=Ycfg(:,:,lineCodes(ii));
     availablePhases=find(any(Zseries));
  
 
  
     Zseries=Zseries(availablePhases,availablePhases); 
     Yshunt=Yshunt(availablePhases,availablePhases); 
    
    % length of the line:
    lineLength=lineLengths(ii);


    
    Yseries=inv(Zseries);
    % Ytilde

    YtildeNMn(availablePhases,availablePhases)=(Yseries/lineLength/Ybase+0.5*Yshunt*lineLength/Ybase);
    YtildeNMm(availablePhases,availablePhases)=(Yseries/lineLength/Ybase);
    YtildeMNn(availablePhases,availablePhases)=(Yseries/lineLength/Ybase);
    YtildeMNm(availablePhases,availablePhases)=(Yseries/lineLength/Ybase+0.5*Yshunt*lineLength/Ybase);
    
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




%% 6. Adding Ycap
Ycap=sparse(zeros(size(Ytilde)));
for ii=1:length(capacitorBuses)
    nIdx=find(strcmp(busNames,capacitorBuses(ii)));
    capIdx=(nIdx-1)*3+1:nIdx*3;
        availablePhases=find(~isnan([capA(ii), capB(ii), capC(ii)]));
   ycap=zeros(3,3);
    ycap=diag([capA(ii), capB(ii), capC(ii)])*1i*1000/Sbase;
    ycap(isnan(ycap))=0; 
    Ycap(capIdx,capIdx)=ycap;
   
end






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




name='IEEE-123';
network=v2struct(name,busNames, nBuses,busSet, JbusSet,nBr,JBranchSet, busNamesWithRegs,connectedPath,...
     lineBuses1, lineBuses2,lineLengths,lineCodes,lineStatuses,...
    CNodes,CIndices,Ytilde, Ybranch, Ycap, regulatorVoltageGains, regulatorCurrentGains,...
    regulatorImpedances,regulatorYNMn, regulatorYNMm,regulatorTypes,epsilon,...
    phaseNodesLinIndices, phaseNodes, N, datapath,...
    Sbase,Vbase, Zbase, Ybase, ...
    availableBusIndices, missingBusIndices, ...
    YbusS, YcapS, Ynet, Y, Y_NS, Y_SS, Ybus);


end

