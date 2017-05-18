function [ network ] = setupYbusIEEE8500( regulatorTypes, epsilon )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if nargin <2 
    epsilon=1e-5;
end
if nargin<1
    regulatorTypes='non-ideal';
end

datapath=[pwd,'/IEEE-8500 feeder data'];

%%
% importing lines:
[lineNames,lineBuses1,linePhases1,lineBuses2,linePhases2,lineLengths,lineUnits,lineCodes,lineStatuses]=importlines([datapath,'/Lines.csv']);
load([datapath,'/lineDataSet2.mat']);
% importing transformers:
[transformerNames,transformerPhases,transformerBuses1,...
    transformerBuses2,transformerkV_pri,transformerkV_sec,transformerMVA,...
    transformerConn_pri,transformerConn_sec,transformerXHL,transformerRHL] = importtransformers([datapath,'/Transformers.csv']);



% importing regulators:
[regulatorNames,regulatorBuses1,regulatorPhases1,regulatorBuses2,regulatorPhases2,~,~,~,~,~,~]=importregulators([datapath,'/Regulators.csv']);
% importing capacitors
[capacitorNames,capacitorBus,capacitorPhase,capacitorkV,capacitorkvar,capacitornumphases,capacitorconnection] = importcapacitors([datapath,'/Capacitors.csv']) ;


%% Some readjustments:
% some simplification in modeling:
% a. removing open switches
nonSwitchLineIndices=find(strcmp(strtrim(lineStatuses),'open')==0);
lineNames=lineNames(nonSwitchLineIndices);
lineBuses1=lineBuses1(nonSwitchLineIndices);
linePhases1=linePhases1(nonSwitchLineIndices);
lineBuses2=lineBuses2(nonSwitchLineIndices);
linePhases2=linePhases2(nonSwitchLineIndices);
lineLengths=lineLengths(nonSwitchLineIndices);
lineUnits=lineUnits(nonSwitchLineIndices);
lineCodes=lineCodes(nonSwitchLineIndices);
lineStatuses=lineStatuses(nonSwitchLineIndices); 


% %b. assuming connectors as ideal connectors
connectorBusHolder={};
cnt=0;
ii=1;
while ii<=length(lineNames)
if (strcmp( strtrim(lower(lineCodes(ii))),'1ph-connector') )
            
            
        cnt=cnt+1;
        connectorBusHolder(cnt,1)=lineBuses1(ii);
        connectorBusHolder(cnt,2)=lineBuses2(ii);
        idx2=find(strcmp(lineBuses1,lineBuses2(ii)));
        lineBuses1(idx2)=lineBuses1(ii);
        
        
        % remove the connector from lines matrix
        lineNames(ii)=[];
        lineBuses1(ii)=[];
        linePhases1(ii)=[];
        lineBuses2(ii)=[];
        linePhases2(ii)=[];
        lineLengths(ii)=[];
        lineUnits(ii)=[];
        lineCodes(ii)=[];
        lineStatuses(ii)=[];
    else
        ii=ii+1;
    end
    
end






% adding transformers to the line data set:
lineNames=[lineNames; transformerNames];
lineBuses1=[lineBuses1;transformerBuses1];
linePhases1=[linePhases1;'ABC'];
lineBuses2=[lineBuses2; transformerBuses2];
linePhases2=[linePhases2; 'ABC'];
lineLengths=[lineLengths; 0];
lineUnits=[lineUnits;'0'];
lineCodes=[lineCodes;'0'];
lineStatuses=[lineStatuses; 'mainTransformer'];


%% 4. Defining base quantities:
Sbase=str2num(transformerMVA)*1000*1000; % from the substation transformer
Vbase=12.47*1000/sqrt(3); % line to neutral conversion (secondary of the substation transformer)
Zbase=(Vbase.^2)/Sbase;
Ybase=1./Zbase;



%% 3.  Organizing bus names
% collecting bus names:
busNames=unique([lineBuses1;lineBuses2;regulatorBuses1]);
% remove secondary of the regulators from the busNames but keeping them in
% the busNamesWithRegs

%% Putting substation index at the very end:
% substation='regxfmr_HVMV_Sub_LSB';
substation=transformerBuses1;
substationIndex=find(strcmp(busNames,substation));

if substationIndex< length(busNames) % if substation is not the end bus
    busNames=[busNames(1:substationIndex-1); busNames(substationIndex+1:end); busNames(substationIndex)];
end



regidx1=find(strcmp(busNames,'_HVMV_Sub_LSB'));
busNames(regidx1)=[];
regLineIdx1=find(strcmp(lineBuses1,'_HVMV_Sub_LSB'));
lineBuses1{regLineIdx1}='regxfmr_HVMV_Sub_LSB';
lineStatuses{regLineIdx1}='reg1';

regidx2=find(strcmp(busNames,'190-8593'));
busNames(regidx2)=[];
regLineIdx2=find(strcmp(lineBuses1,'190-8593'));
lineBuses1{regLineIdx2}='regxfmr_190-8593';
lineStatuses{regLineIdx2}='reg2';




regidx3=find(strcmp(busNames,'190-8581'));
busNames(regidx3)=[];
regLineIdx3=find(strcmp(lineBuses1,'190-8581'));
lineBuses1{regLineIdx3}='regxfmr_190-8581';
lineStatuses{regLineIdx3}='reg3';





regidx4=find(strcmp(busNames,'190-7361'));
busNames(regidx4)=[];
regLineIdx4=find(strcmp(lineBuses1,'190-7361'));
lineBuses1{regLineIdx4}='regxfmr_190-7361';
lineStatuses{regLineIdx4}='reg4';


busNamesWithRegs=[busNames;'_HVMV_Sub_LSB';'190-8593';'190-8581';'190-7361'];
regulatorVoltageGains=cell(4,1); 
regulatorCurrentGains=cell(4,1); 
regulatorImpedances=cell(4,1); 
regulatorYNMn=cell(4,1);
regulatorYNMm=cell(4,1); 










%% 4. Defining base quantities:
Sbase=str2num(transformerMVA)*1000*1000; % from the substation transformer
Vbase=12.47*1000/sqrt(3); % line to neutral conversion (secondary of the substation transformer)
Zbase=(Vbase.^2)/Sbase;
Ybase=1./Zbase;


%% 5. Create Ytilde and Ybranch
% Ybranch is simplified to be without Yshunt (for testing purposes)
% Ytilde considers Yshunt and Ycaps
nBuses=length(busNames);
busSet=1:nBuses;
JbusSet=1:3*(nBuses);
nBr=length(lineNames);
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


for ii=1:length(lineNames)
    jBranchIdx=(ii-1)*3+1:ii*3;
    bus1=lineBuses1(ii);
    bus2=lineBuses2(ii);
    nIdx=find(strcmp(busNames,bus1));
    mIdx=find(strcmp(busNames,bus2));
    jmIdx=(mIdx-1)*3+1:mIdx*3;
     jnIdx=(nIdx-1)*3+1:nIdx*3;
    AR=eye(3);
    
     busInPhases=linePhases1(ii);
    busOutPhases=linePhases2(ii);
    
    if ~isequal(busInPhases,busOutPhases)
        fprintf('the line is not connected properly\n');
        fprintf('different phases at the two end of the line\n');
    end
    
    
    
    switch char(busInPhases)
        case 'ABC'
            availablePhases=[1,2,3];
        case 'AB'
            availablePhases=[1,2];
        case 'BC'
            availablePhases=[2,3];
        case 'AC'
            availablePhases=[1,3];
        case 'CA'
            availablePhases=[1,3];
        case 'CB'
            availablePhases=[2,3];
        case 'BA'
            availablePhases=[1,2];
        case 'A'
            availablePhases=[1];
        case 'B'
            availablePhases=[2];
        case 'C'
            availablePhases=[3];
    end
    
    
    
   
    
    CNodes(ii,nIdx)=1;
    CNodes(ii,mIdx)=-1;
    
    
   
    
    Cmat=zeros(3);
    Cmat(availablePhases,availablePhases)=eye(length(availablePhases));
    CIndices( jBranchIdx, jnIdx) = Cmat;
    CIndices( jBranchIdx, jmIdx) = -Cmat;
    
        YtildeNMn=zeros(3,3);
    YtildeNMm=zeros(3,3);
    YtildeMNn=zeros(3,3);
    YtildeMNm=zeros(3,3);

    
    
    
    switch lineStatuses{ii}
        case 'reg1' % reg 1
            tapA=2;
            tapB=2;
            tapC=1;
            AR(1,1)=1+0.00625*tapA;
            AR(2,2)=1+0.00625*tapB;
            AR(3,3)=1+0.00625*tapC;
             
             
             
             Av=inv(AR);
             Ai=AR;
%              

    ztReg=(Sbase/10000000)*1*0.0001i;
    Zregulator=[ztReg./AR(1,1), 0, 0; 0 ztReg./AR(2,2) 0; 0 0 ztReg./AR(3,3)]; 
    
             
            
                lineCode=strtrim(lower(lineCodes(ii)));
    lineCodeIdx=find(strcmp(lineDataSet(:,1),lower(lineCode)));
    
    
      [Zseries,Yshunt,nphases]=lineDataSet{lineCodeIdx,[3,4,2]};
    if nphases~=length(char(busInPhases))
        fprintf('mismatch in the number of phases per line and the lineCode');
    end
    
    % length of the line:
    lineLength=lineLengths(ii);
    
    
    
    
    % Ytilde

    YtildeNMn(availablePhases,availablePhases)=( inv(Zseries*lineLength)/Ybase+0*0.5*Yshunt/lineLength/Ybase);
    YtildeNMm(availablePhases,availablePhases)=(inv(Zseries*lineLength)/Ybase);
    YtildeMNn(availablePhases,availablePhases)=(inv(Zseries*lineLength)/Ybase);
    YtildeMNm(availablePhases,availablePhases)=(inv(Zseries*lineLength)/Ybase+0*0.5*Yshunt/lineLength/Ybase);
    
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
regulatorYNMn{4}=Ai*inv(eye(3)+ YtildeNMn *inv(Av)*Zregulator)*(YtildeNMn)*inv(Av);
regulatorYNMm{4}=Ai*inv(eye(3)+ YtildeNMn *inv(Av)*Zregulator)*YtildeNMm;
    
    
    
    % Ybranch
    YbranchII=zeros(3,3);
    YbranchII(availablePhases,availablePhases)=(inv(Zseries*lineLength)/Ybase);
    Ybranch(jBranchIdx,jBranchIdx)=YbranchII;
            
            
            
        case 'reg2' % reg 2
            
            tapA=10;
            tapB=5;
            tapC=2;
            
            
            AR(1,1)=1+0.00625*tapA;
            AR(2,2)=1+0.00625*tapB;
            AR(3,3)=1+0.00625*tapC;
             
             Av=inv(AR);
             Ai=AR;
%              

    ztReg=(Sbase/10000000)*1*0.0001i;
    Zregulator=[ztReg./AR(1,1), 0, 0; 0 ztReg./AR(2,2) 0; 0 0 ztReg./AR(3,3)]; 
    
            
            
                lineCode=strtrim(lower(lineCodes(ii)));
    lineCodeIdx=find(strcmp(lineDataSet(:,1),lower(lineCode)));
    
    
    
      [Zseries,Yshunt,nphases]=lineDataSet{lineCodeIdx,[3,4,2]};
    if nphases~=length(char(busInPhases))
        fprintf('mismatch in the number of phases per line and the lineCode');
    end
    
    % length of the line:
    lineLength=lineLengths(ii);
    
    
    % Ytilde

    YtildeNMn(availablePhases,availablePhases)=(inv(Zseries*lineLength)/Ybase+0*0.5*Yshunt/lineLength/Ybase);
    YtildeNMm(availablePhases,availablePhases)=(inv(Zseries*lineLength)/Ybase);
    YtildeMNn(availablePhases,availablePhases)=(inv(Zseries*lineLength)/Ybase);
    YtildeMNm(availablePhases,availablePhases)=(inv(Zseries*lineLength)/Ybase+0*0.5*Yshunt/lineLength/Ybase);
    
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
regulatorVoltageGains{2}=Av;
regulatorCurrentGains{2}=Ai; 
regulatorImpedances{2}=Zregulator; 
regulatorYNMn{2}=Ai*inv(eye(3)+ YtildeNMn *inv(Av)*Zregulator)*(YtildeNMn)*inv(Av);
regulatorYNMm{2}=Ai*inv(eye(3)+ YtildeNMn *inv(Av)*Zregulator)*YtildeNMm;
    
    
    % Ybranch
    YbranchII=zeros(3,3);
    YbranchII(availablePhases,availablePhases)=(inv(Zseries*lineLength)/Ybase);
    Ybranch(jBranchIdx,jBranchIdx)=YbranchII;
            
            
            
        case 'reg3' % reg 3
            tapA=16;
            tapB=11;
            tapC=1;
            
            
            AR(1,1)=1+0.00625*tapA;
            AR(2,2)=1+0.00625*tapB;
            AR(3,3)=1+0.00625*tapC;
            
                   Av=inv(AR);
             Ai=AR;
%              

    ztReg=(Sbase/10000000)*1*0.0001i;
    Zregulator=[ztReg./AR(1,1), 0, 0; 0 ztReg./AR(2,2) 0; 0 0 ztReg./AR(3,3)]; 
            
            
                lineCode=strtrim(lower(lineCodes(ii)));
    lineCodeIdx=find(strcmp(lineDataSet(:,1),lower(lineCode)));
            
    
    
      [Zseries,Yshunt,nphases]=lineDataSet{lineCodeIdx,[3,4,2]};
    if nphases~=length(char(busInPhases))
        fprintf('mismatch in the number of phases per line and the lineCode');
    end
    
    % length of the line:
    lineLength=lineLengths(ii);
    
    
    % Ytilde

    YtildeNMn(availablePhases,availablePhases)=(inv(Zseries*lineLength)/Ybase+0*0.5*Yshunt/lineLength/Ybase);
    YtildeNMm(availablePhases,availablePhases)=(inv(Zseries*lineLength)/Ybase);
    YtildeMNn(availablePhases,availablePhases)=(inv(Zseries*lineLength)/Ybase);
    YtildeMNm(availablePhases,availablePhases)=(inv(Zseries*lineLength)/Ybase+0*0.5*Yshunt/lineLength/Ybase);
    
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
regulatorVoltageGains{3}=Av;
regulatorCurrentGains{3}=Ai; 
regulatorImpedances{3}=Zregulator; 
regulatorYNMn{3}=Ai*inv(eye(3)+ YtildeNMn *inv(Av)*Zregulator)*(YtildeNMn)*inv(Av);
regulatorYNMm{3}=Ai*inv(eye(3)+ YtildeNMn *inv(Av)*Zregulator)*YtildeNMm;
    
    
    % Ybranch
    YbranchII=zeros(3,3);
    YbranchII(availablePhases,availablePhases)=(inv(Zseries*lineLength)/Ybase);
    Ybranch(jBranchIdx,jBranchIdx)=YbranchII;
            
            
        case 'reg4' % reg 4
            tapA=12;
            tapB=12;
            tapC=5;
            AR(1,1)=1+0.00625*tapA;
            AR(2,2)=1+0.00625*tapB;
            AR(3,3)=1+0.00625*tapC;
            
              
                   Av=inv(AR);
             Ai=AR;
            
    ztReg=(Sbase/10000000)*1*0.0001i;
    Zregulator=[ztReg./AR(1,1), 0, 0; 0 ztReg./AR(2,2) 0; 0 0 ztReg./AR(3,3)]; 
            
                lineCode=strtrim(lower(lineCodes(ii)));
    lineCodeIdx=find(strcmp(lineDataSet(:,1),lower(lineCode)));
    
    
      [Zseries,Yshunt,nphases]=lineDataSet{lineCodeIdx,[3,4,2]};
    if nphases~=length(char(busInPhases))
        fprintf('mismatch in the number of phases per line and the lineCode');
    end
    
    % length of the line:
    lineLength=lineLengths(ii);
    
    
    % Ytilde

    YtildeNMn(availablePhases,availablePhases)=(inv(Zseries*lineLength)/Ybase+0*0.5*Yshunt/lineLength/Ybase);
    YtildeNMm(availablePhases,availablePhases)=(inv(Zseries*lineLength)/Ybase);
    YtildeMNn(availablePhases,availablePhases)=(inv(Zseries*lineLength)/Ybase);
    YtildeMNm(availablePhases,availablePhases)=(inv(Zseries*lineLength)/Ybase+0*0.5*Yshunt/lineLength/Ybase);
    
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
            
            
        case 'mainTransformer'
           zt=0.01*(str2num(transformerRHL)+sqrt(-1)*str2num(transformerXHL))*3;
          yt=1./zt;
          Y1=diag([yt; yt; yt]); 
          Y2=(1/3) *[2* yt, -yt, -yt;
                            -yt, 2*yt, -yt; 
                            -yt, -yt, 2*yt]; 
                        
           Y3=(1/sqrt(3))*[-yt, yt, 0; 
                                    0, -yt, yt; 
                                    yt, 0, -yt]; 
                                
                
                                
             Y2hat=Y2+epsilon*abs(yt)*eye(3);                    
          
                 
    
   
    
    
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
            
            
            
            
    
        otherwise
            
            
               lineCode=strtrim(lower(lineCodes(ii)));
    lineCodeIdx=find(strcmp(lineDataSet(:,1),lower(lineCode)));
       [Zseries,Yshunt,nphases]=lineDataSet{lineCodeIdx,[3,4,2]};
    if nphases~=length(char(busInPhases))
        fprintf('mismatch in the number of phases per line and the lineCode');
    end
    
    % length of the line:
    lineLength=lineLengths(ii);
    
    Yseries=inv(Zseries);
    % Ytilde

    YtildeNMn(availablePhases,availablePhases)=(Yseries/lineLength/Ybase+0*0.5*Yshunt*lineLength/Ybase);
    YtildeNMm(availablePhases,availablePhases)=(Yseries/lineLength/Ybase);
    YtildeMNn(availablePhases,availablePhases)=(Yseries/lineLength/Ybase);
    YtildeMNm(availablePhases,availablePhases)=(Yseries/lineLength/Ybase+0*0.5*Yshunt*lineLength/Ybase);
    
    Ytilde(jnIdx, jnIdx)= Ytilde(jnIdx, jnIdx)+YtildeNMn;
    Ytilde(jnIdx, jmIdx)=-YtildeNMm;
    Ytilde(jmIdx,jmIdx)=Ytilde(jmIdx,jmIdx)+YtildeMNm;
    Ytilde(jmIdx,jnIdx)=-YtildeMNn;
    
    
    % Ybranch
    YbranchII=zeros(3,3);
    YbranchII(availablePhases,availablePhases)=(Yseries/lineLength/Ybase);
    Ybranch(jBranchIdx,jBranchIdx)=YbranchII;
            
            
    end
    

    
    
    
    
    
    
    
  
    
    
end

%% create a connected List of Nodes
connectedPath = gen_path(nBuses, CNodes.', sNodes, rNodes); 

%% 6. Adding Ycap
Ycap=sparse(zeros(size(Ytilde)));
for ii=1:length(capacitorNames)
    nIdx=find(strcmp(busNames,capacitorBus(ii)));
    phase=char(capacitorPhase(ii));
    ycap=sqrt(-1)*capacitorkvar(ii)*1000/Sbase;

    switch phase
        case 'A'
            capIdx1=(nIdx-1)*3+1;
            Ycap(capIdx1,capIdx1)=ycap;
        case 'B'
            capIdx2=(nIdx-1)*3+2;
            Ycap(capIdx2,capIdx2)=ycap;
            
        case 'C'
            capIdx3=(nIdx-1)*3+3;
            Ycap(capIdx3,capIdx3)=ycap;
            
        case 'ABC'
            capIdx4=(nIdx-1)*3+1;
            capIdx5=(nIdx-1)*3+2;
            capIdx6=(nIdx-1)*3+3;
          
            
            Ycap(capIdx4,capIdx4)=ycap/3;
            Ycap(capIdx5,capIdx5)=ycap/3;
            Ycap(capIdx6,capIdx6)=ycap/3;
            
            
            
    end
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




name='IEEE-8500';
network=v2struct(name,busNames, nBuses,busSet, JbusSet,nBr,JBranchSet, busNamesWithRegs,connectedPath,...
     lineBuses1, lineBuses2,lineLengths,lineCodes,lineStatuses,...
    CNodes,CIndices,Ytilde, Ybranch, Ycap, regulatorVoltageGains, regulatorCurrentGains,...
    regulatorImpedances,regulatorYNMn, regulatorYNMm,regulatorTypes,epsilon,...
    phaseNodesLinIndices, phaseNodes, N, datapath,...
    Sbase,Vbase, Zbase, Ybase, ...
    availableBusIndices, missingBusIndices, ...
    YbusS, YcapS, Ynet, Y, Y_NS, Y_SS, Ybus);


end

