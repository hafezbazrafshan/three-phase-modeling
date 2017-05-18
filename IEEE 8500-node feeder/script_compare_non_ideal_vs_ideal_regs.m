clear all;
clc;

%% Initial values
epsilon=1e-6; % fixing epsilon
v0mags=[1.048; 1.0486; 1.0492];
 v0phases=[-0.7;-120.7;119.3];

 %% Setting up required networks
 % making network with ideal regs
network_idealRegs=setupYbusIEEE8500('ideal',epsilon);
network_idealRegs=computeNoLoadVoltage(network_idealRegs, v0mags,v0phases);
 network_idealRegs  = setupLoadsIEEE8500( network_idealRegs );
network_idealRegs= performZBus(network_idealRegs); 
network_idealRegs= obtainVoltages(network_idealRegs); 

% making network with non-ideal regs
network_nonIdealRegs=setupYbusIEEE8500('non-ideal',epsilon);
network_nonIdealRegs=computeNoLoadVoltage(network_nonIdealRegs, v0mags,v0phases);
 network_nonIdealRegs  = setupLoadsIEEE8500( network_nonIdealRegs );
network_nonIdealRegs= performZBus(network_nonIdealRegs); 
network_nonIdealRegs= obtainVoltages(network_nonIdealRegs); 


% obtaining dss solution
dssNet=obtainDSSVoltages(network_nonIdealRegs);

%% Retrieving outputs

busNamesWithRegs=dssNet.busNamesWithRegs;
busNames=dssNet.busNames;
N=dssNet.N;
regIdx1=find(strcmp(busNamesWithRegs,'_HVMV_Sub_LSB'));
regIdx2=find(strcmp(busNamesWithRegs,'190-8593'));
regIdx3=find(strcmp(busNamesWithRegs,'190-8581'));
regIdx4=find(strcmp(busNamesWithRegs,'190-7361'));
regIdx=[regIdx1,regIdx2,regIdx3,regIdx4];

% regPrimIdx1=find(strcmp(busNamesWithRegs,'regxfmr_HVMV_Sub_LSB'));
% regPrimIdx2=find(strcmp(busNamesWithRegs,'regxfmr_190-8593')); 
% regPrimIdx3=find(strcmp(busNamesWithRegs,'regxfmr_190-8581'));
% regPrimIdx4=find(strcmp(busNamesWithRegs,'regxfmr_190-7361'));
% capIdx1=find(strcmp(busNamesWithRegs,'R20185')); 
% capIdx2=find(strcmp(busNamesWithRegs,'R42246')); 
% capIdx3=find(strcmp(busNamesWithRegs,'R42247')); 
% capIdx4=find(strcmp(busNamesWithRegs,'R18242')); 
% trIdx1=find(strcmp(busNamesWithRegs, 'HVMV_Sub_HSB')); 
% % trIdx2=find(strcmp(busNamesWithRegs,'610'));

connectedPath=dssNet.connectedPath;

% regPrimIdx1InPath=find(connectedPath==regPrimIdx1); 
% regPrimIdx2InPath=find(connectedPath==regPrimIdx2); 
% regPrimIdx3InPath=find(connectedPath==regPrimIdx3); 
% regPrimIdx4InPath=find(connectedPath==regPrimIdx4); 
% capIdx1InPath=find(connectedPath==capIdx1); 
% capIdx2InPath=find(connectedPath==capIdx2); 
% capIdx3InPath=find(connectedPath==capIdx3); 
% capIdx4InPath=find(connectedPath==capIdx4); 
% trIdx1InPath=find(connectedPath==trIdx1); 
% trIdx2InPath=find(connectedPath==trIdx2); 

% xtickVec=sort([regPrimIdx1InPath; regPrimIdx2InPath; regPrimIdx3InPath;...
%     regPrimIdx4InPath; capIdx1InPath; capIdx2InPath; capIdx3InPath; capIdx4InPath; 
%     trIdx1InPath]); 
xtickVec=[1,100:100:length(connectedPath)];


idealVMag=network_idealRegs.solution.resultsVMag;
idealVPhase=network_idealRegs.solution.resultsVPhase;
nonIdealVMag=network_nonIdealRegs.solution.resultsVMag;
nonIdealVPhase=network_nonIdealRegs.solution.resultsVPhase;
dssVMag=dssNet.DSSsolution.DSSVMag;
dssVPhase=dssNet.DSSsolution.DSSVPhase;


%% comparison tables
if exist('Results')~=7
    mkdir Results
end

cd('Results'); 
fileID = fopen('comparison.txt','w');
fmt = '%10s&%.4f&%.4f&%.4f\n';
fprintf(fileID,fmt,'ideal',max(abs(dssVMag-idealVMag)./abs(dssVMag)));
fprintf(fileID,fmt,'non-ideal',max(abs(dssVMag-nonIdealVMag)./abs(dssVMag)));
fclose(fileID);
cd('..')

%% Plots
x0=2;
y0=2;
width=20;
height=5;
figure1=figure('Units','inches',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');
set(figure1, 'Name', 'Voltages');
phase_names = ['a' 'b' 'c'];


for phase=1:3
		
    subplot(3,1,phase)
%     h1=plot([1:length(busNamesWithRegs)].', idealVMag([connectedPath, regIdx] , phase), 'ks'); 
%     hold on
        h2=plot([1:length(busNamesWithRegs)].', nonIdealVMag([connectedPath, regIdx] , phase), 'bo'); 
hold on
       h3=plot(1:length(busNamesWithRegs), dssVMag([connectedPath,regIdx], phase), 'rx'); 
    hold on

%     set([h1,h2,h3],'markers',10)
set([h2,h3],'markers',10);
set(gca,'box','on');
    set(gca,'fontSize',14); 
set(0,'defaulttextinterpreter','latex')
    xlim([0 N+100])
    ylim([0.9 1.2]);
    grid on;
 xticks=[1, 20:20:length(busNamesWithRegs)].';
           xTickLabels=repmat({''},length(xticks),1);
    xTickLabels([1,5:5:length(xticks)])=cellstr(num2str([1,100:100:length(busNamesWithRegs)].'));
		set(gca,'XTick', xticks );
    set(gca,'XTickLabel',xTickLabels);

%     ax = gca;
%     
%     ax.XAxis.FontSize=10;
    
% 	legendTEXT=legend('Z-Bus (ideal VRs)', 'Z-Bus (non-ideal VRs)',  'Benchmark');
legendTEXT=legend('Z-Bus',  'Benchmark');
    set(legendTEXT,'interpreter','Latex'); 
set(legendTEXT,'fontSize',14); 
set(legendTEXT,'fontWeight','Bold');
set(legend,'orientation','Horizontal'); 
set(legend,'location','North');
    ylabel(sprintf('Phase %c ',phase_names(phase)), 'FontWeight','bold');

    if phase==3
        xlabel('Bus numbers');
    end
end


if exist('Figures')~=7
    mkdir Figures;
end

cd('Figures');
print -dpdf IEEE8500magnitudes.pdf
print -depsc2 IEEE8500magnitudes
cd('..');


%% Plot phases


x0=2;
y0=2;
width=20;
height=5;
figure2=figure('Units','inches',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');
set(figure1, 'Name', 'Phases');
phase_names = ['a' 'b' 'c'];

for phase=1:3
		
    subplot(3,1,phase)
%     h1=plot([1:length(busNamesWithRegs)].', idealVPhase([connectedPath, regIdx] , phase), 'ks'); 
%     hold on
        h2=plot([1:length(busNamesWithRegs)].', nonIdealVPhase([connectedPath, regIdx] , phase), 'bo'); 
hold on
       h3=plot(1:length(busNamesWithRegs), dssVPhase([connectedPath,regIdx], phase), 'rx'); 
    hold on
    
%     set([h1,h2,h3],'markers',10)
set([h2,h3],'markers',10);
set(gca,'box','on');
    set(gca,'fontSize',14); 
set(0,'defaulttextinterpreter','latex')

    xlim([0 N+100])
    if phase==1
    ylim([-60 0]);
    elseif phase==2
    ylim([-180 -120]); 
    else
    ylim([60 120]);
    end
    grid on;

	xticks=[1:20:length(busNamesWithRegs)].';
           xTickLabels=repmat({''},length(xticks),1);
    xTickLabels([1,5:5:length(xticks)])=cellstr(num2str([1,100:100:length(busNamesWithRegs)].'));
		set(gca,'XTick', xticks );
    set(gca,'XTickLabel',xTickLabels);
    
    yRange=ylim;
    
        set(gca,'YTick', [yRange(1):20:yRange(2)]);

%         ax = gca;
%     
%     ax.XAxis.FontSize=10;
% 	legendTEXT=legend('Z-Bus (ideal VRs)', 'Z-Bus (non-ideal VRs)',  'Benchmark');
legendTEXT=legend('Z-Bus',  'Benchmark');
    set(legendTEXT,'interpreter','Latex'); 
set(legendTEXT,'fontSize',14); 
set(legendTEXT,'fontWeight','Bold');
set(legend,'orientation','horizontal'); 
set(legend,'location','North');


    ylabel(sprintf('Phase %c ',phase_names(phase)), 'FontWeight','bold');
    
    if phase==3
        xlabel('Bus numbers');
    end
end

cd('Figures'); 
print -dpdf IEEE8500phases.pdf
print -depsc2 IEEE8500phases
cd('..');
