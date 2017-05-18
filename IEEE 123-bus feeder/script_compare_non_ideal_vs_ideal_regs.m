clear all;
clc;

%% Initial values
epsilon=1e-6; % fixing epsilon
v0mags=[1;1;1];
 v0phases=[0; -120; 120];

 %% Setting up required networks
 % making network with ideal regs
network_idealRegs=setupYbusIEEE123('ideal',epsilon);
network_idealRegs=computeNoLoadVoltage(network_idealRegs, v0mags,v0phases);
 network_idealRegs  = setupLoadsIEEE123( network_idealRegs );
network_idealRegs= performZBus(network_idealRegs); 
network_idealRegs= obtainVoltages(network_idealRegs); 

% making network with non-ideal regs
network_nonIdealRegs=setupYbusIEEE123('non-ideal',epsilon);
network_nonIdealRegs=computeNoLoadVoltage(network_nonIdealRegs, v0mags,v0phases);
 network_nonIdealRegs  = setupLoadsIEEE123( network_nonIdealRegs );
network_nonIdealRegs= performZBus(network_nonIdealRegs); 
network_nonIdealRegs= obtainVoltages(network_nonIdealRegs); 


% obtaining dss solution
kerstingNet=obtainKerstingVoltages(network_nonIdealRegs);

%% Retrieving outputs

busNamesWithRegs=kerstingNet.busNamesWithRegs;
busNames=kerstingNet.busNames;
N=kerstingNet.N;
regIdx1=find(strcmp(busNamesWithRegs,'150r'));
regIdx2=find(strcmp(busNamesWithRegs,'9r')); 
regIdx3=find(strcmp(busNamesWithRegs,'25r'));
regIdx4=find(strcmp(busNamesWithRegs,'160r'));
regIdx=[regIdx1,regIdx2,regIdx3,regIdx4];

% regPrimIdx1=find(strcmp(busNamesWithRegs,'150'));
% regPrimIdx2=find(strcmp(busNamesWithRegs,'9')); 
% regPrimIdx3=find(strcmp(busNamesWithRegs,'25'));
% regPrimIdx4=find(strcmp(busNamesWithRegs,'160'));
% capIdx1=find(strcmp(busNamesWithRegs,'83')); 
% capIdx2=find(strcmp(busNamesWithRegs,'88')); 
% capIdx3=find(strcmp(busNamesWithRegs,'90')); 
% capIdx4=find(strcmp(busNamesWithRegs,'92')); 
% trIdx1=find(strcmp(busNamesWithRegs,'61s')); 
% trIdx2=find(strcmp(busNamesWithRegs,'610'));

connectedPath=kerstingNet.connectedPath;

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
%     trIdx1InPath; trIdx2InPath]); 
xtickVec=[1:3:length(connectedPath)];


idealVMag=network_idealRegs.solution.resultsVMag;
idealVPhase=network_idealRegs.solution.resultsVPhase;
nonIdealVMag=network_nonIdealRegs.solution.resultsVMag;
nonIdealVPhase=network_nonIdealRegs.solution.resultsVPhase;
kerstingVMag=kerstingNet.Kerstingsolution.KerstingVMag;
kerstingVPhase=kerstingNet.Kerstingsolution.KerstingVPhase;


%% comparison tables
if exist('Results')~=7
    mkdir Results
end

cd('Results'); 
fileID = fopen('comparison.txt','w');
fmt = '%10s&%.4f&%.4f&%.4f\n';
fprintf(fileID,fmt,'ideal',max(abs(kerstingVMag-idealVMag)./abs(kerstingVMag)));
fprintf(fileID,fmt,'non-ideal',max(abs(kerstingVMag-nonIdealVMag)./abs(kerstingVMag)));
fclose(fileID);
cd('..')

%% Plots
x0=2;
y0=2;
width=8;
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
       h3=plot(1:length(busNamesWithRegs), kerstingVMag([connectedPath,regIdx], phase), 'rx'); 
    hold on

%     set([h1,h2,h3],'markers',12)
        set([h2,h3],'markers',12)

set(gca,'box','on');
    set(gca,'fontSize',14); 
set(0,'defaulttextinterpreter','latex')
    xlim([-2 N+6]);
    ylim([0.9 1.2]);
    grid on;
    
           xticks=[1:length(busNamesWithRegs)].';
           xTickLabels=repmat({''},length(xticks),1);
    xTickLabels([1,10:10:length(busNamesWithRegs)])=cellstr(num2str([1,10:10:length(busNamesWithRegs)].'));
		set(gca,'XTick', xticks );
    set(gca,'XTickLabel',xTickLabels);

%     ax = gca;
%     
%     ax.XAxis.FontSize=10;
    
% 	legendTEXT=legend('Z-Bus (ideal VRs)', 'Z-Bus (non-ideal VRs)',  'Benchmark');
legendTEXT=legend('Z-Bus', 'Benchmark');
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
print -dpdf IEEE123magnitudes.pdf
print -depsc2 IEEE123magnitudes
cd('..');


%% Plot phases


x0=2;
y0=2;
width=8;
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
       h3=plot(1:length(busNamesWithRegs), kerstingVPhase([connectedPath,regIdx], phase), 'rx'); 
    hold on
    
%     set([h1,h2,h3],'markers',12)
        set([h2,h3],'markers',12)

set(gca,'box','on');
    set(gca,'fontSize',14); 
set(0,'defaulttextinterpreter','latex')

    xlim([-2 N+6])
    if phase==1
    ylim([-10 10]);
    elseif phase==2
    ylim([-130 -110]); 
    else
    ylim([110 130]);
    end
    grid on;

       xticks=[1:length(busNamesWithRegs)].';
           xTickLabels=repmat({''},length(xticks),1);
    xTickLabels([1,10:10:length(busNamesWithRegs)])=cellstr(num2str([1,10:10:length(busNamesWithRegs)].'));
		set(gca,'XTick', xticks );
    set(gca,'XTickLabel',xTickLabels);
%         ax = gca;
%     
%     ax.XAxis.FontSize=10;
% 	legendTEXT=legend('Z-Bus (ideal VRs)', 'Z-Bus (non-ideal VRs)',  'Benchmark');
legendTEXT=legend('Z-Bus', 'Benchmark');
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
print -dpdf IEEE123phases.pdf
print -depsc2 IEEE123phases
cd('..');
