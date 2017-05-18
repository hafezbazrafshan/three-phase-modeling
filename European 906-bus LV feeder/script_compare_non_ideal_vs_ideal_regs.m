clear all;
clc;

%% Initial values
epsilon=1e-6; % fixing epsilon
v0mags=[1.0496;  1.0496;  1.0497];
 v0phases=[0; -120; 120];

 %% Setting up required networks
 % making network with ideal regs
network_idealRegs=setupYbusIEEELV('ideal',epsilon);
network_idealRegs=computeNoLoadVoltage(network_idealRegs, v0mags,v0phases);
 network_idealRegs  = setupLoadsIEEELV( network_idealRegs );
network_idealRegs= performZBus(network_idealRegs); 
network_idealRegs= obtainVoltages(network_idealRegs); 



% obtaining dss solution
dssNet=obtainDSSVoltages(network_idealRegs);

%% Retrieving outputs

busNamesWithRegs=dssNet.busNamesWithRegs;
busNames=dssNet.busNames;
N=dssNet.N;


connectedPath=dssNet.connectedPath;

idealVMag=network_idealRegs.solution.resultsVMag;
idealVPhase=network_idealRegs.solution.resultsVPhase;
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
fprintf(fileID,fmt,'ideal',max(abs(dssVMag([1:8,10:end],:)-idealVMag([1:8,10:end],:))));
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
    h1=plot([1:length(busNamesWithRegs)].', idealVMag(connectedPath , phase), 'bo'); 
    hold on
       h3=plot(1:length(busNamesWithRegs), dssVMag(connectedPath, phase), 'rx'); 
    hold on

    set([h1,h3],'markers',12)
set(gca,'box','on');
    set(gca,'fontSize',14); 
set(0,'defaulttextinterpreter','latex')
    xlim([0 N+30])
    ylim([0.8 1.2]);
    grid on;
 xticks=[1, 6:6:length(busNamesWithRegs)].';
           xTickLabels=repmat({''},length(xticks),1);
    xTickLabels([1,5:5:length(xticks)])=cellstr(num2str([1,30:30:length(busNamesWithRegs)].'));
		set(gca,'XTick', xticks );
    set(gca,'XTickLabel',xTickLabels);
    set(gca,'YTick',[0.85:0.1:1.1])

    
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
print -dpdf IEEELVmagnitudes.pdf
print -depsc2 IEEELVmagnitudes
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
    h1=plot([1:length(busNamesWithRegs)].', idealVPhase(connectedPath , phase), 'bo'); 
    hold on
       h3=plot(1:length(busNamesWithRegs), dssVPhase(connectedPath, phase), 'rx'); 
    hold on
    
    set([h1,h3],'markers',12)
set(gca,'box','on');
    set(gca,'fontSize',14); 
set(0,'defaulttextinterpreter','latex')

   xlim([0 N+30])
    if phase==1
    ylim([-60 0]);
    elseif phase==2
    ylim([-180 -120]); 
    else
    ylim([60 120]);
    end
    
      xticks=[1,6:6:length(busNamesWithRegs)].';
           xTickLabels=repmat({''},length(xticks),1);
    xTickLabels([1,5:5:length(xticks)])=cellstr(num2str([1,30:30:length(busNamesWithRegs)].'));

			set(gca,'XTick', xticks );
    set(gca,'XTickLabel',xTickLabels);
%         ax = gca;
%     
%     ax.XAxis.FontSize=10;
	legendTEXT=legend('Z-Bus ',  'openDSS');
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
print -dpdf IEEELVphases.pdf
print -depsc2 IEEELVphases
cd('..');
