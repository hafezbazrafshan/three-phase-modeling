clear all;
clc;

%% Initial values
epsilon=1e-6; % fixing epsilon
v0mags=[1;1;1];
 v0phases=[0; -120; 120];

 %% Setting up required networks
 % making network with ideal regs
network_idealRegs=setupYbusIEEE37('ideal',epsilon);
network_idealRegs=computeNoLoadVoltage(network_idealRegs, v0mags,v0phases);
 network_idealRegs  = setupLoadsIEEE37( network_idealRegs );
network_idealRegs= performZBus(network_idealRegs); 
network_idealRegs= obtainVoltages(network_idealRegs); 

% making network with non-ideal regs
network_nonIdealRegs=setupYbusIEEE37('non-ideal',epsilon);
network_nonIdealRegs=computeNoLoadVoltage(network_nonIdealRegs, v0mags,v0phases);
 network_nonIdealRegs  = setupLoadsIEEE37( network_nonIdealRegs );
network_nonIdealRegs= performZBus(network_nonIdealRegs); 
network_nonIdealRegs= obtainVoltages(network_nonIdealRegs); 


% obtaining dss solution
dssNet=obtainDSSVoltages(network_nonIdealRegs);

%% Retrieving outputs

busNamesWithRegs=dssNet.busNamesWithRegs;
busNames=dssNet.busNames;
substationIndex=find(strcmp(busNames,'sourcebus'));
busNames{substationIndex}='799s';
N=dssNet.N;
regIdx=find(strcmp(busNamesWithRegs,'799r'));


connectedPath=dssNet.connectedPath;

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
% fprintf(fileID,fmt,'ideal',max(abs(dssVMag([1:8,10:end],:)-idealVMag([1:8,10:end],:))));
% fprintf(fileID,fmt,'non-ideal',max(abs(dssVMag([1:8,10:end],:)-nonIdealVMag([1:8,10:end],:))));
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
    h1=plot([1:length(busNamesWithRegs)].', idealVMag([connectedPath, regIdx] , phase), 'ks'); 
    hold on
        h2=plot([1:length(busNamesWithRegs)].', nonIdealVMag([connectedPath, regIdx] , phase), 'bo'); 
hold on
       h3=plot(1:length(busNamesWithRegs), dssVMag([connectedPath,regIdx], phase), 'rx'); 
    hold on

    set([h1,h2,h3],'markers',12)
set(gca,'box','on');
    set(gca,'fontSize',14); 
set(0,'defaulttextinterpreter','latex')
    xlim([0 N+3]);
    ylim([0.87 1.2]);
    grid on;
    
    xticks=[1:length(busNamesWithRegs)].';
    xTickLabels=cellstr(num2str(xticks))
    xTickLabels(2:2:end)={''};
    set(gca,'XTick', [1:1:length(busNamesWithRegs)],'YTick', [0.9:0.05:1.1] );
    set(gca,'XTickLabel',xTickLabels);
 

    
	legendTEXT=legend('Z-Bus (ideal SVR)', 'Z-Bus (non-ideal SVR)',  'openDSS');
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
print -dpdf IEEE37magnitudes.pdf
print -depsc2 IEEE37magnitudes
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
    h1=plot([1:length(busNamesWithRegs)].', idealVPhase([connectedPath, regIdx] , phase), 'ks'); 
    hold on
        h2=plot([1:length(busNamesWithRegs)].', nonIdealVPhase([connectedPath, regIdx] , phase), 'bo'); 
hold on
       h3=plot(1:length(busNamesWithRegs), dssVPhase([connectedPath,regIdx], phase), 'rx'); 
    hold on
    
    set([h1,h2,h3],'markers',12)
set(gca,'box','on');
    set(gca,'fontSize',14); 
set(0,'defaulttextinterpreter','latex')

    xlim([0 N+3])
    if phase==1
    ylim([-10 10]);
    elseif phase==2
    ylim([-130 -110]); 
    else
    ylim([110 130]);
    end
    grid on;
    
        xticks=[1:length(busNamesWithRegs)].';
    xTickLabels=cellstr(num2str(xticks))
    xTickLabels(2:2:end)={''};

			set(gca,'XTick', xticks );
    set(gca,'XTickLabel',xTickLabels);
%         ax = gca;
%     
%     ax.XAxis.FontSize=10;
	legendTEXT=legend('Z-Bus (ideal SVR)', 'Z-Bus (non-ideal SVR)',  'openDSS');
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
print -dpdf IEEE37phases.pdf
print -depsc2 IEEE37phases
cd('..');
