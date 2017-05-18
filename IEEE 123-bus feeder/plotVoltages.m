function  plotVoltages( network )
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

v2struct(network); 
v2struct(solution); 
% v2struct(DSSsolution);
v2struct(Kerstingsolution);

regIdx1=find(strcmp(busNamesWithRegs,'150r'));
regIdx2=find(strcmp(busNamesWithRegs,'9r')); 
regIdx3=find(strcmp(busNamesWithRegs,'25r'));
regIdx4=find(strcmp(busNamesWithRegs,'160r'));
regIdx=[regIdx1,regIdx2,regIdx3,regIdx4];

x0=2;
y0=2;
width=16;
height=5;
figure1=figure('Units','inches',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');
set(figure1, 'Name', 'Voltages');
phase_names = ['a' 'b' 'c'];


for phase=1:3
		
    subplot(3,1,phase)
    plot([1:length(busNamesWithRegs)].', resultsVMag([connectedPath, regIdx] , phase), 'bo'); 
    hold on
%        plot(1:length(busNamesWithRegs), DSSVMag([connectedPath,regIdx], phase), 'rx'); 
%     hold on
   plot(1:length(busNamesWithRegs), KerstingVMag([connectedPath,regIdx], phase), 'rx'); 
    hold on


    xlim([0 N+5]);
    ylim([0.9 1.2]);
    grid on;
		set(gca,'XTick', [1:10:length(busNamesWithRegs)-1, length(busNamesWithRegs)],'YTick', [0.9:0.05:1.1] );
    set(gca,'XTickLabel',[busNames(connectedPath(1:10:end)); busNamesWithRegs(end)]);
    set(gca,'XMinorTick','on');
	legendTEXT=legend('Z-Bus', 'openDSS');
    set(legendTEXT,'interpreter','Latex'); 
set(legendTEXT,'fontSize',14); 
set(legendTEXT,'fontWeight','Bold');
set(legend,'orientation','Horizontal'); 
set(legend,'location','North');
set(gca,'box','on');
    set(gca,'fontSize',14); 
set(0,'defaulttextinterpreter','latex')
    ylabel(sprintf('Phase %c ',phase_names(phase)), 'FontWeight','bold');

    if phase==3
        xlabel('Nodes');
    end
end


if exist('Figures')~=7
    mkdir Figures;
end

end

