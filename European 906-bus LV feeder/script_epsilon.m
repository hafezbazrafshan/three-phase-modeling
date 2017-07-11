clear all;
clc;

%% Initial values
v0mags=[  1.0496;  1.0496;  1.0497];
 v0phases=[0; -120; 120];
 
 % making network with non-ideal regs
network_nonIdealRegs=setupYbusIEEELV('non-ideal',1e-1);
network_nonIdealRegs=computeNoLoadVoltage(network_nonIdealRegs, v0mags,v0phases);
 network_nonIdealRegs  = setupLoadsIEEELV( network_nonIdealRegs );
network_nonIdealRegs= performZBus(network_nonIdealRegs); 
network_nonIdealRegs= obtainVoltages(network_nonIdealRegs);
firstIterationVoltages=network_nonIdealRegs.solution.resultsVMag;


epsilonLogarithm=[2:10].'; % fixing epsilon
epsilonVec=10.^(-epsilonLogarithm);
a_err=inf(size(epsilonVec)); 
b_err=inf(size(epsilonVec));
c_err=inf(size(epsilonVec));

for ii=1:length(epsilonVec)
    epsilon=epsilonVec(ii);
v0mags=[1;1;1];
 v0phases=[0; -120; 120];
 
 % making network with non-ideal regs
network_nonIdealRegs=setupYbusIEEELV('non-ideal',epsilon);
network_nonIdealRegs=computeNoLoadVoltage(network_nonIdealRegs, v0mags,v0phases);
 network_nonIdealRegs  = setupLoadsIEEELV( network_nonIdealRegs );
network_nonIdealRegs= performZBus(network_nonIdealRegs); 
network_nonIdealRegs= obtainVoltages(network_nonIdealRegs);

% obtaining dss solution
dssNet=obtainDSSVoltages(network_nonIdealRegs);
a_err(ii)=max(abs(network_nonIdealRegs.solution.resultsVMag(:,1)-firstIterationVoltages(:,1)));
b_err(ii)=max(abs(network_nonIdealRegs.solution.resultsVMag(:,2)-firstIterationVoltages(:,2)));
c_err(ii)=max(abs(network_nonIdealRegs.solution.resultsVMag(:,3)-firstIterationVoltages(:,3)));
firstIterationVoltages=network_nonIdealRegs.solution.resultsVMag;

end




%% Plots
x0=2;
y0=2;
width=8;
height=6;
figure1=figure('Units','inches',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');
set(figure1, 'Name', 'Errors');

h1=plot(epsilonLogarithm,a_err,'b--o');
hold on
h2=plot(epsilonLogarithm,b_err,'r--s'); 
hold on
h3=plot(epsilonLogarithm,c_err,'g--x'); 

set([h1,h2,h3],'markerSize',10,'lineWidth',2);

currentAxes=gca;
set(0,'defaulttextinterpreter','latex')
set(currentAxes,'XTick',epsilonLogarithm);
XTickLabelSet=epsilonLogarithm;
set(currentAxes,'XtickLabels',[]);
xlim([epsilonLogarithm(1)-0.2 epsilonLogarithm(end)+0.2]);
ylim([-0.001 max([a_err;b_err;c_err])+0.001]);
VerticalOffset = 0.003;
ax=axis;
xTicks=get(currentAxes,'XTick');
for i = 1:length(XTickLabelSet)
%Create text box and set appropriate properties
     text(xTicks(i),ax(3) - VerticalOffset, ['$\mathbf{10^{' num2str( -XTickLabelSet(i)) '}}$'],...
         'horizontalAlignment','Center','interpreter', 'latex','fontSize',20,'fontName','Times New Roman');   
end
set(currentAxes,'fontSize',16); 
grid on; 
set(currentAxes, 'box','on'); 
t = text(xTicks(5)-0.4,ax(3)-0.007,' $\epsilon/|y_t|$');
s = t.FontSize;
t.FontSize = 20;
% xlabel(' $\epsilon$','FontName','Times New Roman','fontSize',18); 
set(currentAxes,'FontName','Times New Roman');
 ylabel('Max. difference in voltage magnitudes (pu)'); 
legendTEXT=legend([h1, h2, h3], 'phase $a$', 'phase $b$','phase $c$'); 
set(legendTEXT,'interpreter','Latex'); 
set(legendTEXT,'fontSize',20); 
set(legendTEXT,'fontname','Times New Roman');
set(legendTEXT,'fontWeight','Bold');
set(legend,'orientation','Vertical'); 
set(legend,'location','NorthEast');
set(currentAxes,'box','on');
% 
cd('Figures'); 
print -dpdf IEEELVEpsilon.pdf
print -depsc2 IEEELVEpsilon

cd('..'); 


