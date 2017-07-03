
clear all;
clc;

IEEE37 = setupYbusIEEE37('non-ideal',1e-6 );
IEEE37=computeNoLoadVoltage(IEEE37, [1;1;1], [0; -120; 120]);
 IEEE37  = setupLoadsIEEE37( IEEE37 );
IEEE37= performZBus(IEEE37); 
IEEE37= obtainVoltages(IEEE37); 
IEEE37=obtainDSSVoltages(IEEE37);
plotVoltages(IEEE37 )
