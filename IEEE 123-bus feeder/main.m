
clear all;
clc;

IEEE123 = setupYbusIEEE123('non-ideal',1e-6 );
IEEE123=computeNoLoadVoltage(IEEE123, [1;1;1], [0; -120; 120]);
 IEEE123  = setupLoadsIEEE123( IEEE123 );
IEEE123= performZBus(IEEE123); 
IEEE123= obtainVoltages(IEEE123); 
IEEE123=obtainDSSVoltages(IEEE123);
IEEE123=obtainKerstingVoltages(IEEE123);
plotVoltages(IEEE123 )
