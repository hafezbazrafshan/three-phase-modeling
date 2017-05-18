
clear all;
clc;

IEEELV = setupYbusIEEELV('non-ideal',1e-5);
IEEELV=computeNoLoadVoltage(IEEELV, [  1.0496;  1.0496;  1.0497], [0;-120; 120]);
 IEEELV  = setupLoadsIEEELV( IEEELV );
IEEELV= performZBus(IEEELV); 
IEEELV= obtainVoltages(IEEELV); 
IEEELV=obtainDSSVoltages(IEEELV);
plotVoltages(IEEELV )
