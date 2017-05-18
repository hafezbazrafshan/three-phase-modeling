
IEEE8500 = setupYbusIEEE8500('non-ideal',1e-5 );
IEEE8500=computeNoLoadVoltage(IEEE8500, [1.048; 1.0486; 1.0492],[-0.7;-120.7;119.3]);
 IEEE8500  = setupLoadsIEEE8500( IEEE8500 );
IEEE8500= performZBus(IEEE8500); 
IEEE8500= obtainVoltages(IEEE8500); 
IEEE8500=obtainDSSVoltages(IEEE8500);
plotVoltages(IEEE8500 )
