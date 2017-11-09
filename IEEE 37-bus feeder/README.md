# Code description
An explanation is provided here for the most important scripts of this folder to simplify usage.   Some additional scripts are also within the folder that replicate the results of the cited paper. 

The modeling presented here takes into account mixed wye *and* delta ZIP loads, even though the model in the cited paper only assumes wye *or* delta ZIP loads.



## `main.m`
This script is the starting point. It contains all the necessary steps from modeling, to computation of voltages as well as comparisons between solutions. 

## `setupYbusIEEE37.m`
This function uses the data files of the folder [IEEE-37 feeder data](https://github.com/hafezbazrafshan/three-phase-modeling/tree/master/IEEE%2037-bus%20feeder/IEEE-37%20feeder%20data)  provided from [PES](https://ewh.ieee.org/soc/pes/dsacom/testfeeders/) to create the bus admittance matrix.  
### List of outputs
* `network` a structure with the following fields
  * `Sbase` Network power base
  * `Vbase` Network voltage base
  * `Zbase` Network impedance base
  * `Ynet` Bus admittance matrix including the slack bus 
  * `Y` Bus admittance matrix excluding the slack bus
  * `Y_NS` The portion of bus admittance matrix relating the network to the slack bus
  * `Y_SS` The portion of bus admittance matrix relating only to the slack bus

The output of this function should typically be input to the function `computeNoLoadVoltage.m`.

### List of inputs
 * `regulatorTypes` (optional) a string set to either 'ideal' or 'non-ideal' determining whether step-voltage regulators should be modeled as ideal or lossy.  If not provided, 'non-ideal' regulators are assumed.
 * `epsilon` (optional) determines the value of the Epsilon needed to adjust the delta-connected transformers to ensure invertibility. If not provided, a default of 1e-5 is assumed.
 
 ## `computeNoLoadVoltage.m`
 This function computes some required no-load quantities for the network. 
 
 ### List of outputs
 * `network` a structure with the following additional field:
  * `noLoadQuantities` is itself an structure that has the following fields:
   * `v0mag` Vector of magnitude for the slack bus voltage phasor
   * `v0phases` Vector of phases in degrees for the slack bus voltage phasor
   * `v0` Vector of slack bus complex voltage phasor
   * `w` Vector of no-load complex voltage profile for the network
   * `w3phase` The no-load complex voltage profile organized in an (N+1)*3 matrix where the rows determining the bus numbers and the columns determine the phases.
 
 The output of this function should typically be input to the function `setupLoadsIEEE37.m`.
 
 ### List of inputs
 * `network` a network structure created by `setupYbusIEEE37.m`.  
 * `v0mags` (optional) A 3*1 vector of magnitudes for the slack bus voltage phasor.
 * `v0phases` (optional) A 3*1 vector of phases in degrees for the slack bus voltage phasor.
 
 
 ## `setupLoadsIEEE37.m`
 This function updates the `network` structure created previously with information on loads.
 The modeling of loads is based on the load type (constant-power, constant-current, constant-impedance, or any combination of these) as well as the connection (delta, wye, or a combination of these).  The load model is explicitly given in the `calculateIPQII.m` function and is highlighted here first. 
 The general load model is as follows:
```matlab
iL_PQ= c(i,:)*[ fPQ( eVec1.' * v, sL(i,1)); fPQ(eVec2.'*v, sL(i,2))];
iL_I=c(i,:) *[ fI(eVec1.'*v, iL(i,1)); fI(eVec2.'*v, iL(i,2))]; 
iL_Y=c(i,:)*[fY(eVec1.'*v, yL(i,1)); fY(eVec2.'*v, yL(i,2))];
fv(i,1)=g(i,:)*[iL_PQ;iL_I;iL_Y];
```
The vector $c(i,:)$ has two binary entries. The first binary entry corresponds to wye connection. The second binary entry corresponds to delta. If both are set to 1, then the load connection is a mixed wye-delta connection.  The vectors  $eVec1$ and $eVec2$ are needed to decide which phases of the voltage at this node plays a role in the load of the assumed current phase.  The vector $g(i,:)$ determines whether the load is constant-power, constant-current, or constant-impedance.
 
 ### List of outputs
 * `network` a structure with the following additional field:
  * `loadQuantities` a structure with the following fields:
   * `sL_load` Vector of 'nominal power' of constant-power loads  
   * `iL_load` Vector of 'nominal current' of constant-current loads
   * `yL_load` Vector of 'nominal admittance' of constant-impedance loads
   * `ePage`  A page determining which indices of the voltage vector play a role in the considered phase of the current.
   * `gMat` A matrix determining the ZIP load combination.
   * `CMat` A matrix determining the wye, delta, or mixed wye delta connections. (See the above on load modeling)
   
 The output of this function should typically be input to the function `performZBus.m`.
   
   
 ### List of inputs
 * `network` a structure created by `computeNoLoadVoltage.m`
   
## `performZBus.m` 
Computes the Z-Bus load flow.

### List of outputs
* `network` a structure with the following additional field:
 * `ZBusResults` a structure with the following fields:
  * `v` vector of voltages (all phases in vector format)
  * `vsol` solution computed (similar to `v` but does not include the slack voltage)
  * `success` success flag
  * `err` the error at the final iteration
  
  The output of this function should typically be input to the function `obtainVoltages.m`.
 


### List of inputs 
* `network` structure created by 
* `maxIt` (optional) Maximum iterations. If not specified, a default value of 10 is assumed


## `obtainVoltages.m`
This function merely maps the computed voltage vector from the Z-Bus load-flow to appropriate data structure in terms of buses and their phases. 

### List of outputs
*  `network` an structure with the following additional fields:
 * `solution` an structure with the following fields:
  * `resultsVMag` a Matrix (N+1)-by-3 corresponding to computed voltage magnitudes
  * `resultsVPhase` a Matrix (N+1)-by-3 corresponding to computed voltage phases in radians
  * `v3phase` a Matrix (N+1)-by-3 corresponding to computed voltage phasors
  * `v3phaseRegs` a Matrix (N+1+NREgs)-by-3 corresponding to computed voltage phasors including regulator secondary buses
 
 ### List of inputs
 * `network` an structure created by `performZBus.m`
  
