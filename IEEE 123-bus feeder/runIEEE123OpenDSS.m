    clear all;
close all;
clc;

% ****************************************************
% * Initialize OpenDSS
% ****************************************************
% Instantiate the OpenDSS Object
DSSObj = actxserver('OpenDSSEngine.DSS');
% Start up the Solver
if ~DSSObj.Start(0),
disp('Unable to start the OpenDSS Engine')
return
end


% Set up the Text, Circuit, and Solution Interfaces
DSSText = DSSObj.Text;
DSSCircuit = DSSObj.ActiveCircuit;
DSSSolution = DSSCircuit.Solution;


% ****************************************************
% * The DSSText Object
% ****************************************************
% Load in the IEEE-13
DSSText.Command = 'Clear ';
DSSText.Command=' Compile Run_IEEE123Bus';
DSSText.Command='Export Voltages'; 
DSSText.Command='Export Losses'
