function [NodeA,NodeB,Lengthft,Config] = importlines(filename)



%% Import the data
[~, ~, raw] = xlsread(filename,'Sheet1');
raw = raw(4:end,:);
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};

%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells

%% Create output variable
data = reshape([raw{:}],size(raw));

%% Allocate imported array to column variable names
NodeA = data(:,1);
NodeB = data(:,2);
Lengthft = data(:,3);
Config = data(:,4);

%% Clear temporary variables
clearvars data raw R;