% Simple script to find the index of the beggining of mutiple TTL chains.

%% Load the data
D = load('teste.mat'); 

%% Create graphs for visualization
ax1 = axes();
ax2 = axes();

%% Plot complete TTL
hold(ax1, 'on');
gridxy(ax1, D.TTL, 'Color', 'r');

%% Calculate difference in gaps
% Threshold between the difference of the second and first TTL signal plus 5%
threshold = (D.TTL(2) - D.TTL(1));  
gap = diff(D.TTL); % Calculates the difference between all consecutive values
idx = find(gap > threshold); % Locates gaps that exceed the threshold

% The index value found is for the value before the gap, the index of the
% start of the second TTL chain is the idx+1

%% Plot the gap
hold(ax2, 'on');
gridxy(ax2, [D.TTL(idx) D.TTL(idx+1)], 'Color', 'b');



