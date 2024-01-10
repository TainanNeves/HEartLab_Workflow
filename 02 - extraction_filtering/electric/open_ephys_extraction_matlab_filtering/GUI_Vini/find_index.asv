function idx = find_index(TTL)
% Function that locates the index of the beggining of multiple TTL chains

%% Calculate difference in gaps
% Threshold between the difference of the second and first TTL signal plus 5%
threshold = (TTL(2) - D.TTL(1))*1.05;  
gap = diff(TTL); % Calculates the difference between all consecutive values
idx = find(gap > threshold); % Locates gaps that exceed the threshold
idx = idx+1; % Beggining of the new chain

