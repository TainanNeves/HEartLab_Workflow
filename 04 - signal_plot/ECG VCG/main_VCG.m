%% MAIN VCG ANALYSIS

%% Cleaning
clear; clc;


%% Loading Data
load("F:\HEartLab\experiment_analyses\exp28_analysis\data_processed\data_filtered_sync_E28_F01_R09.mat");


%% Organizing data
data = D_SYNC.EL;
ecg_x = data(129,:) - data(145,:);                        % Right to the left
ecg_y = (data(190,:) - data(166,:)) * cos(0.43);    % Top to Botton
ecg_z = data(177,:) - data(161,:);                     % Anterior to Posterior
ecg = [ecg_x; ...
        ecg_y;...
        ecg_z];
ecg = ecg';
clear ecg_x ecg_y ecg_z;



%% 
%%
%%
%% Filtering
low_filt = 0.5;
high_filt = 20;
for i = 1:3
    vcg1(:,i) = fftfilter(ecg(:,i),low_filt,high_filt);
    vcg(:,i) = vcg1(2000:max(length(ecg))-2000,i); % remove edge artifact
end
clear vcg1 low_filt high_filt i;


%% Speed calc
% Speed Calculation
deri(1,:) = derivative_cwt(vcg(:,1)','gaus1',16,1,1);
deri(2,:) = derivative_cwt(vcg(:,2)','gaus1',16,1,1);
deri(3,:) = derivative_cwt(vcg(:,3)','gaus1',16,1,1);
deri = deri';
for i = 1:max(size(deri))
    speed(i) = norm(deri(i,:));
end
maxspd = max(speed);
minspd = min(speed);


%% PLOT - SPEED and VCG

toplim = 4000;

% Plot - Speed and Potentials
figure('Color','w'),
% Plot 1: Speed
subplot(2,1,1)
    plot(speed,'LineWidth',2);
    set(gca,'FontWeight','bold','LineWidth',2,'box','on');
    grid on;
    title('Speed', 'FontWeight','bold');
    ylabel('Velocity (distance / ms)','FontWeight','bold');
    xlabel('(ms)');
    xlim([0, toplim]);

% Plot 2: Potentials
subplot(2,1,2)
    % Filtered signal X
    plot(vcg(:,1),'LineWidth',2);
    hold on;
    % Filtered signal Y
    plot(vcg(:,2),'LineWidth',2);
    hold on;
    % Filtered signal Z
    plot(vcg(:,3),'LineWidth',2); 
    % Configuring the plot
    legend('Vx', 'Vy', 'Vz');
    set(gca,'FontWeight','bold','LineWidth',2,'box','on');
    grid on;
    title('Potentials','FontWeight','bold');
    ylabel('uV','FontWeight','bold');
    xlabel('(Sample)');  % Removed the Position parameter for simplicity
    xlim([0, toplim]);
% Link axes for synchronized zooming
linkaxes(findobj(gcf, 'Type', 'axes'), 'x');


%% Save Variables
save('03_52_37_D2 - vcg_data.mat', 'vcg', 'speed');


%%
%% VCG Plot
figure;
% To have it correctlly as our system of coordinates,
% axis Z is the Vy and axis Y is Vz
scatter3(vcg(:,1), vcg(:,3), vcg(:,2), 5, speed, 'filled');
colormap(jet); 
colorbar; 
clim([min(speed), max(speed)]);
% Add title and axis labels
title('VCG');
xlabel('Vx');
zlabel('Vy'); 
ylabel('Vz'); 
grid on;
view(3);  % 3D view


























