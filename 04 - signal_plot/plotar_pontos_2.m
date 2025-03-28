function f1= plotar_pontos_2 (DATA_O,DATA_E,fs,p17,p18,p19,p20,p21,p22,p23,p24,p25,p26,p27,p28,p29,p30,p31,p32)

%% Figure 01: Picture
I=double(squeeze( DATA_O(:,:,50)));
f=figure('color','white');
imagesc(I);colormap('gray');% axis square
hold on;
text(p17(2), p17(1), num2str(17), 'Color', 'red', 'FontSize', 9);hold on;
text(p18(2), p18(1), num2str(18), 'Color', 'red', 'FontSize', 9);hold on;
text(p19(2), p19(1), num2str(19), 'Color', 'red', 'FontSize', 9);hold on;
text(p20(2), p20(1), num2str(20), 'Color', 'red', 'FontSize', 9);hold on;
text(p21(2), p21(1), num2str(21), 'Color', 'red', 'FontSize', 9);hold on;
text(p22(2), p22(1), num2str(22), 'Color', 'red', 'FontSize', 9);hold on;
text(p23(2), p23(1), num2str(23), 'Color', 'red', 'FontSize', 9);hold on;
text(p24(2), p24(1), num2str(24), 'Color', 'red', 'FontSize', 9);hold on;
text(p25(2), p25(1), num2str(25), 'Color', 'red', 'FontSize', 9);hold on;
text(p26(2), p26(1), num2str(26), 'Color', 'red', 'FontSize', 9);hold on;
text(p27(2), p27(1), num2str(27), 'Color', 'red', 'FontSize', 9);hold on;
text(p28(2), p28(1), num2str(28), 'Color', 'red', 'FontSize', 9);hold on;
text(p29(2), p29(1), num2str(29), 'Color', 'red', 'FontSize', 9);hold on;
text(p30(2), p30(1), num2str(30), 'Color', 'red', 'FontSize', 9);hold on;
text(p31(2), p31(1), num2str(31), 'Color', 'red', 'FontSize', 9);hold on;
text(p32(2), p32(1), num2str(32), 'Color', 'red', 'FontSize', 9);hold on;
set(gca,'fontsize', 14);
title('Cam X');
ylabel('Pixels');xlabel('Pixels');


%% Figure 02: Optic plots
% Selectin Optical signals
channel17=squeeze(DATA_O(p17(1), p17(2),:));
channel18=squeeze(DATA_O(p18(1), p18(2),:));
channel19=squeeze(DATA_O(p19(1), p19(2),:));
channel20=squeeze(DATA_O(p20(1), p20(2),:));
channel21=squeeze(DATA_O(p21(1), p21(2),:));
channel22=squeeze(DATA_O(p22(1), p22(2),:));
channel23=squeeze(DATA_O(p23(1), p23(2),:));
channel24=squeeze(DATA_O(p24(1), p24(2),:));
channel25=squeeze(DATA_O(p25(1), p25(2),:));
channel26=squeeze(DATA_O(p26(1), p26(2),:));
channel27=squeeze(DATA_O(p27(1), p27(2),:));
channel28=squeeze(DATA_O(p28(1), p28(2),:));
channel29=squeeze(DATA_O(p29(1), p29(2),:));
channel30=squeeze(DATA_O(p30(1), p30(2),:));
channel31=squeeze(DATA_O(p31(1), p31(2),:));
channel32=squeeze(DATA_O(p32(1), p32(2),:));

% Plotting Optical Signals
lo=length(channel23);
To=linspace(0,lo/fs,lo);%com a nova frequencia de 4000 que igualamos para o eletrico
f1=figure('color','white','Position', [40 40 900 600]);
subplot(4,4,1);plot(To,channel17);
ylabel('P17', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 lo/fs]);
subplot(4,4,2);plot(To,channel21);
ylabel('P21', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 lo/fs]);
subplot(4,4,3);plot(To,channel25);
ylabel('P25', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 lo/fs]);
subplot(4,4,4);plot(To,channel29);
ylabel('P29', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 lo/fs]);

subplot(4,4,5);plot(To,channel18);
ylabel('P18', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 lo/fs]);
subplot(4,4,6);plot(To,channel22);
ylabel('P22', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 lo/fs]);
subplot(4,4,7);plot(To,channel26);
ylabel('P26', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 lo/fs]);
subplot(4,4,8);plot(To,channel30);
ylabel('P30', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 lo/fs]);

subplot(4,4,9);plot(To,channel19);
ylabel('P19', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 lo/fs]);
subplot(4,4,10);plot(To,channel23);
ylabel('P23', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 lo/fs]);
subplot(4,4,11);plot(To,channel27);
ylabel('P27', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 lo/fs]);
subplot(4,4,12);plot(To,channel31);
ylabel('P31', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 lo/fs]);

subplot(4,4,13);plot(To,channel20);
ylabel('P20', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 lo/fs]);
subplot(4,4,14);plot(To,channel24);
ylabel('P24', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 lo/fs]);
subplot(4,4,15);plot(To,channel28);
ylabel('P28', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 lo/fs]);
subplot(4,4,16);plot(To,channel32);
ylabel('P32', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 lo/fs]);
xlabel('Time [s]');
set(gca,'fontsize', 12); xlim([0 lo/fs]);


%% Figure 03: Electric Plot
% Selectin Optical signals
channel17 = DATA_E(17, :);
channel18 = DATA_E(18, :);
channel19 = DATA_E(19, :);
channel20 = DATA_E(20, :);
channel21 = DATA_E(21, :);
channel22 = DATA_E(22, :);
channel23 = DATA_E(23, :);
channel24 = DATA_E(24, :);
channel25 = DATA_E(25, :);
channel26 = DATA_E(26, :);
channel27 = DATA_E(27, :);
channel28 = DATA_E(28, :);
channel29 = DATA_E(29, :);
channel30 = DATA_E(30, :);
channel31 = DATA_E(31, :);
channel32 = DATA_E(32, :);

% Plotting Optical Signals
lo=length(channel23);
To=linspace(0,lo/fs,lo);%com a nova frequencia de 4000 que igualamos para o eletrico
f1=figure('color','white','Position', [40 40 900 600]);
subplot(4,4,1);plot(To,channel17);
ylabel('el17', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 lo/fs]);
subplot(4,4,2);plot(To,channel21);
ylabel('el21', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 lo/fs]);
subplot(4,4,3);plot(To,channel25);
ylabel('el25', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 lo/fs]);
subplot(4,4,4);plot(To,channel29);
ylabel('el29', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 lo/fs]);

subplot(4,4,5);plot(To,channel18);
ylabel('el18', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 lo/fs]);
subplot(4,4,6);plot(To,channel22);
ylabel('el22', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 lo/fs]);
subplot(4,4,7);plot(To,channel26);
ylabel('el26', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 lo/fs]);
subplot(4,4,8);plot(To,channel30);
ylabel('el30', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 lo/fs]);

subplot(4,4,9);plot(To,channel19);
ylabel('el19', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 lo/fs]);
subplot(4,4,10);plot(To,channel23);
ylabel('el23', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 lo/fs]);
subplot(4,4,11);plot(To,channel27);
ylabel('el27', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 lo/fs]);
subplot(4,4,12);plot(To,channel31);
ylabel('el31', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 lo/fs]);

subplot(4,4,13);plot(To,channel20);
ylabel('el20', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 lo/fs]);
subplot(4,4,14);plot(To,channel24);
ylabel('el24', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 lo/fs]);
subplot(4,4,15);plot(To,channel28);
ylabel('el28', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 lo/fs]);
subplot(4,4,16);plot(To,channel32);
ylabel('el32', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 lo/fs]);
xlabel('Time [s]');
set(gca,'fontsize', 12); xlim([0 lo/fs]);





















