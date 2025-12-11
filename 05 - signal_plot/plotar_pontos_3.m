function f1= plotar_pontos_3 (DATA_O,DATA_E,fs,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16)

%% Figure 01: Picture
I=double(squeeze( DATA_O(:,:,50)));
f=figure('color','white','Position', [40 40 500 500]);
imagesc(I);colormap('gray');% axis square
hold on;
text(p1(2), p1(1), num2str(65), 'Color', 'red', 'FontSize', 9);hold on;
text(p2(2), p2(1), num2str(66), 'Color', 'red', 'FontSize', 9);hold on;
text(p3(2), p3(1), num2str(67), 'Color', 'red', 'FontSize', 9);hold on;
text(p4(2), p4(1), num2str(68), 'Color', 'red', 'FontSize', 9);hold on;
text(p5(2), p5(1), num2str(69), 'Color', 'red', 'FontSize', 9);hold on;
text(p6(2), p6(1), num2str(70), 'Color', 'red', 'FontSize', 9);hold on;
text(p7(2), p7(1), num2str(71), 'Color', 'red', 'FontSize', 9);hold on;
text(p8(2), p8(1), num2str(72), 'Color', 'red', 'FontSize', 9);hold on;
text(p9(2), p9(1), num2str(73), 'Color', 'red', 'FontSize', 9);hold on;
text(p10(2), p10(1), num2str(74), 'Color', 'red', 'FontSize', 9);hold on;
text(p11(2), p11(1), num2str(75), 'Color', 'red', 'FontSize', 9);hold on;
text(p12(2), p12(1), num2str(76), 'Color', 'red', 'FontSize', 9);hold on;
text(p13(2), p13(1), num2str(77), 'Color', 'red', 'FontSize', 9);hold on;
text(p14(2), p14(1), num2str(78), 'Color', 'red', 'FontSize', 9);hold on;
text(p15(2), p15(1), num2str(79), 'Color', 'red', 'FontSize', 9);hold on;
text(p16(2), p16(1), num2str(80), 'Color', 'red', 'FontSize', 9);hold on;
set(gca,'fontsize', 14);
title('Cam 2');
ylabel('Pixels');xlabel('Pixels');


%% Figure 02: Optic plots
% Selectin Optical signals
channel65 = squeeze(DATA_O(p1(1), p1(2),:));
channel66 = squeeze(DATA_O(p2(1), p2(2),:));
channel67 = squeeze(DATA_O(p3(1), p3(2),:));
channel68 = squeeze(DATA_O(p4(1), p4(2),:));
channel69 = squeeze(DATA_O(p5(1), p5(2),:));
channel70 = squeeze(DATA_O(p6(1), p6(2),:));
channel71 = squeeze(DATA_O(p7(1), p7(2),:));
channel72 = squeeze(DATA_O(p8(1), p8(2),:));
channel73 = squeeze(DATA_O(p9(1), p9(2),:));
channel74 = squeeze(DATA_O(p10(1), p10(2),:));
channel75 = squeeze(DATA_O(p11(1), p11(2),:));
channel76 = squeeze(DATA_O(p12(1), p12(2),:));
channel77 = squeeze(DATA_O(p13(1), p13(2),:));
channel78 = squeeze(DATA_O(p14(1), p14(2),:));
channel79 = squeeze(DATA_O(p15(1), p15(2),:));
channel80 = squeeze(DATA_O(p16(1), p16(2),:));

% Plotting Optical Signals
lo=length(channel74);
To=linspace(0,lo/fs,lo);
f1=figure('color','white','Position', [40 40 600 600]);
subplot(4,4,1);plot(To,channel77);
ylabel('P77', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 lo/fs]);
subplot(4,4,2);plot(To,channel78);
ylabel('P78', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 lo/fs]);
subplot(4,4,3);plot(To,channel79);
ylabel('P79', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 lo/fs]);
subplot(4,4,4);plot(To,channel80);
ylabel('P80', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 lo/fs]);


subplot(4,4,5);plot(To,channel73);
ylabel('P73', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 lo/fs]);
subplot(4,4,6);plot(To,channel74);
ylabel('P74', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 lo/fs]);
subplot(4,4,7);plot(To,channel75);
ylabel('P75', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 lo/fs]);
subplot(4,4,8);plot(To,channel76);
ylabel('P76', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 lo/fs]);

subplot(4,4,9);plot(To,channel69);
ylabel('P64', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 lo/fs]);
subplot(4,4,10);plot(To,channel70);
ylabel('P70', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 lo/fs]);
subplot(4,4,11);plot(To,channel71);
ylabel('P71', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 lo/fs]);
subplot(4,4,12);plot(To,channel72);
ylabel('P72', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 lo/fs]);

subplot(4,4,13);plot(To,channel65);
ylabel('P65', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 lo/fs]);
subplot(4,4,14);plot(To,channel66);
ylabel('P66', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 lo/fs]);
subplot(4,4,15);plot(To,channel67);
ylabel('P67', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 lo/fs]);
subplot(4,4,16);plot(To,channel68);
ylabel('P68', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 lo/fs]);
xlabel('Time [s]');


%% Figure 03: Electric plots
% Select Electric signals
channel65 = DATA_E(65,:);
channel66 = DATA_E(66,:);
channel67 = DATA_E(67,:);
channel68 = DATA_E(68,:);
channel69 = DATA_E(69,:);
channel70 = DATA_E(70,:);
channel71 = DATA_E(71,:);
channel72 = DATA_E(72,:);
channel73 = DATA_E(73,:);
channel74 = DATA_E(74,:);
channel75 = DATA_E(75,:);
channel76 = DATA_E(76,:);
channel77 = DATA_E(77,:);
channel78 = DATA_E(78,:);
channel79 = DATA_E(79,:);
channel80 = DATA_E(80,:);

% Plotting Electrical Signals
lo=length(channel74);
To=linspace(0,lo/fs,lo);
f1=figure('color','white','Position', [40 40 600 600]);
subplot(4,4,1);plot(To,channel77);
ylabel('el77', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 lo/fs]);
subplot(4,4,2);plot(To,channel78);
ylabel('el78', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 lo/fs]);
subplot(4,4,3);plot(To,channel79);
ylabel('el79', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 lo/fs]);
subplot(4,4,4);plot(To,channel80);
ylabel('el80', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 lo/fs]);


subplot(4,4,5);plot(To,channel73);
ylabel('el73', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 lo/fs]);
subplot(4,4,6);plot(To,channel74);
ylabel('el74', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 lo/fs]);
subplot(4,4,7);plot(To,channel75);
ylabel('el75', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 lo/fs]);
subplot(4,4,8);plot(To,channel76);
ylabel('el76', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 lo/fs]);

subplot(4,4,9);plot(To,channel69);
ylabel('el64', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 lo/fs]);
subplot(4,4,10);plot(To,channel70);
ylabel('el70', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 lo/fs]);
subplot(4,4,11);plot(To,channel71);
ylabel('el71', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 lo/fs]);
subplot(4,4,12);plot(To,channel72);
ylabel('el72', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 lo/fs]);

subplot(4,4,13);plot(To,channel65);
ylabel('el65', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 lo/fs]);
subplot(4,4,14);plot(To,channel66);
ylabel('el66', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 lo/fs]);
subplot(4,4,15);plot(To,channel67);
ylabel('el67', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 lo/fs]);
subplot(4,4,16);plot(To,channel68);
ylabel('el68', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 lo/fs]);
xlabel('Time [s]');