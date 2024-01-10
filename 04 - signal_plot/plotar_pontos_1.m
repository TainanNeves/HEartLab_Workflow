% this function plots an image of the data with the location of
% the electrodes
% and a plot with th OP after filtering in each electrode location 
% DATAO= original data without filter
% DATA= data filtered 
% pn= pontos dos eletrodos
% fs= sample frequency

function f1= plotar_pontos_1 (DATAO,DATA,fs,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16)
I=double(squeeze( DATAO(:,:,50)));
f=figure('color','white','Position', [40 40 500 500]);
imagesc(I);colormap('gray');% axis square
hold on;
text(p1(2), p1(1), num2str(1), 'Color', 'red', 'FontSize', 9);hold on;
text(p2(2), p2(1), num2str(2), 'Color', 'red', 'FontSize', 9);hold on;
text(p3(2), p3(1), num2str(3), 'Color', 'red', 'FontSize', 9);hold on;
text(p4(2), p4(1), num2str(4), 'Color', 'red', 'FontSize', 9);hold on;
text(p5(2), p5(1), num2str(5), 'Color', 'red', 'FontSize', 9);hold on;
text(p6(2), p6(1), num2str(6), 'Color', 'red', 'FontSize', 9);hold on;
text(p7(2), p7(1), num2str(7), 'Color', 'red', 'FontSize', 9);hold on;
text(p8(2), p8(1), num2str(8), 'Color', 'red', 'FontSize', 9);hold on;
text(p9(2), p9(1), num2str(9), 'Color', 'red', 'FontSize', 9);hold on;
text(p10(2), p10(1), num2str(10), 'Color', 'red', 'FontSize', 9);hold on;
text(p11(2), p11(1), num2str(11), 'Color', 'red', 'FontSize', 9);hold on;
text(p12(2), p12(1), num2str(12), 'Color', 'red', 'FontSize', 9);hold on;
text(p13(2), p13(1), num2str(13), 'Color', 'red', 'FontSize', 9);hold on;
text(p14(2), p14(1), num2str(14), 'Color', 'red', 'FontSize', 9);hold on;
text(p15(2), p15(1), num2str(15), 'Color', 'red', 'FontSize', 9);hold on;
text(p16(2), p16(1), num2str(16), 'Color', 'red', 'FontSize', 9);hold on;
set(gca,'fontsize', 14);
title('Cam 1');
ylabel('Pixels');xlabel('Pixels');

%plotar
channel1=squeeze(DATA(p1(1), p1(2),:));
channel2=squeeze(DATA(p2(1), p2(2),:));
channel3=squeeze(DATA(p3(1), p3(2),:));
channel4=squeeze(DATA(p4(1), p4(2),:));
channel5=squeeze(DATA(p5(1), p5(2),:));
channel6=squeeze(DATA(p6(1), p6(2),:));
channel7=squeeze(DATA(p7(1), p7(2),:));
channel8=squeeze(DATA(p8(1), p8(2),:));
channel9=squeeze(DATA(p9(1), p9(2),:));
channel10=squeeze(DATA(p10(1), p10(2),:));
channel11=squeeze(DATA(p11(1), p11(2),:));
channel12=squeeze(DATA(p12(1), p12(2),:));
channel13=squeeze(DATA(p13(1), p13(2),:));
channel14=squeeze(DATA(p14(1), p14(2),:));
channel15=squeeze(DATA(p15(1), p15(2),:));
channel16=squeeze(DATA(p16(1), p16(2),:));

lo=length(channel10);
To=linspace(0,lo/fs,lo);%com a nova frequencia de 4000 que igualamos para o eletrico
f1=figure('color','white','Position', [40 40 800 600]);
subplot(4,4,1);plot(To,channel13);
ylabel('P13', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 3]);
subplot(4,4,2);plot(To,channel14);
ylabel('P14', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 3]);
subplot(4,4,3);plot(To,channel15);
ylabel('P15', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 3]);
subplot(4,4,4);plot(To,channel16);
ylabel('P16', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 3]);


subplot(4,4,5);plot(To,channel9);
ylabel('P9', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 3]);
subplot(4,4,6);plot(To,channel10);
ylabel('P10', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 3]);
subplot(4,4,7);plot(To,channel11);
ylabel('P11', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 3]);
subplot(4,4,8);plot(To,channel12);
ylabel('P12', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 3]);

subplot(4,4,9);plot(To,channel5);
ylabel('P5', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 3]);
subplot(4,4,10);plot(To,channel6);
ylabel('P6', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 3]);
subplot(4,4,11);plot(To,channel7);
ylabel('P7', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 3]);
subplot(4,4,12);plot(To,channel8);
ylabel('P8', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 3]);

subplot(4,4,13);plot(To,channel1);
ylabel('P1', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 3]);
subplot(4,4,14);plot(To,channel2);
ylabel('P2', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 3]);
subplot(4,4,15);plot(To,channel3);
ylabel('P3', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 3]);
subplot(4,4,16);plot(To,channel4);
ylabel('P4', 'Interpreter','latex');
set(gca,'fontsize', 12); xlim([0 3]);
xlabel('Time [s]');
set(gca,'fontsize', 12); xlim([0 3]);