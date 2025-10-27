function plotar_pontos_mix_1(DATA_O, DATA_E, fs, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16)

%% Figure 01: Picture
I = double(squeeze(DATA_O(:,:,50)));
f = figure('color','white','Position', [40 40 500 500]);
imagesc(I); colormap('gray');
hold on;
text(p1(2), p1(1), num2str(1), 'Color', 'red', 'FontSize', 9);
text(p2(2), p2(1), num2str(2), 'Color', 'red', 'FontSize', 9);
text(p3(2), p3(1), num2str(3), 'Color', 'red', 'FontSize', 9);
text(p4(2), p4(1), num2str(4), 'Color', 'red', 'FontSize', 9);
text(p5(2), p5(1), num2str(5), 'Color', 'red', 'FontSize', 9);
text(p6(2), p6(1), num2str(6), 'Color', 'red', 'FontSize', 9);
text(p7(2), p7(1), num2str(7), 'Color', 'red', 'FontSize', 9);
text(p8(2), p8(1), num2str(8), 'Color', 'red', 'FontSize', 9);
text(p9(2), p9(1), num2str(9), 'Color', 'red', 'FontSize', 9);
text(p10(2), p10(1), num2str(10), 'Color', 'red', 'FontSize', 9);
text(p11(2), p11(1), num2str(11), 'Color', 'red', 'FontSize', 9);
text(p12(2), p12(1), num2str(12), 'Color', 'red', 'FontSize', 9);
text(p13(2), p13(1), num2str(13), 'Color', 'red', 'FontSize', 9);
text(p14(2), p14(1), num2str(14), 'Color', 'red', 'FontSize', 9);
text(p15(2), p15(1), num2str(15), 'Color', 'red', 'FontSize', 9);
text(p16(2), p16(1), num2str(16), 'Color', 'red', 'FontSize', 9);
set(gca,'fontsize', 14);
title('MEA 1 - Electrode Positions');
ylabel('Pixels'); xlabel('Pixels');

%% Figure 02: Combined Optical + Electrical plots
% Optical signals
optic_ch1 = squeeze(DATA_O(p1(1), p1(2),:));
optic_ch2 = squeeze(DATA_O(p2(1), p2(2),:));
optic_ch3 = squeeze(DATA_O(p3(1), p3(2),:));
optic_ch4 = squeeze(DATA_O(p4(1), p4(2),:));
optic_ch5 = squeeze(DATA_O(p5(1), p5(2),:));
optic_ch6 = squeeze(DATA_O(p6(1), p6(2),:));
optic_ch7 = squeeze(DATA_O(p7(1), p7(2),:));
optic_ch8 = squeeze(DATA_O(p8(1), p8(2),:));
optic_ch9 = squeeze(DATA_O(p9(1), p9(2),:));
optic_ch10 = squeeze(DATA_O(p10(1), p10(2),:));
optic_ch11 = squeeze(DATA_O(p11(1), p11(2),:));
optic_ch12 = squeeze(DATA_O(p12(1), p12(2),:));
optic_ch13 = squeeze(DATA_O(p13(1), p13(2),:));
optic_ch14 = squeeze(DATA_O(p14(1), p14(2),:));
optic_ch15 = squeeze(DATA_O(p15(1), p15(2),:));
optic_ch16 = squeeze(DATA_O(p16(1), p16(2),:));

% Electrical signals
elec_ch1 = DATA_E(1, :);
elec_ch2 = DATA_E(2, :);
elec_ch3 = DATA_E(3, :);
elec_ch4 = DATA_E(4, :);
elec_ch5 = DATA_E(5, :);
elec_ch6 = DATA_E(6, :);
elec_ch7 = DATA_E(7, :);
elec_ch8 = DATA_E(8, :);
elec_ch9 = DATA_E(9, :);
elec_ch10 = DATA_E(10, :);
elec_ch11 = DATA_E(11, :);
elec_ch12 = DATA_E(12, :);
elec_ch13 = DATA_E(13, :);
elec_ch14 = DATA_E(14, :);
elec_ch15 = DATA_E(15, :);
elec_ch16 = DATA_E(16, :);

% Plotting Combined Signals
lo = length(optic_ch10);
To = linspace(0, lo/fs, lo);
f1 = figure('color','white','Position', [40 40 1000 800]);
sgtitle('MEA 1: Combined Optical (blue) and Electrical (red) Signals');

% Row 1
subplot(4,4,1); plot_dual_axis(To, optic_ch13, elec_ch13, 'P13', 'el13');
subplot(4,4,2); plot_dual_axis(To, optic_ch14, elec_ch14, 'P14', 'el14');
subplot(4,4,3); plot_dual_axis(To, optic_ch15, elec_ch15, 'P15', 'el15');
subplot(4,4,4); plot_dual_axis(To, optic_ch16, elec_ch16, 'P16', 'el16');

% Row 2
subplot(4,4,5); plot_dual_axis(To, optic_ch9, elec_ch9, 'P9', 'el9');
subplot(4,4,6); plot_dual_axis(To, optic_ch10, elec_ch10, 'P10', 'el10');
subplot(4,4,7); plot_dual_axis(To, optic_ch11, elec_ch11, 'P11', 'el11');
subplot(4,4,8); plot_dual_axis(To, optic_ch12, elec_ch12, 'P12', 'el12');

% Row 3
subplot(4,4,9); plot_dual_axis(To, optic_ch5, elec_ch5, 'P5', 'el5');
subplot(4,4,10); plot_dual_axis(To, optic_ch6, elec_ch6, 'P6', 'el6');
subplot(4,4,11); plot_dual_axis(To, optic_ch7, elec_ch7, 'P7', 'el7');
subplot(4,4,12); plot_dual_axis(To, optic_ch8, elec_ch8, 'P8', 'el8');

% Row 4
subplot(4,4,13); plot_dual_axis(To, optic_ch1, elec_ch1, 'P1', 'el1');
subplot(4,4,14); plot_dual_axis(To, optic_ch2, elec_ch2, 'P2', 'el2');
subplot(4,4,15); plot_dual_axis(To, optic_ch3, elec_ch3, 'P3', 'el3');
subplot(4,4,16); plot_dual_axis(To, optic_ch4, elec_ch4, 'P4', 'el4');

end