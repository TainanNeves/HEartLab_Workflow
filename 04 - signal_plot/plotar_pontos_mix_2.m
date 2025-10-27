function plotar_pontos_mix_2(DATA_O, DATA_E, fs, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, p27, p28, p29, p30, p31, p32)

%% Figure 01: Picture
I = double(squeeze(DATA_O(:,:,50)));
f = figure('color','white','Position', [40 40 500 500]);
imagesc(I); colormap('gray');
hold on;
text(p17(2), p17(1), num2str(17), 'Color', 'red', 'FontSize', 9);
text(p18(2), p18(1), num2str(18), 'Color', 'red', 'FontSize', 9);
text(p19(2), p19(1), num2str(19), 'Color', 'red', 'FontSize', 9);
text(p20(2), p20(1), num2str(20), 'Color', 'red', 'FontSize', 9);
text(p21(2), p21(1), num2str(21), 'Color', 'red', 'FontSize', 9);
text(p22(2), p22(1), num2str(22), 'Color', 'red', 'FontSize', 9);
text(p23(2), p23(1), num2str(23), 'Color', 'red', 'FontSize', 9);
text(p24(2), p24(1), num2str(24), 'Color', 'red', 'FontSize', 9);
text(p25(2), p25(1), num2str(25), 'Color', 'red', 'FontSize', 9);
text(p26(2), p26(1), num2str(26), 'Color', 'red', 'FontSize', 9);
text(p27(2), p27(1), num2str(27), 'Color', 'red', 'FontSize', 9);
text(p28(2), p28(1), num2str(28), 'Color', 'red', 'FontSize', 9);
text(p29(2), p29(1), num2str(29), 'Color', 'red', 'FontSize', 9);
text(p30(2), p30(1), num2str(30), 'Color', 'red', 'FontSize', 9);
text(p31(2), p31(1), num2str(31), 'Color', 'red', 'FontSize', 9);
text(p32(2), p32(1), num2str(32), 'Color', 'red', 'FontSize', 9);
set(gca,'fontsize', 14);
title('MEA 2 - Electrode Positions');
ylabel('Pixels'); xlabel('Pixels');

%% Figure 02: Combined Optical + Electrical plots
% Optical signals
optic_ch17 = squeeze(DATA_O(p17(1), p17(2),:));
optic_ch18 = squeeze(DATA_O(p18(1), p18(2),:));
optic_ch19 = squeeze(DATA_O(p19(1), p19(2),:));
optic_ch20 = squeeze(DATA_O(p20(1), p20(2),:));
optic_ch21 = squeeze(DATA_O(p21(1), p21(2),:));
optic_ch22 = squeeze(DATA_O(p22(1), p22(2),:));
optic_ch23 = squeeze(DATA_O(p23(1), p23(2),:));
optic_ch24 = squeeze(DATA_O(p24(1), p24(2),:));
optic_ch25 = squeeze(DATA_O(p25(1), p25(2),:));
optic_ch26 = squeeze(DATA_O(p26(1), p26(2),:));
optic_ch27 = squeeze(DATA_O(p27(1), p27(2),:));
optic_ch28 = squeeze(DATA_O(p28(1), p28(2),:));
optic_ch29 = squeeze(DATA_O(p29(1), p29(2),:));
optic_ch30 = squeeze(DATA_O(p30(1), p30(2),:));
optic_ch31 = squeeze(DATA_O(p31(1), p31(2),:));
optic_ch32 = squeeze(DATA_O(p32(1), p32(2),:));

% Electrical signals
elec_ch17 = DATA_E(17, :);
elec_ch18 = DATA_E(18, :);
elec_ch19 = DATA_E(19, :);
elec_ch20 = DATA_E(20, :);
elec_ch21 = DATA_E(21, :);
elec_ch22 = DATA_E(22, :);
elec_ch23 = DATA_E(23, :);
elec_ch24 = DATA_E(24, :);
elec_ch25 = DATA_E(25, :);
elec_ch26 = DATA_E(26, :);
elec_ch27 = DATA_E(27, :);
elec_ch28 = DATA_E(28, :);
elec_ch29 = DATA_E(29, :);
elec_ch30 = DATA_E(30, :);
elec_ch31 = DATA_E(31, :);
elec_ch32 = DATA_E(32, :);

% Plotting Combined Signals
lo = length(optic_ch23);
To = linspace(0, lo/fs, lo);
f1 = figure('color','white','Position', [40 40 1000 800]);
sgtitle('MEA 2: Combined Optical (blue) and Electrical (red) Signals');

% Column 1
subplot(4,4,1); plot_dual_axis(To, optic_ch17, elec_ch17, 'P17', 'el17');
subplot(4,4,5); plot_dual_axis(To, optic_ch18, elec_ch18, 'P18', 'el18');
subplot(4,4,9); plot_dual_axis(To, optic_ch19, elec_ch19, 'P19', 'el19');
subplot(4,4,13); plot_dual_axis(To, optic_ch20, elec_ch20, 'P20', 'el20');

% Column 2
subplot(4,4,2); plot_dual_axis(To, optic_ch21, elec_ch21, 'P21', 'el21');
subplot(4,4,6); plot_dual_axis(To, optic_ch22, elec_ch22, 'P22', 'el22');
subplot(4,4,10); plot_dual_axis(To, optic_ch23, elec_ch23, 'P23', 'el23');
subplot(4,4,14); plot_dual_axis(To, optic_ch24, elec_ch24, 'P24', 'el24');

% Column 3
subplot(4,4,3); plot_dual_axis(To, optic_ch25, elec_ch25, 'P25', 'el25');
subplot(4,4,7); plot_dual_axis(To, optic_ch26, elec_ch26, 'P26', 'el26');
subplot(4,4,11); plot_dual_axis(To, optic_ch27, elec_ch27, 'P27', 'el27');
subplot(4,4,15); plot_dual_axis(To, optic_ch28, elec_ch28, 'P28', 'el28');

% Column 4
subplot(4,4,4); plot_dual_axis(To, optic_ch29, elec_ch29, 'P29', 'el29');
subplot(4,4,8); plot_dual_axis(To, optic_ch30, elec_ch30, 'P30', 'el30');
subplot(4,4,12); plot_dual_axis(To, optic_ch31, elec_ch31, 'P31', 'el31');
subplot(4,4,16); plot_dual_axis(To, optic_ch32, elec_ch32, 'P32', 'el32');

end