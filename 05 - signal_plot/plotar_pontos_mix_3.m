function plotar_pontos_mix_3(DATA_O, DATA_E, fs, p65, p66, p67, p68, p69, p70, p71, p72, p73, p74, p75, p76, p77, p78, p79, p80)

%% Figure 01: Picture
I = double(squeeze(DATA_O(:,:,50)));
f = figure('color','white','Position', [40 40 500 500]);
imagesc(I); colormap('gray');
hold on;
text(p65(2), p65(1), num2str(65), 'Color', 'red', 'FontSize', 9);
text(p66(2), p66(1), num2str(66), 'Color', 'red', 'FontSize', 9);
text(p67(2), p67(1), num2str(67), 'Color', 'red', 'FontSize', 9);
text(p68(2), p68(1), num2str(68), 'Color', 'red', 'FontSize', 9);
text(p69(2), p69(1), num2str(69), 'Color', 'red', 'FontSize', 9);
text(p70(2), p70(1), num2str(70), 'Color', 'red', 'FontSize', 9);
text(p71(2), p71(1), num2str(71), 'Color', 'red', 'FontSize', 9);
text(p72(2), p72(1), num2str(72), 'Color', 'red', 'FontSize', 9);
text(p73(2), p73(1), num2str(73), 'Color', 'red', 'FontSize', 9);
text(p74(2), p74(1), num2str(74), 'Color', 'red', 'FontSize', 9);
text(p75(2), p75(1), num2str(75), 'Color', 'red', 'FontSize', 9);
text(p76(2), p76(1), num2str(76), 'Color', 'red', 'FontSize', 9);
text(p77(2), p77(1), num2str(77), 'Color', 'red', 'FontSize', 9);
text(p78(2), p78(1), num2str(78), 'Color', 'red', 'FontSize', 9);
text(p79(2), p79(1), num2str(79), 'Color', 'red', 'FontSize', 9);
text(p80(2), p80(1), num2str(80), 'Color', 'red', 'FontSize', 9);
set(gca,'fontsize', 14);
title('MEA 3 - Electrode Positions');
ylabel('Pixels'); xlabel('Pixels');

%% Figure 02: Combined Optical + Electrical plots
% Optical signals
optic_ch65 = squeeze(DATA_O(p65(1), p65(2),:));
optic_ch66 = squeeze(DATA_O(p66(1), p66(2),:));
optic_ch67 = squeeze(DATA_O(p67(1), p67(2),:));
optic_ch68 = squeeze(DATA_O(p68(1), p68(2),:));
optic_ch69 = squeeze(DATA_O(p69(1), p69(2),:));
optic_ch70 = squeeze(DATA_O(p70(1), p70(2),:));
optic_ch71 = squeeze(DATA_O(p71(1), p71(2),:));
optic_ch72 = squeeze(DATA_O(p72(1), p72(2),:));
optic_ch73 = squeeze(DATA_O(p73(1), p73(2),:));
optic_ch74 = squeeze(DATA_O(p74(1), p74(2),:));
optic_ch75 = squeeze(DATA_O(p75(1), p75(2),:));
optic_ch76 = squeeze(DATA_O(p76(1), p76(2),:));
optic_ch77 = squeeze(DATA_O(p77(1), p77(2),:));
optic_ch78 = squeeze(DATA_O(p78(1), p78(2),:));
optic_ch79 = squeeze(DATA_O(p79(1), p79(2),:));
optic_ch80 = squeeze(DATA_O(p80(1), p80(2),:));

% Electrical signals
elec_ch65 = DATA_E(65, :);
elec_ch66 = DATA_E(66, :);
elec_ch67 = DATA_E(67, :);
elec_ch68 = DATA_E(68, :);
elec_ch69 = DATA_E(69, :);
elec_ch70 = DATA_E(70, :);
elec_ch71 = DATA_E(71, :);
elec_ch72 = DATA_E(72, :);
elec_ch73 = DATA_E(73, :);
elec_ch74 = DATA_E(74, :);
elec_ch75 = DATA_E(75, :);
elec_ch76 = DATA_E(76, :);
elec_ch77 = DATA_E(77, :);
elec_ch78 = DATA_E(78, :);
elec_ch79 = DATA_E(79, :);
elec_ch80 = DATA_E(80, :);

% Plotting Combined Signals
lo = length(optic_ch74);
To = linspace(0, lo/fs, lo);
f1 = figure('color','white','Position', [40 40 1000 800]);
sgtitle('MEA 3: Combined Optical (blue) and Electrical (red) Signals');

% Row 1
subplot(4,4,1); plot_dual_axis(To, optic_ch77, elec_ch77, 'P77', 'el77');
subplot(4,4,2); plot_dual_axis(To, optic_ch78, elec_ch78, 'P78', 'el78');
subplot(4,4,3); plot_dual_axis(To, optic_ch79, elec_ch79, 'P79', 'el79');
subplot(4,4,4); plot_dual_axis(To, optic_ch80, elec_ch80, 'P80', 'el80');

% Row 2
subplot(4,4,5); plot_dual_axis(To, optic_ch73, elec_ch73, 'P73', 'el73');
subplot(4,4,6); plot_dual_axis(To, optic_ch74, elec_ch74, 'P74', 'el74');
subplot(4,4,7); plot_dual_axis(To, optic_ch75, elec_ch75, 'P75', 'el75');
subplot(4,4,8); plot_dual_axis(To, optic_ch76, elec_ch76, 'P76', 'el76');

% Row 3
subplot(4,4,9); plot_dual_axis(To, optic_ch69, elec_ch69, 'P69', 'el69');
subplot(4,4,10); plot_dual_axis(To, optic_ch70, elec_ch70, 'P70', 'el70');
subplot(4,4,11); plot_dual_axis(To, optic_ch71, elec_ch71, 'P71', 'el71');
subplot(4,4,12); plot_dual_axis(To, optic_ch72, elec_ch72, 'P72', 'el72');

% Row 4
subplot(4,4,13); plot_dual_axis(To, optic_ch65, elec_ch65, 'P65', 'el65');
subplot(4,4,14); plot_dual_axis(To, optic_ch66, elec_ch66, 'P66', 'el66');
subplot(4,4,15); plot_dual_axis(To, optic_ch67, elec_ch67, 'P67', 'el67');
subplot(4,4,16); plot_dual_axis(To, optic_ch68, elec_ch68, 'P68', 'el68');

end