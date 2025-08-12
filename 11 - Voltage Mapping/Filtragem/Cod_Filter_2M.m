
%% PARTE 1 - LIMPEZA GERAL

clc; clear; close all;

%% PARTE 2 - CARREGAMENTO DO ARQUIVO

[file, path] = uigetfile('*.oebin', 'Selecione o Arquivo .oebin');
fullfilename = fullfile(path, file); clear file path;

DATA.CONTINUOS = load_open_ephys_binary(fullfilename, 'continuous', 1);
DATA.EVENTS    = load_open_ephys_binary(fullfilename, 'events',     1);

[path, ~, ~] = fileparts(fullfilename);
path_parts = strsplit(path, filesep);
folder_name = path_parts{end};

%% PARTE 3 - CONFIGURAÇÕES INICIAIS

DATA.CONTINUOS.Data = DATA.CONTINUOS.Data * 0.1949999928474426;

DATA.CONTINUOS.Data = (DATA.CONTINUOS.Data/192);
DATA.CONTINUOS.Data = DATA.CONTINUOS.Data * 1000;

Fs = DATA.CONTINUOS.Header.sample_rate;

T = {'MEA 1 - RA', 'MEA 2 - LA'}; MEAs = {1:16, 17:32};

output_folder = folder_name;
if ~isfolder(output_folder)
    mkdir(output_folder);
end

TTL = DATA.EVENTS.Timestamps(1:2:end);
if isempty(TTL)==false
    optictime(1) = (TTL(1) - DATA.CONTINUOS.Timestamps(1));
    optictime(2) = optictime(1) + 10;
end

duracao_total = size(DATA.CONTINUOS.Data, 2) / Fs;
t_meas = (0:size(DATA.CONTINUOS.Data, 2)-1)  / Fs;

%% PARTE 4 - REMOÇÃO DA LINHA DE TENDÊNCIA

w_length = 5; num_segmentos = floor(duracao_total / w_length);

for ee = 1:size(DATA.CONTINUOS.Data, 1)
    sinal = DATA.CONTINUOS.Data(ee, :);
    [sinal_detrend, ~] = detrendSpline(sinal, t_meas, w_length, 0, 'Detrended Signal');
    DATA.CONTINUOS.Data(ee, :) = sinal_detrend;
end

%% PARTE 5 - VISUALIZAÇÃO PRÉ-FILTRAGEM

SF = floor(2 * Fs) + 1; ES = floor(4 * Fs);

for mea_idx = 1:length(MEAs)
    figure(mea_idx); es = MEAs{mea_idx};
    for i = 1:length(es)
        ee = es(i); SS = DATA.CONTINUOS.Data(ee, :);
        subplot(4, 4, i); plot((SF:ES) / Fs, SS(SF:ES));
        xlabel('Tempo (s)'); ylabel('Amplitude (\muV)'); title(['Eletrodo ' num2str(ee)]);
    end
    sgtitle([T{mea_idx} ' - ' folder_name]); savefig(fullfile(output_folder, ['F1_Signal_RAW_' T{mea_idx} '.fig']));
end

%% PARTE 6 - ANÁLISE DO ESPECTRO DE FREQUÊNCIAS

N = length(DATA.CONTINUOS.Data(1, :)); frequencies = (0:N-1) * (Fs / N);

for mea_idx = 1:length(MEAs)
    figure(length(MEAs) + mea_idx); es = MEAs{mea_idx};
    for i = 1:length(es)
        ee = es(i); SS = DATA.CONTINUOS.Data(ee, :);
        Y = fft(SS); magnitude = abs(Y); magnitude = magnitude(1:N/2+1);
        freq = frequencies(1:N/2+1); subplot(4, 4, i); plot(freq, magnitude);
        xlabel('Frequência (Hz)'); ylabel('Magnitude'); title(['Eletrodo ' num2str(ee)]);
        xlim([0 400]);
    end
    sgtitle(['Espectro de Frequências - ' T{mea_idx} ' - ' folder_name]);
    savefig(fullfile(output_folder, ['F2_Signal_Espectro_de_Frequências_' T{mea_idx} '.fig']));
end

%% PARTE 7 - APLICAÇÃO DO FILTRO NOTCH

f0 = [60, 120, 180, 240, 300]; Q = 30; DATA.CONTINUOS.SF = DATA.CONTINUOS.Data;

for freq_idx = 1:length(f0)
    bw = f0(freq_idx) / Q;
    d_notch = designfilt('bandstopiir', ...
        'DesignMethod', 'butter', ...
        'FilterOrder', 2, ...
        'HalfPowerFrequency1', f0(freq_idx) - bw/2, ...
        'HalfPowerFrequency2', f0(freq_idx) + bw/2, ...
        'SampleRate', Fs);
    for ee = 1:size(DATA.CONTINUOS.Data, 1)
        DATA.CONTINUOS.SF(ee, :) = filtfilt(d_notch, DATA.CONTINUOS.SF(ee, :));
    end
end

%% PARTE 8 - VISUALIZAÇÃO PÓS-FILTRAGEM NOTCH

for mea_idx = 1:length(MEAs)
    figure(mea_idx); es = MEAs{mea_idx};
    for i = 1:length(es)
        ee = es(i); SS = DATA.CONTINUOS.SF(ee, :); subplot(4, 4, i);
        plot((SF:ES) / Fs, SS(SF:ES)); xlabel('Tempo (s)'); ylabel('Amplitude (\muV)');
        title(['Eletrodo ' num2str(ee)]);
    end
    sgtitle([T{mea_idx} ' - Filtrado (Notch) - ' folder_name]);
    savefig(fullfile(output_folder, ['F3_Sinal_Filtrado_Pós_NOTCH_' T{mea_idx} '.fig']));
end

%% PARTE 9 - ANÁLISE DO ESPECTRO DE FREQUÊNCIAS 

N = length(DATA.CONTINUOS.SF(1, :)); frequencies = (0:N-1) * (Fs / N);

for mea_idx = 1:length(MEAs)
    figure(length(MEAs) + mea_idx); es = MEAs{mea_idx};
    for i = 1:length(es)
        ee = es(i); SS = DATA.CONTINUOS.SF(ee, :);
        Y = fft(SS); magnitude = abs(Y); magnitude = magnitude(1:N/2+1);
        freq = frequencies(1:N/2+1); subplot(4, 4, i); plot(freq, magnitude);
        xlabel('Frequência (Hz)'); ylabel('Magnitude'); title(['Eletrodo ' num2str(ee)]);
        xlim([0 400]);
    end
    sgtitle(['Espectro de Frequências - ' T{mea_idx} ' - ' folder_name]);
    savefig(fullfile(output_folder, ['F4_Espectro_de_Frequência_Pós_NOTCH' T{mea_idx} '.fig']));
end

%% PARTE 10 - VISUALIZAÇÃO DA GRAVAÇÃO

LEI = [14, 21]; R = {'MEA1', 'MEA2'};

figure;
for i = 1:length(LEI)
    subplot(2, 1, i); SS = DATA.CONTINUOS.SF(LEI(i), :);
    plot((1:length(SS)) / Fs, SS, 'LineWidth', 1);
    xlabel('Tempo (s)'); ylabel('Amplitude (\muV)');
    title(['Sinal Completo - Eletrodo ' num2str(LEI(i)) ' (' R{i} ') - ' folder_name]); grid on;
end
sgtitle(['Visualização de Todo o Sinal - ' folder_name]);
savefig(fullfile(output_folder, 'F5_Sinal_Inteiro.fig'));

%% PARTE 11 - FILTRAGEM PASSA-BANDA

d_highpass = designfilt('highpassiir', ...
    'DesignMethod', 'butter', ...
    'FilterOrder', 2, ...
    'HalfPowerFrequency', 0.5, ...
    'SampleRate', Fs);
d_lowpass = designfilt('lowpassiir', ...
    'DesignMethod', 'butter', ...
    'FilterOrder', 2, ...
    'HalfPowerFrequency', 250, ...
    'SampleRate', Fs);

for ee = 1:size(DATA.CONTINUOS.Data, 1)
    DATA.CONTINUOS.SF(ee, :) = filtfilt(d_highpass, DATA.CONTINUOS.SF(ee, :));
    DATA.CONTINUOS.SF(ee, :) = filtfilt(d_lowpass,  DATA.CONTINUOS.SF(ee, :));
end

%% PARTE 12 - VISUALIZAÇÃO PÓS-FILTRAGEM PASSA-BANDA

for mea_idx = 1:length(MEAs)
    figure(mea_idx); es = MEAs{mea_idx};
    for i = 1:length(es)
        ee = es(i); SS = DATA.CONTINUOS.SF(ee, :); subplot(4, 4, i);
        plot((SF:ES) / Fs, SS(SF:ES)); xlabel('Tempo (s)'); ylabel('Amplitude (\muV)');
        title(['Eletrodo ' num2str(ee)]);
    end
    sgtitle([T{mea_idx} ' - Filtrado (Banda) - ' folder_name]);
    savefig(fullfile(output_folder, ['F6_Sinal_Filtrado_Pós_Filtragem_AltaBaixa_' T{mea_idx} '.fig']));
end

%% PARTE 13 - ANÁLISE DO ESPECTRO DE FREQUÊNCIAS

N = length(DATA.CONTINUOS.SF(1, :)); frequencies = (0:N-1) * (Fs / N);

for mea_idx = 1:length(MEAs)
    figure(6 * length(MEAs) + mea_idx); es = MEAs{mea_idx};
    for i = 1:length(es)
        ee = es(i); SS = DATA.CONTINUOS.SF(ee, :);
        Y = fft(SS); magnitude = abs(Y); magnitude = magnitude(1:N/2+1);
        freq = frequencies(1:N/2+1); subplot(4, 4, i); plot(freq, magnitude);
        xlabel('Frequência (Hz)'); ylabel('Magnitude'); title(['Eletrodo ' num2str(ee)]);
        xlim([0 400]);
    end
    sgtitle(['Espectro de Frequências - ' T{mea_idx} ' (Pós-Passa-Banda) - ' folder_name]);
    savefig(fullfile(output_folder, ['F7_Espectro_de_Frequências_Pós_PassaBaixa_' T{mea_idx} '.fig']));
end

%% PARTE 14 - VISUALIZAÇÃO FINAL

LEI = [14, 21]; R = {'MEA1', 'MEA2'};

figure;
for i = 1:length(LEI)
    subplot(2, 1, i); SS = DATA.CONTINUOS.SF(LEI(i), :);
    plot((1:length(SS)) / Fs, SS, 'LineWidth', 1);
    xlabel('Tempo (s)'); ylabel('Amplitude (\muV)');
    title(['Sinal Completo Filtrado - Eletrodo ' num2str(LEI(i)) ' (' R{i} ') - ' folder_name]); grid on;
end
sgtitle(['Visualização de Todo o Sinal - ' folder_name]);
savefig(fullfile(output_folder, 'F8_Sinal_Inteiro_Final.fig'));

%% PARTE 15 - SALVAMENTO GERAL

FileName = folder_name;

D_EL = struct();
D_EL.Data = DATA.CONTINUOS.SF;
D_EL.Timestamps = DATA.CONTINUOS.Timestamps;
D_EL.Header = DATA.CONTINUOS.Header;
D_EL.Header.channels = size(DATA.CONTINUOS.Data, 1);
D_EL.TTL = TTL;
D_EL.opticalin = TTL(1) - DATA.CONTINUOS.Timestamps(1);
save(fullfile(output_folder, ['electric_data_', FileName, '_filtered.mat']), 'D_EL', '-v7.3');

clear; clc; close all;

%% FIM DO CÓDIGO