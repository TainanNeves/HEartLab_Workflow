
%% PARTE 1 - LIMPEZA GERAL

clc; clear; close all;

%% PARTE 2 - CARREGAMENTO DO ARQUIVO

[file, path] = uigetfile('*_filtered.mat', 'Selecione o arquivo filtrado (*.mat)');
fullfilename = fullfile(path, file);                     

DATA = load(fullfilename); 

[path, ~, ~] = fileparts(fullfilename); 
path_parts = strsplit(path, filesep);
folderName = path_parts{end};

%% PARTE 3 - CONFIGURAÇÕES INICIAIS

T = {'MEA1-RA', 'MEA2-LA'}; 

MEAs = {1:16, 17:32};     

Fs = DATA.D_EL.Header.sample_rate;

percentage_results_m = struct('MEA1', [], 'MEA2', []);
percentage_results_a = struct('MEA1', [], 'MEA2', []);

output_folder = folderName;
if ~isfolder(output_folder)
    mkdir(output_folder);
end

%% PARTE 4 - PLOTAGEM DE VISUALIZAÇÃO

for mea_idx = 1:length(MEAs)
    es = MEAs{mea_idx}; 
    num_eletrodos = length(es); 
    num_figuras = ceil(num_eletrodos / 2);
    
    for fig_idx = 1:num_figuras
        figure;
        for subplot_idx = 1:2
            eletrodo_idx = (fig_idx - 1) * 2 + subplot_idx;
            if eletrodo_idx > num_eletrodos
                break;
            end
            e = es(eletrodo_idx); 
            SS = DATA.D_EL.Data(e, :); 
            subplot(2, 1, subplot_idx);
            plot((1:length(SS)) / Fs, SS / 1000); 
            xlabel('Tempo (s)'); 
            ylabel('Amplitude (mV)');
            title(['Eletrodo ' num2str(e)]);
        end
        sgtitle([T{mea_idx} ' - Fig ' num2str(fig_idx) ' - ' folderName]);
        savefig(fullfile(output_folder, [folderName '_F1_Sign_Filt_' T{mea_idx} '_Fig' num2str(fig_idx) '.fig']));
    end
end

%% PARTE 5 - DETERMINAÇÃO DOS ELETRODOS DE TRABALHO

MEA1 = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 14, 15, 16};
MEA2 = {17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 31, 32};

MEA_all = {MEA1, MEA2}; 

%% PARTE 6 - DETERMINAÇÃO DOS EVENTOS DE ATIVAÇÃO

eletrodos_analise = [8, 22]; janela_tempo = [2 6];

figure('Position', [100 100 900 800], 'Name', 'Marcação Manual');
sgtitle('Análise de Eletrodos', 'FontSize', 12, 'FontWeight', 'normal');
set(gcf, 'DefaultAxesPosition', [0.08 0.08 0.88 0.85]);

for i = 1:2
    subplot(2,1,i); 
    e = eletrodos_analise(i);
    if e <= 16
        mea_idx = 1;
    else
        mea_idx = 2;
    end
    sinal = DATA.D_EL.Data(e, :)/1000; 
    tempo = (0:length(sinal)-1)/Fs; 
    plot(tempo, sinal); 
    xlim(janela_tempo);  
    title(['Eletrodo ' num2str(e) ' - ' T{mea_idx}], 'FontSize', 11);
    xlabel('Tempo (s)'); 
    ylabel('mV'); 
    grid on; 
    hold on;
end

disp('Clique para marcar pontos. ESPAÇO para finalizar.');
while true
    try
        [x, y, button] = ginput(1);
        if button == 32, break; end
        plot(x, y, 'ro', 'MarkerSize', 8);
    catch
        break;
    end
end

savefig(fullfile(output_folder, [folderName '_F2_RA_LA_Marcados_.fig'])); 
disp('Finalizado!');

%% PARTE 7 - DETERMINAÇÃO DA JANELA DE TRABALHO

eletrodo_ref = MEAs{1}(1); 
sinal_ref = DATA.D_EL.Data(eletrodo_ref, :); 
tempo_total = (0:length(sinal_ref)-1)/Fs;

figure('Name', 'Seleção da Janela de Análise', 'NumberTitle', 'off', 'Position', [100 100 900 500]);
plot(tempo_total, sinal_ref / 1000); 
xlabel('Tempo (s)'); 
ylabel('Amplitude (mV)');
title(['Sinal de Referência (Eletrodo ' num2str(eletrodo_ref) ') - Selecione a janela de análise']); 
grid on;

disp('Selecione a janela de análise no gráfico (clique em dois pontos)');
[x, ~] = ginput(2); 
janela_tempo = sort(x);
if diff(janela_tempo) <= 0
    error('Janela inválida selecionada. O segundo ponto deve ser após o primeiro.');
end

hold on;
plot([janela_tempo(1) janela_tempo(1)], ylim, 'r--', 'LineWidth', 1.5);
plot([janela_tempo(2) janela_tempo(2)], ylim, 'r--', 'LineWidth', 1.5);
title(['Janela selecionada: ' num2str(janela_tempo(1)) 's a ' num2str(janela_tempo(2)) 's']);
hold off;

janela_amostras = round(janela_tempo * Fs); 
janela_amostras(1) = max(1, janela_amostras(1));
janela_amostras(2) = min(length(sinal_ref), janela_amostras(2));

janela.tempo_inicio = janela_tempo(1);      
janela.tempo_fim = janela_tempo(2);
janela.amostra_inicio = janela_amostras(1); 
janela.amostra_fim = janela_amostras(2);
janela.duracao = diff(janela_tempo);

disp(['Janela selecionada: ' num2str(janela.tempo_inicio) 's a ' num2str(janela.tempo_fim) 's ('...
      num2str(janela.duracao) 's)']);

figure('Name', 'Janela Selecionada', 'NumberTitle', 'off');
plot(tempo_total(janela_amostras(1):janela_amostras(2)), ...
     sinal_ref(janela_amostras(1):janela_amostras(2)) / 1000);
xlabel('Tempo (s)');
ylabel('Amplitude (mV)');
title(['Sinal na Janela Selecionada (' num2str(janela.tempo_inicio) 's a ' num2str(janela.tempo_fim) 's)']);
grid on;

%% PARTE 8 - MEDIÇÃO MANUAL DA AMPLITUDE PICO-A-PICO

percentage_results_m = struct(...
    'MEA1', struct('Amplitudes', cell(1,16), 'Ativacao', cell(1,16), 'Recuperacao', cell(1,16), ...
                  'Tempos_Ativacao', cell(1,16), 'Tempos_Recuperacao', cell(1,16)), ...
    'MEA2', struct('Amplitudes', cell(1,16), 'Ativacao', cell(1,16), 'Recuperacao', cell(1,16), ...
                  'Tempos_Ativacao', cell(1,16), 'Tempos_Recuperacao', cell(1,16)));

style = struct(...
    'MEA1', struct('color', 'b', 'marker_ativ', 'go', 'marker_rec', 'ro'), ...
    'MEA2', struct('color', 'r', 'marker_ativ', 'g^', 'marker_rec', 'r^'));

mea_list = {'MEA1', 'MEA2'}; 
mea_eletrodos = {MEA1, MEA2};

for mea_idx = 1:2
    current_mea = mea_list{mea_idx};
    eletrodos = mea_eletrodos{mea_idx};
    disp(['Processando ' current_mea '...']); 
    
    for e_idx = 1:length(eletrodos)
        e = eletrodos{e_idx};
        disp(['  Eletrodo ' num2str(e) '...']);
        
        hFig = figure('Name', ['Eletrodo ' num2str(e) ' - ' current_mea], ...
                     'NumberTitle', 'off', 'Position', [100 100 900 500], ...
                     'KeyPressFcn', @(src,evt) setappdata(src, 'keyPressed', evt.Key));
        
        SS = DATA.D_EL.Data(e, janela_amostras(1):janela_amostras(2));
        time_window = tempo_total(janela_amostras(1):janela_amostras(2));
        plot(time_window, SS / 1000, style.(current_mea).color);
        xlabel('Tempo (s)');
        ylabel('Amplitude (mV)');
        title([current_mea ' - Eletrodo ' num2str(e) ' | Janela: ' ...
               num2str(janela.tempo_inicio) 's a ' num2str(janela.tempo_fim) 's']);
        grid on; 
        hold on;
        
        pares = struct('ativacao_val', [], 'ativacao_time', [], ...
                       'recuperacao_val', [], 'recuperacao_time', []);
        proximo_passo = 'ativacao';
        
        done = false;
        while ~done
            figure(hFig);
            waitforbuttonpress;
            key = getappdata(hFig, 'keyPressed');
            
            if strcmp(key, 'space')
                done = true;
                continue;
            
            elseif strcmp(key, 'escape')
                delete(findobj(hFig, 'Type', 'line', 'Marker', 'o'));
                delete(findobj(hFig, 'Type', 'line', 'Marker', '^'));
                delete(findobj(hFig, 'Type', 'text'));
                pares = struct('ativacao_val', [], 'ativacao_time', [], ...
                                'recuperacao_val', [], 'recuperacao_time', []);
                proximo_passo = 'ativacao';
                disp('  Seleção resetada!');
                continue;
            end
            
            pt = get(gca, 'CurrentPoint');
            x = pt(1,1);
            y = pt(1,2) * 1000;
            button = get(gcf, 'SelectionType');
            
            if strcmp(proximo_passo, 'ativacao') && strcmp(button, 'normal')
                pares(end).ativacao_val = y;
                pares(end).ativacao_time = x;
                plot(x, y / 1000, style.(current_mea).marker_ativ, 'MarkerSize', 8, 'LineWidth', 1.5);
                proximo_passo = 'recuperacao';
                disp(['    Ativação: ' num2str(x) ' s, ' num2str(y / 1000) ' mV']);
                
            elseif strcmp(proximo_passo, 'recuperacao') && strcmp(button, 'alt')
                if ~isempty(pares(end).ativacao_val)
                    pares(end).recuperacao_val = y;
                    pares(end).recuperacao_time = x;
                    plot(x, y / 1000, style.(current_mea).marker_rec, 'MarkerSize', 8, 'LineWidth', 1.5);
                    proximo_passo = 'ativacao';
                    pares(end+1) = struct('ativacao_val', [], 'ativacao_time', [], ...
                                          'recuperacao_val', [], 'recuperacao_time', []);
                    disp(['    Recuperação: ' num2str(x) ' s, ' num2str(y / 1000) ' mV']);
                else
                    disp('  AVISO: Marque uma ativação antes da recuperação!');
                end
            end
        end
        
        if isempty(pares(end).ativacao_val) && isempty(pares(end).recuperacao_val)
            pares(end) = [];
        end
        
        amplitudes = [];
        for p = 1:length(pares)
            if ~isempty(pares(p).ativacao_val) && ~isempty(pares(p).recuperacao_val)
                amplitudes = [amplitudes, abs(pares(p).ativacao_val - pares(p).recuperacao_val)];
            end
        end
        
        if ~isempty(amplitudes)
            amp_pico = mean(amplitudes);
            ativacao_vals = [pares.ativacao_val];
            ativacao_times = [pares.ativacao_time];
            recuperacao_vals = [pares.recuperacao_val];
            recuperacao_times = [pares.recuperacao_time];
        else
            amp_pico = NaN;
            ativacao_vals = [];
            ativacao_times = [];
            recuperacao_vals = [];
            recuperacao_times = [];
            warning('  Eletrodo %d: Nenhum par completo foi marcado!', e);
        end
        
        percentage_results_m.(current_mea)(e_idx).Amplitudes = amp_pico;
        percentage_results_m.(current_mea)(e_idx).Ativacao = ativacao_vals;
        percentage_results_m.(current_mea)(e_idx).Recuperacao = recuperacao_vals;
        percentage_results_m.(current_mea)(e_idx).Tempos_Ativacao = ativacao_times;
        percentage_results_m.(current_mea)(e_idx).Tempos_Recuperacao = recuperacao_times;

        savefig(hFig, fullfile(output_folder, [folderName '_FM_Manual_Peak-to-Peak_' current_mea '_E' num2str(e) '.fig']));
        close(hFig);
    end
end

%% PARTE 9 - MEDIÇÃO AUTOMÁTICA DA AMPLITUDE PICO-A-PICO

MaxHR = 300; MPD = (60 / MaxHR) * Fs; MPHF = 0.3;

percentage_results_a = struct(...
    'MEA1', struct('Amplitudes', cell(1,16), 'Ativacao', cell(1,16), 'Recuperacao', cell(1,16), ...
              'Tempos_Ativacao', cell(1,16), 'Tempos_Recuperacao', cell(1,16)), ...
    'MEA2', struct('Amplitudes', cell(1,16), 'Ativacao', cell(1,16), 'Recuperacao', cell(1,16), ...
              'Tempos_Ativacao', cell(1,16), 'Tempos_Recuperacao', cell(1,16)));

for m = 1:length(MEA_all)
    current_mea = ['MEA' num2str(m)]; 
    MEA = [MEA_all{m}{:}];
    
    for e_idx = 1:length(MEA)
        e = MEA(e_idx);
        SS = DATA.D_EL.Data(e, janela.amostra_inicio:janela.amostra_fim);
        time_window = (janela.amostra_inicio:janela.amostra_fim) / Fs;
        
        [pks_a, locs_a] = findpeaks(SS, ...
            'MinPeakHeight', max(SS) * MPHF, ...
            'MinPeakDistance', MPD);
        
        [pks_r, locs_r] = findpeaks(-SS, ...
            'MinPeakHeight', max(-SS) * MPHF, ...
            'MinPeakDistance', MPD);
        pks_r = -pks_r;
        
        if ~isempty(pks_a) && ~isempty(pks_r)
            amp_pico = zeros(min(length(pks_a), length(pks_r)), 1);
            for i = 1:min(length(pks_a), length(pks_r))
                amp_pico(i) = pks_a(i) - pks_r(i);
            end
            mean_amp = mean(abs(amp_pico));
        else
            mean_amp = NaN;
        end
        
        percentage_results_a.(current_mea)(e_idx).Amplitudes = mean_amp;
        percentage_results_a.(current_mea)(e_idx).Ativacao = pks_a;
        percentage_results_a.(current_mea)(e_idx).Recuperacao = pks_r;
        percentage_results_a.(current_mea)(e_idx).Tempos_Ativacao = time_window(locs_a);
        percentage_results_a.(current_mea)(e_idx).Tempos_Recuperacao = time_window(locs_r);
    end
    
    num_eletrodos = length(MEA); 
    num_figuras = ceil(num_eletrodos / 2);
    
    for fig_idx = 1:num_figuras
        figure('Position', [100, 100, 900, 600]);
        for subplot_idx = 1:2
            eletrodo_idx = (fig_idx - 1) * 2 + subplot_idx;
            if eletrodo_idx > num_eletrodos
                break;
            end
            e = MEA(eletrodo_idx);
            e_data = percentage_results_a.(current_mea)(eletrodo_idx);
            subplot(2, 1, subplot_idx);
            hold on;
            SS = DATA.D_EL.Data(e, janela.amostra_inicio:janela.amostra_fim);
            plot(time_window, SS / 1000, 'k', 'DisplayName', 'Sinal');
            
            if ~isempty(e_data.Ativacao)
                plot(e_data.Tempos_Ativacao, e_data.Ativacao / 1000, 'go', ...
                    'MarkerSize', 8, 'LineWidth', 1.5, 'DisplayName', 'Ativação');
            end
            
            if ~isempty(e_data.Recuperacao)
                plot(e_data.Tempos_Recuperacao, e_data.Recuperacao / 1000, 'ro', ...
                    'MarkerSize', 8, 'LineWidth', 1.5, 'DisplayName', 'Recuperação');
            end
            hold off;
            xlabel('Tempo (s)');
            ylabel('Amplitude (mV)');
            title(sprintf('Eletrodo %d - Amp: %.1f mV', e, e_data.Amplitudes / 1000));
            legend('Location', 'best');
            grid on;
            xlim([janela.tempo_inicio, janela.tempo_fim]);
        end
        sgtitle(sprintf('%s - Fig %d - %s', T{m}, fig_idx, folderName));
        savefig(fullfile(output_folder, [folderName '_F3_Auto_Peak-to-Peak_' T{m} '_Fig' num2str(fig_idx) '.fig']));
    end
end

%% PARTE 10 - FREQUÊNCIA DE ATIVAÇÃO

MaxHR = 300; MPD = (60 / MaxHR) * Fs; MPHF = 0.3;

HR_str = struct('MEA1', [], 'MEA2', []);

for m = 1:length(MEA_all)
    current_mea = ['MEA' num2str(m)];
    MEA = [MEA_all{m}{:}];
    num_eletrodos = length(MEA);
    num_figuras = ceil(num_eletrodos / 2);
    HR_values = zeros(1, num_eletrodos);
    
    for fig_idx = 1:num_figuras
        figure('Position', [100, 100, 900, 600]);
        for subplot_idx = 1:2
            eletrodo_idx = (fig_idx - 1) * 2 + subplot_idx;
            if eletrodo_idx > num_eletrodos
                break;
            end
            e = MEA(eletrodo_idx);
            SS = DATA.D_EL.Data(e, janela.amostra_inicio:janela.amostra_fim);
            time_window = (janela.amostra_inicio:janela.amostra_fim) / Fs;
            MPH = MPHF * max(SS);
            [pks, locs] = findpeaks(SS, 'MinPeakHeight', MPH, 'MinPeakDistance', MPD);
            RR_intervals = diff(locs) / Fs;
            
            if length(RR_intervals) >= 1
                HR = 60 / mean(RR_intervals);
                HR_values(eletrodo_idx) = HR;
            else
                HR_values(eletrodo_idx) = NaN;
            end
            
            subplot(2, 1, subplot_idx);
            plot(time_window, SS / 1000, 'b', 'DisplayName', 'Sinal');
            hold on;
            plot(time_window(locs), pks / 1000, 'ro', 'MarkerSize', 8, 'DisplayName', 'Picos');
            hold off;
            xlabel('Tempo (s)');
            ylabel('Amplitude (mV)');
            title(sprintf('Eletrodo %d - HR: %.1f bpm', e, HR_values(eletrodo_idx)));
            legend('Location', 'best');
            grid on;
            xlim([janela.tempo_inicio, janela.tempo_fim]);
        end
        sgtitle(sprintf('%s - Fig %d - %s | Janela: %.1fs a %.1fs', ...
               T{m}, fig_idx, folderName, ...
               janela.tempo_inicio, janela.tempo_fim));
        savefig(fullfile(output_folder, [folderName '_F4_HR_' T{m} '_Fig' num2str(fig_idx) '.fig']));
    end
    HR_str.(current_mea) = HR_values;
end

%% PARTE 11 - ESTATÍSTICA DAS FREQUÊNCIAS DE ATIVAÇÃO E CYCLE LENGTHS

valid_AF = struct();

for m = 1:length(T)
    current_mea = ['MEA' num2str(m)];
    valid_AF.(current_mea) = HR_str.(current_mea)(HR_str.(current_mea) > 0);
end

figure('Name', 'Frequências de Ativação', 'Position', [100, 100, 900, 650], 'Color', 'w');
plot_data = []; 
group_labels = [];
for m = 1:length(T)
    current_data = valid_AF.(['MEA' num2str(m)]);
    plot_data = [plot_data; current_data'];
    group_labels = [group_labels; repmat(T(m), length(current_data), 1)];
end

bp = boxplot(plot_data, group_labels, ...
    'Widths', 0.7, ...
    'Colors', [0 0.4 0.7], ...
    'Symbol', 'k+');

set(gca, 'FontSize', 11, 'LineWidth', 1.2);
title('Distribuição de Frequência de Ativação por MEA', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Frequência de Ativação (Hz)', 'FontSize', 12);
xlabel('Matriz de Eletrodos', 'FontSize', 12);
grid on;
savefig(fullfile(output_folder, [folderName '_F5_AF_Boxplot.fig']));

figure('Name', 'Estatísticas de Ativação', 'Position', [100, 100, 900, 650], 'Color', 'w');
stats = struct();
for m = 1:length(T)
    current_data = valid_AF.(['MEA' num2str(m)]);
    stats.mean(m) = mean(current_data)/60;
    stats.median(m) = median(current_data)/60;
    stats.N(m) = length(current_data);
end

bar_width = 0.8; 
x_pos = 1:length(T);
b = bar(x_pos, stats.mean, bar_width, 'FaceColor', [0.2 0.6 0.8], 'LineWidth', 1.2);
hold on;
b2 = bar(x_pos, stats.median, bar_width/2, 'FaceColor', [0.8 0.4 0.2], 'LineWidth', 1.2);

for m = 1:length(T)
    text(x_pos(m), stats.mean(m)*0.95, sprintf('%.3f Hz', stats.mean(m)), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top', ...
        'Color', 'w', ...
        'FontWeight', 'bold', ...
        'FontSize', 11);    
    text(x_pos(m), stats.median(m)*0.8, sprintf('%.3f Hz', stats.median(m)), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', ...
        'Color', 'w', ...
        'FontWeight', 'bold', ...
        'FontSize', 10);
end

set(gca, 'XTick', x_pos, 'XTickLabel', T, 'FontSize', 11);
ylabel('Frequência de Ativação (Hz)', 'FontSize', 12);
title('Estatística de Frequência de Ativação', 'FontSize', 14, 'FontWeight', 'bold');
legend([b b2], {'Média', 'Mediana'}, 'Location', 'best');

for m = 1:length(T)
    text(x_pos(m), max(ylim)*0.05, sprintf('n = %d', stats.N(m)), ...
        'HorizontalAlignment', 'center', ...
        'FontSize', 10);
end
grid on;
savefig(fullfile(output_folder, [folderName '_F6_AF_Statistics.fig']));

figure('Name', 'Intervalo de Ativação', 'Position', [100, 100, 900, 650], 'Color', 'w');
activation_intervals = (1./stats.median)*1000;
b = bar(x_pos, activation_intervals, bar_width, ...
    'FaceColor', [0.4 0.7 0.4], ...
    'LineWidth', 1.2);
hold on;

for m = 1:length(T)
    text(x_pos(m), activation_intervals(m)*0.95, ...
        sprintf('%.1f ms', activation_intervals(m)), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top', ...
        'Color', 'w', ...
        'FontWeight', 'bold', ...
        'FontSize', 11);
    text(x_pos(m), max(ylim)*0.05, ...
        sprintf('%.3f Hz', stats.median(m)), ...
        'HorizontalAlignment', 'center', ...
        'FontSize', 10);
end

set(gca, 'XTick', x_pos, 'XTickLabel', T, 'FontSize', 11);
ylabel('Intervalo de Ativação (ms)', 'FontSize', 12);
title('Mediana do Intervalo da Frequência de Ativação por MEA', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
ref_interval = 200;
yline(ref_interval, '--', sprintf('Ref: %.0f ms (5 Hz)', ref_interval), ...
    'Color', [0.5 0.5 0.5], 'LineWidth', 1.2, 'LabelVerticalAlignment', 'bottom');
savefig(fullfile(output_folder, [folderName '_F7_Activation_Intervals.fig']));

clc;
disp('===== ANALYSIS COMPLETED =====');
disp('Generated files:');
disp(['1. Boxplot: ' folderName '_F5_AF_Boxplot.fig']);
disp(['2. Statistics: ' folderName '_F6_AF_Statistics.fig']);
disp(['3. Intervals: ' folderName '_F7_Activation_Intervals.fig']);
disp(' ');
disp('Summary statistics:');
for m = 1:length(T)
    disp([T{m} ':']);
    disp(['   Mean frequency: ' num2str(stats.mean(m), '%.3f') ' Hz']);
    disp(['   Median frequency: ' num2str(stats.median(m), '%.3f') ' Hz']);
    disp(['   Activation interval: ' num2str(activation_intervals(m), '%.1f') ' ms']);
    disp(['   Samples: ' num2str(stats.N(m)) ' electrodes']);
    disp(' ');
end

%% PARTE 12 - ESTATÍSTICA DAS AMPLITUDES PICO-A-PICO

colors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980]; 
bar_width = 0.35;

auto_data = {[percentage_results_a.MEA1.Amplitudes]/1000, [percentage_results_a.MEA2.Amplitudes]/1000};
manual_data = {[percentage_results_m.MEA1.Amplitudes]/1000, [percentage_results_m.MEA2.Amplitudes]/1000};

auto_stats = struct(); 
manual_stats = struct();

for m = 1:2
    current_auto = auto_data{m}(~isnan(auto_data{m}));
    current_manual = manual_data{m}(~isnan(manual_data{m}));  
    
    auto_stats(m).MEA = T{m};
    auto_stats(m).Media = mean(current_auto);
    auto_stats(m).Mediana = median(current_auto);
    auto_stats(m).DesvioPadrao = std(current_auto);
    auto_stats(m).Variancia = var(current_auto);
    auto_stats(m).Percentil5 = prctile(current_auto, 5);
    auto_stats(m).N = length(current_auto);

    manual_stats(m).MEA = T{m};
    manual_stats(m).Media = mean(current_manual);
    manual_stats(m).Mediana = median(current_manual);
    manual_stats(m).DesvioPadrao = std(current_manual);
    manual_stats(m).Variancia = var(current_manual);
    manual_stats(m).Percentil5 = prctile(current_manual, 5);
    manual_stats(m).N = length(current_manual);
end

fig1 = figure('Position', [100, 100, 800, 500], 'Name', 'Comparação de Medianas das Amplitudes');
auto_medians = [auto_stats.Mediana];
manual_medians = [manual_stats.Mediana];
x = 1:2;
x_auto = x - bar_width/2;
x_manual = x + bar_width/2;

b_auto = bar(x_auto, auto_medians, bar_width, 'FaceColor', colors(1,:), 'DisplayName', 'Automático');
hold on;
b_manual = bar(x_manual, manual_medians, bar_width, 'FaceColor', colors(2,:), 'DisplayName', 'Manual');
set(gca, 'XTick', x, 'XTickLabel', T);
ylabel('Amplitude (mV)');
legend('show', 'Location', 'best');
grid on;

for i = 1:2
    text(x_auto(i), auto_medians(i), sprintf('%.2f', auto_medians(i)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    text(x_manual(i), manual_medians(i), sprintf('%.2f', manual_medians(i)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end

title('Comparação Automático vs Manual das Medianas das Amplitudes Pico-a-Pico');

saveas(fig1, fullfile(output_folder, [folderName '_F8_Comparacao_Medianas_Amplitudes.fig']));

fig2 = figure('Position', [100, 100, 900, 600], 'Name', 'Distribuição Combinada das Amplitudes');
combined_data = cell(1,4);
for m = 1:2
    combined_data{2*m-1} = auto_data{m}(~isnan(auto_data{m}));
    combined_data{2*m} = manual_data{m}(~isnan(manual_data{m}));
end

groups = [];
labels = {};
for m = 1:2
    groups = [groups (m-1)*2+1 (m-1)*2+2];
    labels = [labels sprintf('%s-Auto', T{m}) sprintf('%s-Manual', T{m})];
end

max_len = max(cellfun(@length, combined_data));
plot_data = NaN(max_len, 4);
for i = 1:4
    plot_data(1:length(combined_data{i}),i) = combined_data{i};
end

boxplot(plot_data, 'positions', groups, 'labels', labels, 'colors', 'rb', 'symbol', 'k+');
hold on;

for i = 1:4
    method_color = colors(mod(i,2)+1,:);
    x = groups(i) + (rand(size(combined_data{i}))-0.5)*0.3;
    scatter(x, combined_data{i}, 30, method_color, 'filled', 'MarkerFaceAlpha', 0.6);
end

ylabel('Amplitude (mV)');
grid on;
set(gca, 'XTickLabelRotation', 45);

h = zeros(2,1);
h(1) = plot(NaN,NaN, 's', 'MarkerFaceColor', colors(1,:), 'MarkerEdgeColor', 'k');
h(2) = plot(NaN,NaN, 's', 'MarkerFaceColor', colors(2,:), 'MarkerEdgeColor', 'k');
legend(h, {'Automático', 'Manual'}, 'Location', 'best');
title('Distribuição das Amplitudes Pico-a-Pico por MEA e Método');

saveas(fig2, fullfile(output_folder, [folderName '_F9_Boxplots_Amplitudes.fig']));

fig3 = figure('Position', [100, 100, 900, 700]);
auto_values = [auto_data{:}];
manual_values = [manual_data{:}];
all_values = [auto_values, manual_values];

percentiles = prctile(all_values, [0 50 75 90 95 97.5 99 100]);
bins = unique(floor(percentiles/10)*10);
bins(end) = inf;

for m = 1:2
    subplot(2, 2, 2*m-1);
    [counts] = histcounts(auto_data{m}, bins, 'Normalization', 'probability');
    bar(1:length(counts), counts*100, 'FaceColor', colors(1,:));
    title([T{m} ' - Automático']);
    xlabel('Amplitude (mV)');
    ylabel('Porcentagem (%)');
    grid on;
    
    xticks(1:length(counts));
    if length(counts) > 10
        xticklabels(compose('%d', bins(1:end-1)));
    else
        labels = arrayfun(@(i) sprintf('%d-%d', bins(i), bins(i+1)), 1:length(counts), 'UniformOutput', false);
        labels{end} = sprintf('>%d', bins(end-1));
        xticklabels(labels);
    end
    
    subplot(2, 2, 2*m);
    [counts] = histcounts(manual_data{m}, bins, 'Normalization', 'probability');
    bar(1:length(counts), counts*100, 'FaceColor', colors(2,:));
    title([T{m} ' - Manual']);
    xlabel('Amplitude (mV)');
    ylabel('Porcentagem (%)');
    grid on;
    xticks(1:length(counts));
    if length(counts) > 10
        xticklabels(compose('%d', bins(1:end-1)));
    else
        labels = arrayfun(@(i) sprintf('%d-%d', bins(i), bins(i+1)), 1:length(counts), 'UniformOutput', false);
        labels{end} = sprintf('>%d', bins(end-1));
        xticklabels(labels);
    end
end
sgtitle('Distribuição das Amplitudes Pico-a-Pico por Intervalo');

saveas(fig3, fullfile(output_folder, [folderName '_F10_Distribuicao_Amplitudes.fig']));

comparison_table = table(...
    [auto_stats.Mediana]', [auto_stats.Percentil5]', [auto_stats.N]', ...
    [manual_stats.Mediana]', [manual_stats.Percentil5]', [manual_stats.N]', ...
    'VariableNames', {'Mediana_Auto', 'P95_Auto', 'N_Auto', ...
                     'Mediana_Manual', 'P95_Manual', 'N_Manual'}, ...
    'RowNames', T);

amplitude_results = struct(...
    'auto_stats', auto_stats, ...
    'manual_stats', manual_stats);

disp(' ');
disp('======================================================================');
disp('ESTATÍSTICAS DAS AMPLITUDES PICO-A-PICO (mV)');
disp('======================================================================');

for m = 1:2
    disp(' ');
    disp(['ESTATÍSTICAS PARA ' T{m} ':']);
    disp('-----------------------------------');
    disp('Método Automático:');
    disp(['  Média:         ' num2str(auto_stats(m).Media, '%.2f')]);
    disp(['  Desvio Padrão: ' num2str(auto_stats(m).DesvioPadrao, '%.2f')]);
    disp(['  Mediana:       ' num2str(auto_stats(m).Mediana, '%.2f')]);
    disp(['  Percentil 95%: ' num2str(auto_stats(m).Percentil5, '%.2f')]);
    disp(['  N:            ' num2str(auto_stats(m).N)]);
    
    disp(' ');
    disp('Método Manual:');
    disp(['  Média:         ' num2str(manual_stats(m).Media, '%.2f')]);
    disp(['  Desvio Padrão: ' num2str(manual_stats(m).DesvioPadrao, '%.2f')]);
    disp(['  Mediana:       ' num2str(manual_stats(m).Mediana, '%.2f')]);
    disp(['  Percentil 95%: ' num2str(manual_stats(m).Percentil5, '%.2f')]);
    disp(['  N:            ' num2str(manual_stats(m).N)]);
    disp('-----------------------------------');
end

figure('Name', 'Estatísticas das Amplitudes Pico-a-Pico', 'Position', [100, 100, 800, 400], 'Color', 'w');
dados_tabela = {};
cabecalho = {'MEA', 'Método', 'Média (mV)', 'Desvio Padrão (mV)', 'Mediana (mV)', 'Percentil 95% (mV)', 'N'};

for m = 1:2
    dados_tabela{end+1, 1} = T{m};
    dados_tabela{end, 2} = 'Automático';
    dados_tabela{end, 3} = auto_stats(m).Media;
    dados_tabela{end, 4} = auto_stats(m).DesvioPadrao;
    dados_tabela{end, 5} = auto_stats(m).Mediana;
    dados_tabela{end, 6} = auto_stats(m).Percentil5;
    dados_tabela{end, 7} = auto_stats(m).N;

    dados_tabela{end+1, 1} = T{m};
    dados_tabela{end, 2} = 'Manual';
    dados_tabela{end, 3} = manual_stats(m).Media;
    dados_tabela{end, 4} = manual_stats(m).DesvioPadrao;
    dados_tabela{end, 5} = manual_stats(m).Mediana;
    dados_tabela{end, 6} = manual_stats(m).Percentil5;
    dados_tabela{end, 7} = manual_stats(m).N;
end

tabela_formatada = cell2table(dados_tabela, 'VariableNames', cabecalho);

t = uitable('Data', dados_tabela, 'ColumnName', cabecalho, ...
            'RowName', [], 'Units', 'Normalized', 'Position', [0.05, 0.1, 0.9, 0.8]);

set(t, 'ColumnWidth', {70, 80, 80, 100, 80, 100, 40});
set(t, 'FontSize', 10);

title('Estatísticas das Amplitudes Pico-a-Pico por MEA e Método', 'FontSize', 12, 'FontWeight', 'bold');

savefig(fullfile(output_folder, [folderName '_F11_Estatisticas_Amplitudes.fig']));

print(fullfile(output_folder, [folderName '_F11_Estatisticas_Amplitudes.png']), '-dpng', '-r300');

%% PARTE 13 - VISÃO GERAL DA MEDIANA DAS AMPLITUDES PICO-A-PICO

cores = [0 0.5 0.8; 0.8 0.3 0.2];

for mea_idx = 1:2
    current_mea = ['MEA' num2str(mea_idx)];
    fig = figure('Position', [100, 100, 1200, 800], 'Color', 'w', 'Name', ['Comparação de Medianas - ' T{mea_idx}], 'NumberTitle', 'off');
    sgtitle(['Comparação de Medianas Manual (M) e Automático (A) - ' T{mea_idx} ' | ' folderName], 'FontSize', 14, 'FontWeight', 'bold');
    
    eletrodos = [MEA_all{mea_idx}{:}]; 
    num_eletrodos = length(eletrodos); 
    num_linhas = ceil(num_eletrodos/4);
    
    for e = 1:num_eletrodos
        subplot(num_linhas, 4, e); 
        hold on; 
        e_idx = find(eletrodos(e) == [MEA_all{mea_idx}{:}]);
        
        if ~isempty(percentage_results_m.(current_mea)(e_idx).Amplitudes)
            med_manual = percentage_results_m.(current_mea)(e_idx).Amplitudes/1000;
            if ~isnan(med_manual)
                line([0 1], [med_manual med_manual], ...
                     'Color', cores(1,:), 'LineWidth', 3, ...
                     'DisplayName', 'Manual');
                text(0.5, med_manual, sprintf('M: %.2f', med_manual), ...
                     'Color', cores(1,:), 'FontWeight', 'bold', ...
                     'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
            end
        end
        
        if ~isempty(percentage_results_a.(current_mea)(e_idx).Amplitudes)
            med_auto = percentage_results_a.(current_mea)(e_idx).Amplitudes/1000;
            if ~isnan(med_auto)
                line([1 2], [med_auto med_auto], ...
                     'Color', cores(2,:), 'LineWidth', 3, ...
                     'DisplayName', 'Auto');
                text(1.5, med_auto, sprintf('A: %.2f', med_auto), ...
                     'Color', cores(2,:), 'FontWeight', 'bold', ...
                     'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
            end
        end
        
        title(['E' num2str(eletrodos(e))], 'FontSize', 10); 
        ylabel('Amplitude (mV)', 'FontSize', 9);
        xlim([0 2]); 
        set(gca, 'XTick', [], 'Box', 'on', 'FontSize', 8); 
        grid on;
     
        vals = [];
        if exist('med_manual', 'var') && ~isnan(med_manual)
            vals = [vals, med_manual];
        end
        if exist('med_auto', 'var') && ~isnan(med_auto)
            vals = [vals, med_auto];
        end   
        
        if ~isempty(vals)
            ylim([0 max(vals)*1.2]);
        else
            ylim([0 1]);
        end
        clear med_manual med_auto
    end

    savefig(fig, fullfile(output_folder, [folderName '_F12_Comparacao_Medianas_' T{mea_idx} '.fig']));
    pos = get(fig, 'Position');
    set(fig, 'Position', [pos(1)+100*(mea_idx-1), pos(2)-100*(mea_idx-1), pos(3), pos(4)]);
end

%% PARTE 14 - SALVAMENTO FINAL

FileName = folderName;

activation_frequency = struct('MEA1', HR_str.MEA1, 'MEA2', HR_str.MEA2);

mean_HR = struct();
mean_HR.MEA1 = mean(HR_str.MEA1(~isnan(HR_str.MEA1)));
mean_HR.MEA2 = mean(HR_str.MEA2(~isnan(HR_str.MEA2)));

mean_activation_frequency = struct('MEA1', mean_HR.MEA1, 'MEA2', mean_HR.MEA2);

peak_to_peak_amplitudes = struct('MEA1', [], 'MEA2', []);

for m = 1:2
    current_mea = ['MEA' num2str(m)];
    MEA = [MEA_all{m}{:}];
    num_eletrodos = length(MEA);
    peak_to_peak_amplitudes.(current_mea) = cell(1, num_eletrodos);
    
    for idx = 1:num_eletrodos
        e = MEA(idx);
        auto_data = percentage_results_a.(current_mea)(idx);
        manual_data = percentage_results_m.(current_mea)(idx);
        
        peak_to_peak_amplitudes.(current_mea){idx} = struct(...
            'Electrode', e, ...
            'Auto_ActivationPks', auto_data.Ativacao, ...
            'Auto_ActivationLocs', auto_data.Tempos_Ativacao, ...
            'Auto_RecoveryPks', auto_data.Recuperacao, ...
            'Auto_RecoveryLocs', auto_data.Tempos_Recuperacao, ...
            'Auto_PeakToPeakAmplitude', auto_data.Amplitudes/1000, ...
            'Manual_ActivationPks', manual_data.Ativacao, ...
            'Manual_ActivationLocs', manual_data.Tempos_Ativacao, ...
            'Manual_RecoveryPks', manual_data.Recuperacao, ...
            'Manual_RecoveryLocs', manual_data.Tempos_Recuperacao, ...
            'Manual_PeakToPeakAmplitude', manual_data.Amplitudes/1000 ...
        );
    end
end

mean_amplitude_per_MEA = struct();
mean_amplitude_per_MEA.MEA1 = mean([percentage_results_m.MEA1.Amplitudes]/1000, 'omitnan');
mean_amplitude_per_MEA.MEA2 = mean([percentage_results_m.MEA2.Amplitudes]/1000, 'omitnan');

mean_peak_to_peak_amplitudes = struct(...
    'MEA1', mean_amplitude_per_MEA.MEA1, ...
    'MEA2', mean_amplitude_per_MEA.MEA2);

D_AVM = struct();
D_AVM.FileName = FileName;
D_AVM.analysis_window = janela;
D_AVM.sample_rate = Fs;
D_AVM.activation_frequency = activation_frequency;
D_AVM.mean_activation_frequency = mean_activation_frequency;
D_AVM.peak_to_peak_amplitudes = peak_to_peak_amplitudes;
D_AVM.mean_peak_to_peak_amplitudes = mean_peak_to_peak_amplitudes;

save(fullfile(output_folder, ['electric_data_' FileName '_analise_tensao.mat']), 'D_AVM', '-v7.3');

clear; clc; close all

%% FIM DO CÓDIGO