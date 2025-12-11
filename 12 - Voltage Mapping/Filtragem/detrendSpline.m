function [s_det, s_pp] = detrendSpline(signal, t, w_length, drawflag, mensaje)
    % Função que remove o ruído de linha de base usando splines cúbicos
    % [s_det, s_pp] = detrendSpline(signal, t, w_length, drawflag, mensaje)
    %
    % Parâmetros:
    % - signal: Sinal de entrada.
    % - t: Vetor de tempo correspondente ao sinal.
    % - w_length: Comprimento da janela para o ajuste do spline (em segundos).
    % - drawflag: Ativa (1) ou desativa (0) a visualização gráfica.
    % - mensaje: Mensagem opcional para o título do gráfico.

    % Verifica se o vetor de tempo começa em zero
    t = t - t(1);

    % Número total de amostras
    L_s = length(signal);

    % Número de segmentos (janelas)
    numSeg = floor(t(end) / w_length);

    % Verifica se há segmentos suficientes
    if numSeg < 2
        error('O comprimento da janela (w_length) é muito grande. Reduza w_length para criar pelo menos 2 segmentos.');
    end

    % Inicializa vetores para os pontos médios
    t_m = zeros(numSeg, 1);
    s_m = zeros(numSeg, 1);

    % Calcula os pontos médios de cada janela
    for m = 1:numSeg
        % Índices da janela atual
        ind_seg = (t >= (m-1)*w_length) & (t < m*w_length);
        % Verifica se há dados na janela
        if sum(ind_seg) == 0
            t_m(m) = NaN;
            s_m(m) = NaN;
            continue;
        end
        % Calcula o ponto médio da janela
        t_aux = t(ind_seg);
        t_m(m) = t_aux(round(length(t_aux)/2));
        s_m(m) = mean(signal(ind_seg));
    end

    % Remove segmentos sem dados
    t_m(isnan(t_m)) = [];
    s_m(isnan(s_m)) = [];

    % Verifica se há pontos suficientes para o spline
    if length(t_m) < 2
        error('Não há pontos suficientes para ajustar o spline. Reduza w_length ou verifique os dados.');
    end

    % Ajusta o spline cúbico
    pp = csaps(t_m, s_m);

    % Avalia o spline no vetor de tempo original
    s_pp = ppval(pp, t);

    % Remove a tendência do sinal
    s_det = signal - s_pp;

    % Visualização gráfica (se drawflag = 1)
    if drawflag
        figure();
        subplot(2, 1, 1);
        plot(t, signal, 'k', 'DisplayName', 'Sinal Original');
        hold on;
        plot(t, s_pp, 'r-.', 'DisplayName', 'Spline Ajustado');
        legend('show');
        xlabel('Tempo (s)');
        ylabel('Amplitude');
        title(mensaje);

        subplot(2, 1, 2);
        plot(t, s_det, 'k', 'DisplayName', 'Sinal Detrended');
        legend('show');
        xlabel('Tempo (s)');
        ylabel('Amplitude');
        title('Sinal após Remoção da Tendência');
    end
end