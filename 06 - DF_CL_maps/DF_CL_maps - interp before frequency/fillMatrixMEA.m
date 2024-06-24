function LAT_matrix = fillMatrixMEA(LAT_E_interp)
    % FILLMATRIXMEA preenche uma matriz 11x11 com os primeiros valores de LAT_E_interp
    %
    %   LAT_matrix = FILLMATRIXMEA(LAT_E_interp) preenche uma matriz 11x11 com os primeiros
    %   121 valores de LAT_E_interp, organizando-os linha por linha da esquerda para a direita.
    %   Se LAT_E_interp contiver menos de 121 valores, a função retornará uma matriz vazia.
    %
    %   Exemplo:
    %       LAT_E_interp = rand(221, 1);  % Exemplo de matriz de entrada
    %       LAT_matrix = fillMatrixMEA(LAT_E_interp);  % Preenche a matriz 11x11
    %
    %   Veja também: OUTRAS_FUNCOES_RELACIONADAS (se houver)
    
    % Verifica se o tamanho de LAT_E_interp é adequado para preencher a matriz 11x11
    if length(LAT_E_interp) < 121
        disp('Erro: A matriz LAT_E_interp não contém dados suficientes para preencher a matriz 11x11.');
        LAT_matrix = [];
        return;
    end
    
    % Inicializa a matriz 11x11 com zeros
    LAT_matrix = zeros(11, 11);
    
    % Preenche a matriz linha a linha da esquerda para a direita com os valores de LAT_E_interp
    for i = 1:11
        for j = 1:11
            LAT_matrix(i, j) = LAT_E_interp((i - 1) * 11 + j);
        end
    end
end
