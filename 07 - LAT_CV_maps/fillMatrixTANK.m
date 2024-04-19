function LAT_matrix = fillMatrixTANK(LAT_E_interp)
    % FILLMATRIXTANK preenche uma matriz 25x25 com os valores de LAT_E_interp
    %
    %   LAT_matrix = FILLMATRIXTANK(LAT_E_interp) preenche uma matriz 25x25 com os valores
    %   de LAT_E_interp, organizando-os linha por linha da esquerda para a direita.
    %   Se LAT_E_interp contiver menos de 625 valores, a função retornará uma matriz vazia.
    %
    %   Exemplo:
    %       LAT_E_interp = rand(1000, 1);  % Exemplo de variável de entrada
    %       LAT_matrix = fillMatrixTANK(LAT_E_interp);  % Preenche a matriz 25x25
    %
    
    % Verifica se o tamanho de LAT_E_interp é adequado para preencher a matriz 25x25
    if length(LAT_E_interp) < 625
        disp('Erro: A matriz LAT_E_interp não contém dados suficientes para preencher a matriz 25x25.');
        LAT_matrix = [];
        return;
    end
    
    % Inicializa a matriz 25x25 com zeros
    LAT_matrix = zeros(25, 25);
    
    % Preenche a matriz linha a linha da esquerda para a direita com os valores de LAT_E_interp
    for i = 1:25
        for j = 1:25
            LAT_matrix(i, j) = LAT_E_interp((i - 1) * 25 + j);
        end
    end
end
