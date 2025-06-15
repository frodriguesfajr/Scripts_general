function x = my_ifft(X, varargin)
    % my_ifft - Implementação manual da IFFT com suporte à opção 'symmetric'
    
    N = length(X);
    X = X(:);  % Garante coluna
    
    % Verifica se foi passada a opção 'symmetric'
    symmetric = nargin > 1 && ischar(varargin{1}) && strcmpi(varargin{1}, 'symmetric');

    % Construir matriz de IFFT
    n = 0:N-1;
    k = n';
    W = exp(1j * 2 * pi * k * n / N);  % Sinal positivo para IFFT

    % Calcular IFFT
    x = (1/N) * W * X;

    % Se for 'symmetric', forçar parte real
    if symmetric
        x = real(x);
    end

    % Retorna resultado no mesmo formato de X (linha ou coluna)
    if isrow(varargin{1}) || isrow(X)
        x = x.';
    end
end
