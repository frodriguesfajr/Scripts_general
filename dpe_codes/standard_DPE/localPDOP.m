function PDOP = localPDOP(UserPosition, SatPosition)
% LOCALPDOP  Calcula o Position Dilution of Precision (PDOP)
%            para uma geometria de satélites GNSS.
%
%   PDOP = localPDOP(UserPosition, SatPosition)
%
%   UserPosition : [1x3] vetor ECEF [m] da posição do usuário.
%   SatPosition  : [Mx3] matriz ECEF [m] com as posições dos M satélites.
%
%   PDOP         : escalar com o índice PDOP. Se M<4 ou geometria singular,
%                  retorna NaN.
%
%   Definição: PDOP = sqrt( trace( (Hᵀ H)⁻¹ (1:3,1:3) ) )
%   onde H = [-u  1], u = vetores unitários LOS (linha-de-visada).

    M = size(SatPosition,1);                  % nº de satélites disponíveis
    if M < 4                                   % pelo menos 4 satélites
        PDOP = NaN;
        return;
    end

    % Vetores dos satélites ao usuário (linha-de-visada)
    r = SatPosition - UserPosition;           % [Mx3] vetor diferença

    % Distâncias (norma Euclidiana de cada vetor LOS)
    d = vecnorm(r,2,2);                        % [Mx1]
    if any(d == 0)                             % evita divisão por zero
        PDOP = NaN;
        return;
    end

    % Vetores unitários de linha-de-visada (da posição do usuário p/ satélite)
    u = r ./ d;                                % [Mx3], cada linha = u_i

    % Matriz de geometria H:
    %   primeira 3 colunas: -u (sinal negativo porque range cresce
    %   quando o usuário se afasta do satélite);
    %   quarta coluna: 1 (derivada em relação ao viés de relógio).
    H = [-u, ones(M,1)];                       % [Mx4]

    % Matriz normal de mínimos quadrados A = Hᵀ H
    A = H.' * H;                               % [4x4]

    % Checa condição numérica para inversão
    if rcond(A) < 1e-12                        % se mal-condicionada, retorna NaN
        PDOP = NaN;
        return;
    end

    % Matriz de covariância teórica dos parâmetros (posição + clock bias)
    C = inv(A);                                % [4x4]

    % PDOP: raiz da soma das variâncias das 3 coordenadas de posição
    PDOP = sqrt( trace( C(1:3,1:3) ) );
end
