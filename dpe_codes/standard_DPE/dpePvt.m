function pos_est_dpe = dpePvt(r, dpe_input)
% DPEPVT  (Direct Position Estimation PVT) – busca aleatória adaptativa em R^3
%         para maximizar a soma das energias de correlação por SV.
%         Retorna a melhor posição estimada (ECEF) ao final das iterações.
%
% ENTRADAS
%   r            : matriz [numSV x N] com |R_i[n]|^2 (energia da correlação) por SV
%                  (N = NsamplesLocal = nº de amostras em 1 ms)
%   dpe_input    : struct com parâmetros e dados (ver extrações abaixo)
%
% SAÍDA
%   pos_est_dpe  : vetor [3x1] (ECEF, metros) com a posição estimada

% ----------------- Extrai parâmetros/constantes do input ---------------------
UserPosition   = dpe_input.UserPosition;   % [1x3] ou [3x1], posição inicial/âncora (ECEF m)
c              = dpe_input.c;              % velocidade da luz (m/s)
numSV          = dpe_input.numSV;          % número de satélites utilizados
fs             = dpe_input.fs;             % taxa de amostragem (Hz)
Tc             = dpe_input.Tc;             % período de chip (s) (ver observações)
SatPosition    = dpe_input.SatPosition;    % [numSV x 3] posições ECEF dos SVs (m)
Niter          = dpe_input.Niter;          % número total de iterações
contraction    = dpe_input.contraction;    % fator de contração do passo (>1)
dmax           = dpe_input.dmax;           % passo espacial máximo (m)
dmin           = dpe_input.dmin;           % passo espacial mínimo (m)
dmax_clk       = dpe_input.dmax_clk;       % passo máx. do viés de relógio (s) (não usado aqui)
dmin_clk       = dpe_input.dmin_clk;       % passo mín. do viés de relógio (s) (não usado aqui)
NormalizaFactor= dpe_input.NormalizaFactor;% fator p/ normalização de amplitude (ver abaixo)

dt             = 1/fs;                     % resolução temporal por amostra (s)
randomDelay    = 0;                         % deslocamento adicional (aqui 0)

% ----------------- Pré-aloca trajetórias/armazenamento -----------------------
gamma_est      = zeros(3, Niter+1);        % trajetória de posição estimada (3D)
amp_est        = zeros(numSV, Niter+1);    % "amplitudes" por SV (métrica derivada de r)
EstRxClkBias   = zeros(1, Niter+1);        % trajetória do viés de relógio (s)

EstRange       = zeros(1, numSV);          % distâncias geométricas SV-receptor (m)
MaxCorr        = zeros(1, numSV);          % melhor valor r_i[..] por SV na iteração

% ----------------- Inicialização --------------------------------------------
% posição inicial aleatória em torno de UserPosition (cubo de lado 200 m):
% (se preferir usar exatamente a posição do usuário, substitua pela linha comentada)
% gamma_est(:,1) = UserPosition + 100*(2*rand(3,1)-1)';       % versão antiga
gamma_est(:,1) = UserPosition(:) + 100*(2*rand(3,1)-1);       % [3x1]

% viés inicial do relógio (em segundos): aqui parte de zero (randomDelay=0)
EstRxClkBias(:,1) = -randomDelay * Tc;

% ----------------- Avalia função objetivo na posição inicial -----------------
for kSV = 1:numSV
    EstRange(kSV) = norm( SatPosition(kSV,:) - gamma_est(:,1).' ); % distância (m)
end

% Estima atraso fracionário por SV a partir da hipótese de posição:
%   EstFracDelay = frac( Range/c + clock_bias + dt, 1 ms )
% OBS: "+ dt" desloca 1 amostra; serve como "nudge" para evitar índice zero.
%      Veja observações ao final sobre o uso de "+dt".
EstFracDelay = mod( EstRange/c + EstRxClkBias(:,1) + dt, 1e-3 );

for kSV = 1:numSV
    % mapeia o atraso fracionário (s) para índice de amostra (1..N), N = 1e-3*fs
    aux = round( EstFracDelay(kSV) * fs );   % índice nominal (0..N)
    % evita índice zero; manda para "N" quando for 0
    aux(aux == 0) = 1e-3 * fs;               % assume N = fs*1ms (inteiro)
    MaxCorr(kSV) = r(kSV, aux);              % pega energia de correlação no índice
end

% Função objetivo J = soma das energias de pico por SV
J_ant        = sum(MaxCorr);
amp_est(:,1) = MaxCorr ./ NormalizaFactor^2; % escala opcional (documente o porquê)

% passos atuais da busca (espacial e de relógio)
d     = dmax;
d_clk = dmax_clk;

% ----------------- Loop de busca aleatória adaptativa ------------------------
for it = 1:Niter-1

    % gera novo ponto candidato em um cubo de aresta 2*d ao redor do atual
    rand_point = gamma_est(:,it) + d * (2*rand(3,1) - 1);

    % relógio: aqui mantido fixo (não há busca conjunta em b neste código)
    rand_clk = 0;  % se quiser buscar b, sorteie em [-d_clk, d_clk] e some ao viés atual

    % avalia J no ponto candidato
    for kSV = 1:numSV
        EstRange(kSV) = norm( SatPosition(kSV,:) - rand_point.' );
    end
    EstFracDelay = mod( EstRange/c + rand_clk + dt, 1e-3 );

    for kSV = 1:numSV
        aux = round( EstFracDelay(kSV) * fs );
        aux(aux == 0) = 1e-3 * fs;
        MaxCorr(kSV)  = r(kSV, aux);
    end
    J = sum(MaxCorr);

    % critério de aceitação: greedy (aceita se melhorou J)
    if J > J_ant
        gamma_est(:, it+1)   = rand_point;
        EstRxClkBias(:,it+1) = rand_clk;
        amp_est(:, it+1)     = MaxCorr ./ NormalizaFactor^2;
        J_ant = J;

        % “reinicia” passo após melhoria
        d     = dmax;
        d_clk = dmax_clk;
    else
        % sem melhoria: mantém ponto anterior e contrai passo
        gamma_est(:, it+1)   = gamma_est(:, it);
        EstRxClkBias(:,it+1) = EstRxClkBias(:, it);
        amp_est(:, it+1)     = amp_est(:, it);

        d     = d     / contraction;
        d_clk = d_clk / contraction;
    end

    % recicla passo se ficou pequeno demais (salto grosseiro/annealing)
    if d     < dmin,     d     = dmax;     end
    if d_clk < dmin_clk, d_clk = dmax_clk; end
end

% ----------------- Saída: última posição aceita --------------------------------
pos_est_dpe = gamma_est(:, it+1);   % [3x1] melhor posição ao final
end
