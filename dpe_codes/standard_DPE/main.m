close all;
clear; 
clc; 
format long;

%% ===================== Cenário / Geometria ================================

% --- Coordenadas aproximadas do CT/UFRJ (ajuste conforme o ponto exato) ---
lat = -22.8596582;      % [graus]
lon = -43.2303236;      % [graus]
h   = 10;               % [m] altura elipsoidal (estimativa)

UserPosition = llh2ecef(lat, lon, h);   % [x y z] ECEF (m)

M = 7;                      % nº de SVs na constelação
numSV = M;                  % redundante, mas usado adiante

% --- Gera geometria "boa" (PDOP baixo) com máscara de 5° e janela alvo -----
[SatPosition, SatPRN, azel_deg, PDOP] = randomConstellationPDOP(UserPosition, M, ...
    'good', 5, [1.5 1.55], 3000);
fprintf('PDOP (good): %.2f\n', PDOP);

% % --- Alternativa: geometria "ruim" (PDOP alto) ---------------------------
% [SatPosition, SatPRN, azel_deg, PDOP] = randomConstellationPDOP(UserPosition, M, ...
%     'good', 5, [5 5.2], 3000);
% fprintf('PDOP (bad): %.2f\n', PDOP);
%% ===================== Parâmetros de Simulação ============================

CNosim = 30:1:50;       % vetor de C/N0 em dB-Hz (cenários)
Nexpe  = 200;             % nº de repetições por C/N0 p/ estatística (RMSE)

% Alocação de memória p/ erros de posição (LS e DPE)
PosErrLS  = zeros(length(CNosim), Nexpe);
PosErrDPE = zeros(length(CNosim), Nexpe);

%% ===================== Parâmetros do Sinal / Front-end ====================

c  = 299792458;         % [m/s] velocidade da luz
f0 = 1575.42e6;         % [Hz] portadora L1 (não utilizada diretamente aqui)
CodePeriod = 1e-3;      % [s] período do C/A (1 ms)
fs  = 50e6;             % [Hz] taxa de amostragem do RX
dt  = 1/fs;             % [s] período de amostragem

% Estrutura de configuração de geração de sinal (réplicas, filtro, etc.)
gnss_input.CodePeriod = CodePeriod;
gnss_input.CoherentIntegrations     = 1;          % 1 ms coerente
gnss_input.NonCoherentIntegrations  = 1;          % uma acumulação não-coerente
gnss_input.m = 1; gnss_input.n = 1;               % parâmetros de escala do chip (usados p/ Tc/Ts)
gnss_input.type = 'BPSK';

% Atenção: Tc e Ts definidos a partir de 1.023 MHz (C/A)
gnss_input.Tc = 1/(gnss_input.n*1.023e6);         % [s] período de chip
gnss_input.Ts = 1/(2*gnss_input.m*1.023e6);       % [s] (não usado aqui diretamente)

% Front-end
gnss_input.fs   = fs; 
gnss_input.dt   = dt;
gnss_input.fn   = 2e6;                            % [Hz] corte do LPF
gnss_input.order= 36;                             % ordem do FIR
gnss_input.numSV= numSV;
gnss_input.SatPosition = SatPosition;
gnss_input.UserPosition = UserPosition;
gnss_input.c   = c;
gnss_input.SatPRN = SatPRN;
gnss_input.Nexpe  = Nexpe;

%% ===================== Parâmetros do LS em 2 etapas =======================

ls_input.UserPosition = UserPosition;             % inicialização
ls_input.SatPosition  = SatPosition;
ls_input.c  = c; 
ls_input.numSV = numSV; 
ls_input.fs = fs;
ls_input.num2stepsIterations = 10;                % iterações máx. do LS 4D

%% ===================== Parâmetros do DPE ==================================

dpe_input.UserPosition = UserPosition;
dpe_input.c  = c; 
dpe_input.numSV = numSV;
dpe_input.fs = fs; 
dpe_input.fn = gnss_input.fn;
dpe_input.Tc = 1/(gnss_input.n*1.023e6);          % [s] período de chip (consistente com geração)
dpe_input.SatPosition = SatPosition;

% Busca aleatória adaptativa
dpe_input.Niter       = 10000;     % nº iterações
dpe_input.contraction = 2;         % fator de contração do passo
dpe_input.dmax        = 10000;     % [m] passo espacial máx.
dpe_input.dmin        = 0.01;      % [m] passo espacial mín.
dpe_input.dmax_clk    = gnss_input.dt/10;   % [s] passo máx. relógio (não usado na sua versão atual)
dpe_input.dmin_clk    = gnss_input.dt/100;  % [s] passo mín. relógio (idem)

% Fator de normalização da métrica (compatível com sua implementação de r)
dpe_input.NormalizaFactor = sqrt(gnss_input.NonCoherentIntegrations) * ...
                            gnss_input.CoherentIntegrations * ...
                            gnss_input.CodePeriod * gnss_input.fs;

dpe_input.CN0_est_ind = 1;   % placeholder (não usado neste script)

%% ===================== Geração de Réplicas / FFT local =====================

sigen = signalGen(gnss_input);  % gera x_local, x_delay, fft_local, etc.

%% ===================== Loop de Simulação (C/N0 x Experimentos) ============

for CNo_idx = 1:length(CNosim)
    for exp_idx = 1:Nexpe
        % vetor C/N0 por SV (aqui todos iguais)
        CNo = CNosim(CNo_idx)*ones(numSV,1);

        % Sinal recebido: soma SVs escalados p/ C/N0 + ruído, e filtrado
        x = receivedSignal(sigen, CNo);

        % Correlação circular via FFT (1 ms) + acumulação não-coerente
        r = correlateSignal(sigen, x);

        % --- Estimador em 2 etapas (LS 4D) ----------------------------------
        pos_est_ls = conv2stepsPVT(r, ls_input);          % função externa
        PosErrLS(CNo_idx, exp_idx) = norm(pos_est_ls - UserPosition);

        % --- DPE (Adaptive Random Search) -----------------------------------
        pos_est_dpe = dpePvt(r, dpe_input);               % função externa
        PosErrDPE(CNo_idx, exp_idx) = norm(pos_est_dpe' - UserPosition);
    end
end

% RMSE por C/N0 (sobre as Nexpe repetições)
RMSE_LS  = sqrt(mean(PosErrLS.^2, 2));
RMSE_DPE = sqrt(mean(PosErrDPE.^2, 2));

%% ===================== CRB (convencional x DPE) ===========================

M    = numSV;                 % nº de SV
Tcoh = CodePeriod;            % [s] tempo de integração coerente (1 ms)

% --- Vetores LOS e matriz de geometria -------------------------------------
u_mat = zeros(M,3); 
rho   = zeros(M,1);
for i = 1:M
    vec       = SatPosition(i,:) - UserPosition;
    rho(i)    = norm(vec);
    u_mat(i,:)= vec ./ rho(i);        % LOS unitário
end
T_geo = [u_mat, ones(M,1)];          % matriz G (equivalente a H), Mx4
Pmat  = (1/c) * u_mat;               % projeção p/ posição (DPE aprox.)

% --- B2: energia da derivada temporal da réplica filtrada -------------------
s  = sigen.x_local(1,:);             % usa a mesma réplica usada na correlação
dt = 1/fs;                           % (sombreamento local de dt anterior — ok)
s  = s ./ sqrt(sum(abs(s).^2) * dt); % normaliza energia: ∫|s|^2 dt = 1
ds = diff(s) / dt;                   % derivada temporal discreta
B_2 = sum(abs(ds).^2) * dt;          % ≈ ∫ |s'(t)|^2 dt

% Alocação dos vetores de bound
fCRB_2SP = zeros(length(CNosim),1);  % bound p/ LS 4D (posição+clock, projeta p/ pos)
fCRB_DPE = zeros(length(CNosim),1);  % bound p/ DPE (apenas posição, aprox.)

for idx = 1:length(CNosim)
    CN0dB   = CNosim(idx);

    % Modelo baseband complexo: SNR = 2 * (C/N0) * Tcoh (linear)
    SNR_lin = 2 * 10^((CN0dB + 10*log10(Tcoh))/10);

    % FIM (Fisher) por satélite para atraso: J_tau = 2 * B2 * SNR
    % (assumindo SVs independentes) -> diagonal MxM
    Jtau = diag( 2 * B_2 * SNR_lin * ones(M,1) );

    % --- CRB 2-steps por propagação (ver sua equação Eq. 41) ---------------
    W    = eye(M);                        % sem ponderação (pode usar WLS com C/N0)
    G    = T_geo;
    Ainv = inv(G.'*W*G);                  % (H'H)^{-1}
    % cov([p; b/c]) >= c^2 * Ainv * G' * W * inv(Jtau) * W * G * Ainv
    C_pdt = c^2 * Ainv * (G.'*W*(Jtau\(W*G))) * Ainv;
    fCRB_2SP(idx) = sqrt(trace(C_pdt(1:3,1:3)));  % bound p/ posição (raiz da soma das variâncias)

    % --- CRB DPE (aprox. pré-correlação, posição apenas) --------------------
    Jp   = Pmat.' * Jtau * Pmat;          % FIM projetada p/ posição
    Cp   = inv(Jp);                       % cov mínima (posição)
    fCRB_DPE(idx) = sqrt(trace(Cp));
end

% Tabela resumo (RMSE vs CRB)
T_SUM = table(CNosim(:), RMSE_LS(:), fCRB_2SP(:), ...
                        RMSE_DPE(:), fCRB_DPE(:), ...
    'VariableNames', {'CN0_dBHz','RMSE_LS_m','CRB_2SP_m', ...
                      'RMSE_DPE_m','CRB_DPE_m'});
disp(T_SUM);

%% ===================== Gráfico (RMSE x Bound) =============================

figure; 
semilogy(CNosim, RMSE_LS,'b-.', ...
         CNosim, fCRB_2SP,'k--', ...
         CNosim, RMSE_DPE,'g-',  ...
         CNosim, fCRB_DPE,'m-.', 'LineWidth',2);
grid on; 
xlabel('C/N_0 [dB-Hz]'); 
ylabel('RMSE / Bound [m]');
legend('RMSE 2-steps','CRB 2-steps','RMSE DPE','CRB DPE','Location','southwest');
title('RMSE vs CRB (convencional x DPE)');
