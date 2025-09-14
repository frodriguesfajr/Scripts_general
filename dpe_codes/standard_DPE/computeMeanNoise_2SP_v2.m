% ===========================================

close all;
clear;
format long;

rng(42)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constantes físicas e portadora de interesse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c   = 299792458;      % [m/s] velocidade da luz
f0  = 1575.42e6;      % [Hz] frequência central L1 GPS (C/A)
%% ---- Parâmetros do sinal (GPS L1 C/A) ---------------------------------
SatPRN       = [12 15 17 19 24 25 32];   % PRNs escolhidos
UserPosition = [3.915394273911475e+06, 2.939638207807819e+05, ...
                5.009529661006817e+06];  % [m] posição ECEF do usuário
numSV = length(SatPRN); % nº de satélites usados efetivamente

CodePeriod = 1e-3; % [s] período do código C/A (1 ms)
CoherentIntegrations = 1; % nº integrações coerentes (coh)
NonCoherentIntegrations = 1; % nº integrações não-coerentes (ncoh)
m = 1; % fator p/ Ts (amostragem do correlator) — aqui simbólico
n = 1; % fator p/ Tc (taxa de chip) — aqui simbólico
type = 'BPSK'; % modulação do código (C/A é BPSK)
Tc = 1/(n*1.023e6);% [s/chip] duração de um chip do C/A (1/1.023 MHz)
Ts = 1/(2*m*1.023e6); % [s] step temporal correlator
fs = 50e6; % [Hz] freq. de amostragem do sinal em IF/baseband
dt = 1/fs; % [s] período de amostragem
fn = 2e6; % [Hz] banda do filtro passa-baixa (p/ front-end simulado)
order = 36; % ordem do FIR usado (quando filtrar sinais)
% Fator de normalização usado no pós-correlator (energia acumulada)
% sqrt(NCoh)*Coh*CodePeriod*fs ≈ escala da energia acumulada nas integrações
NormalizaFactor = sqrt(NonCoherentIntegrations) *...
    CoherentIntegrations * CodePeriod * fs;
%% ---- Parâmetros de cenário (posições ECEF dos satélites) ---------------
% Matriz 10x3 com posições ECEF (em metros). Você usa só as 7 primeiras.
corrSatPosition = 1.0e+07 * [
    2.061934439245598  -0.649651248625989   1.515031823419662
    2.652937217353064   0.225144209072909   0.245090150841460
   -0.036582056912767   1.531548249671478   2.204309909466485
    1.031402024072278   1.611770455054111   1.801481566078742
    1.660579169172477   0.314602454675085   2.047357395811576
    1.986519610236877  -1.663607551240955   0.522461680298307
   -0.265437875221575  -1.573610606389880   2.123719767033498
   -0.948907876055638  -1.438620519096841   2.407317325333366
    1.479532383323564   0.751773843931400   2.450026136317899
    1.777572228593643   2.192914078898261   0.888480393951951
];
numSV = 7;
SatPosition  = corrSatPosition(1:numSV, :);  % Seleciona os satélites ativos)
% OBS: As posições estão na ordem de 10^7 m (compatível com órbitas GNSS).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---- Parâmetros do DPE (busca direta no espaço de estado PVT) ----------
% A ideia da abordagem DPE é: para cada hipótese de estado γ (posição/clock),
% gerar réplicas consistentes (atraso/Doppler) e comparar via correlação.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dpe_param.gamma_guess   = UserPosition;  % chute inicial de posição (aqui: a posição real)
dpe_param.c = c;  % velocidade da luz (para converter atraso <-> distância)
dpe_param.fs = fs; % freq. de amostragem (impacta grade de busca/índices)
dpe_param.numSV = numSV;         % nº de satélites
dpe_param.Tc            = Tc;            % duração do chip (mapeia atraso fracionário <-> índice de chip)
dpe_param.SatPosition   = SatPosition;   % posições ECEF dos satélites
dpe_param.NormalizaFactor = NormalizaFactor;

% Parâmetros do algoritmo de busca (ARS – accelerated random search):
dpe_param.Niter        = 10000;  % nº de iterações de busca
dpe_param.dmax         = 10000;  % [m] passo máximo de perturbação na posição (exploração)
dpe_param.dt           = 1/fs;   % [s] período de amostragem (p/ consistência)
dpe_param.dmax_clk     = dt/10;  % [s] passo máximo de perturbação no clock (se for estimar viés)
dpe_param.dmin_clk     = dt/100; % [s] passo mínimo de perturbação no clock (para evitar estagnar)
dpe_param.contraction  = 2;      % fator de contração do passo (exploration -> exploitation)
dpe_param.dmin         = 0.01;   % [m] passo mínimo de perturbação em posição (reinicia quando muito pequeno)

% Alocação de buffers para o histórico da busca (traçar convergência/diagnóstico):
gamma_est      = zeros(3,  dpe_param.Niter+1); % estimativas sucessivas de posição (x,y,z)
amp_est        = zeros(numSV, dpe_param.Niter+1); % estimativas de amplitude por SV (se você infere a partir da correlação)
EstRxClkBias   = zeros(1,  dpe_param.Niter+1); % estimativas do viés de relógio do receptor (s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                              
%% ---- Simulation parameters"    
CNosim = 30:5:50; % Varredura de C/N0 (dB-Hz)
Nexpe  = 2; % Nº de repetições por ponto (média/RMSE estatístico)
%% Alocação de memória
PosErrLS   = zeros(length(CNosim), Nexpe);  % erros de posição (LS / 2SP) por C/N0 e experimento
PosErrDPE  = zeros(length(CNosim), Nexpe);  % erros de posição (DPE/ARS) por C/N0 e experimento
p = 3;           % colunas de T (3 se sem clock, 4 se incluir clock)
Tmat_vec  = zeros(length(CNosim), Nexpe, numSV, p);    % guarda T (M×p)
Wmat_vec  = zeros(length(CNosim), Nexpe, numSV, numSV);    % guarda W (M×M)
TWTx_vec  = zeros(length(CNosim), Nexpe, p, p);    % opcional: guarda T'WT (p×p)
%% Geração de sinal (sigen struct)
NsamplesLocal = CodePeriod*fs*CoherentIntegrations; % Nº de amostras por janela local (réplica)
NsamplesData  = CodePeriod*fs*CoherentIntegrations*NonCoherentIntegrations; % Nº de amostras do “dado” (coh * ncoh)  
%% Alocação de Memória
Range     = zeros(1,numSV);                 % distâncias geométricas usuário–satélite (cada SV)
x_local   = zeros(numSV,NsamplesLocal);     % réplicas locais dos códigos por SV (amostras)
fft_local = zeros(numSV,NsamplesLocal);     % FFT das réplicas (para correlação rápida no domínio da freq.)
x_delay   = zeros(numSV,NsamplesData);      % versões “atrasadas” dos códigos (sinal recebido ideal sem ruído)
%% Calcula range e fração do delay para cada satélite
for kSV = 1:numSV 
    Range(kSV) = norm(SatPosition(kSV,:) - UserPosition);  % norma Euclidiana -> distância em metros
end

% Atraso fracionário “verdadeiro” (em segundos, reduzido ao período do código)
FracDelay = mod(Range/c, CodePeriod);  
%% Gera replica local e sinais atrasados de acordo com os delays calculados
PrevNCOIndex = -FracDelay/Tc;   % índice fracionário inicial do NCO de código (em chips) para cada SV
randomDelay  = 0;               % (placeholder) jitter aleatório extra (se quiser simular clock do receptor)

for kSV = 1:numSV
    Code  = genCAcode(SatPRN(kSV));               % sequência PRN do SV k
    Tchip = CodePeriod / length(Code);            % duração de um chip “amostral” na construção da réplica
    
    % --- Réplica local (janela NsamplesLocal) ---
    ii = 1:NsamplesLocal;
    % Mapeia cada amostra para um índice de chip da sequência PRN:
    % round((ii-1)/fs/Tchip) -> índice inteiro de chip
    % mod(..., length(Code)) -> wrap-around na sequência (periodicidade do código)
    x_local(kSV,:) = Code(1 + mod(round((ii - 1)/fs/Tchip), length(Code)));

    % --- Código “atrasado” (sinal recebido ideal) para janela NsamplesData ---
    ii = 1:NsamplesData;
    % Aplica o deslocamento fracionário PrevNCOIndex(kSV) + randomDelay
    x_delay(kSV,:) = Code(1 + mod(round(PrevNCOIndex(kSV) + randomDelay + ...
                                        (ii - 1)/fs/Tchip), length(Code)));
end
%% Fitragem do sinal local and gera sua FFT 
% Front-end idealizado: um FIR passa-baixa (fn) aplicado aos sinais
% OBS: o mesmo filtro em x_local e x_delay deixa o “canal” consistente

wn = pi*fn/(pi*fs/2);       % normalização da frequência de corte (formato do fir1)
h  = fir1(order, wn);       % projeto do FIR (janela padrão do MATLAB)

for kSV = 1:numSV
    % filtfilt -> filtro zero-phase (evita atraso de grupo na forma de onda)
    x_delay(kSV,:)  = filtfilt(h, 1, x_delay(kSV,:));
    x_local(kSV,:)  = filtfilt(h, 1, x_local(kSV,:));
    % FFT da réplica local: será usada em correlação rápida via
    %   r = ifft( FFT(x)*conj(FFT(local)) ), depois magnitude^2
    fft_local(kSV,:) = fft(x_local(kSV,:), NsamplesLocal);
end
%% Normaliza Potência do sinal recebido após filtragem
% Normaliza cada “canal” x_delay para ter a mesma energia após o filtro.
% Isso facilita a definição de amplitude/SNR quando somar os SVs e adicionar ruído.
for kSV = 1:numSV
    % Escala tal que sum(x_delay.^2) = NsamplesData  (potência unitarizada)
    x_delay(kSV,:) = x_delay(kSV,:) * sqrt( NsamplesData / sum(x_delay(kSV,:).^2) );
end
%
ls_param.UserPosition = UserPosition;
ls_param.fs = fs;
ls_param.c = c;
ls_param.SatPosition = SatPosition;
sigen.Tc = Tc;
ls_param.num2stepsIterations = 10;
%% gather outputs in a struct
sigen.x_delay = x_delay;
sigen.fft_local = fft_local;
sigen.NsamplesLocal = NsamplesLocal;
sigen.NsamplesData = NsamplesData;
sigen.numSV = numSV;
sigen.NonCoherentIntegrations = NonCoherentIntegrations;
sigen.fn = fn;
sigen.order = order;
sigen.CodePeriod = CodePeriod;
sigen.CoherentIntegrations = CoherentIntegrations;
sigen.fs = fs;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Adiciona ruído AWGN aos sinais transmitidos
for exp_idx=1:Nexpe
    noise = ( sqrt(1/2)*randn(1,NsamplesData) + ...
        1i* sqrt(1/2)*randn(1,NsamplesData) );
    r = correlateSignal(sigen,noise);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Início da Simulação
% Loop externo: percorre os diferentes níveis de C/N0 simulados
for CNo_idx = 1:length(CNosim)
    disp([CNosim(CNo_idx)])   % Exibe no console o C/N0 atual (em dB-Hz)
    % Loop interno: repetições independentes (para calcular estatísticas)
    for exp_idx = 1:Nexpe
        % Cria um vetor com o mesmo valor de C/N0 para todos os satélites
        CNo = CNosim(CNo_idx) * ones(numSV,1);
        
        %% Geração do sinal recebido (sinal + ruído)
        % Função que soma o sinal de todos os satélites, ajusta a amplitude
        % conforme o C/N0 desejado e adiciona ruído AWGN (branco complexo).
        x = receivedSignal(sigen, CNo);
        r = correlateSignal(sigen,x);
        [EstRxPVT, Tmat, Wmat, TWTx]  = conv2stepsPVT(ls_param, r,numSV);
        PosErrLS(CNo_idx,exp_idx) = norm(EstRxPVT(1:3)-UserPosition);
        % Grava nas fatias corretas
        Tmat_vec(CNo_idx, exp_idx, :, :) = Tmat;   % (M×p)
        Wmat_vec(CNo_idx, exp_idx, :, :) = Wmat;   % (M×M)
        TWTx_vec(CNo_idx, exp_idx, :, :) = TWTx;   % (p×p)
        
    end
end
% Para ler de volta uma matriz específica (por exemplo, T do 2º experimento no 1º C/N0):
T_1_2 = squeeze(Tmat_vec(1, 2, :, :));   % vira M×p
W_1_2 = squeeze(Wmat_vec(1, 2, :, :));   % vira M×M
%% Cálculo do RMSE (Root Mean Square Error)
% Para cada nível de C/N0, calcula o erro quadrático médio de posição
% (média dos quadrados dos erros das repetições) e tira a raiz quadrada,
% resultando no RMSE em metros.
RMSE_LS = sqrt(mean(PosErrLS.^2, 2));
%% Cramer Rao Bound computation
SNR = 10.^((CNosim(:) + 10*log10(CodePeriod))/10);  % vetor
M = numSV;           % nº de satélites
u = (SatPosition - UserPosition) ./ vecnorm(SatPosition - UserPosition,2,2); % Mx3
T_true = [-u, ones(M,1)];     % Mx4
W_true = eye(M);              % ajuste se quiser ponderar
% B2: do seu local replica + mesmo filtro
B_2 = sum((diff(x_local(1,:))/dt).^2)/sum(x_local(1,:).^2);

fCRB_LS = zeros(length(CNosim),1);
for k = 1:length(CNosim)
    SNRk = 10^((CNosim(k) + 10*log10(CodePeriod))/10);
    J_tautau_inv = (1/(2*SNRk*B_2))*eye(M);
    G = c * ((T_true.'*W_true*T_true) \ (T_true.'*W_true));  % 4xM
    C_pdt = G * J_tautau_inv * G.';    % 4x4
    C_pos = C_pdt(1:3,1:3);
    fCRB_LS(k) = sqrt(trace(C_pos));
end



%% ===== Tabela CN0 vs RMSE_DPE, CRB_DPE e ZZB_DPE ========================
T_CRB_ZZB = table(CNosim(:), RMSE_LS(:), fCRB_LS(:), ...
                  'VariableNames', {'CN0_dBHz','RMSE_DPE_m','CRB_LS_m'});
disp(T_CRB_ZZB)
% return
%% ===== Plot comparando DPE x CRB-DPE x ZZB-DPE ==========================
figure;
h = semilogy(CNosim, RMSE_LS, 'b-.', ...
             CNosim, fCRB_LS, 'k-');
legend('2SP (RMSE)','CRB 2SP','fontsize',16);
grid on; set(h,'LineWidth',2);
xlabel('C/N_0 [dB-Hz]'); ylabel('Erro [m]');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CAcode = genCAcode(PRN)
% generateCAcode.m generates one of the 32 GPS satellite C/A codes.
%
% CAcode = generateCAcode(PRN)
%
%   Inputs:
%       PRN         - PRN number of the sequence.
%
%   Outputs:
%       CAcode      - a vector containing the desired C/A code sequence 
%                   (chips).  

%--------------------------------------------------------------------------
%                           SoftGNSS v3.0
% 
% Copyright (C) Darius Plausinaitis
% Written by Darius Plausinaitis
% Based on Dennis M. Akos, Peter Rinder and Nicolaj Bertelsen
%--------------------------------------------------------------------------
%This program is free software; you can redistribute it and/or
%modify it under the terms of the GNU General Public License
%as published by the Free Software Foundation; either version 2
%of the License, or (at your option) any later version.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License
%along with this program; if not, write to the Free Software
%Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
%USA.
%--------------------------------------------------------------------------

%CVS record:
%$Id: generateCAcode.m,v 1.1.2.1 2007/01/29 10:22:24 dpl Exp $

%--- Make the code shift array. The shift depends on the PRN number -------
% The g2s vector holds the appropriate shift of the g2 code to generate
% the C/A code (ex. for SV#19 - use a G2 shift of g2s(19) = 471)
g2s = [  5,   6,   7,   8,  17,  18, 139, 140, 141, 251, ...
       252, 254, 255, 256, 257, 258, 469, 470, 471, 472, ...
       473, 474, 509, 512, 513, 514, 515, 516, 859, 860, ...
       861, 862 ... end of shifts for GPS satellites 
       ... Shifts for the ground GPS transmitter are not included
       ... Shifts for EGNOS and WAAS satellites (true_PRN = PRN + 87)
                 145, 175,  52,  21, 237, 235, 886, 657, ...
       634, 762, 355, 1012, 176, 603, 130, 359, 595, 68, ...
       386];

%--- Pick right shift for the given PRN number ----------------------------
g2shift = g2s(PRN);

%--- Generate G1 code -----------------------------------------------------

%--- Initialize g1 output to speed up the function ---
g1 = zeros(1, 1023);
%--- Load shift register ---
reg = -1*ones(1, 10);

%--- Generate all G1 signal chips based on the G1 feedback polynomial -----
for i=1:1023
    g1(i)       = reg(10);
    saveBit     = reg(3)*reg(10);
    reg(2:10)   = reg(1:9);
    reg(1)      = saveBit;
end

%--- Generate G2 code -----------------------------------------------------

%--- Initialize g2 output to speed up the function ---
g2 = zeros(1, 1023);
%--- Load shift register ---
reg = -1*ones(1, 10);

%--- Generate all G2 signal chips based on the G2 feedback polynomial -----
for i=1:1023
    g2(i)       = reg(10);
    saveBit     = reg(2)*reg(3)*reg(6)*reg(8)*reg(9)*reg(10);
    reg(2:10)   = reg(1:9);
    reg(1)      = saveBit;
end

%--- Shift G2 code --------------------------------------------------------
%The idea: g2 = concatenate[ g2_right_part, g2_left_part ];
g2 = [g2(1023-g2shift+1 : 1023), g2(1 : 1023-g2shift)];

%--- Form single sample C/A code by multiplying G1 and G2 -----------------
CAcode = -(g1 .* g2);
end

function r = correlateSignal(sigen, x_delay_noise)

fft_local = sigen.fft_local;
NsamplesLocal = sigen.NsamplesLocal;
numSV = sigen.numSV;
NonCoherentIntegrations = sigen.NonCoherentIntegrations;
% memory allocation
r=zeros(numSV,NsamplesLocal);
% Perform NonCoherentIntegrations times non coherent integrations of
% CoherentIntegrations times coherent integrations.
for kSV=1:numSV
    for idx_nc = 1:NonCoherentIntegrations
        r(kSV,:) = r(kSV,:) + ...
            abs(ifft(fft(x_delay_noise(1,NsamplesLocal*(idx_nc-1) + ...
            1:NsamplesLocal*idx_nc),NsamplesLocal) .* ...
            conj(fft_local(kSV,:)))).^2;
    end
end
end

function x = receivedSignal(sigen,CNo)

x_delay = sigen.x_delay;
NsamplesData = sigen.NsamplesData;
numSV = sigen.numSV;
fs = sigen.fs;
fn = sigen.fn;
order = sigen.order;
%% memory allocation
x_delay_noise=zeros(numSV,NsamplesData);
%% Add AWGN noise to the transmitted signals
noise = ( sqrt(1/2)*randn(1,NsamplesData) + ...
    1i* sqrt(1/2)*randn(1,NsamplesData) );
for kSV=1:numSV
    if CNo(kSV)<100
        % Sets amplitude assuming complex-noise power equal to 1
        %For CNo >=100 no noise is added.
        
        A       = sqrt(10^(CNo(kSV)/10)/fs);
        x_delay_noise(kSV,:) = A * x_delay(kSV,:);
    else
        x_delay_noise(kSV,:)=x_delay(kSV,:);
    end
end
%% Add noice to received signal and filter
x = sum(x_delay_noise, 1)+noise;

wn=pi*fn/(pi*fs/2);
h=fir1(order,wn);
x  = filtfilt(h,1,x);
end

function answer = q(x) 
answer = erfc(x/sqrt(2))/2;
end

function x = wrap( x, x_max )

while( sum( abs(x) > x_max ) ~= 0)
    x(abs(x)>x_max)  =   x(abs(x)>x_max) - sign(x(abs(x)>x_max))*2*x_max;
end

end

function [EstRxPVT, Tmat, Wmat, TWTx] = conv2stepsPVT(ls_param, r, numSV)
% [EstRxPVT, Tmat, Wmat, TWTx]
%  EstRxPVT : vetor estimado [x y z (opcional: delta t)] no fim das iterações
%  Tmat     : matriz de geometria (line-of-sight), última iteração
%  Wmat     : matriz de pesos (aqui identidade)
%  TWTx     : matriz T' * W * T da última iteração

UserPosition = ls_param.UserPosition;
fs          = ls_param.fs;
c           = ls_param.c;
num2stepsIterations = ls_param.num2stepsIterations;
SatPosition = ls_param.SatPosition;

%---------------------------------------------------------
% Pré-alocação de memória
EstRange = zeros(1,numSV);

% Posição inicial com pequeno ruído
RefPos   = UserPosition + 10*(2*rand(3,1)-1)';
EstRxPVT = RefPos;

%--- 1) Estima atrasos fracionários a partir da correlação
[~, maxPos] = max(r,[],2);
maxPos      = maxPos - 1;
EstFracDelay  = maxPos/fs;
EstFracRange  = EstFracDelay * c;

%---------------------------------------------------------
% 2) Iterações LS
for kIter = 1 : num2stepsIterations

    % Monta matriz de geometria H (3 colunas se não incluir clock)
    H = zeros(numSV,3);
    for kSV = 1:numSV
        EstRange(kSV) = norm(SatPosition(kSV,:) - EstRxPVT(1:3));
        numH          = SatPosition(kSV,:) - EstRxPVT(1:3);
        H(kSV,1:3)    = - numH / norm(numH);
        % Se quiser clock: descomente a linha abaixo e
        % inicialize EstRxPVT como [RefPos,0]
        % H(kSV,4)    = 1;
    end

    % Resíduo de pseudodistância (removendo ambiguidade de 1 ms)
    corrP  = (EstFracRange - EstRange') / c;
    corrP_noAmbg = wrap( rem(corrP, 1e-3), 0.5e-3 );
    corrFracPseudorange = corrP_noAmbg * c;

    % LS update
    deltaPVT = ((H' * H) \ H') * corrFracPseudorange;
    EstRxPVT = EstRxPVT + deltaPVT.';
end

%---------------------------------------------------------
% 3) Retornar as matrizes da última iteração
Tmat = H;                     % matriz T (geometria) final
Wmat = eye(numSV);            % pesos identidade (W = I)
TWTx = Tmat' * Wmat * Tmat;   % produto T'WT

end
