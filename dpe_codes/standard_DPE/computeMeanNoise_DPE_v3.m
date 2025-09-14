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

%% gather outputs in a struct
sigen.fs = fs;
sigen.NormalizaFactor = NormalizaFactor;
sigen.x_delay = x_delay;
sigen.x_local = x_local;
sigen.fft_local = fft_local;
sigen.randomDelay = randomDelay;
sigen.NsamplesLocal = NsamplesLocal;
sigen.NsamplesData = NsamplesData;
sigen.numSV = numSV;
sigen.NonCoherentIntegrations = NonCoherentIntegrations;
NsamplesData = sigen.NsamplesData;
sigen.fs = fs;
sigen.fn = fn;
sigen.order = order;
sigen.CodePeriod = CodePeriod;
sigen.CoherentIntegrations = CoherentIntegrations;
sigen.UserPosition = UserPosition;
sigen.fs = fs;
sigen.c = c;
sigen.SatPosition = SatPosition;
sigen.Tc = Tc;
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
        
        %% Correlação (integração coerente e não-coerente)
        % Correlaciona o sinal recebido com as réplicas locais (códigos C/A),
        % realizando as integrações coerentes e não-coerentes configuradas.
        % A saída 'r' contém as energias de correlação para cada SV
        % ao longo dos atrasos de código.
        r = correlateSignal(sigen,x);
        
        %% Estimativa direta de posição (DPE) via ARS
        % Aplica o método DPE (Direct Position Estimation) usando
        % ARS (Accelerated Random Search), que busca diretamente a posição
        % que maximiza a métrica de correlação.
        gamma_est = DPE_est(r, dpe_param);     
        % Calcula o erro de posição em metros (distância Euclidiana entre
        % a posição estimada e a posição real do usuário).
        PosErrDPE(CNo_idx,exp_idx) = norm(gamma_est - UserPosition);
    end
end
%% Cálculo do RMSE (Root Mean Square Error)
% Para cada nível de C/N0, calcula o erro quadrático médio de posição
% (média dos quadrados dos erros das repetições) e tira a raiz quadrada,
% resultando no RMSE em metros.
RMSE_DPE = sqrt(mean(PosErrDPE.^2, 2));

%% Cramer Rao Bound computation
% B_2: “largura efetiva” do espectro do código (razão entre energia da
% derivada temporal e energia do próprio código). Para BPSK no tempo discreto,
% usamos diferenças finitas: diff(x)/dt aproxima dx/dt.
B_2 = sum((diff(x_local(1,:))/dt).^2) / sum(x_local(1,:).^2);

T  = CodePeriod;     % tempo de integração coerente (1 ms para C/A)
D  = T*c;            % distância percorrida pela luz em T (escala para ambiguidades de tempo)
M  = numSV;          % número de satélites/canais

% P: “geometria” (Jacobiano de atraso em relação à posição), em [s/m].
% Para cada SV i, d(tau_i)/d(p) = -(u_i^T)/c. Aqui usamos (User - Sat) / ||User-Sat||
% (o sinal muda com a convenção; no CRB ele entra em P'JP, o que elimina o sinal).
P = 1/c * (UserPosition - SatPosition) ./ sqrt(sum((UserPosition - SatPosition).^2, 2));


varT = T^2/12;

% Converte C/N0 (dB-Hz) para SNR de integração:
% SNR(dB) = C/N0(dB-Hz) + 10log10(T)  => SNR linear = 10^(SNR(dB)/10).
SNRdb = CNosim + 10*log10(T);

% Pré-alocação de vetores para guardar as curvas (se for usar ZZB depois)
fZZLB_DPE = zeros(1, length(CNosim));
fZZB_DPE  = zeros(1, length(CNosim));

%% === CRB DPE ===============================================
% CRB_DPE aqui é obtido propagando o CRB dos atrasos para posição,
% isto é, J_tau = 2*SNR*B_2*I_M e J_p = P' * J_tau * P.
fCRB_DPE  = zeros(1, length(CNosim));   % RMSE 3D teórico (m)





for k = 1:length(SNRdb)
    SNR = 10^(SNRdb(k)/10);

     % --- FIM dos atrasos (cada SV): J_tau(ii) = 2 * SNR * B_2
    Jtau = diag( 2*SNR*B_2 * ones(M,1) );

    % --- Projeção para posição: J_p = P' * J_tau * P
    Jp  = P' * Jtau * P;            % [1/m^2]
    CRB = inv(Jp);                  % [m^2]

    % Métrica escalar (RMSE 3D)
    fCRB_DPE(k) = sqrt( trace(CRB) );              % [m]
end

%% ===== Tabela CN0 vs RMSE_DPE, CRB_DPE e ZZB_DPE ========================
T_CRB_ZZB = table(CNosim(:), RMSE_DPE(:), fCRB_DPE(:), ...
                  'VariableNames', {'CN0_dBHz','RMSE_DPE_m','CRB_DPE_m'});
disp(T_CRB_ZZB)
% return
%% ===== Plot comparando DPE x CRB-DPE x ZZB-DPE ==========================
figure;
h = semilogy(CNosim, RMSE_DPE, 'b-.', ...
             CNosim, fCRB_DPE, 'k-');
legend('DPE (RMSE)','CRB-DPE','fontsize',16);
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

function gamma_est_dpe = DPE_est(r, dpe_param)
% DPE_est
% Estima a posição (ECEF) via DPE usando uma busca aleatória acelerada (ARS)
% maximizando a soma das correlações por satélite.
%
% ENTRADAS
%   r          : matriz numSV x NsamplesLocal com a energia de correlação
%                (já integrada coerente/não-coerente) por SV ao longo do atraso
%   dpe_param  : struct com parâmetros do cenário (fs, Tc, SatPosition, etc.)
%
% SAÍDA
%   gamma_est_dpe : posição ECEF estimada [1x3] (m)

numSV = dpe_param.numSV;
Tc = dpe_param.Tc;
SatPosition = dpe_param.SatPosition;
c = dpe_param.c;
fs = dpe_param.fs;
NormalizaFactor = dpe_param.NormalizaFactor;
gamma_guess = dpe_param.gamma_guess;   % chute inicial para posição
dt = dpe_param.dt;                     % passo de amostragem (1/fs)
dmax = dpe_param.dmax;                 % passo máx. de deslocamento espacial (m) para ARS
dmax_clk = dpe_param.dmax_clk;         % passo máx. para clock (não usado aqui)
Niter = dpe_param.Niter;               % nº de iterações da ARS
dmin = dpe_param.dmin;                 % passo mínimo (quando atinge, volta a dmax)
dmin_clk = dpe_param.dmin_clk;         % idem para clock (não usado)
contraction = dpe_param.contraction;   % fator de contração do passo se não melhora

randomDelay = 0; % “número mágico”: atraso adicional (em períodos de código). Mantido em 0 aqui.

%% Alocação de memória
EstRange=zeros(1,numSV);  % distâncias sat-usuário (m) para uma hipótese de posição
MaxCorr=zeros(1,numSV);   % correlação máxima (pega o valor de r no atraso previsto por cada SV)

% Posição inicial da busca = chute + ruído uniforme (±100 m em cada eixo)
gamma_est(:,1) = gamma_guess + 100*(2*rand(3,1)-1)';

% Estimativa (inicial) do viés de relógio do receptor em períodos de código
% Aqui fixado por randomDelay (0). O algoritmo, mais abaixo, também não atualiza clock.
EstRxClkBias(:,1)=-randomDelay*Tc;

% Calcula as ranges previstas para o chute inicial
for kSV  =   1 : numSV
    EstRange(kSV) = norm(SatPosition(kSV,:) - gamma_est(:,1)'); % ||p_sv - p_user||
end

% Converte ranges em FRAÇÃO de período de código (atraso fracionário)
%   EstFracDelay = (tau_i + clockBias + dt) mod 1ms
% Observação: somar dt aqui é um “offset” que alinha discretização com sua grade de amostras.
EstFracDelay=mod(EstRange/c+EstRxClkBias(:,1)+dt,1e-3);





% Para cada SV, transforma delay fracionário em índice de amostra e lê a correlação r
for kSV = 1:numSV
    aux=round(EstFracDelay(kSV)*fs);  % atraso em segundos -> índice
    aux(aux==0)=1e-3*fs;              % evita índice 0 (usa 1 ms em amostras)
    MaxCorr(kSV)=r(kSV,aux);          % pega a energia de correlação naquele atraso
end

% Custo a maximizar: soma das correlações dos SVs
J_ant=sum(MaxCorr);

% Estimativa de amplitude (não usada na atualização), só armazenada
amp_est(:,1) = MaxCorr./NormalizaFactor^2;

% Inicializa passos de exploração (posição e clock)
d = dmax;
d_clk=dmax_clk;

% ====================== Loop da ARS (accelerated random search) ======================
for it = 1:Niter-1
    % 1) Propõe um novo ponto aleatório em torno da posição corrente
    rand_point = gamma_est(:,it) + d*(2*rand(3,1)-1);

    % 2) (Opcional) Atualizaria clock; aqui foi fixado em zero
    % rand_clk = EstRxClkBias(:,it)+d_clk*(2*rand-1);
    rand_clk = 0;
    

    % 3) Recalcula ranges e atrasos fracionários para o ponto candidato
    for kSV = 1:numSV
        EstRange(kSV) = norm(SatPosition(kSV,:) - rand_point');
    end
    EstFracDelay=mod(EstRange/c+rand_clk+dt,1e-3);

    % 4) Lê a correlação r naquelas hipóteses de atraso para cada SV
    for kSV =   1:numSV
        aux=round(EstFracDelay(kSV)*fs);
        aux(aux==0)=1e-3*fs;
        MaxCorr(kSV)=r(kSV,aux);
    end

    % 5) Avalia o custo total para o candidato
    J = sum(MaxCorr);

    % 6) Critério de aceitação:
    %    - Se melhorou, aceita o ponto e reseta passo para dmax
    %    - Senão, mantém ponto anterior e contrai o passo (exploração mais “fina”)
    if J > J_ant
        gamma_est(:,it+1) = rand_point;
        EstRxClkBias(:,it+1)=rand_clk;
        amp_est(:, it+1) = MaxCorr./NormalizaFactor^2;
        J_ant = J;
        d = dmax;
        d_clk=dmax_clk;
    else
        gamma_est(:,it+1) = gamma_est(:,it);
        EstRxClkBias(:,it+1)=EstRxClkBias(:,it);
        amp_est(:,it+1) = amp_est(:,it);
        d = d/contraction;
        d_clk=d_clk/contraction;
    end

    % 7) Se o passo ficou pequeno demais, “reaquece” (volta ao dmax) para escapar de ótimos locais
    if d < dmin
        d = dmax;
    end
    if d_clk < dmin_clk
        d_clk =dmax_clk;
    end
end

% Saída: última posição aceita (em metros, ECEF)
gamma_est_dpe = gamma_est(:,it+1)';

end
