function [sigen] = signalGen(gnss_input)
% SIGNALGEN  Gera réplicas de código (local e atrasada) + FFT local
%            para um conjunto de SVs, aplicando um LPF FIR e normalizando potência.
%            Útil para pipelines de correlação/FFT em simulação GNSS.

% --------- Parâmetros de entrada (estrutura) ---------------------------------
CodePeriod              = gnss_input.CodePeriod;             % período do código (s), ex.: 1 ms (C/A)
CoherentIntegrations    = gnss_input.CoherentIntegrations;   % nº integrações coerentes (bloco FFT)
NonCoherentIntegrations = gnss_input.NonCoherentIntegrations;% nº integrações não-coerentes (acúmulo)
Tc                      = gnss_input.Tc;                     % período de chip (s) (ver nota nas "melhorias")
fs                      = gnss_input.fs;                     % taxa de amostragem (Hz)
fn                      = gnss_input.fn;                     % freq. de corte do LPF (Hz)
order                   = gnss_input.order;                  % ordem do FIR
numSV                   = gnss_input.numSV;                  % nº de satélites
SatPosition             = gnss_input.SatPosition;            % [numSV x 3] posições ECEF dos SVs (m)
UserPosition            = gnss_input.UserPosition;           % [1 x 3] posição ECEF do usuário (m)
c                       = gnss_input.c;                      % velocidade da luz (m/s)
SatPRN                  = gnss_input.SatPRN;                 % vetor de PRNs (1..32 etc.)
Nexpe                   = gnss_input.Nexpe;                  % nº de "experimentos" para queimar RNG

% --------- Tamanhos dos vetores gerados --------------------------------------
NsamplesLocal = CodePeriod*fs*CoherentIntegrations;           % amostras na réplica local (bloco FFT)
NsamplesData  = CodePeriod*fs*CoherentIntegrations* ...
                NonCoherentIntegrations;                      % amostras na sequência atrasada acumulada

% --------- Pré-alocação ------------------------------------------------------
Range      = zeros(1,numSV);               % distâncias geométricas SV-usuário (m)
x_local    = zeros(numSV,NsamplesLocal);   % réplicas locais (por SV) antes/apos filtro
fft_local  = zeros(numSV,NsamplesLocal);   % FFTs das réplicas locais (por SV)
x_delay    = zeros(numSV,NsamplesData);    % réplicas atrasadas (por SV) antes/depois do LPF

% --------- Distâncias e atrasos fracionários ---------------------------------
for kSV = 1:numSV
    Range(kSV) = norm(SatPosition(kSV,:) - UserPosition); % distância Euclidiana (m)
end
FracDelay   = mod(Range/c, CodePeriod);   % atraso fracionário modulo período de código (s)
PrevNCOIndex = -FracDelay / Tc;           % índice inicial (em chips), sinal negativo = atraso
                                          % (ver nota: se Tc ~= Tchip calculado abaixo, pode haver inconsistência)

randomDelay = 0;                           % deslocamento adicional aleatório (aqui fixo a 0)

% --------- Geração das réplicas (local e atrasada) ----------------------------
for kSV = 1:numSV
    Code  = genCAcode(SatPRN(kSV));                     % sequência C/A do PRN (tip. ±1 ou {0,1})
    Tchip = CodePeriod / length(Code);                  % período de chip calculado da sequência

    % --- Réplica local (bloco de NsamplesLocal) ---
    ii = 1:NsamplesLocal;
    % mapeia amostra -> índice de chip (inteiro), repetindo a sequência por mod
    x_local(kSV,:) = Code( 1 + mod( round((ii-1)/(fs*Tchip)), length(Code) ) );

    % --- Réplica atrasada (bloco de NsamplesData) ---
    ii = 1:NsamplesData;
    % aplica deslocamento em "chips" PrevNCOIndex(kSV)+randomDelay e repete por mod
    x_delay(kSV,:) = Code( 1 + mod( round( PrevNCOIndex(kSV) + randomDelay + ...
                              (ii-1)/(fs*Tchip) ), length(Code) ) );
end

% --------- LPF FIR (passa-baixas) e FFT da réplica local ---------------------
wn = pi*fn/(pi*fs/2);     % normalização p/ fir1: equivale a 2*fn/fs (0..1, onde 1 = Nyquist)
h  = fir1(order, wn);     % FIR linear-fase (janela padrão do fir1, ex.: Hamming)

for kSV = 1:numSV
    % filtfilt -> resposta zero-fase (remove atraso de grupo), duplica ordem efetiva
    x_delay(kSV,:)  = filtfilt(h, 1, x_delay(kSV,:));
    x_local(kSV,:)  = filtfilt(h, 1, x_local(kSV,:));

    % FFT da réplica local (para correlação no domínio da frequência)
    fft_local(kSV,:) = fft(x_local(kSV,:), NsamplesLocal);
end

% --------- Normalização de potência das réplicas atrasadas -------------------
% Escala para que ||x_delay||^2 = NsamplesData  -> potência média ~ 1 por amostra
for kSV = 1:numSV
    x_delay(kSV,:) = x_delay(kSV,:) * sqrt( NsamplesData / sum( x_delay(kSV,:).^2 ) );
end

% --------- "Queima" de amostras de ruído p/ alinhar RNG entre execuções -----
% (não usado adiante; apenas garante reprodutibilidade comparando diferentes runs)
for exp_idx = 1:Nexpe
    noise = ( sqrt(1/2)*randn(1,NsamplesData) + 1i*sqrt(1/2)*randn(1,NsamplesData) ); %#ok<NASGU>
end

% --------- Saída --------------------------------------------------------------
sigen.x_delay                 = x_delay;
sigen.x_local                 = x_local;
sigen.fft_local               = fft_local;
sigen.randomDelay             = randomDelay;
sigen.NsamplesLocal           = NsamplesLocal;
sigen.NsamplesData            = NsamplesData;
sigen.NonCoherentIntegrations = NonCoherentIntegrations;
sigen.fs = fs;  sigen.fn = fn;  sigen.order = order;
sigen.numSV = numSV;  sigen.Nexpe = Nexpe;
end
