close all;
clear; 
clc; 
format long;
% rng(42)

%% Scenario Parameters
% % SatPRN=  [2 12 15 17 19 24 25 28 32];
% M = 15;  % número de SVs
% UserPosition = [3.915394273911475e+06 2.939638207807819e+05 5.009529661006817e+06];
% [SatPosition, SatPRN, azel_deg, PDOP] = randomConstellation(UserPosition, M, 5, [1.5 4.5], 1000);
% 
% 
% fprintf('Constelação gerada: PDOP = %.2f\n', PDOP);
% disp(table(SatPRN(:), azel_deg(:,1), azel_deg(:,2), 'VariableNames', {'PRN','Az_deg','El_deg'}));
% Coordenadas aproximadas do CT/UFRJ (ajuste se tiver o ponto exato do Bloco H)
lat = -22.8596582;      % graus
lon = -43.2303236;      % graus
h   = 10;               % altura em metros (estimativa)

UserPosition = llh2ecef(lat, lon, h);   % [x y z] em metros (ECEF)


% rng(42);
M = 9;
numSV = M;
% % Geometria BOA (PDOP baixo)
% [SatPosition, SatPRN, azel_deg, PDOP] = randomConstellationPDOP(UserPosition, M, 'good', 5, [1.5 3.0], 3000);
% fprintf('PDOP (good): %.2f\n', PDOP);

% Geometria RUIM (PDOP alto)
[SatPosition, SatPRN, azel_deg, PDOP] = randomConstellationPDOP(UserPosition, M, 'bad', 5, [8 20], 3000);
fprintf('PDOP (bad): %.2f\n', PDOP);



% corrSatPosition = genSyntheticGPSPositions(SatPRN, UserPosition, 26560e3, 5); % raio GPS, máscara 5°
% numSV       = length(SatPRN);

% SatPosition=corrSatPosition(1:numSV, :);
% pesos simples por C/N0 (se todos iguais, vira identidade)
w = ones(length(SatPRN),1);  % ou w ~ C/N0 normalizado, ou w ~ sin(el)^2
W = diag(w/max(w));

% dop = computeDOP(SatPosition, UserPosition, W);
% fprintf('HDOP=%.2f  VDOP=%.2f  PDOP=%.2f  GDOP=%.2f  cond(G)=%.1f\n', ...
%          dop.HDOP, dop.VDOP, dop.PDOP, dop.GDOP, dop.condG);

% return


UserPosition=[3.915394273911475e+06 2.939638207807819e+05 5.009529661006817e+06];

%% Simulation parameters
CNosim = 30:5:50;
Nexpe  = 5;

%% Memory
PosErrLS  = zeros(length(CNosim),Nexpe);
PosErrDPE = zeros(length(CNosim),Nexpe);

%% Signal Parameters
c  = 299792458;
f0 = 1575.42e6;
CodePeriod = 1e-3;
fs  = 50e6;
dt  = 1/fs;

gnss_input.CodePeriod = CodePeriod;
gnss_input.CoherentIntegrations= 1;
gnss_input.NonCoherentIntegrations = 1;
gnss_input.m=1; gnss_input.n=1;
gnss_input.type='BPSK';
gnss_input.Tc = 1/(gnss_input.n*1.023e6);
gnss_input.Ts = 1/(2*gnss_input.m*1.023e6);
gnss_input.fs = fs; gnss_input.dt = dt;
gnss_input.fn = 2e6; gnss_input.order = 36;
gnss_input.numSV = numSV;
gnss_input.SatPosition = SatPosition;
gnss_input.UserPosition = UserPosition;
gnss_input.c = c;
gnss_input.SatPRN = SatPRN;
gnss_input.Nexpe = Nexpe;

%% 2-steps parameters (LS 4D agora!)
ls_input.UserPosition = UserPosition;
ls_input.SatPosition  = SatPosition;
ls_input.c = c; ls_input.numSV = numSV; ls_input.fs = fs;
ls_input.num2stepsIterations = 10;

%% DPE parameters (ARS)
dpe_input.UserPosition = UserPosition;
dpe_input.c = c; dpe_input.numSV = numSV;
dpe_input.fs = fs; dpe_input.fn = gnss_input.fn;
dpe_input.Tc = 1/(gnss_input.n*1.023e6);
dpe_input.SatPosition = SatPosition;
dpe_input.Niter=10000; dpe_input.contraction = 2;
dpe_input.dmax = 10000; dpe_input.dmin = 0.01;
dpe_input.dmax_clk=gnss_input.dt/10; dpe_input.dmin_clk=gnss_input.dt/100;
dpe_input.NormalizaFactor = sqrt(gnss_input.NonCoherentIntegrations)*gnss_input.CoherentIntegrations*gnss_input.CodePeriod*gnss_input.fs;
dpe_input.CN0_est_ind = 1;

%% Signal generation
sigen = signalGen(gnss_input);

%% Simulation
for CNo_idx=1:length(CNosim)
    for exp_idx=1:Nexpe
        CNo=CNosim(CNo_idx)*ones(numSV,1);
        x = receivedSignal(sigen,CNo);
        r = correlateSignal(sigen,x);

        % 2-steps (LS 4D)
        pos_est_ls = conv2stepsPVT(r,ls_input);
        PosErrLS(CNo_idx,exp_idx) = norm(pos_est_ls-UserPosition);

        % DPE (ARS)
        pos_est_dpe = DPEarsPVT(r,dpe_input);
        PosErrDPE(CNo_idx,exp_idx) = norm(pos_est_dpe'-UserPosition);
    end
end

RMSE_LS  = sqrt(mean(PosErrLS.^2,2));
RMSE_DPE = sqrt(mean(PosErrDPE.^2,2));

%% ========= CRB (convencional x DPE) =========
M    = numSV;
Tcoh = CodePeriod;

% LOS e geometria
u_mat = zeros(M,3); 
rho = zeros(M,1);
for i = 1:M
    vec      = SatPosition(i,:) - UserPosition;
    rho(i)   = norm(vec);
    u_mat(i,:) = vec ./ rho(i);
end
T_geo = [u_mat, ones(M,1)];          % G do paper (Eq. 11)
Pmat  = (1/c) * u_mat;

% === B2 com normalização de energia e integral discreta ===
s  = sigen.x_local(1,:);             % mesma réplica usada na correlação
dt = 1/fs;
s  = s ./ sqrt(sum(abs(s).^2) * dt); % ||s||^2 = 1 (energia)
ds = diff(s) / dt;                   % derivada temporal
B_2 = sum(abs(ds).^2) * dt;          % ≈ ∫|s'(t)|^2 dt

fCRB_2SP = zeros(length(CNosim),1);
fCRB_DPE = zeros(length(CNosim),1);

for idx = 1:length(CNosim)
    CN0dB   = CNosim(idx);

    % *** MODELO BASEBAND COMPLEXO: Es/N0 = 2 * (C/N0) * Tcoh ***
    SNR_lin = 2 * 10^((CN0dB + 10*log10(Tcoh))/10);
    

    % FIM em atraso por satélite: J_tau = 2 * B2 * SNR  (paper)
    Jtau = diag( 2 * B_2 * SNR_lin * ones(M,1) );

    % --- CRB 2-steps (propagação via LS do paper, Eq. 41) ---
    W    = eye(M);                        % (se quiser, troque por WLS)
    G    = T_geo;
    Ainv = inv(G.'*W*G);                  % (H'H)^{-1}
    C_pdt = c^2 * Ainv * (G.'*W*(Jtau\(W*G))) * Ainv;
    fCRB_2SP(idx) = sqrt(trace(C_pdt(1:3,1:3)));

    % --- CRB DPE (aprox. pré-correlação, p apenas) ---
    Jp   = Pmat.' * Jtau * Pmat;
    Cp   = inv(Jp);
    fCRB_DPE(idx) = sqrt(trace(Cp));
end

T_SUM = table(CNosim(:), RMSE_LS(:), fCRB_2SP(:), ...
                        RMSE_DPE(:), fCRB_DPE(:), ...
    'VariableNames', {'CN0_dBHz','RMSE_LS_m','CRB_2SP_m', ...
                      'RMSE_DPE_m','CRB_DPE_m'});
disp(T_SUM);

%% (Opcional) curvas
figure; semilogy(CNosim, RMSE_LS,'b-.', CNosim, fCRB_2SP,'k--', ...
                 CNosim, RMSE_DPE,'g-',  CNosim, fCRB_DPE,'m-.' , 'LineWidth',2);
grid on; xlabel('C/N_0 [dB-Hz]'); ylabel('RMSE / Bound [m]');
legend('RMSE 2-steps','CRB 2-steps','RMSE DPE','CRB DPE','Location','southwest');
title('RMSE vs CRB (convencional x DPE)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sigen] = signalGen(gnss_input)
CodePeriod = gnss_input.CodePeriod;
CoherentIntegrations = gnss_input.CoherentIntegrations;
NonCoherentIntegrations = gnss_input.NonCoherentIntegrations;
Tc = gnss_input.Tc;
fs = gnss_input.fs;
fn = gnss_input.fn;
order = gnss_input.order;
numSV = gnss_input.numSV;
SatPosition = gnss_input.SatPosition;
UserPosition = gnss_input.UserPosition;
c = gnss_input.c;
SatPRN = gnss_input.SatPRN;
Nexpe = gnss_input.Nexpe;

NsamplesLocal=CodePeriod*fs*CoherentIntegrations;
NsamplesData =CodePeriod*fs*CoherentIntegrations*NonCoherentIntegrations;

Range=zeros(1,numSV);
x_local=zeros(numSV,NsamplesLocal);
fft_local=zeros(numSV,NsamplesLocal);
x_delay=zeros(numSV,NsamplesData);

for kSV=1:numSV 
    Range(kSV) = norm(SatPosition(kSV,:) - UserPosition);    
end
FracDelay=mod(Range/c,CodePeriod);

PrevNCOIndex = -FracDelay/Tc;
randomDelay=0;
for kSV=1:numSV
    Code  = genCAcode(SatPRN(kSV));
    Tchip = CodePeriod / length(Code);
    ii = 1:NsamplesLocal;
    x_local(kSV,:) = Code((1 + mod(round((ii - 1) / fs / Tchip), length(Code))));
    ii = 1:NsamplesData;
    x_delay(kSV,:) = Code((1 + mod(round(PrevNCOIndex(kSV)+randomDelay+(ii - 1)/fs/Tchip), length(Code))));
end

wn=pi*fn/(pi*fs/2);
h=fir1(order,wn);
for kSV=1:numSV
    x_delay(kSV,:)  = filtfilt(h,1,x_delay(kSV,:));
    x_local(kSV,:)  = filtfilt(h,1,x_local(kSV,:));
    fft_local(kSV,:) = fft(x_local(kSV,:),NsamplesLocal);
end

for kSV=1:numSV
    x_delay(kSV,:)  = x_delay(kSV,:)*sqrt((NsamplesData/sum(x_delay(kSV,:).^2)));
end

% dummy noise draws (mantém RNG alinhado se quiser comparar runs)
for exp_idx=1:Nexpe %#ok<NASGU>
    noise = ( sqrt(1/2)*randn(1,NsamplesData) +1i*sqrt(1/2)*randn(1,NsamplesData) ); %#ok<NASGU>
end

sigen.x_delay = x_delay;
sigen.x_local = x_local;
sigen.fft_local = fft_local;
sigen.randomDelay = randomDelay;
sigen.NsamplesLocal = NsamplesLocal;
sigen.NsamplesData = NsamplesData;
sigen.NonCoherentIntegrations = NonCoherentIntegrations;
sigen.fs = fs; sigen.fn = fn; sigen.order = order;
sigen.numSV = numSV; sigen.Nexpe = Nexpe;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = receivedSignal(sigen, CNo)
x_delay = sigen.x_delay;
NsamplesData = sigen.NsamplesData;
numSV = sigen.numSV;
fs = sigen.fs;
fn = sigen.fn;
order = sigen.order;

x_delay_noise=zeros(numSV,NsamplesData);

% ruído complexo unitário (potência 1)
noise = ( sqrt(1/2)*randn(1,NsamplesData) +1i*sqrt(1/2)*randn(1,NsamplesData) );

for kSV=1:numSV
    if CNo(kSV)<100
        % A^2 * fs / 2 = 10^(C/N0/10)  -> A = sqrt(2 * 10^(C/N0/10) / fs)
        % A = sqrt( 2 * 10^(CNo(kSV)/10) / fs );
        A = sqrt(2*10^(CNo(kSV)/10)/fs);

        x_delay_noise(kSV,:) = A * x_delay(kSV,:);
    else
        x_delay_noise(kSV,:)=x_delay(kSV,:);
    end
end

x = sum(x_delay_noise, 1)+noise;

wn=pi*fn/(pi*fs/2);
h=fir1(order,wn);
x  = filtfilt(h,1,x);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r = correlateSignal(sigen, x)
numSV = sigen.numSV;
fft_local = sigen.fft_local;
N  = sigen.NsamplesLocal;
NC = sigen.NonCoherentIntegrations;

r = zeros(numSV, N);
for kSV = 1:numSV
    for idx_nc = 1:NC
        seg = x(1, N*(idx_nc-1)+1 : N*idx_nc);
        Rk  = ifft( fft(seg, N) .* conj(fft_local(kSV,:)) );
        r(kSV,:) = r(kSV,:) + abs(Rk).^2;
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pos_est_ls = conv2stepsPVT(r, ls_input)
    % === Entradas ===
    UserPosition        = ls_input.UserPosition(:);   % 3x1
    c                   = ls_input.c;
    numSV               = ls_input.numSV;
    fs                  = ls_input.fs;
    numIts              = ls_input.num2stepsIterations;
    SatPosition         = ls_input.SatPosition;       % Mx3

    % === 1) Medida de pseudo-distância fracionária (1 ms) em metros ===
    [~, maxPos] = max(r,[],2);        % índice do pico
    maxPos      = maxPos - 1;         % zero-based
    tau_hat     = maxPos / fs;        % s
    rho_hat     = c * tau_hat;        % m (mod 1 ms)

    % === 2) Estado inicial [x;y;z;b] com b (m) ===
    xhat = [UserPosition + 10*(2*rand(3,1)-1) ; 0];   % b inicial = 0 m

    % meia-faixa do wrap em metros (± c*1ms/2)
    halfspan_m = c * 1e-3 / 2;

    % === 3) LS iterativo com normal equations (paper) ===
    for it = 1:numIts
        H        = zeros(numSV, 4);
        rho_pred = zeros(numSV, 1);

        for k = 1:numSV
            vec  = SatPosition(k,:).' - xhat(1:3);   % vetor LOS
            dist = norm(vec);
            if dist < eps, dist = eps; end
            u         = vec / dist;                  % u_i
            H(k,1:3)  = -u.';                        % d rho / d p
            H(k,4)    = 1;                           % d rho / d b
            rho_pred(k) = dist + xhat(4);            % m
        end

        % resíduo em metros, com wrap para ± c*Tc/2
        res = wrap_pm(rho_hat - rho_pred, halfspan_m);

        % === atualização LS exatamente como no paper ===
        % delta = (H'H)^{-1} H' res
        N     = (H.' * H);
        g     = (H.' * res);
        delta = inv(N) * g;

        % atualiza estado
        xhat = xhat + delta;

        % (opcional) critério de parada
        if norm(delta(1:3)) < 1e-3 && abs(delta(4)) < 1e-3
            break;
        end
    end

    % retorna só a posição (linha, para bater com seu código)
    pos_est_ls = xhat(1:3).';
end

function y = wrap_pm(x, halfspan)
    % mapeia x para o intervalo [-halfspan, +halfspan]
    y = mod(x + halfspan, 2*halfspan) - halfspan;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pos_est_dpe = DPEarsPVT(r, dpe_input)
UserPosition = dpe_input.UserPosition;
c = dpe_input.c;
numSV = dpe_input.numSV;
fs = dpe_input.fs;
Tc = dpe_input.Tc;
SatPosition = dpe_input.SatPosition;
Niter = dpe_input.Niter;
contraction = dpe_input.contraction;
dmax = dpe_input.dmax;
dmin = dpe_input.dmin;
dmax_clk = dpe_input.dmax_clk;
dmin_clk = dpe_input.dmin_clk;
NormalizaFactor = dpe_input.NormalizaFactor;
dt=1/fs;
randomDelay = 0;

gamma_est=zeros(3,Niter+1);
amp_est  = zeros(numSV, Niter+1);
EstRxClkBias=zeros(1,Niter+1);

EstRange=zeros(1,numSV); MaxCorr=zeros(1,numSV);

% gamma_est(:,1) = UserPosition+100*(2*rand(3,1)-1)';
gamma_est(:,1) = UserPosition(:) + 100*(2*rand(3,1)-1);

EstRxClkBias(:,1)=-randomDelay*Tc;

for kSV=1:numSV
    EstRange(kSV) = norm(SatPosition(kSV,:) - gamma_est(:,1)');
end
EstFracDelay=mod(EstRange/c+EstRxClkBias(:,1)+dt,1e-3);
for kSV=1:numSV
    aux=round(EstFracDelay(kSV)*fs);
    aux(aux==0)=1e-3*fs;
    MaxCorr(kSV)=r(kSV,aux);
end
J_ant=sum(MaxCorr);
amp_est(:,1) = MaxCorr./NormalizaFactor^2;

d = dmax; d_clk=dmax_clk;

for it = 1:Niter-1
    rand_point = gamma_est(:,it) + d*(2*rand(3,1)-1);
    rand_clk = 0;
    for kSV=1:numSV
        EstRange(kSV) = norm(SatPosition(kSV,:) - rand_point');
    end
    EstFracDelay=mod(EstRange/c+rand_clk+dt,1e-3);
    for kSV=1:numSV
        aux=round(EstFracDelay(kSV)*fs);
        aux(aux==0)=1e-3*fs;
        MaxCorr(kSV)=r(kSV,aux);
    end
    J = sum(MaxCorr);
    if J > J_ant
        gamma_est(:,it+1) = rand_point;
        EstRxClkBias(:,it+1)=rand_clk;
        amp_est(:, it+1) = MaxCorr./NormalizaFactor^2;
        J_ant = J; d = dmax; d_clk=dmax_clk;
    else
        gamma_est(:,it+1) = gamma_est(:,it);
        EstRxClkBias(:,it+1)=EstRxClkBias(:,it);
        amp_est(:,it+1) = amp_est(:,it);
        d = d/contraction; d_clk=d_clk/contraction;
    end
    if d < dmin, d = dmax; end
    if d_clk < dmin_clk, d_clk =dmax_clk; end
end
pos_est_dpe = gamma_est(:,it+1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = wrap_to_pm(x, halfspan)
% mapeia x para o intervalo [-halfspan, +halfspan]
y = mod(x + halfspan, 2*halfspan) - halfspan;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CAcode = genCAcode(PRN)
g2s = [5,6,7,8,17,18,139,140,141,251,252,254,255,256,257,258,469,470,471,472,473,474,509,512,513,514,515,516,859,860,861,862, ...
        145,175,52,21,237,235,886,657,634,762,355,1012,176,603,130,359,595,68,386];
g2shift = g2s(PRN);

g1 = zeros(1,1023); reg = -1*ones(1,10);
for i=1:1023
    g1(i)     = reg(10);
    saveBit   = reg(3)*reg(10);
    reg(2:10) = reg(1:9);
    reg(1)    = saveBit;
end

g2 = zeros(1,1023); reg = -1*ones(1,10);
for i=1:1023
    g2(i)     = reg(10);
    saveBit   = reg(2)*reg(3)*reg(6)*reg(8)*reg(9)*reg(10);
    reg(2:10) = reg(1:9);
    reg(1)    = saveBit;
end

g2 = [g2(1023-g2shift+1 : 1023), g2(1 : 1023-g2shift)];
CAcode = -(g1 .* g2);
end

function SatPosition = genSyntheticGPSPositions(SatPRN, UserPosition, R_orbit, elevMaskDeg)
% Gera posições ECEF (m) para os PRNs informados, coerentes com GPS:
% - Satélites colocados numa casca esférica de raio R_orbit (~26.56e6 m)
% - Azimutes bem espaçados (ângulo áureo) e elevações crescentes acima do mask
% - Interseção da LOS do usuário com a esfera orbital → coordenadas ECEF
%
% Parâmetros:
%   SatPRN        : vetor de PRNs (apenas para tamanho/identidade)
%   UserPosition  : [x y z] do usuário em ECEF (m)
%   R_orbit       : raio orbital (m). Default: 26560e3 (GPS)
%   elevMaskDeg   : máscara de elevação (graus). Default: 5
%
% Regime de operação:
% - Terra esférica (aprox.): lat = atan2(z, sqrt(x^2+y^2)), lon = atan2(y,x)
% - Azimute medido a partir do Norte, sentido horário: 0°=N, 90°=E (convenção ENU)
% - Saída é determinística (não usa aleatoriedade)

    if nargin < 3 || isempty(R_orbit),     R_orbit = 26560e3; end
    if nargin < 4 || isempty(elevMaskDeg), elevMaskDeg = 5;    end

    N = numel(SatPRN);

    % --- 1) Lat/Lon (aprox. esférica) do usuário a partir de ECEF ---
    x = UserPosition(1); y = UserPosition(2); z = UserPosition(3);
    lon = atan2(y, x);
    lat = atan2(z, hypot(x,y));  % geocêntrica (boa o suficiente aqui)

    % Rotação ECEF->ENU no usuário; ENU->ECEF é a transposta
    R_ecef2enu = [ -sin(lon)           cos(lon)            0;
                   -sin(lat)*cos(lon) -sin(lat)*sin(lon)  cos(lat);
                    cos(lat)*cos(lon)  cos(lat)*sin(lon)  sin(lat) ];
    R_enu2ecef = R_ecef2enu.';  % transpose

    % --- 2) Gera (az, el) determinísticos e bem espaçados ---
    % Azimutes via "ângulo áureo" para uniformidade em [0,360)
    az_deg = mod((0:N-1) * 137.508, 360);       % graus
    % Elevações distribuídas entre o mask e ~85° (evita zênite exato)
    el_deg = elevMaskDeg + ((1:N)/(N+1)) * (85 - elevMaskDeg);

    % --- 3) Constrói vetores LOS em ENU e leva para ECEF ---
    SatPosition = zeros(N,3);
    r_u   = UserPosition(:);
    r_u2  = dot(r_u, r_u);

    for k = 1:N
        az = deg2rad(az_deg(k));
        el = deg2rad(el_deg(k));

        % LOS unitário em ENU (conv.: az de Norte→Leste, elevação positiva)
        u_enu  = [cos(el)*sin(az);
                  cos(el)*cos(az);
                  sin(el)];
        u_ecef = R_enu2ecef * u_enu;            % unitário em ECEF

        % --- 4) Interseção r_u + s*u = sat, com ||sat|| = R_orbit ---
        b1   = dot(r_u, u_ecef);                % r_u·u
        disc = b1^2 - (r_u2 - R_orbit^2);
        if disc < 0
            % Numérico muito raro; "empurra" para o horizonte
            disc = 0;
        end
        s = -b1 + sqrt(disc);                   % raiz positiva → sat "acima" do usuário
        SatPosition(k,:) = (r_u + s*u_ecef).';
    end
end

function dop = computeDOP(SatPosition, UserPosition, W)
% Computa HDOP, VDOP, PDOP, TDOP, GDOP e cond(G).
% Regime: ECEF, LOS unitários do usuário para cada SV.

    if nargin < 3 || isempty(W), W = eye(size(SatPosition,1)); end

    M  = size(SatPosition,1);
    u  = zeros(M,3);
    for i=1:M
        v   = SatPosition(i,:) - UserPosition;
        u(i,:)= v / norm(v);
    end
    G = [u, ones(M,1)];            % matriz de geometria (M x 4)

    % DOP com/sem pesos
    Q = inv(G.'*W*G);

    dop.HDOP = sqrt(Q(1,1)+Q(2,2));
    dop.VDOP = sqrt(Q(3,3));
    dop.PDOP = sqrt(Q(1,1)+Q(2,2)+Q(3,3));
    dop.TDOP = sqrt(Q(4,4));
    dop.GDOP = sqrt(trace(Q));
    dop.condG = cond(G);           % número de condição (quanto menor, melhor)
end

function [SatPosition, SatPRN, azel_deg, PDOP] = randomConstellationPDOP( ...
    UserPosition, M, mode, elevMaskDeg, pdopTarget, maxTries)
% Gera M satélites com PDOP aproximado desejado.
% mode: 'good' (baixa DOP) | 'bad' (alta DOP) | 'custom' (intervalo arbitrário)
% pdopTarget: [PDOPmin PDOPmax], ex: good=[1.5 3.5], bad=[8 20]
% elevMaskDeg: máscara de elevação em graus (ex.: 5)
% maxTries: nº máx de tentativas (ex.: 2000)

if nargin < 3 || isempty(mode),        mode = 'good'; end
if nargin < 4 || isempty(elevMaskDeg), elevMaskDeg = 5; end
if nargin < 5 || isempty(pdopTarget)
    if strcmpi(mode,'good'), pdopTarget = [1.5 3.5]; else, pdopTarget = [8 20]; end
end
if nargin < 6 || isempty(maxTries),    maxTries = 2000; end

% --- Conversão ECEF -> base ENU no usuário ---
x = UserPosition(1); y = UserPosition(2); z = UserPosition(3);
lon = atan2(y, x);
hyp = hypot(x, y);
lat = atan2(z, hyp);

e_hat = [-sin(lon),              cos(lon),             0          ];
n_hat = [-sin(lat)*cos(lon), -sin(lat)*sin(lon),  cos(lat)];
u_hat = [ cos(lat)*cos(lon),  cos(lat)*sin(lon),  sin(lat)];
B     = [e_hat(:) n_hat(:) u_hat(:)];  % 3x3 ENU->ECEF

elMask = deg2rad(elevMaskDeg);

bestPDOP = inf;
bestSat  = [];
bestAzEl = [];

for it = 1:maxTries

    switch lower(mode)
        case 'good'
            % Azimutes estratificados para espalhar no céu
            base = (0:M-1)/M*2*pi;
            jitter = (2*pi/M)*(rand(1,M)-0.5);      % pequeno ruído
            az = mod(base + jitter, 2*pi);

            % Elevações variadas (baixa até alta)
            % Amostra uniforme em sin(el) para não enviesar para o horizonte
            sinEl = sin(elMask) + (1 - sin(elMask))*rand(1,M);
            el = asin(sinEl);
            % garante 1-2 satélites bem altos
            nhigh = max(1, round(0.2*M));
            idxh = randperm(M, nhigh);
            el(idxh) = max(el(idxh), deg2rad(70 + 15*rand(1,nhigh)));

        case 'bad'
            % Azimutes aglomerados num setor estreito
            az0 = 2*pi*rand;
            span = deg2rad(50);                      % setor estreito ~50°
            az   = mod(az0 + (rand(1,M)-0.5)*span, 2*pi);

            % Elevações próximas da máscara
            el = elMask + deg2rad(0.5 + 4.5*rand(1,M));  % 0.5°–5° acima da máscara

        otherwise % 'custom'
            % Sorteio amplo; deixa o pdopTarget fazer o filtro
            az  = 2*pi*rand(1,M);
            sinEl = sin(elMask) + (1 - sin(elMask))*rand(1,M);
            el  = asin(sinEl);
    end

    % Vetor LOS em ENU (az a partir do norte, sentido horário)
    ce = cos(el); sa = sin(az); ca = cos(az);
    v_enu = [ (ce .* sa); (ce .* ca); (sin(el)) ];  % 3xM

    % ENU -> ECEF
    v_ecef = B * v_enu;                    % 3xM (unit)

    % Alcances realistas (usuário -> satélite), 20–27 mil km
    R = (2.0e7 + (2.7e7-2.0e7)*rand(1,M));
    SatPosition_try = UserPosition + (v_ecef .* R).';   % Mx3

    % Calcula PDOP
    PDOP_try = localPDOP(UserPosition, SatPosition_try);
    if isnan(PDOP_try), continue; end

    % Teste alvo
    if PDOP_try >= pdopTarget(1) && PDOP_try <= pdopTarget(2)
        SatPosition = SatPosition_try;
        PDOP        = PDOP_try;
        azel_deg    = [rad2deg(az(:)) rad2deg(el(:))];
        SatPRN      = randperm(32, M);
        return
    end

    % Guarda melhor
    if PDOP_try < bestPDOP && strcmpi(mode,'good')
        bestPDOP = PDOP_try; bestSat = SatPosition_try; bestAzEl = [rad2deg(az(:)) rad2deg(el(:))];
    elseif PDOP_try > bestPDOP && strcmpi(mode,'bad')
        % para 'bad', o "melhor" é o de maior PDOP
        bestPDOP = PDOP_try; bestSat = SatPosition_try; bestAzEl = [rad2deg(az(:)) rad2deg(el(:))];
    elseif isempty(bestSat)
        bestPDOP = PDOP_try; bestSat = SatPosition_try; bestAzEl = [rad2deg(az(:)) rad2deg(el(:))];
    end
end

% Se não atingiu a janela alvo, devolve o melhor que conseguiu
SatPosition = bestSat;
PDOP        = bestPDOP;
azel_deg    = bestAzEl;
SatPRN      = randperm(32, M);
end

% ===== Helper: PDOP a partir de posições ECEF =====
function PDOP = localPDOP(UserPosition, SatPosition)
M = size(SatPosition,1);
if M < 4, PDOP = NaN; return; end
r = SatPosition - UserPosition;              % Mx3
d = vecnorm(r,2,2);
if any(d == 0), PDOP = NaN; return; end
u = r ./ d;                                  % Mx3 unit LOS
H = [-u, ones(M,1)];                         % Mx4
A = H.'*H;
if rcond(A) < 1e-12, PDOP = NaN; return; end
C = inv(A);
PDOP = sqrt(trace(C(1:3,1:3)));
end

function ecef = llh2ecef(lat_deg, lon_deg, h_m)
    % WGS-84
    a  = 6378137.0;           % semi-eixo maior [m]
    f  = 1/298.257223563;
    e2 = f*(2-f);

    lat = deg2rad(lat_deg);
    lon = deg2rad(lon_deg);

    sinlat = sin(lat); coslat = cos(lat);
    sinlon = sin(lon); coslon = cos(lon);

    N = a / sqrt(1 - e2*sinlat.^2);

    x = (N + h_m)      * coslat * coslon;
    y = (N + h_m)      * coslat * sinlon;
    z = (N*(1 - e2) + h_m) * sinlat;

    ecef = [x, y, z];
end
