close all;
clear; 
clc; 
format long;
% rng(42)

%% Scenario Parameters
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
numSV=7;
SatPosition=corrSatPosition(1:7, :);
SatPRN=  [12 15 17 19 24 25 32];
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
    % Estima [x,y,z,b] por LS (iterativo) usando atrasos de 1 ms (mod 1e-3)
    % Equações: (8)–(14) do paper, com modelagem em tempo e Jacobiana H=[-u 1]

    % ---- Dados / constantes
    p0  = ls_input.UserPosition(:);      % 3x1
    c   = ls_input.c;
    fs  = ls_input.fs;
    M   = ls_input.numSV;
    S   = ls_input.SatPosition;          % Mx3
    Nit = ls_input.num2stepsIterations;

    % ---- Medição de atraso (em segundos) a partir do pico de correlação
    [~, kmax] = max(r, [], 2);           % índice (1..N)
    tau_hat   = (kmax - 1) / fs;         % s

    % ---- Inicialização do estado [x;y;z;b]
    xhat = [ p0 + 10*(2*rand(3,1)-1) ; 0 ];  % 4x1, b em segundos

    for it = 1:Nit
        % Geometria atual
        rho = zeros(M,1);
        u   = zeros(M,3);
        for i = 1:M
            v       = S(i,:).' - xhat(1:3);    % 3x1
            rho(i)  = norm(v);
            u(i,:)  = (v./rho(i)).';           % 1x3
        end

        % Modelo em tempo: tau_model = rho/c + b
        tau_model = rho/c + xhat(4);           % s

        % Resíduo em tempo com correção de 1 ms (equivale ao wrap do artigo)
        d_tau = tau_hat - tau_model;           % s
        d_tau = d_tau - 1e-3*round(d_tau/1e-3);% mapeia para [-0.5,0.5] ms

        % Jacobiana H = [-u   1]
        H = [-u, ones(M,1)];                   % Mx4

        % Passo de Gauss-Newton na forma (12)/(14)
        dx = (H.'*H) \ (H.' * (c*d_tau));      % dx em metros para xyz e em metros para b*c

        % Atualização (converter b de metros para segundos)
        xhat(1:3) = xhat(1:3) + dx(1:3);
        xhat(4)   = xhat(4)   + dx(4)/c;       % b [s]

        % (opcional) pare se ||dx|| for pequeno
        % if norm(dx) < 1e-3, break; end
    end

    pos_est_ls = xhat(1:3).';  % retorna 1x3 como no resto do seu código
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
