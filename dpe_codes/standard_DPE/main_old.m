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
% Geometria BOA (PDOP baixo)
[SatPosition, SatPRN, azel_deg, PDOP] = randomConstellationPDOP(UserPosition, M, 'good', 5, [1.5 3.0], 3000);
fprintf('PDOP (good): %.2f\n', PDOP);

% % Geometria RUIM (PDOP alto)
% [SatPosition, SatPRN, azel_deg, PDOP] = randomConstellationPDOP(UserPosition, M, 'bad', 5, [8 20], 3000);
% fprintf('PDOP (bad): %.2f\n', PDOP);



% corrSatPosition = genSyntheticGPSPositions(SatPRN, UserPosition, 26560e3, 5); % raio GPS, máscara 5°
% numSV       = length(SatPRN);

% SatPosition=corrSatPosition(1:numSV, :);
% pesos simples por C/N0 (se todos iguais, vira identidade)
w = ones(length(SatPRN),1);  % ou w ~ C/N0 normalizado, ou w ~ sin(el)^2
W = diag(w/max(w));

% dop = computeDOP(SatPosition, UserPosition, W);
% fprintf('HDOP=%.2f  VDOP=%.2f  PDOP=%.2f  GDOP=%.2f  cond(G)=%.1f\n', ...
%          dop.HDOP, dop.VDOP, dop.PDOP, dop.GDOP, dop.condG);


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

%% 2-steps parameters
ls_input.UserPosition = UserPosition;
ls_input.SatPosition  = SatPosition;
ls_input.c = c; ls_input.numSV = numSV; ls_input.fs = fs;
ls_input.num2stepsIterations = 10;

%% DPE parameters
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
        pos_est_dpe = dpePvt(r,dpe_input);
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

%% curvas
figure; semilogy(CNosim, RMSE_LS,'b-.', CNosim, fCRB_2SP,'k--', ...
                 CNosim, RMSE_DPE,'g-',  CNosim, fCRB_DPE,'m-.' , 'LineWidth',2);
grid on; xlabel('C/N_0 [dB-Hz]'); ylabel('RMSE / Bound [m]');
legend('RMSE 2-steps','CRB 2-steps','RMSE DPE','CRB DPE','Location','southwest');
title('RMSE vs CRB (convencional x DPE)');