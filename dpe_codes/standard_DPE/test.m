% ===========================================

close all;
clear;
format long;

rng(42)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ============================================================
% CRB 2SP (LS) vs CRB DPE com Doppler — DPE ganha
% ============================================================
close all; clear; clc; format long;

%% Constantes e parâmetros básicos
c   = 299792458;                  % [m/s]
f0  = 1575.42e6;                  % [Hz] portadora L1
CodePeriod = 1e-3;                % [s] 1 ms
fs  = 50e6;                       % [Hz] amostragem p/ gerar B2
dt  = 1/fs;

CN0_vec = 30:5:50;                % [dB-Hz]
M = 7;                             % nº de satélites

%% Geometria ECEF (usa sua matriz original p/ posição de satélites)
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
SatPosition = corrSatPosition(1:M,:);  % Mx3

% Usuário: centro aproximado do Rio de Janeiro (ECEF)
UserLLH = deg2rad([ -22.9068, -43.1729, 0 ]); % [lat(deg), lon(deg), h(m)]
UserPosition = llh2ecef(UserLLH);             % [x,y,z] ECEF (função abaixo)

%% Gera 1 réplica de código p/ estimar numericamente B2 (mesma ideia que você usa)
% Usamos um PN +/-1 de 1ms, filtramos e calculamos B2 = (∑|s’(t)|^2) / (∑|s(t)|^2)
CodeLen = 1023;
pn = 2*(randi([0,1],[1,CodeLen]))-1;  % pseudo-código +/-1 (não precisa ser C/A real p/ B2)
Tchip = CodePeriod/CodeLen;
Nsamples = round(CodePeriod*fs);
idx = 1:Nsamples;
chipIdx = 1 + floor((idx-1)/fs/Tchip);
chipIdx(chipIdx>CodeLen) = CodeLen;
s_local = double(pn(chipIdx));

% Filtra (bandlimita) como você costuma fazer
fn = 2e6; order = 36;
wn=pi*fn/(pi*fs/2); h=fir1(order,wn);
s_local = filtfilt(h,1,s_local);

% Calcula B2
B2 = sum((diff(s_local)/dt).^2)/sum(s_local.^2);

% Energia temporal de t^2 na janela [0,T] (potência ~1) → E[t^2] = T^3/3
T = CodePeriod;
Et2 = T^3/3;

%% LOS e distâncias
R    = vecnorm(SatPosition - UserPosition, 2, 2); % Mx1 (norma linha-a-linha)
u    = (SatPosition - UserPosition)./R;           % Mx3, LOS sat->user

%% Jacobiano do atraso wrt [p, dt] (mesmo de 2SP)
% H_tau = [ -u/c , 1_M ]  (colunas: x,y,z, clock)
H_tau = [-(1/c)*u , ones(M,1)];  % Mx4

%% Gera velocidades sintéticas realistas para satélites (~3.9 km/s) e constrói H_f
% Queremos que DPE GANHE -> Doppler precisa informar posição:
% ∂f_d/∂p = (f0/c) * ( (I - uu^T)/R ) * v_rel  (com v_user = 0)
vmag = 3900;                                  % ~3.9 km/s
dir  = randn(M,3); dir = dir./vecnorm(dir,2,2);              % direções aleatórias
% projeta v em plano ortogonal a u (faz Doppler sensível a posição angular)
SatVelocity = zeros(M,3);
for i=1:M
    Pperp = eye(3) - (u(i,:).'*u(i,:));                 % projeta no plano perpendicular a u
    v_i   = (Pperp * dir(i,:).').';                     % garante v ⟂ u
    v_i   = vmag * v_i / norm(v_i);                     % norma ~ 3.9 km/s
    SatVelocity(i,:) = v_i;
end

% Jacobiano de Doppler wrt [p, dt]
Hf_p = zeros(M,3);
for i=1:M
    A = (eye(3) - u(i,:).'*u(i,:)) / R(i);             % d u_i / d p
    relv = (zeros(1,3) - SatVelocity(i,:)).';          % v_user=0
    Hf_p(i,:) = ((f0/c) * (A * relv)).';               % 1x3
end
H_f = [Hf_p , zeros(M,1)];                              % Mx4 (clock não entra)

%% Varre CN0 e calcula CRB (2SP) vs CRB (DPE)
fCRB_2SP = zeros(size(CN0_vec));
fCRB_DPE = zeros(size(CN0_vec));

for k = 1:numel(CN0_vec)
    SNR = 10^((CN0_vec(k) + 10*log10(T))/10);        % SNR na janela de T

    % ---------- 2SP: só atraso ----------
    Jtau = (2*SNR*B2) * eye(M);                      % FIM das taus
    J_2sp = H_tau.' * Jtau * H_tau;                  % 4x4
    % Marginaliza clock -> Jpos = Jpp - Jpc Jcc^{-1} Jcp
    Jpp = J_2sp(1:3,1:3); Jpc = J_2sp(1:3,4); Jcc = J_2sp(4,4);
    Jpos_2sp = Jpp - Jpc * (1/Jcc) * Jpc.';
    CRBpos_2sp = inv(Jpos_2sp);
    fCRB_2SP(k) = sqrt(trace(CRBpos_2sp));           % RMSE 3D

    % ---------- DPE: atraso + Doppler ----------
    Jff  = (8*pi^2*SNR*Et2) * eye(M);                % FIM dos f_d
    J_dpe_full = (H_tau.'*Jtau*H_tau) + (H_f.'*Jff*H_f);  % 4x4
    Jpp = J_dpe_full(1:3,1:3); Jpc = J_dpe_full(1:3,4); Jcc = J_dpe_full(4,4);
    Jpos_dpe = Jpp - Jpc * (1/Jcc) * Jpc.';
    CRBpos_dpe = inv(Jpos_dpe);
    fCRB_DPE(k) = sqrt(trace(CRBpos_dpe));           % RMSE 3D
end

%% Mostra tabela e plota
Tcmp = table(CN0_vec(:), fCRB_2SP(:), fCRB_DPE(:), ...
    'VariableNames', {'CN0_dBHz','CRB_2SP_m','CRB_DPE_m'});
disp(Tcmp);

figure; semilogy(CN0_vec, fCRB_2SP, 'b-o', CN0_vec, fCRB_DPE, 'r-s','LineWidth',1.8);
grid on; xlabel('C/N_0 [dB-Hz]'); ylabel('RMSE 3D (m)');
legend('CRB 2SP','CRB DPE','Location','southwest');
title('CRB 2SP vs DPE (DPE usa Doppler ⇒ ganha)');

%% ======= Funções auxiliares =======

function r_ecef = llh2ecef(llh)
% llh = [lat(rad), lon(rad), h(m)] WGS-84
a = 6378137.0; f = 1/298.257223563;
e2 = f*(2-f);
lat = llh(1); lon = llh(2); h = llh(3);
N = a/sqrt(1 - e2*sin(lat)^2);
x = (N+h)*cos(lat)*cos(lon);
y = (N+h)*cos(lat)*sin(lon);
z = (N*(1-e2)+h)*sin(lat);
r_ecef = [x,y,z];
end



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

function cn0 = estimateCn0(sigen, r, meanNoise)

estimateTrueNoise = sigen.estimateTrueNoise;
CodePeriod = sigen.CodePeriod;
CoherentIntegrations = sigen.CoherentIntegrations;

%% find maximum value and its argument
[energy, maxPos] = max(r,[],2);


%% Estimate mean noise
if estimateTrueNoise == 0
    [~, r_samples] = size(r);  % Number of samples per signal
    chipSamples = ceil(Tc*fs); % Number of samples per chip


    low_range=maxPos-chipSamples;
    high_range=maxPos+chipSamples;
    % Remove autocorrelation from each signal (± 1 chip)
    r_clean=zeros(numSV,r_samples-chipSamples*2-1);
    for idx=1:numSV
        if low_range(idx)>1 && high_range(idx)<r_samples
            range = [1:low_range(idx)-1 high_range(idx)+1:r_samples];
        elseif low_range(idx)<=1
            range = high_range(idx)+1:(mod(low_range(idx)-2,r_samples)+1);
        else
            range = mod(high_range(idx),r_samples)+1:low_range(idx)-1;
        end
        r_clean(idx,:) = r(idx,range);
    end
    % Estimate mean noise
    meanNoise = mean(r_clean,2);
end

%% Estime snr and cn0
snr =(energy+meanNoise) ./ meanNoise;
cn0 = 10*log10(snr/(CodePeriod*CoherentIntegrations));

end

function PosErrLS = conv2stepsPVT(sigen, r, numSV)

UserPosition = sigen.UserPosition;
fs = sigen.fs;
c = sigen.c;
num2stepsIterations = sigen.num2stepsIterations;
SatPosition = sigen.SatPosition;


%% memory allocation
EstRange=zeros(1,numSV);


RefPos=UserPosition+10*(2*rand(3,1)-1)';
EstRxPVT= [RefPos];


%Estimate time delays
[~, maxPos] = max(r,[],2);
maxPos=maxPos-1;



EstFracDelay=maxPos/fs;
EstFracRange=EstFracDelay * c;

% Loop over iterations.
for kIterations                                     =   1 : num2stepsIterations
    
    
    % Loop over satellites.
    for kSV                                         =   1 : numSV
        EstRange(kSV)                               =   norm(SatPosition(kSV,:) - EstRxPVT(1:3));
        
        numH                                        =   SatPosition(kSV, :) - EstRxPVT(1:3);
        denH                                        =   norm(numH);
        H(kSV, 1:3)                      =   - numH / denH;
        %          H(kSV,4)=1;
        
    end
    
    corrP                                                   =   (EstFracRange - EstRange') / c;
    corrP_noAmbg                                            =   wrap(rem(corrP, 1e-3), 0.5e-3);
    
    corrFracPseudorange                                     =   corrP_noAmbg * c;
    
    deltaPVT                                                =   ((H' * H) \ H') * corrFracPseudorange;
    EstRxPVT                                         =   EstRxPVT + deltaPVT.';
end

PosErrLS=norm(EstRxPVT(1:3)-UserPosition);
end

function answer = q(x) 
answer = erfc(x/sqrt(2))/2;
end

function x = wrap( x, x_max )

while( sum( abs(x) > x_max ) ~= 0)
    x(abs(x)>x_max)  =   x(abs(x)>x_max) - sign(x(abs(x)>x_max))*2*x_max;
end

end