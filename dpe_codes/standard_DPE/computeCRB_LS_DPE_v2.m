% ===========================================

close all;
clear;
format long;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c   = 299792458;
%% ---- Signal Parameters (GPS E1)

SatPRN       = [12 15 17 19 24 25 32];
UserPosition = [3.915394273911475e+06 2.939638207807819e+05 ...
    5.009529661006817e+06];  
CodePeriod               = 1e-3;
CoherentIntegrations     = 1;
NonCoherentIntegrations  = 1;
n                        = 1;
type                     = 'BPSK';
Tc                       = 1/(n*1.023e6);
fs                       = 50e6;
%% ---- Scenario Parameters
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

numSV        = 3;

SatPosition  = corrSatPosition(1:numSV, :);
                                              
CNosim = 30:5:50;
NsamplesLocal=CodePeriod*fs*CoherentIntegrations;   
NsamplesData=CodePeriod*fs*CoherentIntegrations*NonCoherentIntegrations;    

%% Compute range and fractional delays for each SV
for kSV=1:numSV 
    Range(kSV)                               =   norm(SatPosition(kSV,:) - UserPosition);    
end

FracDelay=mod(Range/c,CodePeriod);
%% Generate local replica and delayed signals according to the computed delays
PrevNCOIndex    =  -  FracDelay/Tc;
randomDelay= 0;
for kSV=1:numSV
    Code                                    =   genCAcode(SatPRN(kSV));
    Tchip                                   =   CodePeriod / length(Code);
    ii                                      =   1 : NsamplesLocal;          % generating LGenBlocks samples
    x_local(kSV,:)                                 =   Code((1 + mod(round((ii - 1) / fs / Tchip), length(Code))));
    ii                                      =   1 : NsamplesData;
    x_delay(kSV,:)                                 =   Code((1 + mod(round(PrevNCOIndex(kSV)+randomDelay+(ii - 1) / fs / Tchip), length(Code))));
end


%% gather outputs in a struct
sigen.x_local = x_local;
sigen.NsamplesData = NsamplesData;
NsamplesData = sigen.NsamplesData;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cramer Rao Bound computation
x_local = sigen.x_local;

dt = 1/fs;
B_2=sum((diff(x_local(1,:))/dt).^2)/sum(x_local(1,:).^2);
T=CodePeriod;
D=T*c;
M=numSV;
% P= 1/c*(UserPosition-SatPosition)./sqrt(sum((UserPosition-SatPosition).^2,2));
normPos = vecnorm(UserPosition-SatPosition,2,2); % 3d norm
P= 1/c*(UserPosition-SatPosition)./normPos;
% === Velocidades sintéticas (se você não tiver efemérides) ===
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fCRB_LS   = zeros(1,length(CNosim));
fCRB_DPE  = zeros(1,length(CNosim));

f0 = 1575.42e6;              % L1
varT = T^2/12;               % momento de 2ª ordem do tempo (0..T) / T (já usado no seu ZZB)

% Jacobiano do atraso wrt posição (somente posição; se quiser clock, acrescente coluna)
% d tau_i / d p = -u_i^T / c, em que u_i é o LOS unitário sat->usuário
u = (SatPosition - UserPosition) ./ normPos;        % Mx3 (LOS sat->user)
H_tau_p = -(1/c) * u;                               % Mx3

% Se você quiser incluir clock (coluna extra), descomente as 2 linhas abaixo
% H_tau = [H_tau_p, ones(M,1)];                      % Mx4  (posição + clock)
% mapIdx = 1:3;                                      % usaremos só a parte de posição no RMSE

% Como o seu DPE atual estima apenas posição (clock foi fixado em 0 no ARS),
% tomamos H_tau = H_tau_p para comparar com o CRB de posição 3D.
H_tau = H_tau_p;
mapIdx = 1:3;
for k = 1:length(CNosim)
    SNR = 10^( (CNosim(k) + 10*log10(T))/10 );

    % ===== LS (2SP): J_tau diagonal por canal; propagar para posição =====
    J_tau = (2*SNR*B_2) * eye(M);              % MxM
    % Propagação: J_p = H_tau^T * J_tau * H_tau
    Jp_LS = H_tau.' * J_tau * H_tau;           % 3x3  (posição)
    CRB_LS = inv(Jp_LS);
    fCRB_LS(k) = sqrt(trace(CRB_LS));          % RMSE 3D [m]

    % ===== DPE: FIM direto em gamma, eliminando amplitudes por Schur =====
    % No modelo atual (apenas atraso), o bloco efetivo em gamma é o mesmo:
    % J_gamma = H_tau^T * J_tau * H_tau
    % (porque d d/d gamma entra só via tau(gamma) e amplitudes somem pelo complemento de Schur)
    Jp_DPE = H_tau.' * J_tau * H_tau;          % 3x3
   
    % --------------------------------------------------------------------

    CRB_DPE = inv(Jp_DPE);
    fCRB_DPE(k) = sqrt(trace(CRB_DPE));        % RMSE 3D [m]
end

%% ====== Mostrar e plotar comparação LS x DPE ============================
T_CRB = table(CNosim(:), fCRB_LS(:), fCRB_DPE(:), ...
    'VariableNames', {'CN0_dBHz','CRB_LS_m','CRB_DPE_m'});
disp(T_CRB);

figure; 
h = semilogy(CNosim, fCRB_LS, 'b-o', CNosim, fCRB_DPE, 'r-s');
grid on; set(h,'LineWidth',1.8);
legend('CRB - LS (2SP)','CRB - DPE','Location','southwest');
xlabel('C/N_0 [dB-Hz]'); ylabel('RMSE 3D [m]');
title('CRB: Two-Step (LS) vs. DPE (hipóteses atuais)');





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
