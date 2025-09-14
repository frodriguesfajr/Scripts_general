% ===========================================

close all;
clear;
format long;

rng(42)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c   = 299792458;
f0  = 1575.42e6;
%% ---- Signal Parameters (GPS E1)
numSV        = 7;
SatPRN       = [12 15 17 19 24 25 32];
UserPosition = [3.915394273911475e+06 2.939638207807819e+05 ...
    5.009529661006817e+06];  
CodePeriod               = 1e-3;
CoherentIntegrations     = 1;
NonCoherentIntegrations  = 1;
m                        = 1;
n                        = 1;
type                     = 'BPSK';
Tc                       = 1/(n*1.023e6);
Ts                       = 1/(2*m*1.023e6);
fs                       = 50e6;
dt                       = 1/fs;
fn                       = 2e6;
order                    = 36;
Niter          = 10000;
gamma_est      = zeros(3, Niter+1);
amp_est        = zeros(numSV, Niter+1);
EstRxClkBias   = zeros(1, Niter+1);
contraction    = 2;
dmax           = 10000;
dmin           = 0.01;
dmax_clk       = dt/10;
dmin_clk       = dt/100;
NormalizaFactor = sqrt(NonCoherentIntegrations)*CoherentIntegrations*CodePeriod*fs;
CN0_est_ind     = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ---- DPE Parameters
dpe_param.gamma_guess = UserPosition;
dpe_param.Niter = 10000;
dpe_param.dmax = 10000;
dpe_param.dt = 0.01;
dpe_param.dmax_clk = dt/10;
dpe_param.dmin_clk = dt/100;
dpe_param.contraction = 2;
dpe_param.dmin = 0.01;
dpe_param.fn = fn;
dpe_param.CN0_est_ind = CN0_est_ind;
gamma_est      = zeros(3, dpe_param.Niter+1);
amp_est        = zeros(numSV,  dpe_param.Niter+1);
EstRxClkBias   = zeros(1,  dpe_param.Niter+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
SatPosition  = corrSatPosition(1:7, :);
% (Opcional) salve para carregar depois
% save('SatPositions.mat','corrSatPosition');
                                              
%% ---- Simulation parameters"    
CNosim = 30:5:50;
Nexpe  = 2;
RIMuse             = 0;
estimateTrueNoise  = 1;
% ---- 2-steps parameters"
num2stepsIterations = 10;
%% Alocação de memória
PosErrLS   = zeros(length(CNosim), Nexpe);
CN0_est    = zeros(length(CNosim), Nexpe);
cn0        = zeros(length(CNosim), numSV, Nexpe);
%% Geração de sinal (sigen struct)
% number of samples calculation
% Number of samples of the Local Replica
NsamplesLocal=CodePeriod*fs*CoherentIntegrations;   
%Number of samples of the Received Signal (Data).
NsamplesData=CodePeriod*fs*CoherentIntegrations*NonCoherentIntegrations;    
%% memory allocation
Range=zeros(1,numSV);
x_local=zeros(numSV,NsamplesLocal);
fft_local=zeros(numSV,NsamplesLocal);
x_delay=zeros(numSV,NsamplesData);
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

%% Filter local signal and generate its FFT 

wn=pi*fn/(pi*fs/2);
h=fir1(order,wn);
for kSV=1:numSV
    x_delay(kSV,:)  = filtfilt(h,1,x_delay(kSV,:));
    x_local(kSV,:)  = filtfilt(h,1,x_local(kSV,:));
    fft_local(kSV,:) = fft(x_local(kSV,:),NsamplesLocal);
end

%% Normalize Received Signal Power after filtering
for kSV=1:numSV
    x_delay(kSV,:)  = x_delay(kSV,:)*sqrt((NsamplesData/sum(x_delay(kSV,:).^2)));
end

%% gather outputs in a struct
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
sigen.estimateTrueNoise = estimateTrueNoise;
sigen.CodePeriod = CodePeriod;
sigen.CoherentIntegrations = CoherentIntegrations;
sigen.UserPosition = UserPosition;
sigen.fs = fs;
sigen.c = c;
sigen.num2stepsIterations = num2stepsIterations;
sigen.SatPosition = SatPosition;
sigen.Tc = Tc;
sigen.dt = dt;
sigen.fs = fs;
sigen.NormalizaFactor = NormalizaFactor;
sigen.dmax = dmax;
sigen.dmax_clk = dmax_clk;
sigen.Niter = Niter;
sigen.dmin = dmin;
sigen.dmin_clk = dmin_clk;
sigen.contraction = contraction;
sigen.fn = fn;
sigen.CN0_est_ind = CN0_est_ind;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if estimateTrueNoise == 0
    meanNoise = nan;
else
    %% Add AWGN noise to the transmitted signals
    mean_noise=zeros(1,Nexpe);
    for exp_idx=1:Nexpe
        noise = ( sqrt(1/2)*randn(1,NsamplesData) + ...
            1i* sqrt(1/2)*randn(1,NsamplesData) );
        r = correlateSignal(sigen,noise);

        mean_noise(exp_idx) = mean(mean(r,2));
    end
    meanNoise=mean(mean_noise);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start Simulation
for CNo_idx = 1:length(CNosim)
    % disp([CNosim(CNo_idx)])
    for exp_idx = 1:Nexpe
        CNo = CNosim(CNo_idx) * ones(numSV,1);
        %% Signal + noise
        x = receivedSignal(sigen, CNo);
        %% Perform coherent/non-coherent integration times
        r = correlateSignal(sigen,x);
        %% Estimate CN0
        %% DPE approach ARS (accelerated random search)
        PosErrDPE(CNo_idx,exp_idx) = DPEarsPVT(r,sigen, dpe_param);
    end
end
% Compute RMSEs
RMSE_DPE   = sqrt(mean(PosErrDPE.^2, 2));


%% Cramer Rao Bound computation
x_local = sigen.x_local;

B_2=sum((diff(x_local(1,:))/dt).^2)/sum(x_local(1,:).^2);
T=CodePeriod;
D=T*c;
M=numSV;
P= 1/c*(UserPosition-SatPosition)./sqrt(sum((UserPosition-SatPosition).^2,2));

Rd=[D^2/12 0 0; 0 D^2/12 0; 0 0 D^2/12 ];
RT=[T^2/12 0 0; 0 T^2/12 0; 0 0 T^2/12 ];
varT=T^2/12;

SNRdb=CNosim+10*log10(T);
fZZLB_DPE = zeros(1,length(CNosim));
fZZB_DPE  = zeros(1,length(CNosim));

%% === CRB DPE ===============================================
SNRdb  = CNosim + 10*log10(T);
fCRB_DPE  = zeros(1,length(CNosim));


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

function PosErrDPE = DPEarsPVT(r, sigen, dpe_param)

numSV = sigen.numSV;
% UserPosition = sigen.UserPosition;
Tc = sigen.Tc;
SatPosition = sigen.SatPosition;
c = sigen.c;
dt = sigen.dt;
fs = sigen.fs;
NormalizaFactor = sigen.NormalizaFactor;
dmax = sigen.dmax;
dmax_clk = sigen.dmax_clk;
Niter = sigen.Niter;
dmin = sigen.dmin;
dmin_clk = sigen.dmin_clk;
contraction = sigen.contraction;
fn = sigen.fn;
CN0_est_ind = sigen.CN0_est_ind;

randomDelay = 0; %%% magic number....

%% memory allocation
EstRange=zeros(1,numSV);
MaxCorr=zeros(1,numSV);

UserPosition = dpe_param.gamma_guess;
Niter = dpe_param.Niter;
dmax = dpe_param.dmax;
% dt = dpe_param.dt;
dmax_clk = dpe_param.dmax_clk;
dmin_clk = dpe_param.dmin_clk;
contraction = dpe_param.contraction;
% dmin = dpe_param.dmin;
% fn = dpe_param.fn;
% CN0_est_ind = dpe_param.CN0_est_ind;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gamma_est(:,1) = UserPosition+100*(2*rand(3,1)-1)';
EstRxClkBias(:,1)=-randomDelay*Tc;


for kSV                                         =   1 : numSV
    % Compute range estimate with corrected satellite position, atmosphere
    % corrections and satellite clock error.
    EstRange(kSV)                               =   norm(SatPosition(kSV,:) - gamma_est(:,1)');
end;

EstFracDelay=mod(EstRange/c+EstRxClkBias(:,1)+dt,1e-3);
% EstFracDelay=EstFracDelay+EstRxClkBias(:,1);
% EstFracDelay(EstFracDelay<0)=EstFracDelay(EstFracDelay<0)+1e-3;
% EstFracDelay(EstFracDelay>1e-3)=EstFracDelay(EstFracDelay>1e-3)-1e-3;

for kSV                                         =   1 : numSV
    
    %     round(FracDelay*Fs)
    aux=round(EstFracDelay(kSV)*fs);
    aux(aux==0)=1e-3*fs;
    MaxCorr(kSV)=r(kSV,aux);
    
end
J_ant=sum(MaxCorr);
amp_est(:,1) = MaxCorr./NormalizaFactor^2;

d = dmax;
d_clk=dmax_clk;

for it = 1:Niter-1        %%% ARS algorithm iterations
    
    % draw a random movement
    rand_point = gamma_est(:,it) + d*(2*rand(3,1)-1);
    % rand_clk = EstRxClkBias(:,it)+d_clk*(2*rand-1);
    rand_clk = 0;
    for kSV                                         =   1 : numSV
        
        EstRange(kSV)                               =   norm(SatPosition(kSV,:) - rand_point');
    end
    
    EstFracDelay=mod(EstRange/c+rand_clk+dt,1e-3);
    
    for kSV                                         =   1 : numSV
        
        aux=round(EstFracDelay(kSV)*fs);
        aux(aux==0)=1e-3*fs;
        MaxCorr(kSV)=r(kSV,aux);
        
    end
    
    
    J = sum(MaxCorr);
    
    % select or discard point
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
    if d < dmin
        d = dmax;
    end
    
    if d_clk < dmin_clk
        d_clk =dmax_clk;
    end
end
% DPE position estimation
PosErrDPE=norm(gamma_est(:,it+1)'-UserPosition);
CN0_est = 10.*log10(amp_est(:, it+1).*(2*fn)); 
CN0_est = CN0_est (CN0_est_ind);

end