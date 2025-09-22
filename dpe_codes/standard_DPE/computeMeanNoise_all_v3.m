% GNSS-DPE tool main script %
% Adrià Gusi, Pau Closas
% v2.0
% This script generates GPS signals and computes the receiver position with
% two-steps and DPE based algorithms.
% Current implementation:
% -7 GPS satellites at a fixed location
% -Code based only
% -Non-coherent and coherent integration methods
% -LS method implemented with Integer-Millisecon Rollover correction
% -DPE algorithm with ARS implementation

% close all;
clear; 
clc; 
format long;

rng(42)

%% Signal generation (sigen struct)
%% Scenario Parameters
% Set GPS SVs and user positions
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
SatPRN=  [12    15    17    19    24    25    32 ];
UserPosition=[3.915394273911475e+06 2.939638207807819e+05 5.009529661006817e+06];
%% Simulation parameters
% Set C/N0 (dB Hz) of each SV signal
CNosim = 30:5:50;
% Number of experiments for each simulated C/N0
Nexpe = 5;
%% Memory allocation.
PosErrLS=zeros(length(CNosim),Nexpe);
PosErrDPE=zeros(length(CNosim),Nexpe);
CN0_est=zeros(length(CNosim),Nexpe);
cn0= zeros(length(CNosim),numSV,Nexpe);
%% Signal Parameters
% GPS E1
c = 299792458; %   Speed of light.
f0 =   1575.42e6; %   Carrier frequency (L1-E1 band).
CodePeriod = 1e-3;
fs = 50e6; % Sampling frequency
dt = 1/fs; % Sampling period
gnss_input.CodePeriod = CodePeriod;
gnss_input.CoherentIntegrations= 1; % Number of Coherent Integrations
gnss_input.NonCoherentIntegrations = 1; % Number o Non Coherent Integrations
gnss_input.m=1; % m
gnss_input.n=1; % n
gnss_input.type='BPSK'; % BPSK, BOCcos, BOCsin
gnss_input.Tc = 1/(gnss_input.n*1.023e6); % Chip time
gnss_input.Ts = 1/(2*gnss_input.m*1.023e6);
gnss_input.fs = fs;
gnss_input.dt = dt; % Sampling period
gnss_input.fn = 2e6;
gnss_input.order = 36;
gnss_input.numSV = numSV;
gnss_input.SatPosition = SatPosition;
gnss_input.UserPosition = UserPosition;
gnss_input.c = c;
gnss_input.SatPRN = SatPRN;
%% 2-steps parameters
% LS iterations
ls_input.UserPosition = UserPosition;
ls_input.SatPosition = SatPosition;
ls_input.c = c;
ls_input.numSV = numSV;
ls_input.fs = gnss_input.fs;
ls_input.num2stepsIterations = 10;
%% DPE parameters (ARS)
dpe_input.UserPosition = UserPosition;
dpe_input.c = c;
dpe_input.numSV = numSV;
dpe_input.fs = gnss_input.fs;
dpe_input.fn = gnss_input.fn;
dpe_input.Tc = 1/(gnss_input.n*1.023e6); % Chip time
dpe_input.SatPosition = SatPosition;
dpe_input.Niter=10000;
dpe_input.contraction = 2;          % contraction parameter
dpe_input.dmax = 10000;
dpe_input.dmin = 0.01;
dpe_input.dmax_clk=gnss_input.dt/10;
dpe_input.dmin_clk=gnss_input.dt/100;
dpe_input.NormalizaFactor = sqrt(gnss_input.NonCoherentIntegrations)*gnss_input.CoherentIntegrations*gnss_input.CodePeriod*gnss_input.fs;
dpe_input.CN0_est_ind = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigen = signalGen(gnss_input);
sigen.estimateTrueNoise = 1;
sigen.Nexpe = Nexpe;
sigen.numSV = numSV;
% meanNoise = computeMeanNoise(sigen);

NsamplesData = sigen.NsamplesData;
% estimateTrueNoise = sigen.estimateTrueNoise;
% Nexpe = sigen.Nexpe;

% if estimateTrueNoise == 0
%     meanNoise = nan;
% else
%% Add AWGN noise to the transmitted signals
mean_noise=zeros(1,Nexpe);
for exp_idx=1:Nexpe
    noise = ( sqrt(1/2)*randn(1,NsamplesData) +1i* ...
        sqrt(1/2)*randn(1,NsamplesData) );
    r = correlateSignal(sigen,noise);
    mean_noise(exp_idx) = mean(mean(r,2));
end
meanNoise=mean(mean_noise);

    %% Start Simulation
for CNo_idx=1:length(CNosim)
    CNosim(CNo_idx)
    for exp_idx=1:Nexpe
        CNo=CNosim(CNo_idx)*ones(numSV,1);
            
        %% Signal + noise 
        x = receivedSignal(sigen,CNo);
            
        %% Perform coherent/non-coherent integration times
        r = correlateSignal(sigen,x);
            
        %% 2-steps: Conventional approach estimation
        PosErrLS(CNo_idx,exp_idx) = conv2stepsPVT(r,ls_input);
                                               
        %% DPE approach ARS (accelerated random search)
        [PosErrDPE(CNo_idx,exp_idx), CN0_est(CNo_idx,exp_idx)] = DPEarsPVT(r,dpe_input);
    end
end
% Compute RMSEs
RMSE_LS=sqrt(mean(PosErrLS.^2,2));
RMSE_DPE=sqrt(mean(PosErrDPE.^2,2));
averageCn0= (mean(mean(cn0,2),3));
    

%% Ziv-Zakai bound computation
computeZZB

T_CRB_ZZB = table(CNosim(:), RMSE_LS(:), fZZLB_2SP(:), RMSE_DPE(:), fZZLB_DPE(:), ...
                  'VariableNames', {'CN0_dBHz','RMSE_LS_m','fZZLB_2SP','RMSE_DPE', 'fZZLB_DPE'});
disp(T_CRB_ZZB)

% return
%% PLOTS
figure,
h=semilogy(CNosim,RMSE_LS,'b-.',CNosim,RMSE_DPE,'b',CNosim,fZZLB_2SP,'r-.',CNosim,fZZLB_DPE,'r');
legend('MLE 2SP','MLE DPE','ZZB 2SP','ZZB DPE', fontsize=16)
grid
set(h,'Linewidth',2)
xlabel('CN0 [dB-Hz]')
ylabel('RMSE [m]')

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

%% number of samples calculation
NsamplesLocal=CodePeriod*fs*CoherentIntegrations;   % Number of samples of the Local Replica
NsamplesData=CodePeriod*fs*CoherentIntegrations*NonCoherentIntegrations;    %Number of samples of the Received Signal (Data).


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
%     fft_local(kSV,:) = fft(x_local(kSV,:),Nsamples);
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
sigen.NonCoherentIntegrations = NonCoherentIntegrations;
sigen.fs = fs;
sigen.fn = fn;
sigen.order = order;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function meanNoise = computeMeanNoise(sigen)


NsamplesData = sigen.NsamplesData;
% estimateTrueNoise = sigen.estimateTrueNoise;
Nexpe = sigen.Nexpe;

% if estimateTrueNoise == 0
%     meanNoise = nan;
% else
    %% Add AWGN noise to the transmitted signals
    mean_noise=zeros(1,Nexpe);
    for exp_idx=1:Nexpe
        noise = ( sqrt(1/2)*randn(1,NsamplesData) +1i* sqrt(1/2)*randn(1,NsamplesData) );
        r = correlateSignal(sigen,noise);
        mean_noise(exp_idx) = mean(mean(r,2));
    end
    meanNoise=mean(mean_noise);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = receivedSignal(sigen, CNo)



x_delay = sigen.x_delay;
NsamplesData = sigen.NsamplesData;
numSV = sigen.numSV;
fs = sigen.fs;
fn = sigen.fn;
order = sigen.order;


%% memory allocation
x_delay_noise=zeros(numSV,NsamplesData);



%% Add AWGN noise to the transmitted signals
noise = ( sqrt(1/2)*randn(1,NsamplesData) +1i* sqrt(1/2)*randn(1,NsamplesData) );
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r = correlateSignal(sigen,x_delay_noise)


% %% load configuration file
% eval(config)
numSV = sigen.numSV;
fft_local = sigen.fft_local;
NsamplesLocal = sigen.NsamplesLocal;
NonCoherentIntegrations = sigen.NonCoherentIntegrations;

%% memory allocation
r=zeros(numSV,NsamplesLocal);
    %% Perform NonCoherentIntegrations times non coherent integrations of CoherentIntegrations times coherent integrations.
for kSV=1:numSV
    for idx_nc=1:NonCoherentIntegrations
        r(kSV,:)=r(kSV,:)+abs(ifft(fft(x_delay_noise(1,NsamplesLocal*(idx_nc-1)+1:NsamplesLocal*idx_nc),NsamplesLocal) .* conj(fft_local(kSV,:)))).^2;
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PosErrLS = conv2stepsPVT(r,ls_input)


UserPosition = ls_input.UserPosition;
c = ls_input.c;
numSV = ls_input.numSV;
fs = ls_input.fs;
num2stepsIterations = ls_input.num2stepsIterations;
SatPosition = ls_input.SatPosition;

%% memory allocation
EstRange=zeros(1,numSV);


RefPos=UserPosition+10*(2*rand(3,1)-1)';
EstRxPVT=[RefPos];


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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [PosErrDPE, CN0_est] = DPEarsPVT(r, dpe_input)

UserPosition = dpe_input.UserPosition;
c = dpe_input.c;
numSV = dpe_input.numSV;
fs = dpe_input.fs;
fn = dpe_input.fn;
Tc = dpe_input.Tc;
SatPosition = dpe_input.SatPosition;
Niter = dpe_input.Niter;
contraction = dpe_input.contraction;
dmax = dpe_input.dmax;
dmin = dpe_input.dmin;
dmax_clk = dpe_input.dmax_clk;
dmin_clk = dpe_input.dmin_clk;
NormalizaFactor = dpe_input.NormalizaFactor;
CN0_est_ind = dpe_input.CN0_est_ind;
dt=1/fs;                            % Sampling period
randomDelay = 0; %%% magic number....
gamma_est=zeros(3,Niter+1);
amp_est = zeros(numSV, Niter+1);
EstRxClkBias=zeros(1,Niter+1);

%% memory allocation
EstRange=zeros(1,numSV);
MaxCorr=zeros(1,numSV);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

