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

close all;
clear all; 
clc; 
format long;


%% load configuration file
config = 'ConfigFile';
eval(config)

%% Memory allocation.
PosErrLS=zeros(length(CNosim),Nexpe);
PosErrDPE=zeros(length(CNosim),Nexpe);
CN0_est=zeros(length(CNosim),Nexpe);
cn0= zeros(length(CNosim),numSV,Nexpe);
%% Signal generation (sigen struct)
sigen = signalGen(config);
meanNoise = computeMeanNoise(config,sigen);

if simulate_mle
    %% Start Simulation
    for CNo_idx=1:length(CNosim)
        CNosim(CNo_idx)
        
        for exp_idx=1:Nexpe
            CNo=CNosim(CNo_idx)*ones(numSV,1);

            %% Signal + noise 
            x = receivedSignal(sigen,config,CNo);
        
         
            %% Apply RIM
            if RIMuse
                x = RIM(x);
            end

            %% Perform coherent/non-coherent integration times
            r = correlateSignal(sigen,config,x);

            %% Estimate CN0
            % cn0(CNo_idx,:,exp_idx) = estimateCn0(r,config,meanNoise);

            %% 2-steps: Conventional approach estimation
            PosErrLS(CNo_idx,exp_idx) = conv2stepsPVT(r,config);

            %% DPE approach ARS (accelerated random search)
            [PosErrDPE(CNo_idx,exp_idx), CN0_est(CNo_idx,exp_idx)] = DPEarsPVT(r,config);

        end
    end
    
    % Compute RMSEs
    RMSE_LS=sqrt(mean(PosErrLS.^2,2));
    RMSE_DPE=sqrt(mean(PosErrDPE.^2,2));
    averageCn0= (mean(mean(cn0,2),3));
    
end 
%% PLOTS
figure,
h=semilogy(CNosim,RMSE_LS,'b-.',CNosim,RMSE_DPE,'b');
legend('MLE 2SP','MLE DPE', fontsize=16)
grid
set(h,'Linewidth',2)
xlabel('CN0 [dB-Hz]')
ylabel('RMSE [m]')
return

if plot_estimated_cn0
    figure
    h=plot(CNosim,averageCn0,CNosim,CNosim);
    legend('Estimated CN0','True CN0', fontsize=16)
    grid
    set(h,'Linewidth',2)
end

function [sigen] = signalGen(config)


%% load configuration file
eval(config)


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



end

function [PosErrDPE, CN0_est] = DPEarsPVT(r,config)

%% load configuration file
eval(config)

randomDelay = 0; %%% magic number....

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

function x = receivedSignal(sigen,config,CNo)


%% load configuration file
eval(config)

x_delay = sigen.x_delay;
NsamplesData = sigen.NsamplesData;


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

function PosErrLS = conv2stepsPVT(r,config)

%% load configuration file
eval(config)

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