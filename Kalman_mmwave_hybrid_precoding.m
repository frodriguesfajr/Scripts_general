%%%%%%%%%%----- Performance of Hybrid Kalman based beamforming in MmWave Channels-----%%%%%%%
% Authors: Anna Vizziello and Pietro Savazzi                                              
% Date: September 21, 2018
% This code calculates the spectral efficiency achieved by the multi-user 
% hybrid precoding algorithms in the paper - A. Vizziello, P. Savazzi and K. Chowdhury, 
% "A Kalman based Hybrid Precoding for Multi-User Millimeter Wave Massive MIMO Systems,"  
% in IEEE Access, 2018.
%
% It has been developed starting from the code, licensed under a Creative Commons
% Attribution-NonCommercial 4.0 International License, that may be found at:
% http://www.profheath.org/research/millimeter-wave-cellular-systems/hybrid-precoding-and-channel-estimation/beamforming-and-precoding/matlab-code-multi-user-hybrid-precoding/
% 
% The orginal code has been modified by adding:
% The hybrid precoders: MMSE, fully-digital MMSE and ZF, the proposed Kalman based one;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%--------------------------------------------------------------------------
clear all;close all;clc;warning('off');
rng(42, 'twister'); % Define a semente com o gerador Mersenne Twister
% ----------------------------- System Parameters -------------------------
Num_users=4; % Number of users
TX_ant=64; %Number of UPA TX antennas
TX_ant_w=sqrt(TX_ant); % width
TX_ant_h=sqrt(TX_ant); % hieght 
ind_TX_w=reshape(repmat([0:1:TX_ant_w-1],TX_ant_h,1),1,TX_ant_w*TX_ant_h);
ind_TX_h=repmat([0:1:TX_ant_h-1],1,TX_ant_w);

RX_ant=4; %Number of UPA RX antennas
RX_ant_w=sqrt(RX_ant); % width 
RX_ant_h=sqrt(RX_ant); % hieght
ind_RX_w=reshape(repmat([0:1:RX_ant_w-1],RX_ant_h,1),1,RX_ant_w*RX_ant_h);
ind_RX_h=repmat([0:1:RX_ant_h-1],1,RX_ant_w);

% ----------------------------- Channel Parameters ------------------------
Num_paths=10; %Number of channel paths
 
% ----------------------------- Simulation Parameters ---------------------
SNR_dB_range=-10:5:20;  % SNR in dB
Rate_SU=zeros(1,length(SNR_dB_range)); % Will carry the single-user MIMO rate (without interference)
Rate_BS=zeros(1,length(SNR_dB_range));% Will carry the rate with analog-only beamsteering
Rate_HP=zeros(1,length(SNR_dB_range)); % Will carry the rate of hybrid ZF precoding 
Rate_HP_MSE=zeros(1,length(SNR_dB_range)); % Will carry the rate of the hybrid MMSE precoder
Rate_HP_Kalman=zeros(1,length(SNR_dB_range)); % Will carry the rate of the proposed hybrid Kalman precoder
Rate_HP_FD_ZF=zeros(1,length(SNR_dB_range)); % Will carry the rate of the fully digital hybrid ZF precoding
Rate_HP_FD_MSE=zeros(1,length(SNR_dB_range)); % Will carry the rate of the fully digital hybrid MMSE precoding

ITER=1; % Number of iterations
% disp([pi*rand(1,Num_paths)-pi/2])
%disp([2*pi*rand(1,Num_paths)])
% disp([sqrt(1/Num_paths)*sqrt(1/2)])
%disp(randn(1, Num_paths))
% parte_real = norminv(rand(1,Num_paths),0,1);
% parte_imag = 1j*norminv(rand(1,Num_paths),0,1);
% disp([randn(1,Num_paths)+1j*randn(1,Num_paths)])
%disp([sqrt(1/Num_paths)*sqrt(1/2)*(randn(1,Num_paths)+1j*randn(1,Num_paths))])
% --------------- Simulation starts ---------------------------------------
for iter=1:1:ITER
    
    % Generate user channels 
    % [H,a_TX,a_RX]=generate_channels(Num_users,TX_ant_w,TX_ant_h,RX_ant_w,RX_ant_h,Num_paths); 
    % H is a 3-dimensional matrix, with Num_users,RX_ant,TX_ant dimensions
    H=zeros(Num_users,RX_ant_w*RX_ant_h,TX_ant_w*TX_ant_h);  % One user channel
    a_TX=zeros(TX_ant_w*TX_ant_h,Num_users); % TX steering vector
    a_RX=zeros(RX_ant_w*RX_ant_h,Num_users); % RX steering vector
    ind_TX_w=reshape(repmat([0:1:TX_ant_w-1],TX_ant_h,1),1,TX_ant_w*TX_ant_h);
    ind_TX_h=repmat([0:1:TX_ant_h-1],1,TX_ant_w);
    ind_RX_w=reshape(repmat([0:1:RX_ant_w-1],RX_ant_h,1),1,RX_ant_w*RX_ant_h);
    ind_RX_h=repmat([0:1:RX_ant_h-1],1,RX_ant_w);
    % % Constructing the channels
    for u=1:1:Num_users
         AoD_el(u,:)=pi*rand(1,Num_paths)-pi/2;
         AoD_az(u,:)=2*pi*rand(1,Num_paths);
         AoA_el(u,:)=pi*rand(1,Num_paths)-pi/2;
         AoA_az(u,:)=2*pi*rand(1,Num_paths);
         alpha(u,:)= sqrt(1/Num_paths)*(norminv(rand(1,Num_paths),0,1) + ...
             1j*norminv(rand(1,Num_paths),0,1));
         Temp_Channel=zeros(RX_ant_w*RX_ant_h,TX_ant_w*TX_ant_h);
         for l=1:1:Num_paths
             parte1 = 1/(TX_ant_w*TX_ant_h);
             parte2 = ind_TX_w*sin(AoD_az(u,l)) * sin(AoD_el(u,l));
             parte3 = ind_TX_h*cos(AoD_el(u,l));
             parte4 = exp(1j*pi*(parte2+parte3));
             a_TX(:,u)=transpose(sqrt(parte1)*parte4);
             parte5 = 1/(RX_ant_w*RX_ant_h);
             parte6 = ind_RX_w*sin(AoA_az(u,l))*sin(AoA_el(u,l));
             parte7 = ind_RX_h*cos(AoA_el(u,l));
             parte8 = exp(1j*pi*(parte6+parte7));
             a_RX(:,u)=transpose(sqrt(parte5)*parte8);
             parte9 = (TX_ant_w*TX_ant_h)*(RX_ant_w*RX_ant_h);
             parte10 = alpha(u,l)*a_RX(:,u);
             Temp_Channel=Temp_Channel+sqrt(parte9)*parte10*a_TX(:,u)';
         end
         H(u,:,:)=Temp_Channel;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Stage 1 of the proposed algorithm (Analog precoding)
    Frf=zeros(TX_ant,Num_users); % BS RF precoders 
    Wrf=zeros(RX_ant,Num_users); % MS RF precoders 

    for u=1:1:Num_users
        Frf(:,u)=a_TX(:,u);
        Wrf(:,u)=a_RX(:,u);
    end      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Constructin the effective channels
    for u=1:1:Num_users
        Channel=zeros(RX_ant,TX_ant);
        Channel(:,:)= H(u,:,:);
        He(u,:)=Wrf(:,u)'*Channel*Frf ;    % Effective channels
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % effective channel for fully digital precoding
    for u=1:1:Num_users
         Channel=zeros(RX_ant,TX_ant);
         Channel(:,:)= H(u,:,:);
         He_fd(u,:)=Wrf(:,u)'*Channel;    
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Baseband zero-forcing precoding
    Fbb=He'*(He*He')^(-1); 
    for u=1:1:Num_users % Normalization of the hybrid precoders
        Fbb(:,u)=Fbb(:,u)/sqrt((Frf*Fbb(:,u))'*(Frf*Fbb(:,u)));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Fully-digital zero-forcing precoding
    Ffd=He_fd'*pinv(He_fd*He_fd');
    for u=1:1:Num_users % Normalization of the hybrid precoders
         Ffd(:,u)=Ffd(:,u)/sqrt((Ffd(:,u))'*(Ffd(:,u)));
    end 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Spectral efficiency calculations
    count=0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for SNR_dB=SNR_dB_range
        count=count+1;
        SNR=10^(.1*SNR_dB)/Num_users; % SNR value
        sigma2=1/SNR;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % MMSE baseband precoder
        FbbMSE=inv(He'*He+Num_users*sigma2*Frf'*Frf)*He';
        for u=1:1:Num_users % Normalization of the hybrid precoders
            FbbMSE(:,u)=FbbMSE(:,u)/sqrt((Frf*FbbMSE(:,u))'*(Frf*FbbMSE(:,u)));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Fully-Digital MMSE Precoding
        FfdMSE=inv(He_fd'*He_fd+Num_users*sigma2*eye(TX_ant))*He_fd';
        for u=1:1:Num_users % Normalization of the hybrid precoders
            FfdMSE(:,u)=FfdMSE(:,u)/sqrt((FfdMSE(:,u))'*(FfdMSE(:,u)));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Kalman baseband precoder
        Fbbk=eye(Num_users,Num_users);
        RN1=Fbbk*Fbbk';
        Qm=eye(Num_users)*sigma2;
        ITERK=10; % Number of Kalman iterations
        for ii=1:ITERK
            Hk=He;
            K=RN1*Hk'*pinv(Hk*RN1*Hk'+Qm);
            errk=(eye(Num_users)-Hk*Fbbk);
            errk=errk/norm(errk);
            Fbbk=Fbbk+K*errk; 
            RN=RN1-K*Hk*RN1;
            RN1=RN; 
        end

        for u=1:1:Num_users % Normalization of the hybrid precoders
            Fbbk(:,u)=Fbbk(:,u)/sqrt((Frf*Fbbk(:,u))'*(Frf*Fbbk(:,u)));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for u=1:1:Num_users
            Int_set=[]; % interference index
            for i=1:1:Num_users
                if(i~=u)
                    Int_set=[Int_set i]; 
                end
            end
            Channel=zeros(RX_ant,TX_ant);
            Channel(:,:)= H(u,:,:);
            [U_channel S_channel V_channel]=svd(Channel);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Single-user rate
            % Rate_SU(count)=Rate_SU(count);
            % disp([log2(1+SNR)])
            % disp([S_channel(1,1)^2])
            % disp([log2(1+SNR*S_channel(1,1)^2)/(Num_users*ITER)])
            Rate_SU(count)=Rate_SU(count)+log2(1+SNR*S_channel(1,1)^2)/(Num_users*ITER);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Analog-only beamforming
            SINR_BS=(SNR*(abs(Wrf(:,u)'*Channel*Frf(:,u)).^2))/(SNR*sum((abs(Wrf(:,u)'*Channel*Frf(:,Int_set)).^2))+1);
            % disp([log2(1+SINR_BS)/(Num_users*ITER)])
            Rate_BS(count)=Rate_BS(count)+log2(1+SINR_BS)/(Num_users*ITER);            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Hybrid Precoding
        termo_a = eye(Num_users);
        termo_b = Fbb*Fbb';
        termo_c = He*(termo_b)*He';
        termo_d = log2(det(termo_a+SNR*(termo_c)));
        Rate_HP(count)=Rate_HP(count)+termo_d/(Num_users*ITER);
        % Hybrid Precoding MMSE
        Rate_HP_MSE(count)=Rate_HP_MSE(count)+log2(det(eye(Num_users)+SNR*(He*(FbbMSE*FbbMSE')*He')))/(Num_users*ITER);
        % Hybrid Precoding Kalman
        Rate_HP_Kalman(count)=Rate_HP_Kalman(count)+log2(det(eye(Num_users)+SNR*(He*(Fbbk*Fbbk')*He')))/(Num_users*ITER);
        % ZF fully digital precoding
        Rate_HP_FD_ZF(count)=Rate_HP_FD_ZF(count)+log2(det(eye(Num_users)+SNR*(He_fd*(Ffd*Ffd')*He_fd')))/(Num_users*ITER);
        % MSE fully digital precoding
        Rate_HP_FD_MSE(count)=Rate_HP_FD_MSE(count)+log2(det(eye(Num_users)+SNR*(He_fd*(FfdMSE*FfdMSE')*He_fd')))/(Num_users*ITER);
    end % End of SNR loop
end % End of ITER loop

%Plotting the spectral efficiency
    figure(1),plot(SNR_dB_range,Rate_SU,'-mv','LineWidth',1.5); hold on;
    plot(SNR_dB_range,Rate_HP,'-rs','LineWidth',1.5);
    hold on; plot(SNR_dB_range,Rate_HP_MSE,'-b*','LineWidth',1.5);
    hold on; plot(SNR_dB_range,Rate_HP_Kalman,'-go','LineWidth',1.5);
    hold on; plot(SNR_dB_range,Rate_HP_FD_MSE,'--k','LineWidth',2);
if Num_paths==1
    hold on; plot(SNR_dB_range,Rate_BS,'-ro','LineWidth',1.5);
    legend('Single-user (No Interference)','ZF Hybrid Precoding','MMSE Hybrid Precoding','Kalman Hybrid Precoding','MSE Fully-Digital Precoding','Analog-only Beamsteering');

else
    hold on; plot(SNR_dB_range,Rate_BS,'-ko','LineWidth',1.5);
    legend('Single-user (No Interference)','ZF Hybrid Precoding','MMSE Hybrid Precoding','Kalman Hybrid Precoding','MSE Fully-Digital Precoding','Analog-only Beamsteering');
end
xlabel('SNR (dB)','FontSize',12);
ylabel('Spectral Efficiency (bps/ Hz)','FontSize',12);
grid;
