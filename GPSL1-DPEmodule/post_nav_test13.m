%% Limpeza do ambiente
clear;           % Remove todas as variáveis da memória
clc;             % Limpa a janela de comandos
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize constants, settings =========================================
settings = initSettings();
fileName = 'C:\Repository\Scripts_general\GPSL1-DPEmodule\IF_Data_Set\Medium Urban in TST with one NLOS.dat';
% fileName = 'C:\Repository\Scripts_general\GPSL1-DPEmodule\IF_Data_Set\Open Sky GPS L1.dat';
% fileName = 'C:\Repository\Scripts_general\GPSL1-DPEmodule\IF_Data_Set\Open Sky Beidou B1I.dat';
[fid, message] = fopen(fileName, 'rb');
data_gps = prepareGNSSProcessing(fileName, 0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[acqResults, data_aq] = acquireSatellites(settings, fid, 1);
[trackResults, channel, outfil] = trackSatellites(fileName, acqResults, ...
    fid, settings, 1);
settings.outfile_root = outfil; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate navigation solutions =====================================
disp('   Calculating navigation solutions...');
[navSolutions, eph_ch, activeChnList] = ...
    verifyAndInitializeGNSS(settings, trackResults);

% Inicialização de vetores auxiliares
[subFrameStart, TOW] = deal(inf(1, settings.numberOfChannels));
%% Loop sobre canais ativos
for channelNr = activeChnList
    PRN = trackResults(channelNr).PRN;
    fprintf('Decoding NAV for PRN %02d --------------------------\n', PRN);

    % Decodifica e valida efemérides
    [success, eph, sfStart, tow] = ...
        decodeAndValidateNav(trackResults(channelNr), settings);

    if success
        eph_ch(PRN) = eph;
        subFrameStart(channelNr) = sfStart;
        TOW(channelNr) = tow;
        fprintf('  [OK] Ephemeris completa para PRN %02d!\n', PRN);
    else
        activeChnList = setdiff(activeChnList, channelNr);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (isempty(activeChnList) || (size(activeChnList, 2) < 4))
    disp(['Too few satellites with ephemeris data for postion ' ...
        'calculations. Exiting!']);
    return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcula os limites comuns de amostra para os canais ativos
[sampleStart, sampleEnd] = getSampleWindow(trackResults, ...,
    subFrameStart, activeChnList);
% Ajusta os limites para janela comum entre canais
sampleStart = max(sampleStart) + 1;
sampleEnd = min(sampleEnd) - 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loop pelos canais GNSS prontos para uso (decodificados corretamente)
for channelNr = activeChnList
    % Obtém o índice da primeira amostra do subquadro detectado
    sampleStart(channelNr) = ...
        trackResults(channelNr).absoluteSample(subFrameStart(channelNr));
    % Obtém o índice da última amostra disponível no rastreamento
    sampleEnd(channelNr) = ...
        trackResults(channelNr).absoluteSample(end);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Determina a janela de amostras comum entre todos os canais ativos
% Define o início como a última posição de início entre os canais
sampleStart = max(sampleStart) + 1;
% Define o fim como a primeira posição de término entre os canais
sampleEnd = min(sampleEnd) - 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cálculo do passo de amostragem entre soluções de navegação
% Ex: se samplingFreq = 5 MHz e navSolPeriod = 1000 ms,
% então 5e6 * 1 = 5e6 → uma solução a cada 5 milhões de amostras
measSampleStep = fix(settings.samplingFreq * settings.navSolPeriod / 1000);
%% Número total de épocas de navegação possíveis dentro da janela
measNrSum = fix((sampleEnd - sampleStart) / measSampleStep);
LS_error = zeros(measNrSum,1);
%% Inicializações para cálculo de posição
% Elevações dos satélites
satElev = inf(1, settings.numberOfChannels);  
% Canais que seguem aptos para posicionamento
readyChnList = activeChnList;
% Tempo local ainda indefinido
localTime = inf;                              
%% Mensagem para o usuário
fprintf('Positions are being computed. Please wait...\n');
el_arr = zeros(channelNr, measNrSum);
az_arr = zeros(channelNr, measNrSum);
rawP_arr = zeros(length(TOW), measNrSum);
DOP_arr = zeros(5, measNrSum);
dt_arr = zeros(1, measNrSum);
latitude_arr = zeros(1, measNrSum);
longitude_arr = zeros(1, measNrSum);
height_arr = zeros(1, measNrSum);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loop principal para calcular a posição em cada época de medição
navSolutions = [];
measNrSum = 4;
for currMeasNr = 1:measNrSum
    fprintf('Fix: Processing %02d of %02d \n', currMeasNr, measNrSum);
    %% Filtra canais com elevação acima do limiar (em graus)
    activeChnList = intersect(find(satElev >= settings.elevationMask),...
        readyChnList);
    %% Determina o índice da amostra atual para esta época
    currMeasSample = sampleStart + measSampleStep * (currMeasNr - 1);
    %% Calcula pseudodistâncias brutas e tempo de transmissão dos satélites
    [rawP_arr(:, currMeasNr), transmitTime, localTime, codePhase] = ...
        calculatePseudoranges(trackResults, subFrameStart, TOW, ...
                              currMeasSample, localTime, activeChnList,...
                              settings);
    %% Calcula as posições dos satélites e correções de relógio
    % tempo de transmissão para cada canal do receptor
    tT_meas = transmitTime(activeChnList);
    % lista dos PRNs de cada satélite
    prn_list_meas = [trackResults(activeChnList).PRN];
    [satPos, satClkCorr] = satpos(tT_meas, prn_list_meas,eph_ch);

    % Solução de posição é possível se houver pelo menos 4 
    % satélites visíveis
    if numel(activeChnList) > 3

        %--- Calcular pseudodistâncias e posições de satélites ------------
        rawP = rawP_arr(activeChnList,currMeasNr)' + satClkCorr*settings.c;

        % Calcule o sigma_pr médio para todos os satélites ativos
        sigma_pr_all = zeros(length(activeChnList),1);

        
        % Corrige as pseudodistâncias com correção de relógio dos satélites
        rawP_arr_sat = rawP_arr(activeChnList, currMeasNr)';
        clkCorrRawP = rawP_arr_sat + satClkCorr * settings.c;

        %% Aplica algoritmo de mínimos quadrados para estimar posição + dt
         [xyzdt, el_ls, az_ls, DOP_ls, satPositions_ls] = ...
             leastSquarePos(satPos, clkCorrRawP, settings);
         rec_pos = xyzdt(1:3)';      % 3×1

         % --- 4) Prepara variável para σ_pr e Jacobiano H_i --------------
         M = numel(activeChnList);
         sigma_pr_all = zeros(M,1);
         J_tautau     = zeros(4,4);
         T            = zeros(M,4);

         for i = 1:M
             ch      = activeChnList(i);
             sat_p   = satPos(:,i);
             %--- modelo realista de σ_pr ---------------------------------
             CNo_db  = mean(trackResults(ch).CNo.VSMValue);
             CNo_lin = 10^(CNo_db/10);
             sigma_i = settings.c / (2*sqrt(2)*pi * settings.codeFreqBasis * sqrt(CNo_lin));
             sigma_pr_all(i) = sigma_i;
             %--- Jacobiano H_i = dτ_i/d[x y z dt] ------------------------
             diff = sat_p - rec_pos;      
             ri   = norm(diff);
             H_i  = [ -diff'/ri, 1 ];    % 1×4
             J_tautau = J_tautau + (H_i'*H_i)/(sigma_i^2);  % soma FIM de τ
             T(i,:)    = H_i;                             % reutiliza H_i em T
         end

         %--- 5) CRB _DPE_ (direto) --------------------------------------------
         Cov_DPE = inv(J_tautau);
         std_DPE = sqrt(diag(Cov_DPE));  % [σx σy σz σdt]
         navSolutions.CRLB_DPE_2D(currMeasNr) = norm(std_DPE(1:2));
         navSolutions.CRLB_DPE_3D(currMeasNr) = norm(std_DPE(1:3));
         
         %--- 6) CRB _LS_ (via Eq.41) ------------------------------------------
         W        = diag(1./(sigma_pr_all.^2));
         A        = (T' * W * T)\(T' * W);               % (T^H W T)^(-1) T^H W
         Cov_LS   = A * inv(J_tautau) * A';
         std_LS   = sqrt(diag(Cov_LS));
         navSolutions.CRLB_LS_2D(currMeasNr) = norm(std_LS(1:2));
         navSolutions.CRLB_LS_3D(currMeasNr) = norm(std_LS(1:3));
         % 
         % % Inicialização da matriz FIM
         % FIM = zeros(4, 4);
         % % Posição estimada atual
         % rec_pos = xyzdt(1:3);
         % J_tautau = zeros(4,4);
         % % --- bloco DPE: γ = [x y z dt] para caso estático (1 antena) ---
         % J_gammaGamma = zeros(4,4);
         % % rec_params   = [xyzdt'; 0;0;0; dt_arr(currMeasNr)]; % [x y z vx vy vz dt]
         % 
         % for idx = 1:length(activeChnList)
         %    ch = activeChnList(idx);
         %    % Pegue o C/N0 para cada canal (pode usar média, última amostra etc)
         %    CNo_db = mean(trackResults(ch).CNo.VSMValue);
         %    CNo_lin  = 10^(CNo_db/10);       % converte para linear
         %    % desvio padrão [m]
         %    sigma_i = settings.c / (2*sqrt(2)*pi * settings.codeFreqBasis * sqrt(CNo_lin));
         %    sigma_pr_all(idx) = sigma_i;
         %    % fprintf('Sigma_pr (modelo realista): %.3f m\n', sigma_i)
         %    % posição do satélite
         %    sat_pos = satPos(:, idx);
         %    diff_vec = sat_pos - rec_pos';
         %    r_i = norm(diff_vec);
         %    % Jacobiano de tau_i em relação a [x y z dt]
         %    H_i = [diff_vec'/r_i, 1];
         %    J_tautau = J_tautau + (H_i'*H_i)/(sigma_i^2);
         %    % soma na FIM
         %    J_gammaGamma = J_gammaGamma + (H_i' * H_i) / sigma_i^2;
         %    FIM = FIM + (H_i' * H_i) / sigma_i^2;    % soma na FIM          
         % end
         % 
         % 
         % 
         % M = length(activeChnList);
         % T = zeros(M,4);
         % W = diag(1./(sigma_pr_all.^2));
         % for i=1:M
         %     diff   = satPos(:,i) - rec_pos';
         %     ri     = norm(diff);
         %     T(i,:) = [ -diff'/ri, 1 ];    % mesmo H_i acima
         % end
         % % Cov_tautau = inv(J_tautau);    % Cov( [x y z dt] ) em unidades de tempo (s^2)
         % % 
         % A = (T' * W * T)\(T' * W);      % = (T^H W T)^{-1} T^H W
         % % C_LS = settings.c^2 * (A * Cov_tautau * A');  % Cov([x y z dt]) em m^2
         % 
         % 
         % Cov_DPE = settings.c^2 * inv(J_gammaGamma);  
         % CRLB_DPE_std = sqrt(diag(Cov_DPE));  % [σ_x σ_y σ_z σ_dt]
         % navSolutions.CRLB_DPE_xyzdt(currMeasNr, :) = CRLB_DPE_std;
         % 
         % Cov_LS     = inv(J_tautau);
         % C_LS       = settings.c^2*(A * Cov_LS * A');   % eq. (41)
         % CRLB_LS_std= sqrt(diag(C_LS));
         % navSolutions.CRLB_LS_2D(currMeasNr) = norm(CRLB_LS_std(1:2));
         % navSolutions.CRLB_LS_3D(currMeasNr) = norm(CRLB_LS_std(1:3));
         % 
         % 
         % CRLB_LS_std = sqrt(diag(C_LS));               % desvios-padrão [m, m, m, s]
         % CRLB_LS_2D = norm(CRLB_LS_std(1:2));
         % CRLB_LS_3D = norm(CRLB_LS_std(1:3));
         % 
         % CRLB_matrix = inv(FIM);
         % % Armazenamento dos desvios padrão teóricos (CRLB)
         % CRLB_std_xyzdt = sqrt(diag(CRLB_matrix));
         % navSolutions.CRLB_xyzdt(currMeasNr, :) = CRLB_std_xyzdt;
         % % Opcional: Calcular erro 2D e 3D mínimo teórico
         % navSolutions.CRLB_2D(currMeasNr) = norm(CRLB_std_xyzdt(1:2));
         % navSolutions.CRLB_3D(currMeasNr) = norm(CRLB_std_xyzdt(1:3));
         % % ERROS TEÓRICOS 2D/3D
         % navSolutions.CRLB_DPE_2D(currMeasNr) = norm(CRLB_DPE_std(1:2));
         % navSolutions.CRLB_DPE_3D(currMeasNr) = norm(CRLB_DPE_std(1:3));
         


         el_arr(activeChnList, currMeasNr) = el_ls;
         az_arr(activeChnList, currMeasNr) = az_ls;
         DOP_arr(:, currMeasNr) = DOP_ls;

        % Desvio de relógio do receptor
        dt_arr(currMeasNr) = (currMeasNr == 1) * 0 + ...
                                      (currMeasNr > 1) * xyzdt(4);

        % Atualiza elevação dos satélites para o próximo ciclo
        satElev = el_arr(:, currMeasNr)';

        % Calcula pseudodistância corrigida final
        rawP_arr_sat = rawP_arr_sat + satClkCorr' * settings.c - xyzdt(4);

        %% Conversão de coordenadas ECEF para latitude/longitude/altura
        % 5 iterações para convergência
        [lat, long, height] = cart2geo(xyzdt(1), xyzdt(2), xyzdt(3), 5);  
        latitude_arr(currMeasNr) = lat;
        longitude_arr(currMeasNr) = long;
        height_arr(currMeasNr) = height;

        m2lat = 1/110734;
        m2lon = 1/103043;

        error_norm = norm(([lat, long] - ...
            settings.gt_llh(1:2))./[m2lat m2lon]);
        LS_error(currMeasNr,1) = error_norm;
        % disp(error_norm)
        % return

    else
        %% Caso não haja satélites suficientes para solução
        warning(['   Epoch ', num2str(currMeasNr), ...
                 ': Not enough satellites for position fix.']);
    end

    navSolutions.az = az_arr;
    navSolutions.el = el_arr;
    navSolutions.PRN = trackResults.PRN;


    navSolutions.latitude(currMeasNr) = lat;
    navSolutions.longitude(currMeasNr) = long;
    navSolutions.height(currMeasNr) = height;
    navSolutions.X(currMeasNr) = xyzdt(1);
    navSolutions.Y(currMeasNr) = xyzdt(2);
    navSolutions.Z(currMeasNr) = xyzdt(3);
    navSolutions.rawP = rawP_arr;
    %% DPE análise
    tic

    [navSolutions] = DPE_module_mod...
        (currMeasNr,navSolutions,activeChnList,...
        trackResults,currMeasSample,satPositions_ls,...
        tT_meas,localTime,...
        settings,satElev,fid,xyzdt(4),satClkCorr);
    % return
    % return
    navSolutions.DPE_processingtime(currMeasNr) = toc;




    % %%%%%%%%%%%%%%%%%%%
    % % localTime = localTime - xyzdt(4)/settings.c;       
    % % navSolutions.localTime(currMeasNr) = localTime;
    % % localTime = localTime + measSampleStep/settings.samplingFreq ;
    % %%%%%%%%%%%%%%%%%%%%%%%
    navSolutions.LLH_error(currMeasNr,2)=...
    norm(([navSolutions.DPE_latitude(currMeasNr),...
        navSolutions.DPE_longitude(currMeasNr)]...
        -settings.gt_llh(1:2))./[m2lat m2lon]);    
    %%
    navSolutions.LLH_error(currMeasNr,1)=...
        norm(([navSolutions.latitude(currMeasNr),...
        navSolutions.longitude(currMeasNr)]...
        -settings.gt_llh(1:2))./[m2lat m2lon]);
    GT_ECEF =  llh2xyz(settings.gt_llh.*[pi/180 pi/180 1]);
    navSolutions.LLH_error(currMeasNr,3)= ...
        sqrt(sum(([navSolutions.X(currMeasNr)...
        navSolutions.Y(currMeasNr)...
        navSolutions.Z(currMeasNr)]-GT_ECEF).^2));

    pos_xyz = ...
        llh2xyz([navSolutions.DPE_latitude(currMeasNr)/180*pi,...
        navSolutions.DPE_longitude(currMeasNr)/180*pi,...
        navSolutions.DPE_height(currMeasNr)]);

    navSolutions.LLH_error(currMeasNr,4)=  ...
        sqrt(sum(([pos_xyz]-GT_ECEF).^2));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   

    %% Impressão do erro de posição, se existir
    fprintf('\nCurrent 2D Error of DPE: %1f\n',...
         (navSolutions.LLH_error(currMeasNr,2)));
    fprintf('Current 2D Error of LS: %.6f m\n', LS_error(currMeasNr, 1));

    % === Prints the 3D errors of both Least Squares and DPE ==============

    fprintf('\nCurrent 3D Error of DPE     : %1f\n',...
        (navSolutions.LLH_error(currMeasNr,4)));
    fprintf('Current 3D Error of LS      : %1f\n',...
        (navSolutions.LLH_error(currMeasNr,3)));
    
    % return
end

% return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % ls_2d_list_new = navSolutions.LLH_error(:,1)';
% ls_2d_list_new = LS_error(:,1)';
% dpe_2d_list_new = navSolutions.LLH_error(:,2)';
% navSolPeriod = 500;
% 
% % === Plot 2D positioning error ===========================================
% 
% figure;
% 
% % Plotagem das duas curvas com mesmo estilo
% plot(ls_2d_list_new, 'Color', [80, 29, 138]/255, 'LineWidth', 3);
% hold on;
% plot(dpe_2d_list_new, 'Color', [230, 57, 70]/255, 'LineWidth', 3);  % Usei outra cor para diferenciar
% 
% % Legenda para as duas curvas
% legend('Least Squares (Scalar Tracking)', 'DPE Method', 'FontSize', 12);
% 
% % Eixos e título
% xlabel(sprintf('Epoch (%d ms)', navSolPeriod), 'FontSize', 12);
% ylabel('2D error (m)', 'FontSize', 12);
% xlim([1, max(length(ls_2d_list_new), length(dpe_2d_list_new))]);
% title('2D Positioning Error');
% 
% % Grade e formatação
% grid on;
% set(gca, 'FontSize', 12);
% 
% % disp('   Processing is complete for this data block');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save(['C:\Repository\Scripts_general\GPSL1-DPEmodule\navSolutions_execution.mat'], 'navSolutions');
% 'C:\Repository\Scripts_general\GPSL1-DPEmodule\navSolutions_execution.mat'
ls_2d_list = navSolutions.LLH_error(:,1)';
dpe_2d_list = navSolutions.LLH_error(:,2)';
ls_3d_list = navSolutions.LLH_error(:,3)';
dpe_3d_list = navSolutions.LLH_error(:,4)';
navSolPeriod = 500;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1);
% % Plotagem das duas curvas com mesmo estilo
% plot(ls_2d_list, 'Color', [80, 29, 138]/255, 'LineWidth', 3);
% hold on;
% plot(dpe_2d_list, 'Color', [230, 57, 70]/255, 'LineWidth', 3);
% % Legenda para as duas curvas
% legend('Least Squares (Scalar Tracking)', 'DPE Method', 'FontSize', 12);
% % Eixos e título
% xlabel(sprintf('Epoch (%d ms)', navSolPeriod), 'FontSize', 12);
% ylabel('2D error (m)', 'FontSize', 12);
% xlim([1, max(length(ls_2d_list), length(dpe_2d_list))]);
% title('2D Positioning Error');
% % Grade e formatação
% grid on;
% set(gca, 'FontSize', 12);
figure(1);
% Plot das curvas de erro LS e DPE
plot(ls_2d_list,  'Color', [80, 29, 138]/255, 'LineWidth', 3); hold on;
plot(dpe_2d_list, 'Color', [230, 57, 70]/255, 'LineWidth', 3);

% % Plot do CRLB 2D
% plot(navSolutions.CRLB_2D, '--k', 'LineWidth', 3);
% Limites CRLB
plot(navSolutions.CRLB_LS_2D,  '--', 'Color', [80,29,138]/255, 'LineWidth', 2);
plot(navSolutions.CRLB_DPE_2D, '--', 'Color', [230,57,70]/255, 'LineWidth', 2);
% Legend e labels
legend( ...
  'Erro LS 2D', ...
  'Erro DPE 2D', ...
  'CRLB LS 2D', ...
  'CRLB DPE 2D', ...
  'Location','Best' ...
);

% Eixos e título
xlabel(sprintf('Epoch (%d ms)', navSolPeriod), 'FontSize', 12);
ylabel('2D error (m)',                 'FontSize', 12);
xlim([1, max(length(ls_2d_list), length(dpe_2d_list))]);
title('2D Positioning Error vs. CRLB', 'FontSize', 14);

% Grade e formatação
grid on;
set(gca, 'FontSize', 12);
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(2);
% % Plotagem das duas curvas com mesmo estilo
% plot(ls_3d_list, 'Color', [80, 29, 138]/255, 'LineWidth', 3);
% hold on;
% plot(dpe_3d_list, 'Color', [230, 57, 70]/255, 'LineWidth', 3);
% % Legenda para as duas curvas
% legend('Least Squares (Scalar Tracking)', 'DPE Method', 'FontSize', 12);
% % Eixos e título
% xlabel(sprintf('Epoch (%d ms)', navSolPeriod), 'FontSize', 12);
% ylabel('3D error (m)', 'FontSize', 12);
% xlim([1, max(length(ls_3d_list), length(dpe_3d_list))]);
% title('3D Positioning Error');
% % Grade e formatação
% grid on;
% set(gca, 'FontSize', 12);
% % disp('   Processing is complete for this data block');
figure(2);
% Plot das curvas de erro LS e DPE
plot(ls_3d_list,  'Color', [80, 29, 138]/255, 'LineWidth', 3); hold on;
plot(dpe_3d_list, 'Color', [230, 57, 70]/255, 'LineWidth', 3);

% Plot do CRLB 3D
% plot(navSolutions.CRLB_3D, '--k', 'LineWidth', 3);
plot(navSolutions.CRLB_LS_3D, '--', 'Color', [80,29,138]/255, 'LineWidth', 2);
plot(navSolutions.CRLB_DPE_3D,'--', 'Color', [230,57,70]/255, 'LineWidth', 2);
% Legenda
legend( ...
  'Erro LS 3D', ...
  'Erro DPE 3D', ...
  'CRLB LS 3D', ...
  'CRLB DPE 3D', ...
  'Location','Best' ...
);

% Eixos e título
xlabel(sprintf('Epoch (%d ms)', navSolPeriod), 'FontSize', 12);
ylabel('3D error (m)',                      'FontSize', 12);
xlim([1, max(length(ls_3d_list), length(dpe_3d_list))]);
title('3D Positioning Error vs. CRLB',       'FontSize', 14);

% Grade e formatação
grid on;
set(gca, 'FontSize', 12);
hold off;
