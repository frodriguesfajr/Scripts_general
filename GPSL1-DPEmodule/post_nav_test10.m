%% Limpeza do ambiente
clear;           % Remove todas as variáveis da memória
clc;             % Limpa a janela de comandos
% addpath ('.\include')  
% addpath ('.\common') 
%% Initialize constants, settings =========================================
settings = initSettings();
[fid, message] = fopen(settings.fileName, 'rb');
data_gps = prepareGNSSProcessing(0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[acqResults, data_aq] = acquireSatellites(settings, fid, 0);
[trackResults, channel, outfil] = trackSatellites(acqResults, ...
    fid, settings, 0);
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

        % Calcule o sigma_pr médio para todos os satélites ativos
        sigma_pr_all = zeros(length(activeChnList),1);
        for idx = 1:length(activeChnList)
            ch = activeChnList(idx);
            % Pegue o C/N0 para cada canal (pode usar média, última amostra etc)
            CNo = mean(trackResults(ch).CNo.VSMValue);
            sigma_pr_all(idx) = settings.c / (10^(CNo/10)); % modelo simplificado
        end

        % Você pode usar média ou valor típico para sigma_pr:
        sigma_pr = mean(sigma_pr_all);


        % Corrige as pseudodistâncias com correção de relógio dos satélites
        rawP_arr_sat = rawP_arr(activeChnList, currMeasNr)';
        clkCorrRawP = rawP_arr_sat + satClkCorr * settings.c;

        %% Aplica algoritmo de mínimos quadrados para estimar posição + dt
         [xyzdt, el_ls, az_ls, DOP_ls, satPositions_ls] = ...
             leastSquarePos(satPos, clkCorrRawP, settings);


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

    %% Impressão do erro de posição, se existir
    fprintf('\nCurrent 2D Error of DPE: %1f\n',...
         (navSolutions.LLH_error(currMeasNr,2)));
    fprintf('Current 2D Error of LS: %.6f m\n', LS_error(currMeasNr, 1));
    
    % return
end

% return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ls_2d_list_new = navSolutions.LLH_error(:,1)';
ls_2d_list_new = LS_error(:,1)';
dpe_2d_list_new = navSolutions.LLH_error(:,2)';
navSolPeriod = 500;

% === Plot 2D positioning error ===========================================

figure;

% Plotagem das duas curvas com mesmo estilo
plot(ls_2d_list_new, 'Color', [80, 29, 138]/255, 'LineWidth', 3);
hold on;
plot(dpe_2d_list_new, 'Color', [230, 57, 70]/255, 'LineWidth', 3);  % Usei outra cor para diferenciar

% Legenda para as duas curvas
legend('Least Squares (Scalar Tracking)', 'DPE Method', 'FontSize', 12);

% Eixos e título
xlabel(sprintf('Epoch (%d ms)', navSolPeriod), 'FontSize', 12);
ylabel('2D error (m)', 'FontSize', 12);
xlim([1, max(length(ls_2d_list_new), length(dpe_2d_list_new))]);
title('2D Positioning Error');

% Grade e formatação
grid on;
set(gca, 'FontSize', 12);

% disp('   Processing is complete for this data block');