%% Limpeza do ambiente
clear;           % Remove todas as variáveis da memória
clc;             % Limpa a janela de comandos

%% Carregamento dos arquivos .mat com structs de dados de simulação GNSS
load('settings_dpe_sim.mat');       % Carrega configurações da simulação (struct 'settings')
load('trackResults_dpe_sim.mat');   % Carrega resultados de rastreamento de satélites (struct 'trackResults')
load('navSolutions_dpe_sim.mat');   % Carrega soluções de navegação anteriores (struct 'navSolutions')

%%
%% Inicialização de identificador de arquivo
fid = 4;         % Pode ser usado posteriormente para logging em arquivo ou output formatado
%% Inclusão de diretórios com funções auxiliares (bibliotecas internas)
addpath ('C:\Repository\GPSL1-DPEmodule\include')   % Adiciona diretório com cabeçalhos/funções auxiliares
addpath ('C:\Repository\GPSL1-DPEmodule\common')    % Adiciona diretório com funções comuns usadas no processamento
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Verificação da duração mínima do sinal processado
if (settings.msToProcess < 36000)  % Verifica se há menos de 36 segundos de dados (36000 ms)
    % Exibe mensagem de erro e encerra o script
    disp('Record is too short. Exiting!');
    navSolutions = [];  % Zera estrutura de soluções de navegação
    eph          = [];  % Zera estrutura de efemérides (se esperada adiante)
    return                % Encerra o script imediatamente
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inicialização de vetores auxiliares
[subFrameStart, TOW] = deal(inf(1, settings.numberOfChannels));

% Identificação de canais ativos
activeChnList = find([trackResults.status] ~= '-');

% Loop sobre canais ativos
for channelNr = activeChnList
    PRN = trackResults(channelNr).PRN;
    fprintf('Decoding NAV for PRN %02d ----------------------------\n', PRN);

    % Decodifica e valida efemérides
    [success, eph, sfStart, tow] = decodeAndValidateNav(trackResults(channelNr), settings);

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
    disp('Too few satellites with ephemeris data for postion calculations. Exiting!');
    navSolutions = [];
    eph          = [];
    return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcula os limites comuns de amostra para os canais ativos
[sampleStart, sampleEnd] = getSampleWindow(trackResults, subFrameStart, activeChnList);
% Ajusta os limites para janela comum entre canais
sampleStart = max(sampleStart) + 1;
sampleEnd = min(sampleEnd) - 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
measSampleStep = calcSampleStep(settings);
measNrSum = calcNumEpochs(sampleStart, sampleEnd, measSampleStep);

% Inicializações de variáveis e mensagens
[satElev, readyChnList, localTime] = initPositionVars(settings, activeChnList);

fprintf('Positions are being computed. Please wait...\n');

% Inicialização das matrizes de soluções de navegação
[P, M, N, navSolutions] = initNavSolutions(channelNr, settings, measNrSum);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Positions are being computed. Please wait... \n');
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
satElev = inf(1, settings.numberOfChannels);  % Elevações dos satélites (placeholder)
readyChnList = activeChnList;                 % Canais que seguem aptos para posicionamento
localTime = inf;                              % Tempo local ainda indefinido
%% Mensagem para o usuário
fprintf('Positions are being computed. Please wait...\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loop principal para calcular a posição em cada época de medição
for currMeasNr = 1:measNrSum
    fprintf('Fix: Processing %02d of %02d \n', currMeasNr, measNrSum);
    %% Filtra canais com elevação acima do limiar (em graus)
    activeChnList = intersect(find(satElev >= settings.elevationMask), readyChnList);
    %% Armazena os PRNs dos canais válidos para esta época
    navSolutions.PRN(activeChnList, currMeasNr) = [trackResults(activeChnList).PRN];
    %% Inicializa variáveis com NaNs para garantir valores conhecidos
    navSolutions.el(:, currMeasNr)           = NaN;
    navSolutions.az(:, currMeasNr)           = NaN;
    navSolutions.transmitTime(:, currMeasNr) = NaN;
    navSolutions.satClkCorr(:, currMeasNr)   = NaN;
    %% Determina o índice da amostra atual para esta época
    currMeasSample = sampleStart + measSampleStep * (currMeasNr - 1);
    %% Calcula pseudodistâncias brutas e tempo de transmissão dos satélites
    [navSolutions.rawP(:, currMeasNr), transmitTime, localTime, codePhase] = ...
        calculatePseudoranges(trackResults, subFrameStart, TOW, ...
                              currMeasSample, localTime, activeChnList, settings);

    %% Armazena o tempo de transmissão para os canais ativos
    navSolutions.transmitTime(activeChnList, currMeasNr) = transmitTime(activeChnList);

    %% Calcula as posições dos satélites e correções de relógio
    [satPositions, satClkCorr] = satpos(...
        transmitTime(activeChnList), ...
        [trackResults(activeChnList).PRN], ...
        eph_ch);

    navSolutions.satClkCorr(activeChnList, currMeasNr) = satClkCorr;

    %% Solução de posição é possível se houver pelo menos 4 satélites visíveis
    if numel(activeChnList) > 3
        % Corrige as pseudodistâncias com correção de relógio dos satélites
        clkCorrRawP = navSolutions.rawP(activeChnList, currMeasNr)' + ...
                      satClkCorr * settings.c;

        %% Aplica algoritmo de mínimos quadrados para estimar posição + dt
        [xyzdt, navSolutions.el(activeChnList, currMeasNr), ...
         navSolutions.az(activeChnList, currMeasNr), ...
         navSolutions.DOP(:, currMeasNr), ~] = ...
            leastSquarePos(satPositions, clkCorrRawP, settings);

        % Armazena as coordenadas ECEF
        navSolutions.X(currMeasNr) = xyzdt(1);
        navSolutions.Y(currMeasNr) = xyzdt(2);
        navSolutions.Z(currMeasNr) = xyzdt(3);

        % Desvio de relógio do receptor
        navSolutions.dt(currMeasNr) = (currMeasNr == 1) * 0 + ...
                                      (currMeasNr > 1) * xyzdt(4);

        % Guarda o índice de amostra correspondente à medição
        navSolutions.currMeasSample(currMeasNr) = currMeasSample;

        % Atualiza elevação dos satélites para o próximo ciclo
        satElev = navSolutions.el(:, currMeasNr)';

        % Calcula pseudodistância corrigida final
        navSolutions.correctedP(activeChnList, currMeasNr) = ...
            navSolutions.rawP(activeChnList, currMeasNr) + ...
            satClkCorr' * settings.c - xyzdt(4);

        %% Conversão de coordenadas ECEF para latitude/longitude/altura
        [navSolutions.latitude(currMeasNr), ...
         navSolutions.longitude(currMeasNr), ...
         navSolutions.height(currMeasNr)] = ...
            cart2geo(xyzdt(1), xyzdt(2), xyzdt(3), 5);  % 5 iterações para convergência

        %% Determina a zona UTM correspondente
        navSolutions_utmZone = findUtmZone(...
            navSolutions.latitude(currMeasNr), ...
            navSolutions.longitude(currMeasNr));

        %% Conversão para coordenadas locais UTM (ENU)
        [navSolutions.E(currMeasNr), ...
         navSolutions.N(currMeasNr), ...
         navSolutions.U(currMeasNr)] = ...
            cart2utm(xyzdt(1), xyzdt(2), xyzdt(3), navSolutions_utmZone);
        m2lat = 1/110734;
        m2lon = 1/103043;
        navSolutions.LLH_error(currMeasNr,1)=...
            norm(([navSolutions.latitude(currMeasNr),...
            navSolutions.longitude(currMeasNr)]...
            -settings.gt_llh(1:2))./[m2lat m2lon]);
        LS_error(currMeasNr,1)=...
            norm(([navSolutions.latitude(currMeasNr),...
            navSolutions.longitude(currMeasNr)]...
            -settings.gt_llh(1:2))./[m2lat m2lon]);

    else
        %% Caso não haja satélites suficientes para solução
        warning(['   Epoch ', num2str(currMeasNr), ...
                 ': Not enough satellites for position fix.']);
        % Marca os valores como inválidos (NaN ou 0)
        navSolutions.X(currMeasNr)         = NaN;
        navSolutions.Y(currMeasNr)         = NaN;
        navSolutions.Z(currMeasNr)         = NaN;
        navSolutions.dt(currMeasNr)        = NaN;
        navSolutions.DOP(:, currMeasNr)    = zeros(5, 1);
        navSolutions.latitude(currMeasNr)  = NaN;
        navSolutions.longitude(currMeasNr) = NaN;
        navSolutions.height(currMeasNr)    = NaN;
        navSolutions.E(currMeasNr)         = NaN;
        navSolutions.N(currMeasNr)         = NaN;
        navSolutions.U(currMeasNr)         = NaN;
        navSolutions.az(activeChnList, currMeasNr) = NaN;
        navSolutions.el(activeChnList, currMeasNr) = NaN;
    end
    % [navSolutions] = DPE_module...
    %     (currMeasNr,navSolutions,activeChnList,...
    %     trackResults,currMeasSample,satPositions,...
    %     transmitTime(activeChnList),localTime,...
    %     settings,satElev,fid,xyzdt(4),satClkCorr);
    % return
    %% Impressão do erro de posição, se existir
    if isfield(navSolutions, 'LLH_error')
        % fprintf('Current 2D Error of LS: %.6f m\n', navSolutions.LLH_error(currMeasNr, 1));
        fprintf('Current 2D Error of LS: %.6f m\n', LS_error(currMeasNr, 1));
    end
end

% return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ls_2d_list_new = navSolutions.LLH_error(:,1)';
ls_2d_list_new = LS_error(:,1)';
navSolPeriod = 500;

% === Plot 2D positioning error ===========================================

figure;

plot(ls_2d_list_new, 'Color', [80, 29, 138]/255, 'LineWidth', 3);
hold on;

legend('Least Squares (Scalar Tracking)', 'FontSize', 12);
xlabel(sprintf('Epoch (%d ms)', navSolPeriod), 'FontSize', 12);
ylabel('2D error (m)', 'FontSize', 12);
xlim([1, length(ls_2d_list_new)]);
title('2D Positioning Error');
grid on;
set(gca, 'FontSize', 12);

return