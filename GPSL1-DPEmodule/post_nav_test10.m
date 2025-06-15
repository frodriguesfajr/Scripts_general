%% Limpeza do ambiente
clear;           % Remove todas as variáveis da memória
clc;             % Limpa a janela de comandos
%% Initialize constants, settings =========================================
settings = initSettings();
%% Generate plot of raw data and ask if ready to start processing =========
fprintf('Probing data (%s)...\n', settings.fileName)
%% Check the number of arguments ==========================================
fileNameStr = settings.fileName;
[fid, message] = fopen(fileNameStr, 'rb');
if (fid > 0)
    % Move the starting point of processing. Can be used to start 
    % the signal processing at any point in the data record (e.g. for long
    % records).
    fseek(fid, settings.skipNumberOfBytes, 'bof');    
    % Find number of samples per spreading code
    samplesPerCode = round(settings.samplingFreq / ...
        (settings.codeFreqBasis / settings.codeLength));
    % Read 10ms of signal
    [data, count] = fread(fid, [1, 10*samplesPerCode], settings.dataType);
    fclose(fid);
    if (count < 10*samplesPerCode)
        % The file is to short
        error('Could not read enough data from the data file.');
    end
    % %--- Initialization -------------------------------------------------
    % figure(100);
    % clf(100);
    % timeScale = 0 : 1/settings.samplingFreq : 5e-3;    
    % %--- Time domain plot -----------------------------------------------
    % subplot(2, 2, 1);
    % plot(1000 * timeScale(1:round(samplesPerCode/50)), ...
    %     data(1:round(samplesPerCode/50)));
    % axis tight;
    % grid on;
    % title ('Time domain plot');
    % xlabel('Time (ms)'); ylabel('Amplitude');    
    % %--- Frequency domain plot ------------------------------------------
    % subplot(2,2,2);
    % pwelch(data-mean(data), 16384, 1024, 2048, settings.samplingFreq/1e6)
    % axis tight;
    % grid on;
    % title ('Frequency domain plot');
    % xlabel('Frequency (MHz)'); ylabel('Magnitude');    
    % %--- Histogram ------------------------------------------------------
    % subplot(2, 2, 3.5);
    % histogram(data, -128:128)
    % dmax = max(abs(data)) + 1;
    % axis tight;
    % adata = axis;
    % axis([-dmax dmax adata(3) adata(4)]);
    % grid on;
    % title ('Histogram'); 
    % xlabel('Bin'); ylabel('Number in bin');
else
    %=== Error while opening the data file ================================
    error('Unable to read file %s: %s.', fileNameStr, message);
end % if (fid > 0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 1.1) Open the data file for the processing and seek to desired point.
% %
% % 2.1) Acquire satellites
% %
% % 3.1) Initialize channels (preRun.m).
% % 3.2) Pass the channel structure and the file identifier to the tracking
% % function. It will read and process the data. The tracking results are
% % stored in the trackResults structure. The results can be accessed this
% % way (the results are stored each millisecond):
% % trackResults(channelNumber).XXX(fromMillisecond : toMillisecond), where
% % XXX is a field name of the result (e.g. I_P, codePhase etc.)
% %
% % 4) Pass tracking results to the navigation solution function. It will
% % decode navigation messages, find satellite positions, measure
% % pseudoranges and find receiver position.
% %
% % 5) Plot the results.
% 
%% Initialization =========================================================
addpath ('.\include')  
addpath ('.\common') 
disp ('Starting processing...'); 
settings = initSettings();
[fid, message] = fopen(settings.fileName, 'rb');
%Initialize the multiplier to adjust for the data type
if (settings.fileType==1) 
    dataAdaptCoeff=1;
else
    dataAdaptCoeff=2;
end 
if (fid > 0)    
    % Move the starting point of processing. Can be used to start the
    % signal processing at any point in the data record (e.g. good for long
    % records or for signal processing in blocks).
    fseek(fid, dataAdaptCoeff*settings.skipNumberOfBytes, 'bof'); 
    %% Acquisition ========================================================
    % Do acquisition if it is not disabled in settings or if the variable
    % acqResults does not exist.
    if ((settings.skipAcquisition == 0) || ~exist('acqResults', 'var'))
        % Find number of samples per spreading code
        samplesPerCode = round(settings.samplingFreq / ...
            (settings.codeFreqBasis / settings.codeLength));
        % Read data for acquisition. 11ms of signal are needed for the fine
        % frequency estimation
        data  = fread(fid, dataAdaptCoeff*11*samplesPerCode, ...
            settings.dataType)'; 
        if (dataAdaptCoeff==2)    
            data1=data(1:2:end);    
            data2=data(2:2:end);    
            data=data1 + 1i .* data2;    
        end
        %--- Do the acquisition -------------------------------------------
        disp ('   Acquiring satellites...');
        acqResults = acquisition(data, settings);
        % plotAcquisition(acqResults);
    end 
    % if ((settings.skipAcquisition == 0) || ~exist('acqResults', 'var'))
    %% Initialize channels and prepare for the run ========================
    % Start further processing only if a GNSS signal was acquired (the
    % field FREQUENCY will be set to 0 for all not acquired signals)
    if (any(acqResults.carrFreq))
        channel = preRun(acqResults, settings);
        % showChannelStatus(channel, settings);
    else
        % No satellites to track, exit
        disp('No GNSS signals detected, signal processing finished.');
        trackResults = [];
        return;
    end
    %% Make folder for every IF file ======================================
    outfile_root=strfind(settings.fileName,'\');
    outfile_root=settings.fileName(outfile_root(end)+1:end-4); 
    if exist(outfile_root, 'dir') ~= 7
        mkdir(outfile_root);
    end
    settings.outfile_root=outfile_root; 
    %% Track the signal ===================================================
    % Checks for tracking file with specific DLL correlator spacing and PDI
    if settings.MMT ~= 1
        if ~exist([settings.outfile_root, '\trackResults_gps_', ...,
                settings.outfile_root...,
                '_',num2str(settings.dllCorrelatorSpacing),'_',...
                num2str(settings.DPE_cohInt),'.mat'])
            startTime = now;
            disp (['   Tracking started at ', datestr(startTime)]);
            % Process all channels for given data block
            [trackResults, channel] = tracking(fid, channel, settings);
            disp(['   Tracking is over (elapsed time ', ...
                datestr(now - startTime, 13), ')'])
            % Auto save the acquisition & tracking results to a file 
            % to allow running the positioning solution afterwards.
            disp(['   Saving Acq & Tracking results to ' ...
                'file "trackingResults.mat"'])
            % === Saves tracking result to folder =========================
            % tracking file with specific DLL correlator spacing and PDI
         save([settings.outfile_root, '\trackResults_gps_', ...
                    settings.outfile_root,'_',...
                    num2str(settings.dllCorrelatorSpacing),'_',...
                    num2str(settings.DPE_cohInt),'.mat'],...
                    'trackResults', 'settings', 'acqResults',...
                    'channel');
        else
            % In case of new settings, save new settings to tracking file
            new_settings = settings;
            % === Loads prev tracking result =============================
            load([settings.outfile_root, '\trackResults_gps_',...
                settings.outfile_root,'_',...
                num2str(settings.dllCorrelatorSpacing),'_',...
                num2str(settings.DPE_cohInt),'.mat']);
            settings=new_settings;
            if ~exist('trackResults','var')
                trackResults = trackResults_gps;
            end
            save([settings.outfile_root, '\trackResults_gps_', ...
                settings.outfile_root,'_',...
                num2str(settings.dllCorrelatorSpacing),'_',...
                num2str(settings.DPE_cohInt),'.mat'],...
                'trackResults', 'settings', 'acqResults',...
                'channel');
        end % if ~exist([settings.outfile_root, 
    else 
        % Proceed tracking with MMT
        if ~exist([settings.outfile_root, '\trackResults_gps_MMT_',...
                settings.outfile_root,'_',...
                num2str(settings.DPE_cohInt),'.mat'])
            startTime = now;
            disp (['   Tracking started at ', datestr(startTime)]);
            % Process all channels for given data block
            [trackResults, channel] = tracking_MMT(fid, channel, settings);
            disp(['   Tracking is over (elapsed time ', ...
                datestr(now - startTime, 13), ')'])
            % Auto save the acquisition & tracking results to a file to 
            % allow running the positioning solution afterwards.
            disp(['   Saving Acq & Tracking results to ' ...
                'file "trackingResults.mat"'])

            % === Saves tracking result to folder =========================
            save([settings.outfile_root, '\trackResults_gps_MMT_', ...
                settings.outfile_root,'_',...
                    num2str(settings.DPE_cohInt),'.mat'],...
                    'trackResults', 'settings', 'acqResults',...
                'channel');
        else
        % === In case of new settings, save new settings to tracking file =====
        new_settings = settings;
        % === Loads prev tracking result ======================================
        load([settings.outfile_root, '\trackResults_gps_MMT_',...
            settings.outfile_root,'_',...
                    num2str(settings.DPE_cohInt),'.mat']);
        settings=new_settings;
            if ~exist('trackResults','var')
                trackResults = trackResults_gps;
            end
        save([settings.outfile_root, '\trackResults_gps_MMT_', ...
            settings.outfile_root,'_',...
                    num2str(settings.DPE_cohInt),'.mat'],...
                    'trackResults', 'settings', 'acqResults',...
            'channel');

        end

    end
    %% Calculate navigation solutions =====================================
    disp('   Calculating navigation solutions...');
    % "fid" added as input for postNavigation for loading IF data for DPE
    % module
    % [navSolutions, eph] = postNavigation(trackResults, settings, fid);
    % 
    % disp('   Processing is complete for this data block');
    %% Plot all results =======================================================
    % disp ('   Ploting results...');
    % if settings.plotTracking
    % plotTracking(1:settings.numberOfChannels, trackResults, settings);
    % end
    % %     plotNavigation(navSolutions, settings);
    % %     disp('Post processing of the signal is over.');
else
    % Error while opening the data file.
    error('Unable to read file %s: %s.', settings.fileName, message);
end % if (fid > 0)

%% Calculate navigation solutions =======================================
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

%% Inicialização estrutura contendo as efemérides (parâmetros orbitais 
% Keplerianos e correções)
campos = {
    'C_ic', 'omega_0', 'C_is', 'i_0', 'C_rc', 'omega', 'omegaDot',...
    'IODE_sf3', 'iDot', 'idValid', 'weekNumber', 'accuracy',...,
    'health', 'T_GD', 'IODC', 't_oc', 'a_f2', 'a_f1', 'a_f0',...,
    'IODE_sf2', 'C_rs', 'deltan', 'M_0', 'C_uc', 'e', 'C_us',...,
    'sqrtA', 't_oe', 'TOW'
};

% Valores iniciais (todos zeros, exceto idValid = false)
valores = num2cell(zeros(1, numel(campos)));
valores{strcmp(campos, 'idValid')} = false;

% Criar struct
eph_ch = cell2struct(valores, campos, 2);

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
% measSampleStep = calcSampleStep(settings);
% measNrSum = calcNumEpochs(sampleStart, sampleEnd, measSampleStep);

% fprintf('Positions are being computed. Please wait...\n');
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
    navSolutions = [];
    navSolutions.latitude(currMeasNr) = lat;
    % 22.2999120726705
    navSolutions.longitude(currMeasNr) = long;
    % 114.179901819529
    navSolutions.height(currMeasNr) = height;
    % -37.346079277806
    navSolutions.X(currMeasNr) = xyzdt(1);
    % % navSolutions.X(currMeasNr) = -2418267.60022463
    navSolutions.Y(currMeasNr) = xyzdt(2);
    % 5385941.55281702
    navSolutions.Z(currMeasNr) = xyzdt(3);
    % 2405157.75721639
    navSolutions.rawP = rawP_arr;
    % navSolutions.rawP(activeChnList,currMeasNr)
    % navSolutions.rawP(activeChnList,currMeasNr)
    %       22499455.0384896
    %       23881607.0207043
    %       21626163.4200389
    %       20626320.6936134
    disp('ok48')
    % return
    tic

    [navSolutions] = DPE_module...
        (currMeasNr,navSolutions,activeChnList,...
        trackResults,currMeasSample,satPositions_ls,...
        tT_meas,localTime,...
        settings,satElev,fid,xyzdt(4),satClkCorr);
    navSolutions.DPE_processingtime(currMeasNr) = toc;
    %%%%%%%%%%%%%%%%%%%
    % localTime = localTime - xyzdt(4)/settings.c;       
    % navSolutions.localTime(currMeasNr) = localTime;
    % localTime = localTime + measSampleStep/settings.samplingFreq ;
    %%%%%%%%%%%%%%%%%%%%%%%
    navSolutions.LLH_error(currMeasNr,2)=...
    norm(([navSolutions.DPE_latitude(currMeasNr),...
        navSolutions.DPE_longitude(currMeasNr)]...
        -settings.gt_llh(1:2))./[m2lat m2lon]);    

    navSolutions.LLH_error(currMeasNr,1)=...
        norm(([navSolutions.latitude(currMeasNr),...
        navSolutions.longitude(currMeasNr)]...
        -settings.gt_llh(1:2))./[m2lat m2lon]);

    %% Impressão do erro de posição, se existir
    fprintf('\nCurrent 2D Error of DPE: %1f\n',...
        (navSolutions.LLH_error(currMeasNr,2)));
    fprintf('Current 2D Error of LS: %.6f m\n', LS_error(currMeasNr, 1));
    
    return
end

return
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
% disp('   Processing is complete for this data block');