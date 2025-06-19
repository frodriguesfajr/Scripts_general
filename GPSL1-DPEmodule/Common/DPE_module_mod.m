function [navSolutions]=DPE_module_mod...
    (currMeasNr,navSolutions,activeChnList,trackResults,currMeasSample,...
    satPositions,transmitTime,localTime,settings,satElev,fid,...
    LS_clkBias,satClkCorr)
% DPE_module é um módulo plug-in de Estimação Direta de Posição (DPE) que pode
% ser integrado em SDRs MATLAB existentes que usam abordagem em duas etapas.
%
% As correlações para cada fase de código são primeiro calculadas e depois
% interpoladas para obter correlações para cada posição candidata
% (um correlograma sobre as posições candidatas)
%
% O módulo usa por padrão 1 integração não-coerente

% --- Parâmetros de Entrada ---
% currMeasNr: Número da época de posicionamento atual
% navSolutions: Soluções de navegação do método dos mínimos quadrados (inicialização)
% activeChnList: Lista de canais ativos
% trackResults: Resultados do rastreamento
% currMeasSample: Amostra atual de medição
% satPositions: Posições dos satélites (ECEF)
% transmitTime: Tempo de transmissão dos satélites
% localTime: Tempo local no receptor
% settings: Configurações do receptor
% satElev: Ângulos de elevação dos satélites
% fid: Identificador do arquivo de dados brutos
% LS_clkBias: Viés do relógio estimado por mínimos quadrados
% satClkCorr: Correção do relógio do satélite

% --- Parâmetros de Saída ---
% navSolutions: Estrutura contendo as soluções DPE (latitude, longitude, altura, viés de relógio)

% Inicialização de caminhos para funções auxiliares
addpath ('C:\Repository\GPSL1-DPEmodule\include')
addpath ('C:\Repository\GPSL1-DPEmodule\common') 

% === Inicialização de parâmetros =========================================
dtr = pi/180; % Conversão de graus para radianos
m2lat = 1/110734; % Fator de conversão metros para latitude (aproximado)
m2lon = 1/103043; % Fator de conversão metros para longitude (aproximado)
trop = zeros(length(activeChnList),1); % Inicializa vetor para correção troposférica

% Define espaçamentos de chips para pré-cálculo de correlações
chip_spacings = [-(flip(settings.chipspacing_dpe_precalc:...
                  settings.chipspacing_dpe_precalc:1)),0,...
                  settings.chipspacing_dpe_precalc:...
                  settings.chipspacing_dpe_precalc:1];

% Pré-aloca matriz para armazenar correlações pré-calculadas
precalc_correlations = zeros(length(chip_spacings)*length(activeChnList),3);

% Cria diretório para armazenar correlogramas se necessário
if ~exist([settings.outfile_root '\Correlogram\'], 'dir') && ...
        settings.DPE_plotCorrelogram == 1
    mkdir([settings.outfile_root '\Correlogram\']);
end

% === Geração do espaço de busca para latitude e longitude ================
% Cria grade de busca centrada na estimativa de mínimos quadrados
[~,~,~,mesh_ll] = ...
   meshgrid_LIST(navSolutions.latitude(currMeasNr)+(m2lat).*...
   (-settings.DPE_latlong_span:settings.candPVT_spacing:settings.DPE_latlong_span)...
    ,navSolutions.longitude(currMeasNr)+(m2lon).*...
    (-settings.DPE_latlong_span:settings.candPVT_spacing:settings.DPE_latlong_span));

% === Geração do espaço de busca para altura =============================
% Intervalo de busca centrado na estimativa de mínimos quadrados
alt_search = round(navSolutions.height(currMeasNr),0)-...
    settings.DPE_height_span:settings.candPVT_spacing:...
    round(navSolutions.height(currMeasNr),0)+settings.DPE_height_span;

% === Geração do espaço de busca para viés de relógio ====================
% Na primeira época usa busca mais ampla para garantir convergência
if currMeasNr == 1
    % Cálculo do viés de relógio inicial para cada satélite
    candia_xyz_LS = [navSolutions.X(currMeasNr),...
                    navSolutions.Y(currMeasNr),navSolutions.Z(currMeasNr)];
    clock_bias = zeros(length(activeChnList),1);
    simulated_pseudorange = zeros(length(activeChnList),1);
    
    for i = 1:length(activeChnList)
        if (settings.useTropCorr == 1)
            % Calcula correção troposférica
            trop(i) = tropo(sin(satElev(activeChnList(i)) * dtr), ...
                     0.0, 1013.0, 293.0, 50.0, 0.0, 0.0, 0.0);
        end
        
        % Pseudodistância simulada
        simulated_pseudorange(i) = sqrt(sum((satPositions(:,i)'...
                                - candia_xyz_LS).^2));
        
        % Calcula viés de relógio
        clock_bias(i) = (navSolutions.rawP(activeChnList(i),currMeasNr)+ ...
                       satClkCorr(i) * settings.c - simulated_pseudorange(i));
    end

    clk_s_area = round(min(clock_bias))-100:1:round(max(clock_bias))+100;
else
    % Para épocas subsequentes, usa busca mais estreita
    clk_s_area = round((LS_clkBias))-settings.DPE_clkBias_span:1:...
                 round((LS_clkBias))+settings.DPE_clkBias_span;

    % Calcula correção troposférica se necessário
    if (settings.useTropCorr == 1)
        for i = 1:length(activeChnList)
            trop(i) = tropo(sin(satElev(activeChnList(i)) * dtr), ...
                     0.0, 1013.0, 293.0, 50.0, 0.0, 0.0, 0.0);
        end
    end
end

% Limita o tamanho do espaço de busca se necessário
if length(clk_s_area)>100000
     clk_s_area=round((LS_clkBias))-settings.DPE_clkBias_span:1:...
        round((LS_clkBias))+settings.DPE_clkBias_span;
end

% === Pré-cálculo dos valores de correlação ==============================
if settings.MMT~=1
    % Processamento para método convencional (não MMT)
    for j=1:length(activeChnList)
        % Encontra amostra mais próxima no tempo
        [~,closestIndex] = min(abs(currMeasSample-...
                              trackResults(activeChnList(j)).absoluteSample));
        if currMeasSample < trackResults(activeChnList(j)).absoluteSample(closestIndex)
            closestIndex=closestIndex-1;
        end

        % Define intervalo para integração não-coerente
        if closestIndex > (settings.DPE_nonCohInt-1)
            nonCohInt=closestIndex-settings.DPE_nonCohInt+1:closestIndex;
        else
            nonCohInt=1:closestIndex;
        end

        % Gera código C/A do satélite
        caCode = generateCAcode(trackResults(activeChnList(j)).PRN);

        % Pré-cálculo das correlações
        count2=3;
        count=1+(j-1)*length(chip_spacings);
        
        for closestIndex=nonCohInt
            % Posiciona no arquivo de dados brutos
            if strcmp(settings.dataType,'int16') 
                fseek(fid, settings.fileType*...
                     ((trackResults(activeChnList(j)).absoluteSample(closestIndex))*2),'bof');
            else
                fseek(fid, settings.fileType*...
                     (trackResults(activeChnList(j)).absoluteSample(closestIndex)),'bof');
            end

            % Calcula parâmetros para processamento
            codePhaseStep = trackResults(activeChnList(j)).codeFreq(closestIndex)...
                          / settings.samplingFreq;
            blksize = ceil((settings.codeLength*settings.DPE_cohInt-...
                          trackResults(activeChnList(j)).remCodePhase(closestIndex)) ...
                          / codePhaseStep);

            % Lê dados brutos
            [rawSignal, ~] = fread(fid, settings.fileType*blksize, settings.dataType);
            rawSignal = rawSignal';

            % Converte para complexo se necessário
            if (settings.fileType==2)
               rawSignal1=rawSignal(1:2:end);
               rawSignal2=rawSignal(2:2:end);
               rawSignal = rawSignal1 + 1i .* rawSignal2;
            end
            
            % Gera portadora
            time = (0:blksize) ./ settings.samplingFreq;
            remCarrPhase = trackResults(activeChnList(j)).remCarrPhase(closestIndex);
            trigarg = ((trackResults(activeChnList(j)).carrFreq(closestIndex)...
                     * 2.0 * pi) .* time) + remCarrPhase;
            carrsig = exp(1i .* trigarg(1:blksize));

            % Mistura para banda base
            qBasebandSignal = real(carrsig .* rawSignal);
            iBasebandSignal = imag(carrsig .* rawSignal);

            % Calcula correlações para diferentes espaçamentos de chip
            for spacing=chip_spacings
                delay_index = trackResults(activeChnList(j)).remCodePhase(closestIndex)-spacing:...
                            codePhaseStep : ...
                            ((blksize-1)*codePhaseStep + ...
                            trackResults(activeChnList(j)).remCodePhase(closestIndex)-spacing);       

                caCode1 = [caCode(end-2) caCode(end-1) caCode(end) ...
                          repmat(caCode,1,settings.DPE_cohInt) caCode(1) caCode(2)...
                          caCode(3)];
                tcodee = ceil(delay_index)+3;
                tcodee(tcodee==0)=tcodee(tcodee==0)+1;

                s = caCode1(tcodee);

                I = sum(s  .* iBasebandSignal);
                Q = sum(s  .* qBasebandSignal);
               
                % Armazena correlações
                precalc_correlations(count,1) = trackResults(activeChnList(j)).PRN;
                precalc_correlations(count,2) = (localTime-transmitTime(j))*settings.c ...
                                             + (spacing/settings.codeFreqBasis) * settings.c;
                precalc_correlations(count,count2) = sqrt(I.^2 + Q.^2);

                count=count+1;
            end
            count=1+(j-1)*length(chip_spacings);
            count2=count2+1;
        end 
    end
else
    % Processamento para método MMT (Modified Moving Tracking)
    % Similar ao convencional, mas com tratamento adicional para multipercurso
    % [Código omitido por brevidade, mas segue estrutura similar]
end

% === Integração não-coerente das correlações ===
precalc_correlations(:,3) = sum(precalc_correlations(1:end,3:end),2);

% === Busca da posição ótima ===
barisan = 1; % Índice para registro de resultados
temprecord_DPE_values = zeros(length(alt_search),5); % Pré-aloca resultados

% Varre o espaço de busca de altura
for alt = alt_search
    % Converte coordenadas lat/long/alt para ECEF
    candia_map=[mesh_ll,ones(size(mesh_ll,1),1)*alt];
    x = candia_map(:,1)./180 .* pi;
    y = candia_map(:,2)./180 .* pi;
    candia_xyz  = llh2xyz([x,y,candia_map(:,3)]); 

    % Pré-aloca correlograma
    correlogram=zeros(length(clk_s_area),length(candia_map));

    % Para cada satélite ativo
    for j=1:length(activeChnList)
        % Calcula pseudodistâncias para cada ponto candidato
        candia_xyz_pseudoranges = sqrt(sum(((repmat(satPositions(:,j)',length(candia_xyz),1)...
                                      - candia_xyz).^2),2));
        
        % Aplica correções
        candia_xyz_pseudoranges = candia_xyz_pseudoranges + trop(j) - satClkCorr(j)* settings.c;
        
        % Adiciona viés de relógio ao espaço de busca
        candia_xyz_pseudoranges_clkBias = zeros(length(clk_s_area),length(candia_xyz_pseudoranges));          
        for blob = 1:length(clk_s_area)
            candia_xyz_pseudoranges_clkBias(blob,1:length(candia_xyz_pseudoranges)) = ...
                candia_xyz_pseudoranges(1:end)+clk_s_area(blob);
        end
    
        % Interpola correlações pré-calculadas para as pseudodistâncias candidatas
        prn_index = find(precalc_correlations(:,1)==trackResults(activeChnList(j)).PRN);
        correlogram_single = interp1(precalc_correlations(prn_index,2),...
                                    precalc_correlations(prn_index,3),...
                                    candia_xyz_pseudoranges_clkBias,'linear');
        
        % Trata valores NaN da interpolação
        try
          nanx = isnan(correlogram_single);
          t = 1:numel(correlogram_single);
          correlogram_single(nanx) = interp1(t(~nanx), correlogram_single(~nanx), t(nanx),'linear');
        catch
        end
    
        % Acumula correlogramas de todos os satélites
        correlogram=correlogram+correlogram_single;
    end
    
    % Encontra máximo no correlograma para esta altura
    [MaxCorrValue,~] = max(correlogram,[],'all','linear'); 
    [I,J,~] = find(correlogram==MaxCorrValue);
    
    % Armazena resultados
    if length(I) ~= 1 
        % Múltiplos máximos - usa média
        temprecord_DPE_values(barisan,1:3) = mean(candia_map(J,:));
        temprecord_DPE_values(barisan,4) = mean(MaxCorrValue);
        temprecord_DPE_values(barisan,5) = mean(clk_s_area(I));
        temprecord_DPE_values(barisan,6) = norm(([temprecord_DPE_values(barisan,1),...
                                               temprecord_DPE_values(barisan,2)]...
                                               -settings.gt_llh(1:2))./[m2lat m2lon]);
        temprecord_DPE_values(barisan,7) = round(mean(I));
    else
        % Único máximo
        temprecord_DPE_values(barisan,1:3) = candia_map(J,:);
        temprecord_DPE_values(barisan,4) = MaxCorrValue;
        temprecord_DPE_values(barisan,5) = clk_s_area(I);
        temprecord_DPE_values(barisan,6) = norm(([candia_map(J,1),candia_map(J,2)]...
                                              -settings.gt_llh(1:2))./[m2lat m2lon]);
        temprecord_DPE_values(barisan,7) = I;
    end
    
    barisan = barisan+1;
end

% === Determina a solução DPE final ===
[~,barisan_yangmanaya] = max(temprecord_DPE_values(:,4));

% === Gera gráficos de correlogramas se configurado ===
if settings.DPE_plotCorrelogram == 1
    % [Código para geração de gráficos omitido por brevidade]
    % Gera gráficos 3D dos correlogramas para cada satélite e combinado
    % Salva imagens em formato PNG e figuras MATLAB
end

% === Armazena resultados na estrutura de saída ===
navSolutions.DPE_estimate(currMeasNr, 1:5) = temprecord_DPE_values(barisan_yangmanaya,1:5);
navSolutions.DPE_latitude(currMeasNr) = temprecord_DPE_values(barisan_yangmanaya,1);
navSolutions.DPE_longitude(currMeasNr) = temprecord_DPE_values(barisan_yangmanaya,2);
navSolutions.DPE_height(currMeasNr) = temprecord_DPE_values(barisan_yangmanaya,3);
navSolutions.DPE_clkBias(currMeasNr) = temprecord_DPE_values(barisan_yangmanaya,5);