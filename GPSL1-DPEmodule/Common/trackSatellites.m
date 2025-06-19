function [trackResults, channel, outfile_root] = ...
    trackSatellites(fileName, acqResults, fid, settings, show_fig)
% Função principal para rastreamento de satélites GNSS

%       show_fig    - Flag para exibir gráficos diagnósticos 
%       (1 = exibir, 0 = não exibir)

% Inicializa variáveis de saída
trackResults = [];

% Verifica se existem sinais adquiridos para rastrear
if (any(acqResults.carrFreq))
    % Pré-configura os canais com base nos resultados de aquisição
    channel = preRun(acqResults, settings);
    % Opção para mostrar status dos canais (comentado)
    % showChannelStatus(channel, settings);
else
    % Caso não haja sinais detectados, exibe mensagem e retorna
    disp('No GNSS signals detected, signal processing finished.');
    return;
end

if show_fig == 1
    showChannelStatus(channel, settings)
end
%% Cria diretório para armazenar os arquivos de saída
% Extrai o nome base do arquivo de dados
outfile_root = strfind(fileName,'\');
outfile_root = fileName(outfile_root(end)+1:end-4); 

% Verifica se o diretório já existe, caso contrário cria
outfile_root = ['.\GPSL1-DPEmodule\' outfile_root];
if exist(outfile_root, 'dir') ~= 7
    mkdir(outfile_root);
end

% Armazena o caminho raiz nas configurações

settings.outfile_root = outfile_root;   
%% Processo de rastreamento principal
% Verifica o modo de operação (normal ou MMT - Modified Moving Tracking)
if settings.MMT ~= 1
    % Modo de rastreamento convencional
    
    % Define o nome do arquivo de resultados com parâmetros específicos
    trackFile = [settings.outfile_root, '\trackResults_gps_', ...
                settings.outfile_root, '_', ...
                num2str(settings.dllCorrelatorSpacing), '_', ...
                num2str(settings.DPE_cohInt), '.mat'];
    
    % Verifica se o arquivo de resultados já existe
    if ~exist(trackFile, 'file')
        % Caso não exista, inicia novo processo de rastreamento
        
        trackFile = strrep(trackFile, '.\GPSL1-DPEmodule', '');
        % return
        % Registra hora de início
        startTime = now;
        disp(['   Tracking started at ', datestr(startTime)]);
        
        % Executa o algoritmo de rastreamento convencional
        [trackResults, channel] = tracking(fid, channel, settings);
        
        % Exibe tempo decorrido
        disp(['   Tracking is over (elapsed time ', ...
             datestr(now - startTime, 13), ')']);
        
        % Salva os resultados no arquivo
        disp(['   Saving Acq & Tracking results to file "trackingResults.mat"']);
        save(trackFile, 'trackResults', 'settings', 'acqResults', 'channel');
    else
        % Caso o arquivo já exista, carrega os resultados prévios
        
        % Preserva as novas configurações
        new_settings = settings;
        disp(trackFile)
        % Carrega resultados anteriores
        load(trackFile);
        
        % Restaura as novas configurações
        settings = new_settings;
        
        % Compatibilidade com versões anteriores
        if ~exist('trackResults','var')
            trackResults = trackResults_gps;
        end
        % disp(trackFile)
        % disp(outfile_root)
        % disp(pwd)
        % Caminho completo usando fullfile (compatível com qualquer SO)
        [~, filename1, ext] = fileparts(trackFile);
        filename1 = [filename1 ext];
        trackFile = fullfile(pwd, outfile_root, filename1);
        % disp(trackFile)
        % if ~exist(trackFile, 'file')
        %     save(trackFile, 'trackResults', 'settings',...
        %         'acqResults', 'channel');
        % end
        % 
        % return
        % % Salva com as configurações atualizadas
        % save(trackFile, 'trackResults', 'settings', 'acqResults', 'channel');
    end
else
    % Modo MMT (Modified Moving Tracking)
    
    % Define nome do arquivo para modo MMT
    trackFile = [settings.outfile_root, '\trackResults_gps_MMT_', ...
                settings.outfile_root, '_', ...
                num2str(settings.DPE_cohInt), '.mat'];
    
    % Verifica se o arquivo já existe
    if ~exist(trackFile, 'file')
        % Novo rastreamento MMT
        
        % Registra hora de início
        startTime = now;
        disp(['   Tracking started at ', datestr(startTime)]);
        
        % Executa rastreamento MMT
        [trackResults, channel] = tracking_MMT(fid, channel, settings);
        
        % Exibe tempo decorrido
        disp(['   Tracking is over (elapsed time ', ...
             datestr(now - startTime, 13), ')']);
        
        % Salva resultados
        disp(['   Saving Acq & Tracking results to file "trackingResults.mat"']);
        save(trackFile, 'trackResults', 'settings', 'acqResults', 'channel');
    else
        % Carrega resultados MMT existentes
        
        % Preserva novas configurações
        new_settings = settings;
        
        % Carrega resultados anteriores
        load(trackFile);
        
        % Restaura novas configurações
        settings = new_settings;
        
        % Compatibilidade com versões anteriores
        if ~exist('trackResults','var')
            trackResults = trackResults_gps;
        end
        
        % Salva com configurações atualizadas
        save(trackFile, 'trackResults', 'settings', 'acqResults', 'channel');
    end
end