function [acqResults, data] = acquireSatellites(settings, fid, show_fig)
% ACQUIRESATELLITES Realiza a aquisição de sinais de satélites GNSS
%
%   [acqResults, data] = acquireSatellites(settings, fid)
%
%   Inputs:
%       settings    - Estrutura com parâmetros de configuração
%       fid         - Identificador do arquivo de dados (opcional)
%       show_fig    - Flag para exibir gráficos diagnósticos 
%       (1 = exibir, 0 = não exibir)
%
%   Outputs:
%       acqResults  - Estrutura com resultados da aquisição
%       data        - Dados brutos usados na aquisição
%
%   Esta função realiza:
%   1. Configuração inicial baseada nas settings
%   2. Leitura dos dados brutos do arquivo
%   3. Processamento de aquisição para todos os satélites visíveis
%   4. Retorno dos resultados estruturados

%% Inicialização e verificação de parâmetros
disp('Iniciando processo de aquisição de satélites...');

% Verifica se precisa abrir o arquivo
if nargin < 2 || isempty(fid)
    % Abre o arquivo se não foi fornecido
    [fid, message] = fopen(settings.fileName, 'rb');
    if fid <= 0
        error('Não foi possível abrir o arquivo %s: %s', settings.fileName, message);
    end
    closeFile = true; % Flag para fechar o arquivo após uso
else
    closeFile = false;
end

% Configura coeficiente de adaptação para tipo de dados
if (settings.fileType == 1) 
    dataAdaptCoeff = 1; % Dados reais
else
    dataAdaptCoeff = 2; % Dados complexos (I/Q)
end

%% Preparação dos dados para aquisição
% Calcula amostras por código C/A
samplesPerCode = round(settings.samplingFreq / ...
                 (settings.codeFreqBasis / settings.codeLength));

% Posiciona no arquivo considerando offset e tipo de dados
fseek(fid, dataAdaptCoeff*settings.skipNumberOfBytes, 'bof');

% Lê 11ms de dados (necessários para estimativa fina de frequência)
data = fread(fid, dataAdaptCoeff*11*samplesPerCode, settings.dataType)';

% Fecha o arquivo se foi aberto nesta função
if closeFile
    fclose(fid);
end

% Converte dados complexos se necessário
if (dataAdaptCoeff == 2)
    data = data(1:2:end) + 1i * data(2:2:end);
end

%% Processamento de Aquisição
disp('   Realizando aquisição de satélites...');

% Chama função principal de aquisição
acqResults = acquisition(data, settings);

% % Exibe resultados (opcional)
% if settings.showAcquisitionResults
%     plotAcquisition(acqResults);
% end

if show_fig == 1
    plotAcquisition(acqResults);
end

% %% Verificação dos resultados
% detectedSatellites = find(acqResults.carrFreq > 0);
% 
% if isempty(detectedSatellites)
%     warning('Nenhum satélite foi detectado na aquisição!');
% else
%     fprintf('   Satélites detectados: %s\n', num2str(detectedSatellites'));
% end

end