function data = prepareGNSSProcessing(show_fig)
% PREPAREGNSSPROCESSING Prepara dados GNSS brutos para processamento
%
%   data = prepareGNSSProcessing(show_fig)
%
%   Input:
%       show_fig - Flag para exibir gráficos diagnósticos 
%       (1 = exibir, 0 = não exibir)
%
%   Output:
%       data - Vetor contendo os dados brutos lidos do arquivo
%
%   Esta função realiza:
%   1. Carrega configurações básicas usando initSettings()
%   2. Abre e lê o arquivo de dados GNSS
%   3. Calcula parâmetros derivados importantes
%   4. Opcionalmente exibe gráficos diagnósticos dos dados

% Inicializa configurações padrão do processamento
settings = initSettings();

% 1. Abertura do arquivo e verificação inicial
fileName = settings.fileName;
fprintf('Preparando processamento para arquivo: %s\n', fileName);

% Tenta abrir o arquivo no modo binário para leitura
[fid, message] = fopen(fileName, 'rb');

% Verifica se o arquivo foi aberto com sucesso
if fid <= 0
    error('Não foi possível abrir o arquivo %s: %s', fileName, message);
end

% 2. Cálculo de parâmetros derivados importantes
% Calcula o número de amostras por código C/A (Code Acquisition)
samplesPerCode = round(settings.samplingFreq / ...
                      (settings.codeFreqBasis / settings.codeLength));

% 3. Leitura dos dados brutos
% Posiciona o ponteiro no início do arquivo (com possível offset)
fseek(fid, settings.skipNumberOfBytes, 'bof');

% Lê 10ms de dados (equivalente a 10 períodos do código C/A)
[data, count] = fread(fid, [1, 10*samplesPerCode], settings.dataType);

% Fecha o arquivo após a leitura
fclose(fid);

% Verifica se a quantidade de dados lida é suficiente
if count < 10*samplesPerCode
    warning('Arquivo contém menos dados que o esperado.');
end

% Atualiza a estrutura de configurações com parâmetros calculados
settings.samplesPerCode = samplesPerCode;
settings.fileName = fileName;

% 4. Exibição dos gráficos diagnósticos (se solicitado)
if show_fig == 1
    showDataDiagnostics(data, settings, samplesPerCode);
end

end

%% Função auxiliar para visualização diagnóstica
function showDataDiagnostics(data, settings, samplesPerCode)
% SHOWDATADIAGNOSTICS Exibe gráficos diagnósticos dos dados GNSS
%
%   Inputs:
%       data - Dados brutos a serem analisados
%       settings - Estrutura com parâmetros de configuração
%       samplesPerCode - Número de amostras por código C/A

% Cria figura para os gráficos diagnósticos
figure('Name', 'Diagnóstico dos Dados GNSS', 'NumberTitle', 'off');

% Subplot 1: Domínio do tempo (primeiros 5ms)
subplot(2, 2, 1);
% Cria escala de tempo para os primeiros 5ms
timeScale = 0 : 1/settings.samplingFreq : 5e-3;
% Calcula número de amostras para plotar (1/20 do código C/A ou todos os 
% dados disponíveis)
samplesToPlot = min(round(samplesPerCode/20), length(data));
% Plota o sinal no domínio do tempo
plot(1000 * timeScale(1:samplesToPlot), data(1:samplesToPlot));
axis tight; grid on;
title('Sinal no Domínio do Tempo');
xlabel('Tempo (ms)'); ylabel('Amplitude');

% Subplot 2: Domínio da frequência
subplot(2, 2, 2);
% Calcula e plota a densidade espectral de potência usando Welch's method
pwelch(data-mean(data), 16384, 1024, 2048, settings.samplingFreq/1e6);
axis tight; grid on;
title('Espectro de Frequência');
xlabel('Frequência (MHz)'); ylabel('Magnitude (dB)');

% Subplot 3: Histograma da distribuição de amplitudes
subplot(2, 2, 3.5);
% Plota histograma com bins específicos para dados de 8 bits (-128 a 127)
histogram(data, -128:128);
% Ajusta os limites do eixo X para cobrir todo o range dos dados
dmax = max(abs(data)) + 1;
axis tight;
adata = axis;
axis([-dmax dmax adata(3) adata(4)]);
grid on;
title('Histograma de Amplitudes'); 
xlabel('Valor Digital'); ylabel('Número de Ocorrências');
end