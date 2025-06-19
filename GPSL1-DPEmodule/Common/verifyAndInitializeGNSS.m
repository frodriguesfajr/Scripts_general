function [navSolutions, eph, activeChnList] = verifyAndInitializeGNSS(settings, trackResults)
% VERIFYANDINITIALIZEGNSS Verifica condições iniciais e inicializa estruturas para processamento GNSS
%
%   [navSolutions, eph, activeChnList] = verifyAndInitializeGNSS(settings, trackResults)
%
%   Inputs:
%       settings    - Estrutura com parâmetros de configuração
%       trackResults - Resultados do rastreamento dos satélites
%
%   Outputs:
%       navSolutions - Estrutura para soluções de navegação (inicializada vazia)
%       eph          - Estrutura para efemérides (inicializada)
%       activeChnList - Lista de canais ativos
%
%   Esta função:
%   1. Verifica se a duração do sinal é suficiente
%   2. Inicializa vetores auxiliares
%   3. Identifica canais ativos
%   4. Prepara estrutura para armazenar efemérides

%% 1. Verificação da duração mínima do sinal
if settings.msToProcess < 36000 % Menos de 36 segundos (36000 ms)
    disp('Registro de dados muito curto. Processamento abortado!');
    navSolutions = [];
    eph = [];
    activeChnList = [];
    return;
end

%% 2. Inicialização de vetores auxiliares
% Cria vetores infinitos para marcação de subframes e TOW (Time Of Week)
[subFrameStart, TOW] = deal(inf(1, settings.numberOfChannels));

%% 3. Identificação de canais ativos
% Encontra todos os canais com status diferente de '-' (inativo)
activeChnList = find([trackResults.status] ~= '-');

%% 4. Inicialização da estrutura de efemérides
% Lista de campos necessários para armazenar parâmetros orbitais
ephFields = {
    'C_ic', 'omega_0', 'C_is', 'i_0', 'C_rc', 'omega', 'omegaDot',...
    'IODE_sf3', 'iDot', 'idValid', 'weekNumber', 'accuracy',...
    'health', 'T_GD', 'IODC', 't_oc', 'a_f2', 'a_f1', 'a_f0',...
    'IODE_sf2', 'C_rs', 'deltan', 'M_0', 'C_uc', 'e', 'C_us',...
    'sqrtA', 't_oe', 'TOW'
};

% Cria célula de valores iniciais (zeros para numéricos, false para idValid)
initValues = num2cell(zeros(1, numel(ephFields)));
initValues{strcmp(ephFields, 'idValid')} = false;

% Converte para estrutura
eph = cell2struct(initValues, ephFields, 2);

%% 5. Inicialização da estrutura de soluções de navegação
navSolutions = struct();
end