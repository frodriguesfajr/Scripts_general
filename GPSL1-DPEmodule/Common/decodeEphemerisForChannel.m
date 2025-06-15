function [eph_ch, subFrameStart, TOW] = ...
    decodeEphemerisForChannel(trackResults, channelNr, settings, activeChnList)
% Função para decodificar e validar efemérides para um único canal
%
% Entradas:
%   - trackResults: struct com resultados de rastreamento
%   - channelNr: índice do canal a ser processado
%   - settings: configurações do sistema GNSS
%   - eph_ch: struct de efemérides por PRN (entrada e saída)
%   - subFrameStart: vetor com inícios de subquadros
%   - TOW: vetor com Time of Week por canal
%   - activeChnList: lista atual de canais ativos
%
% Saídas:
%   - eph_ch: struct atualizado com efemérides decodificadas
%   - subFrameStart: vetor atualizado
%   - TOW: vetor atualizado
%   - activeChnList: lista atualizada com canais válidos

    PRN = trackResults(channelNr).PRN;

    fprintf('Decoding NAV for PRN %02d ----------------------------\n', PRN);

    %% ======================= NAV Decoding =============================
    try
        [eph_ch(PRN), subFrameStart(channelNr), TOW(channelNr)] = ...
            NAVdecoding(trackResults(channelNr).I_P, settings);
    catch ME
        activeChnList = setdiff(activeChnList, channelNr);
        fprintf('  [ERRO] Falha ao decodificar NAV do PRN %02d: %s\n', PRN, ME.message);
        return;
    end

    %% =================== Validação das Efemérides =====================
    try
        if isempty(eph_ch(PRN).IODC) || ...
           isempty(eph_ch(PRN).IODE_sf2) || ...
           isempty(eph_ch(PRN).IODE_sf3)

            activeChnList = setdiff(activeChnList, channelNr);
            fprintf('  [FALHA] Ephemeris incompleta para PRN %02d.\n', PRN);
        else
            fprintf('  [OK] Ephemeris completa para PRN %02d!\n', PRN);
        end
    catch ME
        activeChnList = setdiff(activeChnList, channelNr);
        fprintf('  [ERRO] Validação das efemérides falhou para PRN %02d: %s\n', PRN, ME.message);
    end
end
