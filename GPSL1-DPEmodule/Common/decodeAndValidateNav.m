function [success, eph, subFrameStart, TOW] = decodeAndValidateNav(track, settings)
    success = false;
    eph = struct(); subFrameStart = inf; TOW = inf;

    try
        [eph, subFrameStart, TOW] = NAVdecoding(track.I_P, settings);

        % Validação de campos essenciais das efemérides
        if isempty(eph.IODC) || isempty(eph.IODE_sf2) || isempty(eph.IODE_sf3)
            fprintf('  [FALHA] Ephemeris incompleta para PRN %02d.\n', track.PRN);
            return;
        end

        success = true;

    catch ME
        fprintf('  [ERRO] Falha ao decodificar/validar NAV do PRN %02d: %s\n', ...
                track.PRN, ME.message);
    end
end
