function step = calcSampleStep(settings)
    % Calcula o passo de amostragem (samples entre soluções de navegação)
    step = fix(settings.samplingFreq * settings.navSolPeriod / 1000);
end