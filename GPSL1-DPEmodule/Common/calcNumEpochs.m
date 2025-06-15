function N = calcNumEpochs(sampleStart, sampleEnd, sampleStep)
    % Calcula número total de épocas de navegação possíveis
    N = fix((sampleEnd - sampleStart) / sampleStep);
end