function [startSamples, endSamples] = getSampleWindow(trackResults, subFrameStart, activeChnList)
    nChannels = length(trackResults);
    startSamples = zeros(1, nChannels);
    endSamples = inf(1, nChannels);
    
    for ch = activeChnList
        startSamples(ch) = trackResults(ch).absoluteSample(subFrameStart(ch));
        endSamples(ch) = trackResults(ch).absoluteSample(end);
    end
end
