function [satElev, readyChnList, localTime] = initPositionVars(settings, activeChnList)
    satElev = inf(1, settings.numberOfChannels);
    readyChnList = activeChnList;
    localTime = inf;
end