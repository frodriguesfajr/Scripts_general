clear;
clc;
% function [navSolutions, eph] = postNavigation(trackResults, settings, fid)
load('settings_dpe_sim.mat');  % Carrega a struct 'settings'
load('trackResults_dpe_sim.mat');  % Carrega a struct 'trackResults'
load('navSolutions_dpe_sim.mat');  % Carrega a struct 'navSolutions'
fid = 4;
addpath ('C:\Repository\GPSL1-DPEmodule\include')
addpath ('C:\Repository\GPSL1-DPEmodule\common') 
if (settings.msToProcess < 36000) 
    % Show the error message and exit
    disp('Record is to short. Exiting!');
    navSolutions = [];
    eph          = [];
    return
end
subFrameStart  = inf(1, settings.numberOfChannels);
TOW  = inf(1, settings.numberOfChannels);
activeChnList = find([trackResults.status] ~= '-');
for channelNr = activeChnList
    PRN = trackResults(channelNr).PRN;
    fprintf('Decoding NAV for PRN %02d -------------------- \n', PRN);
    try
    [eph_test(PRN), subFrameStart(channelNr), TOW(channelNr)] = ...
                                  NAVdecoding(trackResults(channelNr).I_P,...
                                  settings);  
    catch
    end

    try
    if (isempty(eph_test(PRN).IODC) || isempty(eph_test(PRN).IODE_sf2) || ...
        isempty(eph_test(PRN).IODE_sf3))
        activeChnList = setdiff(activeChnList, channelNr);
        fprintf('    Ephemeris decoding fails for PRN %02d !\n', PRN);
    else
        fprintf('    Three requisite messages for PRN %02d all decoded!\n', PRN);
    end
    catch
        activeChnList = setdiff(activeChnList, channelNr);
        fprintf('    Ephemeris decoding fails for PRN %02d !\n', PRN);
    end
end
if (isempty(activeChnList) || (size(activeChnList, 2) < 4))
    disp('Too few satellites with ephemeris data for postion calculations. Exiting!');
    navSolutions = [];
    eph          = [];
    return
end
sampleStart = zeros(1, settings.numberOfChannels);
sampleEnd = inf(1, settings.numberOfChannels);

for channelNr = activeChnList
    sampleStart(channelNr) = ...
          trackResults(channelNr).absoluteSample(subFrameStart(channelNr));
    
    sampleEnd(channelNr) = trackResults(channelNr).absoluteSample(end);
end
sampleStart = max(sampleStart) + 1;  
sampleEnd = min(sampleEnd) - 1;
measSampleStep = fix(settings.samplingFreq * settings.navSolPeriod/1000);
measNrSum = fix((sampleEnd-sampleStart)/measSampleStep);
satElev  = inf(1, settings.numberOfChannels);
readyChnList = activeChnList;
localTime = inf;
fprintf('Positions are being computed. Please wait... \n');
measNrSum = 4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
navSolutions_PRN = zeros(channelNr,measNrSum);
navSolutions_el = zeros(settings.numberOfChannels,measNrSum);
navSolutions_az = zeros(settings.numberOfChannels,measNrSum);
navSolutions_transmitTime = zeros(settings.numberOfChannels,measNrSum);
navSolutions_satClkCorr = zeros(settings.numberOfChannels,measNrSum);
navSolutions_rawP = zeros(settings.numberOfChannels,measNrSum);
navSolutions_DOP = zeros(channelNr+1,measNrSum);
navSolutions_X = zeros(1,measNrSum);
navSolutions_Y = zeros(1,measNrSum);
navSolutions_Z = zeros(1,measNrSum);
navSolutions_dt = zeros(1,measNrSum);
navSolutions_currMeasSample = zeros(1,measNrSum);
navSolutions_correctedP = zeros(channelNr,measNrSum);
navSolutions_latitude = zeros(1,measNrSum);
navSolutions_longitude = zeros(1,measNrSum);
navSolutions_height = zeros(1,measNrSum);
navSolutions_utmZone  = 0;
navSolutions_E = zeros(1,measNrSum);
navSolutions_N = zeros(1,measNrSum);
navSolutions_U = zeros(1,measNrSum);


for currMeasNr = 1:measNrSum
    fprintf('Fix: Processing %02d of %02d \n', currMeasNr,measNrSum);
    activeChnList = intersect(find(satElev >= settings.elevationMask), ...
                              readyChnList);
    navSolutions_PRN(activeChnList, currMeasNr) = [trackResults(activeChnList).PRN];
    navSolutions_el(:, currMeasNr) = NaN(settings.numberOfChannels, 1);
    navSolutions_az(:, currMeasNr) = NaN(settings.numberOfChannels, 1);
    navSolutions_transmitTime(:, currMeasNr) = NaN(settings.numberOfChannels, 1);
    navSolutions_satClkCorr(:, currMeasNr) = NaN(settings.numberOfChannels, 1);                                                                  
    currMeasSample = sampleStart + measSampleStep*(currMeasNr-1);
    [navSolutions_rawP(:, currMeasNr),transmitTime,localTime,codePhase]=  ...
                     calculatePseudoranges(trackResults,subFrameStart,TOW, ...
                     currMeasSample,localTime,activeChnList, settings);  
    navSolutions_transmitTime(activeChnList, currMeasNr) = transmitTime(activeChnList);
    [satPositions, satClkCorr] = satpos(transmitTime(activeChnList), ...
                                 [trackResults(activeChnList).PRN], eph_test);
    navSolutions_satClkCorr(activeChnList, currMeasNr) = satClkCorr;
    if size(activeChnList, 2) > 3
        clkCorrRawP = navSolutions_rawP(activeChnList, currMeasNr)' + satClkCorr * settings.c;
        [xyzdt,navSolutions_el(activeChnList, currMeasNr), ...
               navSolutions_az(activeChnList, currMeasNr), ...
               navSolutions_DOP(:, currMeasNr),satPositions] =...
                       leastSquarePos(satPositions, clkCorrRawP, settings);
        % disp('ok11')
        % return
        navSolutions_X(currMeasNr)  = xyzdt(1);
        navSolutions_Y(currMeasNr)  = xyzdt(2);
        navSolutions_Z(currMeasNr)  = xyzdt(3);       
        if (currMeasNr == 1)
        navSolutions_dt(currMeasNr) = 0;  % in unit of (m)
        else
            navSolutions_dt(currMeasNr) = xyzdt(4);  
        end
        navSolutions_currMeasSample(currMeasNr) = currMeasSample;
        satElev = navSolutions_el(:, currMeasNr)';
        navSolutions_correctedP(activeChnList, currMeasNr) = ...
                navSolutions_rawP(activeChnList, currMeasNr) + ...
                satClkCorr' * settings.c - xyzdt(4);
        [navSolutions_latitude(currMeasNr), ...
         navSolutions_longitude(currMeasNr), ...
         navSolutions_height(currMeasNr)] = cart2geo(...
                                            navSolutions_X(currMeasNr), ...
                                            navSolutions_Y(currMeasNr), ...
                                            navSolutions_Z(currMeasNr), ...
                                            5);
        
        navSolutions_utmZone = findUtmZone(navSolutions_latitude(currMeasNr), ...
                                       navSolutions_longitude(currMeasNr));        
        [navSolutions_E(currMeasNr), ...
         navSolutions_N(currMeasNr), ...
         navSolutions_U(currMeasNr)] = cart2utm(xyzdt(1), xyzdt(2), ...
                                                xyzdt(3), ...
                                                navSolutions_utmZone);
    else
        disp(['   Measurement No. ', num2str(currMeasNr), ...
                       ': Not enough information for position solution.']);
        navSolutions_X(currMeasNr)           = NaN;
        navSolutions_Y(currMeasNr)           = NaN;
        navSolutions_Z(currMeasNr)           = NaN;
        navSolutions_dt(currMeasNr)          = NaN;
        navSolutions_DOP(:, currMeasNr)      = zeros(5, 1);
        navSolutions_latitude(currMeasNr)    = NaN;
        navSolutions_longitude(currMeasNr)   = NaN;
        navSolutions_height(currMeasNr)      = NaN;
        navSolutions_E(currMeasNr)           = NaN;
        navSolutions_N(currMeasNr)           = NaN;
        navSolutions_U(currMeasNr)           = NaN;
        navSolutions_az(activeChnList, currMeasNr) = ...
                                             NaN(1, length(activeChnList));
        navSolutions_el(activeChnList, currMeasNr) = ...
                                             NaN(1, length(activeChnList));
    end 
    % tic
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % [navSolutions] = DPE_module...
    %     (currMeasNr,navSolutions,activeChnList,...
    %     trackResults,currMeasSample,satPositions,...
    %     transmitTime(activeChnList),localTime,...
    %     settings,satElev,fid,xyzdt(4),satClkCorr);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('ok45')
    
    % function [navSolutions]=DPE_module...
    % (currMeasNr,navSolutions,activeChnList,trackResults,currMeasSample,...
    % satPositions,transmitTime,localTime,settings,satElev,fid,...

    [navSolutions] = DPE_module...
        (currMeasNr,navSolutions,activeChnList,...
        trackResults,currMeasSample,satPositions,...
        transmitTime(activeChnList),localTime,...
        settings,satElev,fid,xyzdt(4),satClkCorr);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    navSolutions.DPE_processingtime(currMeasNr) = toc;
    localTime = localTime - xyzdt(4)/settings.c;       
    navSolutions.localTime(currMeasNr) = localTime;
    localTime = localTime + measSampleStep/settings.samplingFreq ;
    m2lat = 1/110734;
    m2lon = 1/103043;

    navSolutions.LLH_error(currMeasNr,1)=...
        norm(([navSolutions.latitude(currMeasNr),...
        navSolutions.longitude(currMeasNr)]...
        -settings.gt_llh(1:2))./[m2lat m2lon]);

    navSolutions.LLH_error(currMeasNr,2)=...
    norm(([navSolutions.DPE_latitude(currMeasNr),...
        navSolutions.DPE_longitude(currMeasNr)]...
        -settings.gt_llh(1:2))./[m2lat m2lon]);
    fprintf('\nCurrent 2D Error of DPE     : %1f\n',...
        (navSolutions.LLH_error(currMeasNr,2)));
    fprintf('Current 2D Error of LS      : %1f\n',...
        (navSolutions.LLH_error(currMeasNr,1)));
    return
    GT_ECEF =  llh2xyz(settings.gt_llh.*[pi/180 pi/180 1]);
    navSolutions.LLH_error(currMeasNr,3)= ...
        sqrt(sum(([navSolutions.X(currMeasNr)...
        navSolutions.Y(currMeasNr)...
        navSolutions.Z(currMeasNr)]-GT_ECEF).^2));

    pos_xyz = ...
        llh2xyz([navSolutions.DPE_latitude(currMeasNr)/180*pi,...
        navSolutions.DPE_longitude(currMeasNr)/180*pi,...
        navSolutions.DPE_height(currMeasNr)]);

    navSolutions.LLH_error(currMeasNr,4)=  ...
        sqrt(sum(([pos_xyz]-GT_ECEF).^2));

    % === Prints the 3D errors of both Least Squares and DPE ==============

    fprintf('\nCurrent 3D Error of DPE     : %1f\n',...
        (navSolutions.LLH_error(currMeasNr,4)));
    fprintf('Current 3D Error of LS      : %1f\n',...
        (navSolutions.LLH_error(currMeasNr,3)));

% end
end