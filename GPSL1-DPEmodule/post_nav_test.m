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
    return
    % function [navSolutions]=DPE_module...
    % (currMeasNr,navSolutions,activeChnList,trackResults,currMeasSample,...
    % satPositions,transmitTime,localTime,settings,satElev,fid,...
    % LS_clkBias,satClkCorr)
    % ------------------------ DPE_module v1.0 ----------------------------
    LS_clkBias = xyzdt(4);
    % === Initialize parameters ===========================================
    dtr = pi/180;
    m2lat = 1/110734;
    m2lon = 1/103043;
    trop = zeros(length(activeChnList),1);
    chip_spacings = [-(flip(settings.chipspacing_dpe_precalc:...
        settings.chipspacing_dpe_precalc:1)),0,...
        settings.chipspacing_dpe_precalc:...
        settings.chipspacing_dpe_precalc:1];
    precalc_correlations    = ...
        zeros(length(chip_spacings)...
        *length(activeChnList),3);
    % === Make folder to store correlograms ===============================
    % if ~exist([settings.outfile_root '\Correlogram\'], 'dir') && ...
    %         settings.DPE_plotCorrelogram == 1
    %     mkdir([settings.outfile_root '\Correlogram\']);
    % end
    % === Generate Latitude and Longitude search space ====================
    % Lat-long search space is of plus-minus "settings.DPE_latlong_span" 
    % meters
    % Centered at Least Squares' lat-long estimate
    [~,~,~,mesh_ll] = ...
        meshgrid_LIST(navSolutions.latitude(currMeasNr)+(m2lat).*...
        (-settings.DPE_latlong_span:settings.candPVT_spacing:settings.DPE_latlong_span)...
        ,navSolutions.longitude(currMeasNr)+(m2lon).*...
        (-settings.DPE_latlong_span:settings.candPVT_spacing:settings.DPE_latlong_span));
    % === Generate Height search space ====================================
    % Height search space is of plus-minus "settings.DPE_height_span" meters
    % Centered at Least Squares' height estimate
    alt_search = round(navSolutions_height(currMeasNr),0)-...
        settings.DPE_height_span:settings.candPVT_spacing:...
        round(navSolutions_height(currMeasNr),0)+settings.DPE_height_span;
    % === Generate Clock Bias search space ================================
    % Use a wider search space for first epoch to ensure clock bias
    % convergence 
    if currMeasNr == 1
        candia_xyz_LS = [navSolutions_X(currMeasNr),...
            navSolutions_Y(currMeasNr),navSolutions_Z(currMeasNr)];
        clock_bias=zeros(length(activeChnList),1);
        simulated_pseudorange = zeros(length(activeChnList),1);
        for i = 1:length(activeChnList)
            if (settings.useTropCorr == 1)
                % === Estimate tropospheric error for each satellite ======
                trop(i) = tropo(sin(satElev(activeChnList(i)) * dtr), ...
                    0.0, 1013.0, 293.0, 50.0, 0.0, 0.0, 0.0);
            end
            simulated_pseudorange(i)=sqrt(sum((satPositions(:,i)'...
                - candia_xyz_LS).^2));
            clock_bias(i) = ...
                (navSolutions_rawP(activeChnList(i),currMeasNr)+ ...
                satClkCorr(i) * settings.c - simulated_pseudorange(i));
        end
        clk_s_area=...
            round(min(clock_bias))-100:1:round(max(clock_bias))+100;
    else
    % Since Least Squares' clock bias is already stable and relatively accurate
    % after the first epoch, use narrow search space in consecutive positioning 
    % epochs
        clk_s_area=round((LS_clkBias))-settings.DPE_clkBias_span:1:...
            round((LS_clkBias))+settings.DPE_clkBias_span;
        if (settings.useTropCorr == 1)
            for i = 1:length(activeChnList)
            % === Estimate tropospheric error for each satellite ==========
                trop(i) = tropo(sin(satElev(activeChnList(i)) * dtr), ...
                    0.0, 1013.0, 293.0, 50.0, 0.0, 0.0, 0.0);
            end
        end
    end
    % If clock bias search space is too large
    if length(clk_s_area)>100000
        clk_s_area=round((LS_clkBias))-settings.DPE_clkBias_span:1:...
            round((LS_clkBias))+settings.DPE_clkBias_span;
    end
    % === Pre-calculate the correlation values ================================
    % A pure, traditional DPE implementation would require iterative
    % correlations for every candidate position, which is highly
    % computationally extensive. This section pre-calculates the correlations
    % in order of "settings.chipspacing_dpe_precalc" chip spacing
    % to save computational time
    if settings.MMT~=1        
        for j=1:length(activeChnList)
            % === Obtain closest tracking measurements to the current =====
            % === measurement sample ======================================
            [~,closestIndex] = min(abs(currMeasSample-...
                trackResults(activeChnList(j)).absoluteSample));
            if currMeasSample < ...
                    trackResults(activeChnList(j)).absoluteSample(closestIndex)
                closestIndex=closestIndex-1;
            end
            % === Move IF data to the current measurement sample ==============
            % fread would move the measurement sample automatically for every 
            % iteration of j. This code prevents the current sample to move.
            if closestIndex > (settings.DPE_nonCohInt-1)
                nonCohInt=closestIndex-settings.DPE_nonCohInt+1:closestIndex;
            else
                nonCohInt=1:closestIndex;
            end
            % Generate indexes to store the pre-calculated correlations in an
            % array
            count2=3;
            count=1+(j-1)*length(chip_spacings);
            % === Get the C/A code ============================================
            caCode = generateCAcode(trackResults(activeChnList(j)).PRN);
            % === Pre-calculate the correlations =================================
         for closestIndex=nonCohInt
            % Loop repeats depending on the non-coherent integration time
            % chosen
            % === Record sample number (based on 8-bit and 16-bit samples) ====
            if strcmp(settings.dataType,'int16') 
                fseek(fid, ...
                    settings.fileType*...
                    ((trackResults(activeChnList(j)).absoluteSample(closestIndex))*2),'bof');
            else
                fseek(fid, ...
                    settings.fileType*...
                    (trackResults(activeChnList(j)).absoluteSample(closestIndex)),'bof');
            end

            % === Update the phasestep based on code freq (variable) and ======
            % === sampling frequency (fixed) ==================================
            codePhaseStep = ...
                trackResults(activeChnList(j)).codeFreq(closestIndex)...
                / settings.samplingFreq;

            % === Find the size of a "block" or code period in whole samples ==
            blksize = ceil((settings.codeLength*settings.DPE_cohInt-...
                trackResults(activeChnList(j)).remCodePhase(closestIndex)) ...
                / codePhaseStep);

            % === Read in the appropriate number of samples to process =======
            [rawSignal, ~] = fread(fid, ...
                settings.fileType*blksize, settings.dataType);
            
            rawSignal = rawSignal';

            if (settings.fileType==2)
                rawSignal1=rawSignal(1:2:end);
                rawSignal2=rawSignal(2:2:end);
                rawSignal = rawSignal1 + 1i .* rawSignal2;%transpose vector
            end

            % === Get the argument to sin/cos functions =======================
            time    = (0:blksize) ./ settings.samplingFreq;
            remCarrPhase = ...
                trackResults(activeChnList(j)).remCarrPhase(closestIndex);
            trigarg = ...
                ((trackResults(activeChnList(j)).carrFreq(closestIndex)...
                * 2.0 * pi) .* time) + remCarrPhase;

            % === Compute the carrier signal ==================================
            carrsig = exp(1i .* trigarg(1:blksize));

            % === Mix the carrier replica to baseband =========================
            qBasebandSignal = real(carrsig .* rawSignal);
            iBasebandSignal = imag(carrsig .* rawSignal);

            for spacing=chip_spacings

            % === Define index into the code vector ===========================
            delay_index = ...
                trackResults(activeChnList(j)).remCodePhase(closestIndex)-spacing: ...
                codePhaseStep : ...
                ((blksize-1)*codePhaseStep + ...
                trackResults(activeChnList(j)).remCodePhase(closestIndex)-spacing);       
 
                caCode1 = [caCode(end-2) caCode(end-1) caCode(end) ...
                    repmat(caCode,1,settings.DPE_cohInt) caCode(1) caCode(2)...
                    caCode(3)];
                tcodee = ceil(delay_index)+3;
                tcodee(tcodee==0)=tcodee(tcodee==0)+1;

                s = caCode1(tcodee);

                I = sum(s  .* iBasebandSignal);
                Q = sum(s  .* qBasebandSignal);

            % === Store the correlations and its corresponding code delay =====
              precalc_correlations(count,1) = ...
                  trackResults(activeChnList(j)).PRN;
              precalc_correlations(count,2) = ...
                  (localTime-transmitTime(j))*settings.c ...
                  + (spacing/settings.codeFreqBasis) * settings.c;
              precalc_correlations(count,count2) = sqrt(I.^2 + Q.^2);
              
              count=count+1;
            end
              count=1+(j-1)*length(chip_spacings);
              count2=count2+1;
         end
         % End for pre-calculating the corr values for a single satellite
%
        end % End for pre-calculating the corr values for all satellites
% 
    else % Use MMT-DPE

         for j=1:length(activeChnList)
             % === Obtain closest tracking measurements to the current =========
             % === measurement sample ==========================================
             [~,closestIndex] = min(abs(currMeasSample-...
                 trackResults(activeChnList(j)).absoluteSample));
             if currMeasSample < ...
                     trackResults(activeChnList(j)).absoluteSample(closestIndex)
                 closestIndex=closestIndex-1;
             end
             
             % === Move IF data to the current measurement sample ==============
             % fread would move the measurement sample automatically for every 
             % iteration of j. This code prevents the current sample to move.
             if closestIndex > (settings.DPE_nonCohInt-1)
                 nonCohInt=closestIndex-settings.DPE_nonCohInt+1:closestIndex;
             else
                 nonCohInt=1:closestIndex;
             end

             % Generate indexes to store the pre-calculated correlations in an
             % array
             count2=3;
             count=1+(j-1)*length(chip_spacings);

             % === Get the C/A code ============================================
             caCode = generateCAcode(trackResults(activeChnList(j)).PRN);
             
        % === Pre-calculate the correlations =================================
        for closestIndex=nonCohInt
            % Loop repeats depending on the non-coherent integration time
            % chosen
            
            % === Record sample number (based on 8-bit and 16-bit samples) ====
            if strcmp(settings.dataType,'int16') 
                fseek(fid, ...
                    settings.fileType*...
                    ((trackResults(activeChnList(j)).absoluteSample(closestIndex))*2),'bof');
            else
                fseek(fid, ...
                    settings.fileType*...
                    (trackResults(activeChnList(j)).absoluteSample(closestIndex)),'bof');
            end
            % === Update the phasestep based on code freq (variable) and ======
            % === sampling frequency (fixed) ==================================
            codePhaseStep = ...
                trackResults(activeChnList(j)).codeFreq(closestIndex)...
                / settings.samplingFreq;
            
            % === Find the size of a "block" or code period in whole samples ==
            blksize = ceil((settings.codeLength*settings.DPE_cohInt-...
                trackResults(activeChnList(j)).remCodePhase(closestIndex)) ...
                / codePhaseStep);

            % === Read in the appropriate number of samples to process =======
            [rawSignal, ~] = fread(fid, ...
                settings.fileType*blksize, settings.dataType);

            rawSignal = rawSignal';

            if (settings.fileType==2)
                rawSignal1=rawSignal(1:2:end);
                rawSignal2=rawSignal(2:2:end);
                rawSignal = rawSignal1 + 1i .* rawSignal2;%transpose vector
            end

            % === Get the argument to sin/cos functions =======================
            time    = (0:blksize) ./ settings.samplingFreq;
            remCarrPhase = ...
                trackResults(activeChnList(j)).remCarrPhase(closestIndex);
            trigarg = ...
                ((trackResults(activeChnList(j)).carrFreq(closestIndex)...
                * 2.0 * pi) .* time) + remCarrPhase;

            % === Compute the carrier signal ==================================
            carrsig = exp(1i .* trigarg(1:blksize));

            % === Mix the carrier replica to baseband =========================
            qBasebandSignal = real(carrsig .* rawSignal);
            iBasebandSignal = imag(carrsig .* rawSignal);


            % === Get MMT estimates from tracking =============================
            LOS_delay = trackResults(activeChnList(j)).MMT_codeDelay_LOS(closestIndex);

            NLOS_delay = trackResults(activeChnList(j)).MMT_codeDelay_NLOS(closestIndex);

            LOS_carr_phase=trackResults(activeChnList(j)).MMT_carrPhase_LOS(closestIndex);

            NLOS_carr_phase=trackResults(activeChnList(j)).MMT_carrPhase_NLOS(closestIndex);

            LOS_ampl= trackResults(activeChnList(j)).MMT_amp_LOS(closestIndex);

            NLOS_ampl = trackResults(activeChnList(j)).MMT_amp_NLOS(closestIndex);

            % Compute complex amplitudes
            A = LOS_ampl*cos(LOS_carr_phase);
            C = LOS_ampl*sin(LOS_carr_phase);
            B = NLOS_ampl*cos(NLOS_carr_phase);
            D = NLOS_ampl*sin(NLOS_carr_phase);

            for spacing=chip_spacings
                delay1=spacing+LOS_delay;

                delay_index_0 = ...
                trackResults(activeChnList(j)).remCodePhase(closestIndex)  : ...
                    codePhaseStep : ...
                    ((blksize-1)*codePhaseStep + ...
                    trackResults(activeChnList(j)).remCodePhase(closestIndex) );


                delay_index_1 = ...
                trackResults(activeChnList(j)).remCodePhase(closestIndex) -delay1 : ...
                    codePhaseStep : ...
                    ((blksize-1)*codePhaseStep + ...
                    trackResults(activeChnList(j)).remCodePhase(closestIndex) -delay1);

                delay_index_2 = ...
                trackResults(activeChnList(j)).remCodePhase(closestIndex) -NLOS_delay : ...
                    codePhaseStep : ...
                    ((blksize-1)*codePhaseStep + ...
                    trackResults(activeChnList(j)).remCodePhase(closestIndex) -NLOS_delay);

                delay_index_12 = ...
                trackResults(activeChnList(j)).remCodePhase(closestIndex) -(NLOS_delay-delay1)  : ...
                    codePhaseStep : ...
                    ((blksize-1)*codePhaseStep + ...
                    trackResults(activeChnList(j)).remCodePhase(closestIndex) -(NLOS_delay-delay1) );

                caCode1 = [caCode(1023-2) caCode(1023-1) caCode(1023)...
                    repmat(caCode,1,settings.DPE_cohInt)...
                    caCode(1) caCode(2) caCode(3)];

                tcodee_0 = ceil(delay_index_0)+3;
                tcodee_1 = ceil(delay_index_1)+3;
                tcodee_2 = ceil(delay_index_2)+3;
                tcodee_12 = ceil(delay_index_12)+3;


                s_0 = caCode1(tcodee_0);
                s_1 = caCode1(tcodee_1);
                s_2 = caCode1(tcodee_2);
                s_12 = caCode1(tcodee_12);

                R_mm_12 = sum(s_0.*s_12);
                R_mm_0  = sum(s_0.*s_0); 


            % === Store the correlations and its corresponding code delay =====
              precalc_correlations(count,1) = ...
                  trackResults(activeChnList(j)).PRN;
              precalc_correlations(count,2) = ...
                  (localTime-transmitTime(j))*settings.c ...
                  + (spacing/settings.codeFreqBasis) * settings.c;
              precalc_correlations(count,count2) =  (2*real((A-1i*C)*...
                  (sum(iBasebandSignal.*s_1) + 1i*sum(qBasebandSignal.*s_1)))...
                  +2*real((B-1i*D)*...
                  (sum(iBasebandSignal.*s_2) + 1i*sum(qBasebandSignal.*s_2))) - ...
                  2*(A*B + C*D)*R_mm_12 - (A^2 + B^2 + C^2 +D^2)*R_mm_0);
              
              count=count+1;
            end
                count=1+(j-1)*length(chip_spacings);
                count2=count2+1;
        end 
        % End for pre-calculating the corr values for a single satellite
 
         end % End for pre-calculating the corr values for all satellites

    end

    % === Non-coherent integration of each satellite's autocorrelation function
    precalc_correlations(:,3) = sum(precalc_correlations(1:end,3:end),2);

    % === Initialize index for recording estimates across different heights ===
    barisan = 1;

    % === Pre-allocate space to store lat-long-clk bias estimates============== 
    % === for every altitude ==================================================
    temprecord_DPE_values = zeros(length(alt_search),5);

    % === Repeat lat-long-clkbias search space over height search space =======
    for alt = alt_search

        candia_map=[mesh_ll,ones(size(mesh_ll,1),1)*alt];

        x = candia_map(:,1)./180 .* pi;
        y = candia_map(:,2)./180 .* pi;
        candia_xyz  = llh2xyz([x,y,candia_map(:,3)]); 
        % Use llh2xyz by Todd Walter, 2001

        % === Pre-allocate space for correlogram ==============================
        correlogram=zeros(length(clk_s_area),length(candia_map));


        for j=1:length(activeChnList)

            % === Obtain ranges from satellite j to each candidate position ===
            candia_xyz_pseudoranges =  ...
                sqrt(sum(((repmat(satPositions(:,j)',length(candia_xyz),1)...
                - candia_xyz).^2),2));

            % === Compensate for the tropospheric and SV clock error ==========
            candia_xyz_pseudoranges = candia_xyz_pseudoranges...
                + trop(j) - satClkCorr(j)* settings.c;

            % === Estimate clock bias =========================================
            candia_xyz_pseudoranges_clkBias = zeros(length(clk_s_area),...
                length(candia_xyz_pseudoranges));          
            for blob = 1:length(clk_s_area)
                candia_xyz_pseudoranges_clkBias(blob,1:...
                    length(candia_xyz_pseudoranges)) = ...
                    candia_xyz_pseudoranges(1:end)+clk_s_area(blob);
            end

            prn_index = find(precalc_correlations(:,1)==...
                trackResults(activeChnList(j)).PRN);

            % === Obtain correlogram from single satellite ====================
            correlogram_single = ...
                interp1(precalc_correlations(prn_index,2),...
                precalc_correlations(prn_index,3),...
                candia_xyz_pseudoranges_clkBias,'linear');

            % === Removes NaN values from the interpolation and replaces it ===
            % === with another interpolation ==================================
            try
              nanx = isnan(correlogram_single);
              t    = 1:numel(correlogram_single);
              correlogram_single(nanx) = interp1(t(~nanx),...
              correlogram_single(~nanx), t(nanx),'linear');
            catch
            end

            % === Sum up the correlograms from each satellite =================
            correlogram=correlogram+correlogram_single;

        end % for j=1:length(activeChnList)

        % === Obtain the lat-long-clock bias estimate for a specific height ===
        [MaxCorrValue,~] = max(correlogram,[],'all','linear'); 
        % MaxCorrValue is the max value

        [I,J,~] =find(correlogram==MaxCorrValue); 
        % I holds the index for clock bias
        % J holds the index for latitude and longitude

        % === Record the lat-long-clock bias estimate for a specific height ===
        % In the case of more than one global maxima over the navigation
        % domain due to thermal noise, the mean of coordinates with max
        % correlation value will be used
        if length(I) ~= 1 
            temprecord_DPE_values(barisan,1:3)  = mean(candia_map(J,:));
            % Columns 1 to 3 holds the lat-long-height
            temprecord_DPE_values(barisan,4)    = mean(MaxCorrValue);
            % Column 4 holds the correlation value
            temprecord_DPE_values(barisan,5)    = mean(clk_s_area(I));
            % Column 5 holds the clock bias
            temprecord_DPE_values(barisan,6) = ...
                norm(([temprecord_DPE_values(barisan,1),...
                temprecord_DPE_values(barisan,2)]...
                -settings.gt_llh(1:2))./[m2lat m2lon]);
            % Column 6 holds the 2D error of the lat-long-height (for
            % evaluation purposes)
            temprecord_DPE_values(barisan,7) = round(mean(I));
            % Column 7 holds the clock bias row to be used to plot the
            % correlograms
        else
            temprecord_DPE_values(barisan,1:3)  = candia_map(J,:);
            % Columns 1 to 3 holds the lat-long-height
            temprecord_DPE_values(barisan,4)    = MaxCorrValue;
            % Column 4 holds the correlation value
            temprecord_DPE_values(barisan,5)    = clk_s_area(I);
            % Column 5 holds the clock bias
            temprecord_DPE_values(barisan,6) = norm(([candia_map(J,1),...
                candia_map(J,2)]...
                -settings.gt_llh(1:2))./[m2lat m2lon]);
            % Column 6 holds the 2D error of the lat-long-height (for
            % evaluation purposes)
            temprecord_DPE_values(barisan,7) = I;
            % Column 7 holds the clock bias row to be used to plot the
            % correlograms
        end
 
        barisan = barisan+1;
        clear I J MaxCorrValue
 
    end % for alt = alt_search 

    % === Obtain DPE estimate =================================================
    [~,barisan_yangmanaya] = max(temprecord_DPE_values(:,4));
    % barisan_yangmanaya holds the index of lat-long-height estimate of DPE

    % === Plot correlogram at DPE's estimated altitude ========================
    % Added for evaluation of performance - can be disabled in initSettings.m
    % by settings.DPE_plotCorrelogram == 0

    % The following basically repeates the above process at DPE's estimated
    % altitude
    if settings.DPE_plotCorrelogram == 1
        alt = temprecord_DPE_values(barisan_yangmanaya,3);

        candia_map=[mesh_ll,ones(size(mesh_ll,1),1)*alt];

        x = candia_map(:,1)./180 .* pi;
        y = candia_map(:,2)./180 .* pi;
        candia_xyz  = llh2xyz([x,y,candia_map(:,3)]); 
        % Use llh2xyz by Todd Walter, 2001

        correlogram=zeros(length(clk_s_area),length(candia_map));

        for j=1:length(activeChnList)

            candia_xyz_pseudoranges =  ...
                sqrt(sum(((repmat(satPositions(:,j)',length(candia_xyz),1)...
                - candia_xyz).^2),2));

            candia_xyz_pseudoranges = candia_xyz_pseudoranges...
                + trop(j) - satClkCorr(j)* settings.c;

            candia_xyz_pseudoranges_clkBias = zeros(length(clk_s_area),...
                length(candia_xyz_pseudoranges));          
            for blob = 1:length(clk_s_area)
                candia_xyz_pseudoranges_clkBias(blob,1:...
                    length(candia_xyz_pseudoranges)) = ...
                    candia_xyz_pseudoranges(1:end)+clk_s_area(blob);
            end

            dissstt=(localTime-transmitTime(j))*settings.c;

            codePhase3=...
                dissstt-(repmat(dissstt,...
                size(candia_xyz_pseudoranges_clkBias,1),...
                size(candia_xyz_pseudoranges_clkBias,2))-...
                candia_xyz_pseudoranges_clkBias)...
                ;

            prn_index = find(precalc_correlations(:,1)==...
                trackResults(activeChnList(j)).PRN);

            correlogram_single = ...
                interp1(precalc_correlations(prn_index,2),...
                precalc_correlations(prn_index,3),codePhase3,'linear');

            try
              nanx = isnan(correlogram_single);
              t    = 1:numel(correlogram_single);
              correlogram_single(nanx) = interp1(t(~nanx),...
                  correlogram_single(~nanx), t(nanx),'linear');
            catch
            end

        % === Plot correlogram for single satellite ===========================
            figure;
            W = reshape(correlogram_single((temprecord_DPE_values(barisan_yangmanaya,7)),:),...
                sqrt(length(candia_map)),sqrt(length(candia_map)));
            reshape_cand_llh_1 = reshape(candia_map(:,1),...
                [sqrt(length(candia_map)),sqrt(length(candia_map))]);
            reshape_cand_llh_2 = reshape(candia_map(:,2),...
                [sqrt(length(candia_map)),sqrt(length(candia_map))]);
            surf(reshape_cand_llh_1,reshape_cand_llh_2,W); % plot the result
            hold on;
            gtt = scatter3(settings.gt_llh(1),settings.gt_llh(2),...
                max(max(W)),50,'filled');
            gtt.MarkerFaceColor = "#000000";
            ylabel( 'Longitude', 'Interpreter', 'none');
            xlabel( 'Latitude', 'Interpreter', 'none');
            zlabel('Correlation');
            set(gca,'FontSize',10);
            set(gca, 'FontName', 'Arial');
            grid on;
            grid minor;
            view(0,90);
            xlim([min(candia_map(:,1)),max(candia_map(:,1))]);
            ylim([min(candia_map(:,2)),max(candia_map(:,2))]);
            colorbar;
            legend('','Ground Truth');
            ax = gca;
            ax.XAxis.LineWidth = 3;
            ax.YAxis.LineWidth = 3;
            ax.ZAxis.LineWidth = 3;
            ax.GridAlpha = 1;
            ax.MinorGridAlpha = 1;
            title(['PRN ', num2str(trackResults(activeChnList(j)).PRN)],...
                ['Epoch ', num2str(currMeasNr)]);
            hold off;
            % Save figure in PNG
            saveas(gcf, [pwd ,'\',settings.outfile_root,'\Correlogram\',...
                settings.outfile_root,'_Correlogram_Epoch',num2str(currMeasNr),...
                '_PRN',num2str(trackResults(activeChnList(j)).PRN),...
                '.png']);
            % Save figure in MATLAB figure
            saveas(gcf, [pwd ,'\',settings.outfile_root,'\Correlogram\',...
                settings.outfile_root,'_Correlogram_Epoch',num2str(currMeasNr),...
                '_PRN',num2str(trackResults(activeChnList(j)).PRN)]);

        % === Sum up the correlograms from each satellite =====================
        correlogram=correlogram+correlogram_single;

        end % for j=1:length(activeChnList)


    % === Plot correlogram from all satellites ================================
        figure;
        W = reshape(correlogram(round(temprecord_DPE_values(barisan_yangmanaya,7)),:),...
            sqrt(length(candia_map)),sqrt(length(candia_map)));
        reshape_cand_llh_1 = reshape(candia_map(:,1),...
            [sqrt(length(candia_map)),sqrt(length(candia_map))]);
        reshape_cand_llh_2 = reshape(candia_map(:,2),...
            [sqrt(length(candia_map)),sqrt(length(candia_map))]);
        surf(reshape_cand_llh_1,reshape_cand_llh_2,W); % plot the result
        hold on;
        gtt = scatter3(settings.gt_llh(1),settings.gt_llh(2),...
            max(max(W)),50,'filled');
        gtt.MarkerFaceColor = "#000000";
        dpelatlong = scatter3(temprecord_DPE_values(barisan_yangmanaya,1),...
            temprecord_DPE_values(barisan_yangmanaya,2),...
            max(max(W)),50,'filled');
        dpelatlong.MarkerFaceColor = "#A2142F";
        stllatlong = scatter3(navSolutions.latitude(currMeasNr),...
            navSolutions.longitude(currMeasNr),...
            max(max(W)),50,'filled');
        stllatlong.MarkerFaceColor = "#00FF00";
        ylabel( 'Longitude', 'Interpreter', 'none');
        xlabel( 'Latitude', 'Interpreter', 'none');
        zlabel('Correlation');
        set(gca,'FontSize',10);
        set(gca, 'FontName', 'Arial');
        grid on;
        grid minor;
        view(0,90);
        xlim([min(candia_map(:,1)),max(candia_map(:,1))]);
        ylim([min(candia_map(:,2)),max(candia_map(:,2))]);
        colorbar;
        legend('','Ground Truth','DPE position','Least Squares (STL)');
        ax = gca;
        ax.XAxis.LineWidth = 3;
        ax.YAxis.LineWidth = 3;
        ax.ZAxis.LineWidth = 3;
        ax.GridAlpha = 1;
        ax.MinorGridAlpha = 1;
        title(['Epoch ', num2str(currMeasNr)], ['DPE 2D Positioning Error = ', ...
            num2str(temprecord_DPE_values(barisan_yangmanaya,6)),' meters']);
        hold off;
        % Save figure in PNG
        saveas(gcf, [pwd ,'\',settings.outfile_root,'\Correlogram\',...
        settings.outfile_root,'_Correlogram_Epoch',num2str(currMeasNr),...
        '.png']);
        % Save figure in MATLAB figure
        saveas(gcf, [pwd ,'\',settings.outfile_root,'\Correlogram\',...
        settings.outfile_root,'_Correlogram_Epoch',num2str(currMeasNr)]);

        close all
    end

    % === Record DPE estimates ================================================
    navSolutions.DPE_estimate(currMeasNr, 1:5) =  ...
        temprecord_DPE_values(barisan_yangmanaya,1:5);

    navSolutions.DPE_latitude(currMeasNr) =  ...
        temprecord_DPE_values(barisan_yangmanaya,1);

    navSolutions.DPE_longitude(currMeasNr) =  ...
        temprecord_DPE_values(barisan_yangmanaya,2);

    navSolutions.DPE_height(currMeasNr) =  ...
        temprecord_DPE_values(barisan_yangmanaya,3);

    navSolutions.DPE_clkBias(currMeasNr) = ...
        temprecord_DPE_values(barisan_yangmanaya,5);


    [navSolutions2] = DPE_module...
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