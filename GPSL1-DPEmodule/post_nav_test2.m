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
% eph_test = struct();
for channelNr = activeChnList
    PRN = trackResults(channelNr).PRN;
    fprintf('Decoding NAV for PRN %02d -------------------- \n', PRN);
    try
    [eph_ch(PRN), subFrameStart(channelNr), TOW(channelNr)] = ...
                                  NAVdecoding(trackResults(channelNr).I_P,...
                                  settings);
    catch
    end
    try
    if (isempty(eph_ch(PRN).IODC) || isempty(eph_ch(PRN).IODE_sf2) || ...
        isempty(eph_ch(PRN).IODE_sf3))
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
    % disp([channelNr, subFrameStart(channelNr), trackResults(channelNr).absoluteSample(subFrameStart(channelNr))])
    %disp(sampleStart(channelNr))
    sampleEnd(channelNr) = trackResults(channelNr).absoluteSample(end);
    % disp(sampleEnd(channelNr))
end
sampleStart = max(sampleStart) + 1;  
sampleEnd = min(sampleEnd) - 1;
% disp(sampleEnd)

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
    % disp([trackResults(activeChnList).PRN])
    navSolutions_el(:, currMeasNr) = NaN(settings.numberOfChannels, 1);
    navSolutions_az(:, currMeasNr) = NaN(settings.numberOfChannels, 1);
    navSolutions_transmitTime(:, currMeasNr) = NaN(settings.numberOfChannels, 1);
    navSolutions_satClkCorr(:, currMeasNr) = NaN(settings.numberOfChannels, 1);                                                                  
    currMeasSample = sampleStart + measSampleStep*(currMeasNr-1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate pseudorange
    transmitTime = inf(1, settings.numberOfChannels);
    codePhase = zeros(length(activeChnList),1);
    
    for channelNr = activeChnList
    % Find index of I_P stream whose integration contains current 
    % measurment point location 
        for index = 1: length(trackResults(channelNr).absoluteSample)
            % disp([index, trackResults(channelNr).absoluteSample(index)])
            if(trackResults(channelNr).absoluteSample(index) > currMeasSample )
                break
            end 
        end
        index = index - 1; 
        % Update the phasestep based on code freq and sampling frequency
        codePhaseStep = trackResults(channelNr).codeFreq(index) / settings.samplingFreq; 
        % Code phase from start of a PRN code to current measurement sample location 
        codePhase(channelNr)  = (trackResults(channelNr).remCodePhase(index) +  ...
            codePhaseStep * (currMeasSample - ...
            trackResults(channelNr).absoluteSample(index) ));
        % Transmitting Time (in unite of s)at current measurement sample location
        % codePhase/settings.codeLength: fraction part of a PRN code
        % index - subFrameStart(channelNr): integer number of PRN code
        parte1 = (codePhase(channelNr)/settings.codeLength/settings.DPE_cohInt ...
                              + index  - ...
                              subFrameStart(channelNr)) ;
        parte2 = parte1*settings.codeLength * settings.DPE_cohInt /...
                              settings.codeFreqBasis + TOW(channelNr);
        transmitTime(channelNr) =  ...
            (codePhase(channelNr)/settings.codeLength/settings.DPE_cohInt ...
                              + index  - ...
                              subFrameStart(channelNr)) * ...
                              settings.codeLength * settings.DPE_cohInt /...
                              settings.codeFreqBasis + TOW(channelNr);
    end
    
    % At first time of fix, local time is initialized by transmitTime and 
    % settings.startOffset

    
    if (localTime == inf)
         maxTime   = max(transmitTime(activeChnList));
         localTime = maxTime + settings.startOffset/1000;  
    end
    
    %--- Convert travel time to a distance ------------------------------------
    % The speed of light must be converted from meters per second to meters
    % per millisecond. 
    pseudoranges    = (localTime - transmitTime) * settings.c;
    navSolutions_rawP(:, currMeasNr) = pseudoranges;
    navSolutions_transmitTime(activeChnList, currMeasNr) = ...
        transmitTime(activeChnList);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    transmitTime = transmitTime(activeChnList);
    prnList = [trackResults(activeChnList).PRN];
    eph = eph_ch;
    %% Initialize constants ===================================================
    numOfSatellites = size(prnList, 2);
    % GPS constatns
    gpsPi          = 3.1415926535898;  % Pi used in the GPS coordinate system
    %--- Constants for satellite position calculation -------------------------
    Omegae_dot     = 7.2921151467e-5;  % Earth rotation rate, [rad/s]
    GM             = 3.986005e14;      % Earth's universal gravitational constant,
                                   % [m^3/s^2]
    F              = -4.442807633e-10; % Constant, [sec/(meter)^(1/2)]
    %% Initialize results =====================================================
    satClkCorr   = zeros(1, numOfSatellites);
    satPositions = zeros(3, numOfSatellites);
    %% Process each satellite =================================================
    for satNr = 1 : numOfSatellites 
        prn = prnList(satNr);
        %% Find initial satellite clock correction --------------------------------
        %--- Find time difference ---------------------------------------------
        % disp(transmitTime(satNr))
        % disp(eph(prn).t_oc)
        dt = check_t(transmitTime(satNr) - eph(prn).t_oc);
        % disp(eph(prn).a_f2)
        % disp(eph(prn).a_f1)
        % disp(eph(prn).a_f0)
        % disp(eph(prn).T_GD)
        %--- Calculate clock correction ---------------------------------------
        satClkCorr(satNr) = (eph(prn).a_f2 * dt + eph(prn).a_f1) * dt + ...
            eph(prn).a_f0 - ...
            eph(prn).T_GD;
        time = transmitTime(satNr) - satClkCorr(satNr);
        %% Find satellite's position ----------------------------------------------
        % Restore semi-major axis
        a   = eph(prn).sqrtA * eph(prn).sqrtA;
        % Time correction
        tk  = check_t(time - eph(prn).t_oe);
        % disp(tk)
        % Initial mean motion
        n0  = sqrt(GM / a^3);
        % disp(n0)
        % Mean motion
        n   = n0 + eph(prn).deltan;
        % Mean anomaly
        M   = eph(prn).M_0 + n * tk;        
        % Reduce mean anomaly to between 0 and 360 deg
        M   = rem(M + 2*gpsPi, 2*gpsPi);
        %Initial guess of eccentric anomaly
        E   = M; 
        %--- Iteratively compute eccentric anomaly ----------------------------
        % disp(eph(prn).e)
        for ii = 1:10
            E_old   = E;
            E       = M + eph(prn).e * sin(E);
            dE      = rem(E - E_old, 2*gpsPi);
            if abs(dE) < 1.e-12
                % Necessary precision is reached, exit from the loop
                break;
            end
        end
        % Reduce eccentric anomaly to between 0 and 360 deg
        E   = rem(E + 2*gpsPi, 2*gpsPi);      

        % Relativistic correction
        dtr = F * eph(prn).e * eph(prn).sqrtA * sin(E);

        %Calculate the true anomaly
        nu   = atan2(sqrt(1 - eph(prn).e^2) * sin(E), cos(E)-eph(prn).e);

        %Compute angle phi
        phi = nu + eph(prn).omega;
        
        %Reduce phi to between 0 and 360 deg
        phi = rem(phi, 2*gpsPi);

        %Correct argument of latitude
        u = phi + ...
            eph(prn).C_uc * cos(2*phi) + ...
            eph(prn).C_us * sin(2*phi);	
    	
        % Correct radius
        r = a * (1 - eph(prn).e*cos(E)) + ...
            eph(prn).C_rc * cos(2*phi) + ...
            eph(prn).C_rs * sin(2*phi);
        % Correct inclination
        i = eph(prn).i_0 + eph(prn).iDot * tk + ...
            eph(prn).C_ic * cos(2*phi) + ...
            eph(prn).C_is * sin(2*phi);
    
        % 2.9 SV position in orbital plane
        xk1 = cos(u)*r;        
        yk1 = sin(u)*r;
        

        %Compute the angle between the ascending node and the Greenwich meridian
        Omega = eph(prn).omega_0 + (eph(prn).omegaDot - Omegae_dot)*tk - ...
                Omegae_dot * eph(prn).t_oe;
        %Reduce to between 0 and 360 deg
        % disp(Omega)
        Omega = rem(Omega + 2*gpsPi, 2*gpsPi);
        %--- Compute satellite coordinates ------------------------------------
        xk = xk1 * cos(Omega) - yk1 * cos(i)*sin(Omega);
        yk = xk1 * sin(Omega) + yk1 * cos(i)*cos(Omega);
        zk = yk1 * sin(i);
        satPositions(1, satNr) = xk;
        satPositions(2, satNr) = yk;
        satPositions(3, satNr) = zk;
        %% Include relativistic correction in clock correction --------------------
        satClkCorr(satNr) = (eph(prn).a_f2 * dt + eph(prn).a_f1) * dt + ...
                            eph(prn).a_f0 - ...
                            eph(prn).T_GD + dtr;       
        % disp(satClkCorr(satNr))
    end % for satNr = 1 : numOfSatellites
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    navSolutions_satClkCorr(activeChnList, currMeasNr) = satClkCorr;
    if size(activeChnList, 2) > 3
        clkCorrRawP = navSolutions_rawP(activeChnList, currMeasNr)' + satClkCorr * settings.c;
        %% leastSquarePos
        %=== Initialization =======================================================
        satpos = satPositions;
        obs = clkCorrRawP;
        nmbOfIterations = 10;
        dtr     = pi/180;
        pos     = zeros(4, 1);   % center of earth
        X       = satpos;
        nmbOfSatellites = size(satpos, 2);
        A       = zeros(nmbOfSatellites, 4);
        omc     = zeros(nmbOfSatellites, 1);
        az      = zeros(1, nmbOfSatellites);
        el      = az;
        Rot_X_up = zeros(3, nmbOfSatellites);
        %=== Iteratively find receiver position ===================================
        for iter = 1:nmbOfIterations
             for i = 1:nmbOfSatellites
                 if iter == 1
                     %--- Initialize variables at the first iteration --------------
                     Rot_X = X(:, i);
                     % disp(Rot_X)
                     trop = 2;
                 else
                     % disp(iter)
                     % disp(pos)
                     % return
                     %--- Update equations -----------------------------------------
                     rho2 = (X(1, i) - pos(1))^2 + (X(2, i) - pos(2))^2 + ...
                         (X(3, i) - pos(3))^2;
                     traveltime = sqrt(rho2) / settings.c ;
                     %--- Correct satellite position (do to earth rotation) --------
                     % Convert SV position at signal transmitting time to position 
                     % at signal receiving time. ECEF always changes with time as 
                     % earth rotates.
                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                     % Returns rotated satellite ECEF coordinates due to 
                     % Earth rotation during signal travel time
                     % Rot_X = e_r_corr(traveltime, X(:, i));
                     X_sat = X(:, i);
                     
                     Omegae_dot = 7.292115147e-5;           %  rad/sec
                     %--- Find rotation angle -----------------------------
                     omegatau   = Omegae_dot * traveltime;
                     %--- Make a rotation matrix --------------------------
                     R3 = [ cos(omegatau)    sin(omegatau)   0;
                         -sin(omegatau)    cos(omegatau)   0;
                         0                0               1];
                     %--- Do the rotation ---------------------------------
                     X_sat_rot = R3 * X_sat;
                     Rot_X = X_sat_rot;
                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                     Rot_X_up(:,i) = Rot_X;
                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                     %--- Find the elevation angle of the satellite ----------------
                     % [az(i), el(i), ~] = topocent(pos(1:3, :), Rot_X - pos(1:3, :));
                     X_topocent = pos(1:3, :);
                     dx_topocent = Rot_X - pos(1:3, :);
                     % function [Az, El, D] = topocent(X_topocent, dx_topocent)
                     %TOPOCENT  Transformation of vector dx into topocentric coordinate
                     dtr = pi/180;
                     % [phi, lambda, h] = togeod(6378137, 298.257223563,...
                     %      X_topocent(1), X_topocent(2), X_topocent(3));
                     %TOGEOD   Subroutine to calculate geodetic 
                     % coordinates latitude, longitude,
                     % height given Cartesian coordinates X,Y,Z, and reference ellipsoid
                     % values semi-major axis (a) and the inverse of flattening (finv).
                     % % function [dphi, dlambda, h] = togeod(a, finv, X_togeod, Y_togeod, Z_togeod)
                     a = 6378137; % semi-major axis of the reference ellipsoid
                     finv = 298.257223563; % inverse of flattening of the reference ellipsoid
                     X_togeod = X_topocent(1);
                     Y_togeod = X_topocent(2);
                     Z_togeod = X_topocent(3);
                     %TOGEOD   Subroutine to calculate geodetic 
                     % coordinates latitude, longitude,
                     % height given Cartesian coordinates X,Y,Z, 
                     % and reference ellipsoid
                     % values semi-major axis (a) and the inverse of flattening (finv).
                     h_trop       = 0;
                     tolsq   = 1.e-10;
                     maxit   = 10;
                     % compute radians-to-degree factor
                     rtd     = 180/pi;
                     % disp([X_togeod, Y_togeod, Z_togeod])
                     
                     % compute square of eccentricity
                     if finv < 1.e-20
                         esq = 0;
                     else
                         esq = (2 - 1/finv) / finv;
                     end

                     
                     oneesq  = 1 - esq;
                     % first guess
                     % P is distance from spin axis
                     P = sqrt(X_togeod^2+Y_togeod^2);
                     % disp(P)
                     % direct calculation of longitude
                     if P > 1.e-20
                         dlambda = atan2(Y_togeod,X_togeod) * rtd;
                     else
                         dlambda = 0;
                     end

                     if (dlambda < 0)
                         dlambda = dlambda + 360;
                     end

                     % r is distance from origin (0,0,0)
                     r = sqrt(P^2 + Z_togeod^2);
                     if r > 1.e-20
                         sinphi = Z_togeod/r;
                     else
                         sinphi = 0;
                     end
                     dphi = asin(sinphi);
                     % initial value of height  =  distance from origin minus
                     % approximate distance from origin to surface of ellipsoid
                     if r < 1.e-20
                         h_trop = 0;
                         return
                     end
                     h_trop = r - a*(1-sinphi*sinphi/finv);
                     
                     % iterate
                     for iii = 1:maxit
                         sinphi  = sin(dphi);
                         cosphi  = cos(dphi);
                         % compute radius of curvature in prime vertical direction
                         N_phi   = a/sqrt(1-esq*sinphi*sinphi);
                         % compute residuals in P and Z
                         dP      = P - (N_phi + h_trop) * cosphi;
                         dZ      = Z_togeod - (N_phi*oneesq + h_trop) * sinphi;
                         % update height and latitude
                         h_trop       = h_trop + (sinphi*dZ + cosphi*dP);
                         dphi    = dphi + (cosphi*dZ - sinphi*dP)/(N_phi + h_trop);
                         % test for convergence
                         if (dP*dP + dZ*dZ < tolsq)
                             break;
                         end
                         % Not Converged--Warn user
                         if iii == maxit
                             fprintf([' Problem in TOGEOD, ...' ...
                                 'did not converge in %2.0f',...
                                 ' iterations\n'], iii);
                         end
                     end % for iii = 1:maxit
                     dphi = dphi * rtd;
                     phi = dphi;
                     lambda = dlambda;
                     h = h_trop;
                     
                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                     cl  = cos(lambda * dtr);
                     sl  = sin(lambda * dtr);
                     cb  = cos(phi * dtr); 
                     sb  = sin(phi * dtr);
                     % disp([cl, sl, cb, sb])
                     
                     F   = [-sl -sb*cl cb*cl;
                         cl -sb*sl cb*sl;
                         0    cb   sb];
                     local_vector = F' * dx_topocent;
                     E_topocent   = local_vector(1);
                     N_topocent   = local_vector(2);
                     U_topocent   = local_vector(3);
                     hor_dis = sqrt(E_topocent^2 + N_topocent^2);
                     if hor_dis < 1.e-20
                         Az_topocent = 0;
                         El_topocent = 90;
                     else
                         Az_topocent = atan2(E_topocent, N_topocent)/dtr;
                         El_topocent = atan2(U_topocent, hor_dis)/dtr;
                     end
                     if Az_topocent < 0
                         Az_topocent = Az_topocent + 360;
                     end
                     D   = sqrt(dx_topocent(1)^2 + dx_topocent(2)^2 + dx_topocent(3)^2);
                     az(i) = Az_topocent;
                     el(i) = El_topocent;
                     % disp([D, az(i), el(i)])
                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                     if (settings.useTropCorr == 1)
                         % trop = tropo(sin(el(i) * dtr), ...
                         %     0.0, 1013.0, 293.0, 50.0, 0.0, 0.0, 0.0);
                         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                         % function ddr = tropo(sinel, hsta, p, tkel, hum, hp, htkel, hhum)
                         % disp(el(i))
                         
                         sinel = sin(el(i) * dtr);
                         % disp(sinel)   
                         
                         hsta_1 = 0.0;
                         p = 1013.0;
                         tkel = 293.0;
                         hum = 50.0;
                         hp = 0.0;
                         htkel = 0.0;
                         hhum = 0.0;
                         a_e    = 6378.137;     % semi-major axis of earth ellipsoid
                         b0     = 7.839257e-5;
                         tlapse = -6.5;
                         tkhum  = tkel + tlapse*(hhum-htkel);
                         atkel  = 7.5*(tkhum-273.15) / (237.3+tkhum-273.15);
                         e0     = 0.0611 * hum * 10^atkel;
                         tksea  = tkel - tlapse*htkel;
                         em     = -978.77 / (2.8704e6*tlapse*1.0e-5);
                         tkelh  = tksea + tlapse*hhum;
                         e0sea  = e0 * (tksea/tkelh)^(4*em);
                         tkelp  = tksea + tlapse*hp;
                         psea   = p * (tksea/tkelp)^em;

                         
                         if sinel < 0
                             sinel = 0;
                         end                         
                         tropo_1   = 0;
                         done    = 'FALSE';
                         refsea_1  = 77.624e-6 / tksea;
                         htop_1    = 1.1385e-5 / refsea_1;
                         refsea_1  = refsea_1 * psea;
                         ref_1     = refsea_1 * ((htop_1-hsta_1)/htop_1)^4;
                         % disp(sinel)
                         % disp(ref_1)
                         
                         while 1
                             rtop_1 = (a_e+htop_1)^2 - (a_e+hsta_1)^2*(1-sinel^2);
                             % check to see if geometry is crazy
                             if rtop_1 < 0
                                 rtop_1 = 0; 
                             end 
                             rtop_1 = sqrt(rtop_1) - (a_e+hsta_1)*sinel;
                             a    = -sinel/(htop_1-hsta_1);
                             b    = -b0*(1-sinel^2) / (htop_1-hsta_1);
                             % disp([rtop_1, a , b])
                             rn   = zeros(8,1);

                             for iiii = 1:8
                                 rn(iiii) = rtop_1^(iiii+1);
                             end
                             % disp(rn)

                             alpha = [2*a, 2*a^2+4*b/3, a*(a^2+3*b),...
                                 a^4/5+2.4*a^2*b+1.2*b^2, 2*a*b*(a^2+3*b)/3,...
                                 b^2*(6*a^2+4*b)*1.428571e-1, 0, 0];

                             if b^2 > 1.0e-35
                                 alpha(7) = a*b^3/2; 
                                 alpha(8) = b^4/9; 
                             end

                             dr = rtop_1;
                             dr = dr + alpha*rn;
                             tropo_1 = tropo_1 + dr*ref_1*1000;

                             if done == 'TRUE '
                                 ddr = tropo_1; 
                                 break; 
                             end

                             done    = 'TRUE ';
                             refsea  = (371900.0e-6/tksea-12.92e-6)/tksea;
                             htop    = 1.1385e-5 * (1255/tksea+0.05)/refsea;
                             ref_1     = refsea * e0sea * ((htop_1-hsta_1)/htop)^4;
                         end
                         %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                         trop = ddr;
                         % disp(ddr)
                         % return
                     else
                         % Do not calculate or apply the tropospheric corrections
                         trop = 0;
                     end
                     % disp(trop)
                     % return
                 end % if iter == 1 ... ... else 
                 
                 % disp(Rot_X)
                 % disp(pos(1:3))
                 
                 %--- Apply the corrections ----------------------------------------
                 % disp(obs(i))
                 diff_rot_x = Rot_X - pos(1:3);
                 % disp(diff_rot_x)
                 % disp(pos(4))
                 % disp(norm(diff_rot_x, 'fro'))
                 norm_value = norm(diff_rot_x, 'fro');
                 
                 % disp(norm_value)
                 omc(i) = ( obs(i) - norm(Rot_X - pos(1:3), 'fro') - pos(4) - trop );
                 % disp(omc(i))
                 %disp(omc(i))
                 % disp(obs(i))
                 % disp(norm(Rot_X - pos(1:3), 'fro') -  )
                 % return
                 % disp(pos(1:3))
                 %--- Construct the A matrix ---------------------------------------
                 A(i, :) =  [ (-(Rot_X(1) - pos(1))) / norm(Rot_X - pos(1:3), 'fro') ...
                     (-(Rot_X(2) - pos(2))) / norm(Rot_X - pos(1:3), 'fro') ...
                     (-(Rot_X(3) - pos(3))) / norm(Rot_X - pos(1:3), 'fro') ...
                     1 ];
                 % return
                 % if iter==2
                 %     % disp(iter)
                 %     disp([iter,i])
                 %     % disp(obs(i))
                 %     disp(omc(i))
                 % 
                 % 
                 % 
                 %     return
                 % end
                 % disp([i, A(i, :)])
                 % disp(A(i, :))
                 % A_2 =  [ -(Rot_X(1) - pos(1))/norm_value ...
                 %     -(Rot_X(2) - pos(2))/norm_value ...
                 %     -(Rot_X(3) - pos(3))/norm_value ...
                 %     1 ];
                 % 
                 % disp(A(i, :))
                 % disp(A_2)
                 % % disp()
                 % return
                 
             
            end % for i = 1:nmbOfSatellites
            % These lines allow the code to exit gracefully in case of any errors
            % disp(A)
            if rank(A) ~= 4
                pos     = zeros(1, 4);
                dop     = inf(1, 5);
                fprintf('Cannot get a converged solution! \n');
                return
            end
            % disp(rank(A))
            % return
            %--- Find position update (in the least squares sense)-----------------
            x   = A \ omc;
            %--- Apply position update --------------------------------------------
            pos = pos + x;
            % disp(iter)
            % disp(pos)           
        end % for iter = 1:nmbOfIterations
        %--- Fixing result --------------------------------------------------------
        pos = pos';
        % Saves satellite position for more accurate satellite position for DPE 
        satpos=Rot_X_up; 
        %=== Calculate Dilution Of Precision ======================================
        % Commented to always calculate DOP
        % if nargout  == 4
        %--- Initialize output ------------------------------------------------
        dop     = zeros(1, 5);

        %--- Calculate DOP ----------------------------------------------------
        Q       = inv(A'*A);

        dop(1)  = sqrt(trace(Q));                       % GDOP    
        dop(2)  = sqrt(Q(1,1) + Q(2,2) + Q(3,3));       % PDOP
        dop(3)  = sqrt(Q(1,1) + Q(2,2));                % HDOP
        dop(4)  = sqrt(Q(3,3));                         % VDOP
        dop(5)  = sqrt(Q(4,4));                         % TDOP
        % % end  % if nargout  == 4
        xyzdt = pos;
        navSolutions_el(activeChnList, currMeasNr) = el;
        navSolutions_az(activeChnList, currMeasNr) = az;
        navSolutions_DOP(:, currMeasNr) = dop;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        
        X = navSolutions_X(currMeasNr);
        Y = navSolutions_Y(currMeasNr);
        Z = navSolutions_Z(currMeasNr);
        i = 5;
        % disp([X,Y,Z,i])
        %CART2GEO Conversion of Cartesian coordinates (X,Y,Z) to geographical
        %coordinates (phi, lambda, h) on a selected reference ellipsoid.
        %
        %==========================================================================
        a = [6378388 6378160 6378135 6378137 6378137];
        f = [1/297 1/298.247 1/298.26 1/298.257222101 1/298.257223563];
        lambda = atan2(Y,X);
        ex2 = (2-f(i))*f(i)/((1-f(i))^2);
        c = a(i)*sqrt(1+ex2);
        phi = atan(Z/((sqrt(X^2+Y^2)*(1-(2-f(i)))*f(i))));
        h = 0.1; oldh = 0;
        iterations = 0;
        while abs(h-oldh) > 1.e-12
            oldh = h;
            N = c/sqrt(1+ex2*cos(phi)^2);
            phi = atan(Z/((sqrt(X^2+Y^2)*(1-(2-f(i))*f(i)*N/(N+h)))));
            h = sqrt(X^2+Y^2)/cos(phi)-N;
            iterations = iterations + 1;
            if iterations > 100
                fprintf('Failed to approximate h with desired precision. h-oldh: %e.\n', h-oldh);
                 break;
            end
        end
        phi = phi*180/pi;
        lambda = lambda*180/pi;
        % disp([phi, lambda, h])
        % disp([navSolutions_latitude(currMeasNr), ...
        %  navSolutions_longitude(currMeasNr), ...
        %  navSolutions_height(currMeasNr)])
        %%%%%%%%%%%%%% end cart2geo.m %%%%%%%%%%%%%%%%%%%

        
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
    
    fprintf('Current 2D Error of LS      : %1f\n',...
        (navSolutions.LLH_error(currMeasNr,1)));

end