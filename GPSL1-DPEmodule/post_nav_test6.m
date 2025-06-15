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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % function [eph, subFrameStart,TOW] = NAVdecoding(I_P_InputBits,settings)
    I_P_InputBits = trackResults(channelNr).I_P;
    % disp(I_P_InputBits(1))
    % return
    
    % findPreambles finds the first preamble occurrence in the bit stream of
    % each channel. The preamble is verified by check of the spacing between
    % preambles (6sec) and parity checking of the first two words in a
    % subframe. At the same time function returns list of channels, that are in
    % tracking state and with valid preambles in the nav data stream.
    %
    %[eph, subFrameStart,TOW] = CNAVdecoding(I_P_InputBits)
    %
    %   Inputs:
    %       I_P_InputBits   - output from the tracking function
    %
    %   Outputs:
    %       subFrameStart   - Starting positions of the first message in the 
    %                       input bit stream I_P_InputBits in each channel. 
    %                       The position is CNAV bit(20ms before convolutional decoding) 
    %                       count since start of tracking. Corresponding value will
    %                       be set to inf if no valid preambles were detected in
    %                       the channel.
    %       TOW             - Time Of Week (TOW) of the first message(in seconds).
    %                       Corresponding value will be set to inf if no valid preambles
    %                       were detected in the channel.
    %       eph             - SV ephemeris. 
    %--- Initialize ephemeris structute  --------------------------------------
    % This is in order to make sure variable 'eph' for each SV has a similar 
    % structure when only one or even none of the three requisite messages
    % is decoded for a given PRN.
    eph_new = eph_structure_init();


    % Starting positions of the first message in the input bit stream
    subFrameStart_new = inf;
    

    % TOW of the first message
    TOW_new = inf;

    %% Bit and frame synchronization ====================================
    % Preamble search can be delayed to a later point in the tracking results
    % to avoid noise due to tracking loop transients
    searchStartOffset = 1000; 
    %--- Generate the preamble pattern ----------------------------------------
    preamble_bits = [1 -1 -1 -1 1 -1 1 1];

    % "Upsample" the preamble - make 20 vales per one bit. The preamble must be
    % found with precision of a sample.

    preamble_ms = kron(preamble_bits, ones(1, 20/settings.DPE_cohInt));
    
    % Correlate tracking output with preamble =================================
    % Read output from tracking. It contains the navigation bits. The start
    % of record is skiped here to avoid tracking loop transients.
    bits = I_P_InputBits(1 + searchStartOffset : end);
    % Now threshold the output and convert it to -1 and +1
    bits(bits > 0)  =  1;
    bits(bits <= 0) = -1;
    % Correlate tracking output with the preamble
    % tlmXcorrResult = xcorr(bits, preamble_ms);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Cross-correlation of two column vectors.
    x = bits;
    y_mod = preamble_ms(:);% y was validated to be a vector, make it a column.
    maxlagDefault_mod = max(size(x,1),size(y_mod,1)) - 1;
    maxlag_mod = maxlagDefault_mod;
    % Compute cross-correlation for vector inputs. Output is clipped based on
    % maxlag but not padded if maxlag >= max(size(x,1),size(y,1)).
    nx_mod = numel(x);
    ny_mod = numel(y_mod);
    m_mod = max(nx_mod,ny_mod);
    maxlagDefault_mod = m_mod-1;
    mxl_mod = maxlagDefault_mod;
    m_mod = 2*m_mod;
    while true
        r_mod = m_mod;
        for p_mod = [2 3 5 7]
            while (r_mod > 1) && (mod(r_mod, p_mod) == 0)
                r_mod = r_mod / p_mod;
            end
        end
        if r_mod == 1
            break;
        end
        m_mod = m_mod + 1;
    end
    m2_mod = m_mod;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x = x';
    X_mod = fft(x,m2_mod,1);
    Y_mod = fft(y_mod,m2_mod,1);
    % X_mod = X_mod(1:20);
    % Y_mod = Y_mod(1:20);
    produto = X_mod.*conj(Y_mod);
    if isreal(x) && isreal(y_mod)
        c1_mod = ifft(produto,[],1,'symmetric');
    else
        c1_mod = ifft(produto,[],1);
    end
    %   IFFT(X,[],DIM) or IFFT(X,N,DIM) is the inverse discrete Fourier
    %   transform of X across the dimension DIM.
    %
    %   IFFT(..., 'symmetric') causes IFFT to treat X as conjugate symmetric
    %   along the active dimension.  This option is useful when X is not exactly
    %   conjugate symmetric merely because of round-off error.  See the
    %   reference page for the specific mathematical definition of this
    %   symmetry.
    % Keep only the lags we want and move negative lags before positive
    % lags.
    c_mod2 = [c1_mod(m2_mod - mxl_mod + (1:mxl_mod)); c1_mod(1:mxl_mod+1)];
    tlmXcorrResult_new = c_mod2';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % disp(isequaln(tlmXcorrResult, tlmXcorrResult_new))
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tlmXcorrResult = tlmXcorrResult_new;
    % Find all starting points off all preamble like patterns =================
    clear index
    clear index2
    
    xcorrLength = (length(tlmXcorrResult) +  1) /2;
    
    % list_1 = abs(tlmXcorrResult(4500:8999));
    % index_test = find(list_1 > 6)';
    % list_2 = list_1(1:50);
    % 
    % for i = 1:length(list_1)
    %     if list_1(i) > 6
    %         disp([i, list_1(i)])
    %     end
    % 
    % end

    

    % index_test2 = find(list_2 > 6)';
    
    %--- Find at what index/ms the preambles start ------------------------
     list_preamble = abs(tlmXcorrResult(xcorrLength : xcorrLength * 2 - 1));
    if settings.DPE_cohInt == 1

        index_preamble = find(list_preamble > 153)' + searchStartOffset;

    else % Coherent integration is higher than 1 ms

        % index_preamble = find(...
        % abs(tlmXcorrResult(xcorrLength : xcorrLength * 2 - 1)) > 6)' + ...
        % searchStartOffset;
        index_preamble = find(list_preamble > 6)' + searchStartOffset;

        % index_preamble = find(abs(tlmXcorrResult(xcorrLength : xcorrLength * 2 - 1)) > 6)' + searchStartOffset;
        
     

    end
    % Analyze detected preamble like patterns ================================
    for i = 1:size(index_preamble) % For each occurrence
    
        %--- Find distances in time between this occurrence and the rest of
        %preambles like patterns. If the distance is 6000 milliseconds (one
        %subframe), the do further verifications by validating the parities
        %of two GPS words


        % disp(i)
        index_preamble2 = index_preamble - index_preamble(i);
        % disp(index_preamble2)
        % 
        if (~isempty(find(index_preamble2 == 6000/settings.DPE_cohInt, 1)))
            %=== Re-read bit values for preamble verification ==============
            % Preamble occurrence is verified by checking the parity of
            % the first two words in the subframe. Now it is assumed that
            % bit boundaries a known. Therefore the bit values over 20ms are
            % combined to increase receiver performance for noisy signals.
            % in Total 62 bits mast be read :
            % 2 bits from previous subframe are needed for parity checking;
            % 60 bits for the first two 30bit words (TLM and HOW words).
            % The index is pointing at the start of TLM word.
            int_1 = index_preamble(i)-40/settings.DPE_cohInt;
            int_2 = (index_preamble(i) + (20 * 60) / settings.DPE_cohInt);
            int_3 = int_1: int_2-1;
            bits_100 = I_P_InputBits(int_3)';

            bits = I_P_InputBits(index_preamble(i)-40/settings.DPE_cohInt...
                : (index_preamble(i) + (20 * 60) / settings.DPE_cohInt)-1)';
            %--- Combine the 20 values of each bit ------------------------
            if settings.DPE_cohInt < 20
                 bits = reshape(bits, 20 / settings.DPE_cohInt, ...
                     ((size(bits, 1) / 20)* settings.DPE_cohInt));
                 bits = sum(bits,1); 
            else % Coherent integration is 20 ms
                 bits=bits';
            end
            % Now threshold and make it -1 and +1
            bits(bits > 0)  = 1;
            bits(bits <= 0) = -1;
            
            %--- Check the parity of the TLM and HOW words ----------------
            if (navPartyChk(bits(1:32)) ~= 0) && ...
                    (navPartyChk(bits(31:62)) ~= 0)
                % Parity was OK. Record the preamble start position. Skip
                % the rest of preamble pattern checking for this channel
                % and process next channel.

                subFrameStart_new = index(i);
                break;
            end % if parity is OK ...
            % disp('aqui 56')
            % return
            

        end % if (~isempty(find(index2 == 6000)))
    end % for i = 1:size(index)
    % return
    % disp('aqui')
    % return
    % Exclude channel from the active channel list if no valid preamble was
    % detected
    if subFrameStart_new == inf
        disp('Could not find valid preambles in channel! ');
        return
    end
    
    %% Decode ephemerides ===============================================
    %=== Convert tracking output to navigation bits =======================
    %--- Copy 5 sub-frames long record from tracking output ---------------
    navBitsSamples = I_P_InputBits(subFrameStart_new - 20/settings.DPE_cohInt: ...
        subFrameStart_new + (1500 * 20)/settings.DPE_cohInt -1)';
    %--- Group every 20 values of bits into columns ------------------------
    
    if settings.DPE_cohInt<20
        navBitsSamples = reshape(navBitsSamples, ...
            20/settings.DPE_cohInt, ...
            (size(navBitsSamples, 1) / 20)* settings.DPE_cohInt);
    else
        % If PDI = 20 ms, the bits do not need to be grouped again
        navBitsSamples=navBitsSamples';
    end
    %--- Sum all samples in the bits to get the best estimate -------------
    navBits = sum(navBitsSamples,1); % Does nothing when PDI = 20 ms
    %--- Now threshold and make 1 and 0 -----------------------------------
    % The expression (navBits > 0) returns an array with elements set to 1
    % if the condition is met and set to 0 if it is not met.
    navBits = (navBits > 0);
    
    %--- Convert from decimal to binary -----------------------------------
    % The function ephemeris expects input in binary form. In Matlab it is
    % a string array containing only "0" and "1" characters.
    navBitsBin = dec2bin(navBits);
    %=== Decode ephemerides and TOW of the first sub-frame ================
    [eph_new, TOW_new] = ephemeris(navBitsBin(2:1501)', navBitsBin(1));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % function [eph_new_mod, TOW_new_mod] = ephemeris_mod(nav_bits1, nav_bits2)
    % Check if there is enough data ==========================================
    nav_bits1 = navBitsBin(2:1501)';
    nav_bits2 = navBitsBin(1);
    if length(nav_bits1) < 1500
        error('The parameter BITS must contain 1500 bits!');  
    end
    % Pi used in the GPS coordinate system
    gpsPi = 3.1415926535898; 
    %% Decode all 5 sub-frames ================================================
    for i = 1:5
        %--- "Cut" one sub-frame's bits ---------------------------------------
        subframe = nav_bits1(300*(i-1)+1 : 300*i);
        %--- Correct polarity of the data bits in all 10 words ----------------
        for j = 1:10
            [subframe(30*(j-1)+1 : 30*j)] = ...
                checkPhase(subframe(30*(j-1)+1 : 30*j), nav_bits2);

            nav_bits2 = subframe(30*j);
        end

        %--- Decode the sub-frame id ------------------------------------------
        % For more details on sub-frame contents please refer to GPS IS.
        subframeID = bin2dec(subframe(50:52));

        %--- Decode sub-frame based on the sub-frames id ----------------------
        % The task is to select the necessary bits and convert them to decimal
        % numbers. For more details on sub-frame contents please refer to GPS
        % ICD (IS-GPS-200D).
        switch subframeID
            case 1  %--- It is subframe 1 -------------------------------------
                % It contains WN, SV clock corrections, health and accuracy
                eph_new_mod.weekNumber  = bin2dec(subframe(61:70)) + 1024;
                eph_new_mod.accuracy    = bin2dec(subframe(73:76));
                eph_new_mod.health      = bin2dec(subframe(77:82));
                eph_new_mod.T_GD        = twosComp2dec(subframe(197:204)) * 2^(-31);
                eph_new_mod.IODC        = bin2dec([subframe(83:84) subframe(197:204)]);
                eph_new_mod.t_oc        = bin2dec(subframe(219:234)) * 2^4;
                eph_new_mod.a_f2        = twosComp2dec(subframe(241:248)) * 2^(-55);
                eph_new_mod.a_f1        = twosComp2dec(subframe(249:264)) * 2^(-43);
                eph_new_mod.a_f0        = twosComp2dec(subframe(271:292)) * 2^(-31);
                eph_new_mod.idValid(1)  = 1;

            case 2  %--- It is subframe 2 -------------------------------------
                % It contains first part of ephemeris parameters
                eph_new_mod.IODE_sf2    = bin2dec(subframe(61:68));
                eph_new_mod.C_rs        = twosComp2dec(subframe(69: 84)) * 2^(-5);
                eph_new_mod.deltan      = ...
                    twosComp2dec(subframe(91:106)) * 2^(-43) * gpsPi;
                eph_new_mod.M_0         = ...
                    twosComp2dec([subframe(107:114) subframe(121:144)]) ...
                    * 2^(-31) * gpsPi;
                eph_new_mod.C_uc        = twosComp2dec(subframe(151:166)) * 2^(-29);
                eph_new_mod.e           = ...
                    bin2dec([subframe(167:174) subframe(181:204)]) ...
                    * 2^(-33);
                eph_new_mod.C_us        = twosComp2dec(subframe(211:226)) * 2^(-29);
                eph_new_mod.sqrtA       = ...
                    bin2dec([subframe(227:234) subframe(241:264)]) ...
                    * 2^(-19);
                eph_new_mod.t_oe        = bin2dec(subframe(271:286)) * 2^4;
                eph_new_mod.idValid(1)  = 2;

            case 3  %--- It is subframe 3 -------------------------------------
                % It contains second part of ephemeris parameters
                eph_new_mod.C_ic        = twosComp2dec(subframe(61:76)) * 2^(-29);
                eph_new_mod.omega_0     = ...
                    twosComp2dec([subframe(77:84) subframe(91:114)]) ...
                    * 2^(-31) * gpsPi;
                eph_new_mod.C_is        = twosComp2dec(subframe(121:136)) * 2^(-29);
                eph_new_mod.i_0         = ...
                    twosComp2dec([subframe(137:144) subframe(151:174)]) ...
                    * 2^(-31) * gpsPi;
                eph_new_mod.C_rc        = twosComp2dec(subframe(181:196)) * 2^(-5);
                eph_new_mod.omega       = ...
                    twosComp2dec([subframe(197:204) subframe(211:234)]) ...
                    * 2^(-31) * gpsPi;
                eph_new_mod.omegaDot    = twosComp2dec(subframe(241:264)) * 2^(-43) * gpsPi;
                eph_new_mod.IODE_sf3    = bin2dec(subframe(271:278));
                eph_new_mod.iDot        = twosComp2dec(subframe(279:292)) * 2^(-43) * gpsPi;
                eph_new_mod.idValid(3)  = 3;

            case 4  %--- It is subframe 4 -------------------------------------
                % Almanac, ionospheric model, UTC parameters.
                % SV health (PRN: 25-32).
                % Not decoded at the moment.

            case 5  %--- It is subframe 5 -------------------------------------
                % SV almanac and health (PRN: 1-24).
                % Almanac reference week number and time.
                % Not decoded at the moment.

        end % switch subframeID ...

    end % for all 5 sub-frames ...

    %% Compute the time of week (TOW) of the first sub-frames in the array ====
    % Also correct the TOW. The transmitted TOW is actual TOW of the next
    % subframe and we need the TOW of the first subframe in this data block
    % (the variable subframe at this point contains bits of the last subframe). 
    % The duration of 5 subframe is 30 s, thus the start time of the first 
    % subframe
    TOW_new_mod = bin2dec(subframe(31:47)) * 6 - 30;
    eph_new_mod.TOW = TOW_new_mod;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




    disp(isequaln(eph_new, eph_new_mod))
    disp(isequaln(TOW_new, TOW_new_mod))
    

    return
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
    eph_new          = [];
    return
end