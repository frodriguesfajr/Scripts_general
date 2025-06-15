function [eph_new_mod, TOW_new_mod] = ephemeris_mod(nav_bits1, nav_bits2)
%% Check if there is enough data ==========================================
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
