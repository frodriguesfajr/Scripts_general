function [SatPosition, SatPRN, azel_deg, PDOP] = randomConstellationPDOP( ...
    UserPosition, M, mode, elevMaskDeg, pdopTarget, maxTries)
% Gera M satélites com PDOP aproximado desejado.
% mode: 'good' (baixa DOP) | 'bad' (alta DOP) | 'custom' (intervalo arbitrário)
% pdopTarget: [PDOPmin PDOPmax], ex: good=[1.5 3.5], bad=[8 20]
% elevMaskDeg: máscara de elevação em graus (ex.: 5)
% maxTries: nº máx de tentativas (ex.: 2000)

if nargin < 3 || isempty(mode),        mode = 'good'; end
if nargin < 4 || isempty(elevMaskDeg), elevMaskDeg = 5; end
if nargin < 5 || isempty(pdopTarget)
    if strcmpi(mode,'good'), pdopTarget = [1.5 3.5]; else, pdopTarget = [8 20]; end
end
if nargin < 6 || isempty(maxTries),    maxTries = 2000; end

% --- Conversão ECEF -> base ENU no usuário ---
x = UserPosition(1); y = UserPosition(2); z = UserPosition(3);
lon = atan2(y, x);
hyp = hypot(x, y);
lat = atan2(z, hyp);

e_hat = [-sin(lon),              cos(lon),             0          ];
n_hat = [-sin(lat)*cos(lon), -sin(lat)*sin(lon),  cos(lat)];
u_hat = [ cos(lat)*cos(lon),  cos(lat)*sin(lon),  sin(lat)];
B     = [e_hat(:) n_hat(:) u_hat(:)];  % 3x3 ENU->ECEF

elMask = deg2rad(elevMaskDeg);

bestPDOP = inf;
bestSat  = [];
bestAzEl = [];

for it = 1:maxTries

    switch lower(mode)
        case 'good'
            % Azimutes estratificados para espalhar no céu
            base = (0:M-1)/M*2*pi;
            jitter = (2*pi/M)*(rand(1,M)-0.5);      % pequeno ruído
            az = mod(base + jitter, 2*pi);

            % Elevações variadas (baixa até alta)
            % Amostra uniforme em sin(el) para não enviesar para o horizonte
            sinEl = sin(elMask) + (1 - sin(elMask))*rand(1,M);
            el = asin(sinEl);
            % garante 1-2 satélites bem altos
            nhigh = max(1, round(0.2*M));
            idxh = randperm(M, nhigh);
            el(idxh) = max(el(idxh), deg2rad(70 + 15*rand(1,nhigh)));

        case 'bad'
            % Azimutes aglomerados num setor estreito
            az0 = 2*pi*rand;
            span = deg2rad(50);                      % setor estreito ~50°
            az   = mod(az0 + (rand(1,M)-0.5)*span, 2*pi);

            % Elevações próximas da máscara
            el = elMask + deg2rad(0.5 + 4.5*rand(1,M));  % 0.5°–5° acima da máscara

        otherwise % 'custom'
            % Sorteio amplo; deixa o pdopTarget fazer o filtro
            az  = 2*pi*rand(1,M);
            sinEl = sin(elMask) + (1 - sin(elMask))*rand(1,M);
            el  = asin(sinEl);
    end

    % Vetor LOS em ENU (az a partir do norte, sentido horário)
    ce = cos(el); sa = sin(az); ca = cos(az);
    v_enu = [ (ce .* sa); (ce .* ca); (sin(el)) ];  % 3xM

    % ENU -> ECEF
    v_ecef = B * v_enu;                    % 3xM (unit)

    % Alcances realistas (usuário -> satélite), 20–27 mil km
    R = (2.0e7 + (2.7e7-2.0e7)*rand(1,M));
    SatPosition_try = UserPosition + (v_ecef .* R).';   % Mx3

    % Calcula PDOP
    PDOP_try = localPDOP(UserPosition, SatPosition_try);
    if isnan(PDOP_try), continue; end

    % Teste alvo
    if PDOP_try >= pdopTarget(1) && PDOP_try <= pdopTarget(2)
        SatPosition = SatPosition_try;
        PDOP        = PDOP_try;
        azel_deg    = [rad2deg(az(:)) rad2deg(el(:))];
        SatPRN      = randperm(32, M);
        return
    end

    % Guarda melhor
    if PDOP_try < bestPDOP && strcmpi(mode,'good')
        bestPDOP = PDOP_try; bestSat = SatPosition_try; bestAzEl = [rad2deg(az(:)) rad2deg(el(:))];
    elseif PDOP_try > bestPDOP && strcmpi(mode,'bad')
        % para 'bad', o "melhor" é o de maior PDOP
        bestPDOP = PDOP_try; bestSat = SatPosition_try; bestAzEl = [rad2deg(az(:)) rad2deg(el(:))];
    elseif isempty(bestSat)
        bestPDOP = PDOP_try; bestSat = SatPosition_try; bestAzEl = [rad2deg(az(:)) rad2deg(el(:))];
    end
end

% Se não atingiu a janela alvo, devolve o melhor que conseguiu
SatPosition = bestSat;
PDOP        = bestPDOP;
azel_deg    = bestAzEl;
SatPRN      = randperm(32, M);
end
