function SatPosition = genSyntheticGPSPositions(SatPRN, UserPosition, R_orbit, elevMaskDeg)
% Gera posições ECEF (m) para os PRNs informados, coerentes com GPS:
% - Satélites colocados numa casca esférica de raio R_orbit (~26.56e6 m)
% - Azimutes bem espaçados (ângulo áureo) e elevações crescentes acima do mask
% - Interseção da LOS do usuário com a esfera orbital → coordenadas ECEF
%
% Parâmetros:
%   SatPRN        : vetor de PRNs (apenas para tamanho/identidade)
%   UserPosition  : [x y z] do usuário em ECEF (m)
%   R_orbit       : raio orbital (m). Default: 26560e3 (GPS)
%   elevMaskDeg   : máscara de elevação (graus). Default: 5
%
% Regime de operação:
% - Terra esférica (aprox.): lat = atan2(z, sqrt(x^2+y^2)), lon = atan2(y,x)
% - Azimute medido a partir do Norte, sentido horário: 0°=N, 90°=E (convenção ENU)
% - Saída é determinística (não usa aleatoriedade)

    if nargin < 3 || isempty(R_orbit),     R_orbit = 26560e3; end
    if nargin < 4 || isempty(elevMaskDeg), elevMaskDeg = 5;    end

    N = numel(SatPRN);

    % --- 1) Lat/Lon (aprox. esférica) do usuário a partir de ECEF ---
    x = UserPosition(1); y = UserPosition(2); z = UserPosition(3);
    lon = atan2(y, x);
    lat = atan2(z, hypot(x,y));  % geocêntrica (boa o suficiente aqui)

    % Rotação ECEF->ENU no usuário; ENU->ECEF é a transposta
    R_ecef2enu = [ -sin(lon)           cos(lon)            0;
                   -sin(lat)*cos(lon) -sin(lat)*sin(lon)  cos(lat);
                    cos(lat)*cos(lon)  cos(lat)*sin(lon)  sin(lat) ];
    R_enu2ecef = R_ecef2enu.';  % transpose

    % --- 2) Gera (az, el) determinísticos e bem espaçados ---
    % Azimutes via "ângulo áureo" para uniformidade em [0,360)
    az_deg = mod((0:N-1) * 137.508, 360);       % graus
    % Elevações distribuídas entre o mask e ~85° (evita zênite exato)
    el_deg = elevMaskDeg + ((1:N)/(N+1)) * (85 - elevMaskDeg);

    % --- 3) Constrói vetores LOS em ENU e leva para ECEF ---
    SatPosition = zeros(N,3);
    r_u   = UserPosition(:);
    r_u2  = dot(r_u, r_u);

    for k = 1:N
        az = deg2rad(az_deg(k));
        el = deg2rad(el_deg(k));

        % LOS unitário em ENU (conv.: az de Norte→Leste, elevação positiva)
        u_enu  = [cos(el)*sin(az);
                  cos(el)*cos(az);
                  sin(el)];
        u_ecef = R_enu2ecef * u_enu;            % unitário em ECEF

        % --- 4) Interseção r_u + s*u = sat, com ||sat|| = R_orbit ---
        b1   = dot(r_u, u_ecef);                % r_u·u
        disc = b1^2 - (r_u2 - R_orbit^2);
        if disc < 0
            % Numérico muito raro; "empurra" para o horizonte
            disc = 0;
        end
        s = -b1 + sqrt(disc);                   % raiz positiva → sat "acima" do usuário
        SatPosition(k,:) = (r_u + s*u_ecef).';
    end
end

