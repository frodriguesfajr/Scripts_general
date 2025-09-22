function ecef = llh2ecef(lat_deg, lon_deg, h_m)
% LLH2ECEF  Converte coordenadas geodésicas (latitude, longitude, altura)
%           no elipsóide WGS-84 para coordenadas cartesianas ECEF [x y z].
%
%   ecef = llh2ecef(lat_deg, lon_deg, h_m)
%
%   lat_deg : latitude geodésica   [graus]  (Norte positivo)
%   lon_deg : longitude geodésica  [graus]  (Leste positivo)
%   h_m     : altura elipsoidal    [m]
%
%   ecef    : vetor 1x3 [x y z] em metros no sistema Earth-Centered Earth-Fixed (ECEF)

    % --- Parâmetros do elipsóide WGS-84 --------------------------------------
    a  = 6378137.0;           % semi-eixo maior (raio equatorial) em metros
    f  = 1/298.257223563;     % achatamento (flattening)
    e2 = f*(2-f);             % primeira excentricidade ao quadrado: e² = f(2−f)

    % --- Conversão de graus para radianos ------------------------------------
    lat = deg2rad(lat_deg);
    lon = deg2rad(lon_deg);

    % --- Cálculo de senos e cossenos -----------------------------------------
    sinlat = sin(lat);  coslat = cos(lat);
    sinlon = sin(lon);  coslon = cos(lon);

    % --- Raio de curvatura na direção do meridiano prime vertical ------------
    % N = a / sqrt(1 - e² sin²(lat))
    % distância do centro ao plano da latitude em função da excentricidade
    N = a / sqrt(1 - e2 * sinlat.^2);

    % --- Conversão final para coordenadas ECEF -------------------------------
    % Fórmulas clássicas:
    %   x = (N + h) * cos(lat) * cos(lon)
    %   y = (N + h) * cos(lat) * sin(lon)
    %   z = (N * (1 - e²) + h) * sin(lat)
    x = (N + h_m)        * coslat * coslon;
    y = (N + h_m)        * coslat * sinlon;
    z = (N * (1 - e2) + h_m) * sinlat;

    % --- Saída: vetor 1x3 ----------------------------------------------------
    ecef = [x, y, z];
end
