function [satPositions, satClkCorr] = satpos(transmitTime, prnList,eph)
% SATPOS Calcula as coordenadas X, Y, Z ECEF (Earth-Centered, Earth-Fixed) 
% e a correção dos relógios dos satélites GPS no tempo de transmissão 
% TRANSMITTIME, com base nos efemérides EPH. As coordenadas são calculadas
% para cada satélite na lista PRNLIST. com base nas efemérides
%
% Entradas:
%   transmitTime  - tempo de transmissão para cada canal do receptor
%   prnList       - lista dos PRNs (identificadores dos satélites)
%   eph           - estrutura contendo as efemérides (parâmetros orbitais 
%                   Keplerianos e correções).
%
% Saídas:
%   satPositions  - posições dos satélites no sistema 
%                   ECEF [X; Y; Z] (em metros )
%   satClkCorr    - correções dos relógios dos satélites (em segundos)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           SoftGNSS v3.0
%--------------------------------------------------------------------------
% Based on Kai Borre 04-09-96
% Copyright (c) by Kai Borre
% Updated by Darius Plausinaitis, Peter Rinder and Nicolaj Bertelsen
%
% CVS record:
% $Id: satpos.m,v 1.1.2.15 2006/08/22 13:45:59 dpl Exp $
%% Inicialização de constantes e variáveis --------------------------------
numOfSatellites = size(prnList, 2); % número de satélites a processar

% Constantes do GPS
gpsPi      = 3.1415926535898;  % Valor de pi usado no sistema GPS
Omegae_dot = 7.2921151467e-5;  % Velocidade angular da Terra [rad/s]
GM         = 3.986005e14;      % Constante gravitacional da Terra [m^3/s^2]
F          = -4.442807633e-10; % Constante para correção relativística

% Inicializa matrizes de resultados =======================================
satClkCorr   = zeros(1, numOfSatellites);
satPositions = zeros(3, numOfSatellites);

%% Processa cada satélite individualmente ---------------------------------
for satNr = 1 : numOfSatellites
    
    prn = prnList(satNr);
    %% Correção inicial do relógio do satélite ----------------------------
    
    %--- Diferença de tempo desde a referência ----------------------------
    dt = check_t(transmitTime(satNr) - eph(prn).t_oc);
    % Calcula o erro do relógio do satélite usando um modelo quadrático
    satClkCorr(satNr) = (eph(prn).a_f2 * dt + eph(prn).a_f1) * dt + ...
                         eph(prn).a_f0 - ...
                         eph(prn).T_GD;
    % (a_f0, a_f1, a_f2).T_GD é a correção de atraso de 
    % grupo (importante para dual-frequency).

    time = transmitTime(satNr) - satClkCorr(satNr);

    %% Cálculo da posição orbital do satélite -----------------------------
    % semi-eixo maior restaurado
    a   = eph(prn).sqrtA * eph(prn).sqrtA;
    % diferença de tempo ao instante de efemérides
    tk  = check_t(time - eph(prn).t_oe);
    % Movimento médio inicial (rad/s)
    n0  = sqrt(GM / a^3);
    % Movimento médio corrigido
    n   = n0 + eph(prn).deltan;
    % Anomalia média atualizada
    M   = eph(prn).M_0 + n * tk;
    % Limitada a [0, 2pi]
    M   = rem(M + 2*gpsPi, 2*gpsPi);

    % Estimativa inicial para anomalia excêntrica
    E   = M;
    
    % --- Anomalia Excêntrica (E) - Solução Iterativa ---------------------
    for ii = 1:10
        E_old   = E;
        E       = M + eph(prn).e * sin(E); % Equação de Kepler
        dE      = rem(E - E_old, 2*gpsPi);

        if abs(dE) < 1.e-12 % Convergência
            break;
        end
    end
    
    % Limitada a [0, 2pi]
    E   = rem(E + 2*gpsPi, 2*gpsPi);

    % Correção relativística
    dtr = F * eph(prn).e * eph(prn).sqrtA * sin(E);

    % anomalia verdadeira e argumentos corrigidos
    nu = atan2(sqrt(1 - eph(prn).e^2) * sin(E), cos(E) - eph(prn).e);
    phi = rem(nu + eph(prn).omega, 2*gpsPi);
    
    %% Correções Harmônicas
    % latitude argumental, raio e inclinação corrigidos
    u = phi + eph(prn).C_uc * cos(2*phi) + eph(prn).C_us * sin(2*phi);
    r = a * (1 - eph(prn).e*cos(E)) + eph(prn).C_rc * cos(2*phi) + ...
        eph(prn).C_rs * sin(2*phi);
    i = eph(prn).i_0 + eph(prn).iDot * tk + ...
        eph(prn).C_ic * cos(2*phi) + eph(prn).C_is * sin(2*phi);
    
    %% Posição no Plano Orbital (xk1, yk1)
    % coordenadas no plano orbital
    xk1 = cos(u)*r;
    yk1 = sin(u)*r;
    
    %% Rotação para ECEF (Considerando a Rotação da Terra)
    % cálculo do ângulo do nó ascendente
    Omega = eph(prn).omega_0 + (eph(prn).omegaDot - Omegae_dot)*tk - ...
        Omegae_dot * eph(prn).t_oe;
    Omega = rem(Omega + 2*gpsPi, 2*gpsPi);

    % transforma para coordenadas ECEF
    xk = xk1 * cos(Omega) - yk1 * cos(i)*sin(Omega);
    yk = xk1 * sin(Omega) + yk1 * cos(i)*cos(Omega);
    zk = yk1 * sin(i);

    satPositions(:, satNr) = [xk; yk; zk];

    %% Corrige novamente o relógio com efeito relativístico incluído
    satClkCorr(satNr) = (eph(prn).a_f2 * dt + eph(prn).a_f1) * dt + ...
                         eph(prn).a_f0 - eph(prn).T_GD + dtr;    
                     
%% Cálculo da velocidade do satélite em ECEF (atualmente não utilizado)
% Esse bloco calcula a velocidade vetorial dos satélites no referencial 
% fixo à Terra (ECEF), o que é útil para algoritmos que requerem doppler,
% dinâmica orbital ou estimativas de velocidade do usuário. 
% Embora não esteja ativo na versão padrão, pode ser útil for necessário
% obter também a velocidade dos satélites além da posição.
% % 2.12 Derivada da anomalia excêntrica (E)
% % A taxa de variação de E é usada para calcular as velocidades angulares
% dE = n / (1 - eph(prn).e * cos(E));
% 
% % 2.13 Derivada do argumento da latitude verdadeira (f)
% % Calcula a taxa de variação do ângulo verdadeiro do satélite na órbita
% dphi = sqrt(1 - eph(prn).e^2) * dE / (1 - eph(prn).e * cos(E));
% 
% % 2.14 e 2.15 Derivadas de u, r e i com correções harmônicas
% 
% % Derivada do argumento da latitude corrigido
% du = dphi + ...
%      2 * dphi * (-eph(prn).C_uc * sin(2*phi) + ...
%                   eph(prn).C_us * cos(2*phi));
% 
% % Derivada do raio corrigido (distância do centro da Terra ao satélite)
% dr = a * eph(prn).e * dE * sin(E) + ...
%      2 * dphi * (-eph(prn).C_rc * sin(2*phi) + ...
%                   eph(prn).C_rs * cos(2*phi));
% 
% % Derivada da inclinação corrigida
% di = eph(prn).iDot + ...
%      2 * dphi * (-eph(prn).C_ic * sin(2*phi) + ...
%                   eph(prn).C_is * cos(2*phi));
% 
% % Derivada do ângulo do nó ascendente corrigido
% dOmega = eph(prn).omegaDot - Omegae_dot;
% 
% % 2.16 Cálculo da velocidade no plano orbital
% 
% % Velocidade ao longo dos eixos do plano orbital
% dxk1 = dr * cos(u) - r * du * sin(u);
% dyk1 = dr * sin(u) + r * du * cos(u);
% 
% % 2.17 Conversão para velocidade no sistema ECEF (Terra fixo)
% 
% % Componente X da velocidade em ECEF
% satVolocity(1, satNr) = -yk * dOmega - ...
%                         (dyk1 * cos(i) - zk * di) * sin(Omega) + ...
%                         dxk1 * cos(Omega);
% 
% % Componente Y da velocidade em ECEF
% satVolocity(2, satNr) = xk * dOmega + ...
%                         (dyk1 * cos(i) - zk * di) * cos(Omega) + ...
%                         dxk1 * sin(Omega);
% 
% % Componente Z da velocidade em ECEF
% satVolocity(3, satNr) = dyk1 * sin(i) + yk1 * di * cos(i);
% 
% %% Correção relativística na taxa de variação do relógio do satélite
% 
% % Correção relativística derivada (taxa de variação)
% dtrRat = F * eph(prn).e * eph(prn).sqrtA * cos(E) * dE;
% 
% % Estimativa da taxa de variação do relógio do satélite
% % (normalmente pequena e pode ser negligenciada)
% satClkCorrRat(satNr) = 2 * eph(prn).a_f2 * dt + ...
%                         eph(prn).a_f1 + dtrRat;
                  
   
end % for satNr = 1 : numOfSatellites
