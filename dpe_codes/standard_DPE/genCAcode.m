function CAcode = genCAcode(PRN)
% GEN CACODE  Gera a sequência de 1023 chips do código C/A (Coarse/Acquisition)
%             do sistema GPS para um dado satélite (PRN).
%
%   CAcode = genCAcode(PRN)
%
% ENTRADA
%   PRN : número de identificação do satélite (1–32 GPS clássico)
%
% SAÍDA
%   CAcode : vetor 1x1023 com os chips do código C/A, em valores ±1
%
%   Implementa o gerador de código C/A padrão IS-GPS-200:
%     - Dois registradores de 10 bits (G1 e G2) com realimentações
%     - Sequências de 1023 bits
%     - O código do satélite PRN é o produto G1 .* G2 deslocado de acordo
%       com uma dupla “tap sequence” (g2shift) específica para cada PRN.

    % --------- Tabela de deslocamentos G2 (em número de chips) -------------
    % Para cada PRN, define quais taps da sequência G2 formam a combinação
    % que, junto com G1, gera o C/A code do satélite correspondente.
    % Esses números seguem a especificação IS-GPS-200.
    g2s = [5,6,7,8,17,18,139,140,141,251,252,254,255,256,257,258,...
           469,470,471,472,473,474,509,512,513,514,515,516,...
           859,860,861,862, ...    % PRNs 1–32
           145,175,52,21,237,235,886,657,634,762,...
           355,1012,176,603,130,359,595,68,386]; % extensões para PRNs >32

    % Deslocamento para o PRN solicitado
    g2shift = g2s(PRN);

    % --------- Geração da sequência G1 (registrador LFSR G1) ---------------
    % Registrador de 10 bits, inicializado com -1 (=chip "1" em notação ±1)
    g1  = zeros(1,1023); 
    reg = -1*ones(1,10);
    for i = 1:1023
        g1(i)   = reg(10);              % saída é o último bit
        saveBit = reg(3)*reg(10);       % realimentação taps 3 e 10
        reg(2:10) = reg(1:9);           % shift para a direita
        reg(1)    = saveBit;            % novo bit na posição 1
    end

    % --------- Geração da sequência G2 (registrador LFSR G2) ---------------
    g2  = zeros(1,1023);
    reg = -1*ones(1,10);
    for i = 1:1023
        g2(i)   = reg(10);              % saída é o último bit
        % realimentação taps 2,3,6,8,9,10 (padrão C/A)
        saveBit = reg(2)*reg(3)*reg(6)*reg(8)*reg(9)*reg(10);
        reg(2:10) = reg(1:9);
        reg(1)    = saveBit;
    end

    % --------- Combinação e deslocamento ------------------------------------
    % Para o PRN desejado, a sequência G2 é deslocada g2shift chips
    % (circularmente) antes de multiplicar por G1.
    g2 = [g2(1023 - g2shift + 1 : 1023), g2(1 : 1023 - g2shift)];

    % Código C/A resultante: produto G1 .* G2, em ±1.
    % Sinal negativo (-) mantém convenção onde '1' lógico = -1 (chip de fase 0).
    CAcode = -(g1 .* g2);
end
