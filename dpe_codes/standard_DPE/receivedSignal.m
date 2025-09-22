function x = receivedSignal(sigen, CNo)
% RECEIVEDSIGNAL  Gera o sinal GNSS recebido como soma dos SVs
%                 escalados pelo C/N0 desejado e somados a ruído AWGN.
%
%   x = receivedSignal(sigen, CNo)
%
%   sigen : estrutura produzida por signalGen, com as réplicas filtradas
%   CNo   : vetor [1 x numSV] com C/N0 de cada SV em dB-Hz
%
%   x     : vetor complexo 1 x NsamplesData, sinal composto + ruído + LPF

% --------- Extrai parâmetros da estrutura sigen ------------------------------
x_delay      = sigen.x_delay;       % réplicas atrasadas por SV, já filtradas
NsamplesData = sigen.NsamplesData;  % número total de amostras do bloco
numSV        = sigen.numSV;         % número de SVs simulados
fs           = sigen.fs;            % taxa de amostragem (Hz)
fn           = sigen.fn;            % freq. de corte do LPF (Hz)
order        = sigen.order;         % ordem do FIR

% Pré-alocação para as réplicas escaladas + ruído
x_delay_noise = zeros(numSV, NsamplesData);

% --------- Geração de ruído complexo de potência 1 ---------------------------
% randn(1,N) -> N(0,1); multiplicação por sqrt(1/2) dá variância 0.5
% => parte real e imaginária com variância 0.5 -> ruído complexo de potência 1
noise = ( sqrt(1/2)*randn(1,NsamplesData) + ...
          1i*sqrt(1/2)*randn(1,NsamplesData) );

% --------- Escala cada SV conforme o C/N0 desejado --------------------------
for kSV = 1:numSV
    if CNo(kSV) < 100  % limiar arbitrário: se <100 dB-Hz, aplica escala
        % Relação usada:
        %   A^2 * fs / 2 = 10^(C/N0/10)
        % -> A = sqrt( 2 * 10^(C/N0/10) / fs )
        A = sqrt( 2 * 10^( CNo(kSV)/10 ) / fs );

        % Multiplica a réplica atrasada pelo fator de amplitude
        x_delay_noise(kSV,:) = A * x_delay(kSV,:);
    else
        % Se C/N0 >= 100 dB-Hz (caso especial), mantém amplitude original
        x_delay_noise(kSV,:) = x_delay(kSV,:);
    end
end

% --------- Soma de todos os SVs + ruído -------------------------------------
% Sinal recebido é a soma das réplicas de todos os satélites + ruído AWGN
x = sum( x_delay_noise, 1 ) + noise;

% --------- Filtro passa-baixas de front-end ---------------------------------
% Mesmo filtro FIR usado em signalGen, para simular front-end RX
wn = pi*fn/(pi*fs/2);       % normalização para fir1 (equivale a 2*fn/fs)
h  = fir1(order, wn);       % FIR linear de ordem 'order'
x  = filtfilt(h, 1, x);     % filtfilt -> zero-fase (sem atraso de grupo)
end
