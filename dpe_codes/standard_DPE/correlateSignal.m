function r = correlateSignal(sigen, x)
% CORRELATESIGNAL  Correlação circular via FFT com integração não-coerente.
%   r = correlateSignal(sigen, x)
%
% ENTRADAS
%   sigen.numSV                -> número de SVs
%   sigen.fft_local(kSV, :)    -> FFT da réplica local (tamanho N) do SV k
%   sigen.NsamplesLocal (N)    -> nº de amostras por bloco coerente (1 ms)
%   sigen.NonCoherentIntegrations (NC) -> nº de blocos não-coerentes a acumular
%   x (1 x (N*NC))             -> sinal recebido (soma dos SVs + ruído), contínuo em tempo
%
% SAÍDA
%   r (numSV x N)              -> energia de correlação |R_k[n]|^2 acumulada em NC blocos
%
% OBS: Implementa correlação circular:
%      R_k = IFFT( FFT(seg) .* conj(FFT(codigo_local)) ),
%      e acumula não-coerentemente: r += |R_k|^2.

numSV      = sigen.numSV;
fft_local  = sigen.fft_local;        % [numSV x N]  (pré-computado em signalGen)
N          = sigen.NsamplesLocal;    % amostras por bloco coerente (ex.: 1 ms)
NC         = sigen.NonCoherentIntegrations;

% Pré-aloca saída: uma “curva de correlação” por SV (acumulada em NC blocos)
r = zeros(numSV, N);

for kSV = 1:numSV
    % Para cada SV, somamos a energia de correlação de NC blocos consecutivos
    for idx_nc = 1:NC
        % Seleciona o bloco coerente idx_nc (tamanho N) dentro de x
        % ATENÇÃO: assume que length(x) == N*NC e que os blocos são contíguos, sem overlap.
        seg = x(1, N*(idx_nc-1)+1 : N*idx_nc);

        % Correlação circular no domínio da frequência:
        % - FFT(seg, N): espectro do bloco recebido
        % - conj(fft_local(kSV,:)): espectro conjugado da réplica local do SV k
        % - IFFT(...) retorna a correlação circular tempo-discreto
        Rk  = ifft( fft(seg, N) .* conj(fft_local(kSV,:)) );

        % Integração não-coerente: acumula energia do pico (e arredores)
        % |Rk|^2 remove a fase; soma em 'r' agrega os NC blocos
        r(kSV,:) = r(kSV,:) + abs(Rk).^2;
    end
end
end
