function dop = computeDOP(SatPosition, UserPosition, W)
% Computa HDOP, VDOP, PDOP, TDOP, GDOP e cond(G).
% Regime: ECEF, LOS unitários do usuário para cada SV.

    if nargin < 3 || isempty(W), W = eye(size(SatPosition,1)); end

    M  = size(SatPosition,1);
    u  = zeros(M,3);
    for i=1:M
        v   = SatPosition(i,:) - UserPosition;
        u(i,:)= v / norm(v);
    end
    G = [u, ones(M,1)];            % matriz de geometria (M x 4)

    % DOP com/sem pesos
    Q = inv(G.'*W*G);

    dop.HDOP = sqrt(Q(1,1)+Q(2,2));
    dop.VDOP = sqrt(Q(3,3));
    dop.PDOP = sqrt(Q(1,1)+Q(2,2)+Q(3,3));
    dop.TDOP = sqrt(Q(4,4));
    dop.GDOP = sqrt(trace(Q));
    dop.condG = cond(G);           % número de condição (quanto menor, melhor)
end


