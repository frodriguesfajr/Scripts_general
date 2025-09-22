function pos_est_ls = conv2stepsPVT(r, ls_input)
    % === Entradas ===
    UserPosition        = ls_input.UserPosition(:);   % 3x1
    c                   = ls_input.c;
    numSV               = ls_input.numSV;
    fs                  = ls_input.fs;
    numIts              = ls_input.num2stepsIterations;
    SatPosition         = ls_input.SatPosition;       % Mx3

    % === 1) Medida de pseudo-distância fracionária (1 ms) em metros ===
    [~, maxPos] = max(r,[],2);        % índice do pico
    maxPos      = maxPos - 1;         % zero-based
    tau_hat     = maxPos / fs;        % s
    rho_hat     = c * tau_hat;        % m (mod 1 ms)

    % === 2) Estado inicial [x;y;z;b] com b (m) ===
    xhat = [UserPosition + 10*(2*rand(3,1)-1) ; 0];   % b inicial = 0 m

    % meia-faixa do wrap em metros (± c*1ms/2)
    halfspan_m = c * 1e-3 / 2;

    % === 3) LS iterativo com normal equations (paper) ===
    for it = 1:numIts
        H        = zeros(numSV, 4);
        rho_pred = zeros(numSV, 1);

        for k = 1:numSV
            vec  = SatPosition(k,:).' - xhat(1:3);   % vetor LOS
            dist = norm(vec);
            if dist < eps, dist = eps; end
            u         = vec / dist;                  % u_i
            H(k,1:3)  = -u.';                        % d rho / d p
            H(k,4)    = 1;                           % d rho / d b
            rho_pred(k) = dist + xhat(4);            % m
        end

        % resíduo em metros, com wrap para ± c*Tc/2
        res = wrap_pm(rho_hat - rho_pred, halfspan_m);

        % === atualização LS exatamente como no paper ===
        % delta = (H'H)^{-1} H' res
        N     = (H.' * H);
        g     = (H.' * res);
        delta = inv(N) * g;

        % atualiza estado
        xhat = xhat + delta;

        % (opcional) critério de parada
        if norm(delta(1:3)) < 1e-3 && abs(delta(4)) < 1e-3
            break;
        end
    end

    % retorna só a posição (linha, para bater com seu código)
    pos_est_ls = xhat(1:3).';
end

function y = wrap_pm(x, halfspan)
    % mapeia x para o intervalo [-halfspan, +halfspan]
    y = mod(x + halfspan, 2*halfspan) - halfspan;
end