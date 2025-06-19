function [CRLB_xyz, sigma_xyz] = computeCRLB(satPos, recPos, sigma_pseudorange)
    % Transpõe satPos se necessário (espera Nx3)
    if size(satPos, 1) == 3 && size(satPos, 2) ~= 3
        satPos = satPos';
    end

    N = size(satPos, 1);

    % Verifica número mínimo de satélites
    if N < 4
        warning('Poucos satélites (%d) para cálculo estável do CRLB. Retornando NaN.', N);
        CRLB_xyz = NaN(3);
        sigma_xyz = NaN(3,1);
        return
    end

    H = zeros(N, 3);
    recPos = recPos(:)';

    for i = 1:N
        diff = satPos(i,:) - recPos;
        H(i,:) = diff / norm(diff);
    end

    Q = sigma_pseudorange^2 * eye(N);
    F = H' / Q * H;

    % Verifica condicionamento da matriz Fisher
    if rcond(F) < 1e-12
        warning('Matriz Fisher mal condicionada. Usando pseudo-inversa.');
        CRLB_xyz = pinv(F);
    else
        CRLB_xyz = inv(F);
    end

    sigma_xyz = sqrt(diag(CRLB_xyz));
end
