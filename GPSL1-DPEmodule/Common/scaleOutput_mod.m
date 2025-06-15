function c = scaleOutput_mod(scale,c,x,y)
% Scale correlation as specified.
if strcmp(scale,'none')
    return
end
ySupplied = nargin == 4;
m = size(x,1);
if ySupplied && (m ~= size(y,1))
    error(message('MATLAB:xcorr:NoScale'));
end

if strcmp(scale,'biased')
    % Scales the raw cross-correlation by 1/M.
    c = c./m;
elseif strcmp(scale,'unbiased')
    % Scales the raw correlation by 1/(M-abs(lags)).
    L = (size(c,1) - 1)/2;
    scaleUnbiased = (m - abs(-L:L)).';
    scaleUnbiased(scaleUnbiased <= 0) = 1;
    c = c./scaleUnbiased;
else % 'normalized'/'coeff'
    % Normalizes the sequence so that the auto-correlations
    % at zero lag are identically 1.0.
    if ySupplied
        % Compute autocorrelations at zero lag.
        % scale = norm(x)*norm(y) is numerically superior but slower.
        cxx0 = sum(abs(x).^2);
        cyy0 = sum(abs(y).^2);
        scaleCoeffCross = sqrt(cxx0*cyy0);
        c = c./scaleCoeffCross;
    elseif size(c,2) == 1
        % Autocorrelation of a vector. Normalize by c[0].
        mid = (size(c,1) + 1)/2; % row corresponding to zero lag.
        c = c./c(mid);
    else
        % Compute the indices corresponding to the columns that are
        % autocorrelations.
        [~,n] = size(x);
        % Note that size(c,2) = n^2.
        kvec = 1:n+1:n*n; % a 1-by-n row vector
        % kvec is an index vector such that for an n-by-n matrix A,
        % A(kvec(j)) = A(j,j).
        mid = (size(c,1) + 1)/2; % row index corresponding to zero lag
        trow = sqrt(c(mid,kvec)); % a 1-by-n row vector
        tmat = trow.'*trow; % an n-by-n matrix, tmat(i,j) = trow(i)*trow(j)
        scaleCoeffAuto = tmat(:).'; % a 1-by-n^2 row vector
        % The autocorrelations at zero-lag are normalized to one.
        c = c./scaleCoeffAuto;
    end
end
