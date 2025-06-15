function c = crosscorr_mod(x,y,maxlag)
% Compute cross-correlation for vector inputs. Output is clipped based on
% maxlag but not padded if maxlag >= max(size(x,1),size(y,1)).
nx = numel(x);
ny = numel(y);
m = max(nx,ny);
maxlagDefault = m-1;
mxl = min(maxlag,maxlagDefault);
if ~allfinite(x) || ~allfinite(y)
    c1 = conv(x,conj(flip(y)));
    if mxl <= maxlagDefault
        % Clip if maxlag is small. This would be a no-op for larger maxlag.
        % Account for 0 lag corresponding to c1(ny).
        c1 = c1(max(1,ny-mxl):(ny+min(mxl,nx-1)));
    end
    % Pad the head or tail if nx >= ny or nx < ny, respectively.
    % Note that we may need to both clip and pad: xcorr(1:5,1:3,3).
    c = zeros(2*mxl+1,1,'like',c1);
    offset = max(0,mxl-ny+1);
    c(offset+(1:numel(c1))) = c1;
else
    m2 = findTransformLength_mod(m);
    X = fft(x,m2,1);
    Y = fft(y,m2,1);
    if isreal(x) && isreal(y)
        c1 = ifft(X.*conj(Y),[],1,'symmetric');
    else
        c1 = ifft(X.*conj(Y),[],1);
    end
    % Keep only the lags we want and move negative lags before positive
    % lags.
    c = [c1(m2 - mxl + (1:mxl)); c1(1:mxl+1)];
end

%--------------------------------------------------------------------------