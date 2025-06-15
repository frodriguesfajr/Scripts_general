function [c,lags] = xcorr_mod(x,varargin)
%XCORR Cross-correlation function estimates.
%   C = XCORR(A,B), where A and B are length M vectors (M>1), returns
%   the length 2*M-1 cross-correlation sequence C. If A and B are of
%   different length, the shortest one is zero-padded. C will be a
%   row vector if A is a row vector, and a column vector if A is a
%   column vector.
%
%   XCORR produces an estimate of the correlation between two random
%   (jointly stationary) sequences:
%          C(m) = E[A(n+m)*conj(B(n))] = E[A(n)*conj(B(n-m))]
%   It is also the deterministic correlation between two deterministic
%   signals.
%
%   C = XCORR(A), where A is a length M vector, returns the length 2*M-1
%   auto-correlation sequence C. The zeroth lag of the output correlation
%   is in the middle of the sequence, at element M.
%
%   C = XCORR(A), where A is an M-by-N matrix (M>1), returns a large matrix
%   with 2*M-1 rows and N^2 columns containing the cross-correlation
%   sequences for all combinations of the columns of A; the first N columns
%   of C contain the delays and cross correlations using the first column
%   of A as the reference, the next N columns of C contain the delays and
%   cross correlations using the second column of A as the reference, and
%   so on.
%
%   C = XCORR(...,MAXLAG) computes the (auto/cross) correlation over the
%   range of lags:  -MAXLAG to MAXLAG, i.e., 2*MAXLAG+1 lags.
%   If missing, default is MAXLAG = M-1.
%
%   [C,LAGS] = XCORR(...)  returns a vector of lag indices (LAGS).
%
%   XCORR(...,SCALEOPT), normalizes the correlation according to SCALEOPT:
%     'biased' - scales the raw cross-correlation by 1/M.
%     'unbiased' - scales the raw correlation by 1/(M-abs(lags)).
%     'normalized' or 'coeff' - normalizes the sequence so that the 
%                               auto-correlations at zero lag are 
%                               identically 1.0.
%     'none' - no scaling (this is the default).
%
%   % Example:
%   % Compute and plot the cross-correlation of two 16-sample
%   % exponential sequences
%
%   N = 16;
%   n = 0:N-1;
%   a = 0.84;
%   b = 0.92;
%   xa = a.^n;
%   xb = b.^n;
%   [r,lags] = xcorr(xa,xb);
%   stem(lags,r)
%
%   See also XCOV, CORRCOEF, CONV, COV.

%   Copyright 1988-2021 The MathWorks, Inc.

%   References:
%     S.J. Orfanidis, "Optimum Signal Processing. An Introduction"
%     2nd Ed. Macmillan, 1988.

narginchk(1,4);
if isnumeric(x)
    if isa(x,'uint64') || isa(x,'int64')
        error(message('MATLAB:xcorr:InvalidInputType'));
    end
elseif ~islogical(x)
    error(message('MATLAB:xcorr:InvalidInputType'));
end

[ySupplied,maxlagInput,scale] = ...
    matlab.internal.math.parseXcorrOptions(varargin{:});

% Transform the input so that computations can be performed on columns.
if size(x,1) == 1 && ~isscalar(x)
    % Shift the leading non-singleton dimension to the fore and make a
    % recursive call. Row vector input becomes a column vector, N-D
    % inputs have leading ones shifted out.
    [x,nshift] = shiftdim(x);
    if nargout == 2
        [c,lags] = xcorr_mod(x,varargin{:});
    else
        c = xcorr_mod(x,varargin{:});
    end
    c = shiftdim(c,-nshift);
    return
end

% Calculate the cross-correlation.
if ySupplied
    % Cross-correlation of two column vectors.
    if ~iscolumn(x)
        error(message('MATLAB:xcorr:MismatchedAB'));
    end
    % disp('aqui 3')
    y = varargin{1}(:); % y was validated to be a vector, make it a column.
    % disp(y)
    maxlagDefault = max(size(x,1),size(y,1)) - 1;
    if isempty(maxlagInput)
        maxlag = maxlagDefault;
    else
        maxlag = maxlagInput;
    end
    
    if isempty(x) || isempty(y)
        if isa(x,'single') || isa(y,'single')
            c = zeros(2*maxlag+1,1,'single');
        else
            c = zeros(2*maxlag+1,1);
        end
        lags = -maxlag:maxlag;
        return;
    else
        % Perform the cross-correlation.
        % c = crosscorr(x,y,maxlag);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % % function c = crosscorr(x,y,maxlag)
        % Compute cross-correlation for vector inputs. Output is clipped based on
        % maxlag but not padded if maxlag >= max(size(x,1),size(y,1)).
        nx = numel(x);
        ny = numel(y);
        m_mod = max(nx,ny);
        maxlagDefault_mod = m_mod-1;
        mxl_mod = min(maxlag,maxlagDefault_mod);
        if ~allfinite(x) || ~allfinite(y)
            c1_mod = conv(x,conj(flip(y)));
            if mxl <= maxlagDefault
                % Clip if maxlag is small. This would be a no-op for larger maxlag.
                % Account for 0 lag corresponding to c1(ny).
                c1_mod = c1_mod(max(1,ny-mxl_mod):(ny+min(mxl_mod,nx-1)));
            end
            % Pad the head or tail if nx >= ny or nx < ny, respectively.
            % Note that we may need to both clip and pad: xcorr(1:5,1:3,3).
            c_mod = zeros(2*mxl+1,1,'like',c1);
            offset = max(0,mxl_mod-ny+1);
            c_mod(offset+(1:numel(c1_mod))) = c1_mod;
        else
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            m = m_mod;
            m = 2*m;
            while true
                r = m;
                for p = [2 3 5 7]
                    while (r > 1) && (mod(r, p) == 0)
                        r = r / p;
                    end
                end
                if r == 1
                    break;
                end
                m = m + 1;
            end
            m2_mod = m;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            X_mod = fft(x,m2_mod,1);
            Y_mod = fft(y,m2_mod,1);
            % disp(Y_mod)
            if isreal(x) && isreal(y)
                c1_mod = ifft(X_mod.*conj(Y_mod),[],1,'symmetric');
            else
                c1_mod = ifft(X_mod.*conj(Y_mod),[],1);
            end
            % Keep only the lags we want and move negative lags before positive
            % lags.
            c_mod = [c1_mod(m2_mod - mxl_mod + (1:mxl_mod)); c1_mod(1:mxl_mod+1)];
        end
        c = c_mod;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end
% Pad the output with zeros.
if maxlag > maxlagDefault
    zeropad = zeros(maxlag - maxlagDefault,size(c,2),'like',c);
    c = [zeropad; c; zeropad];
end

if nargout == 2
    lags = -maxlag:maxlag;
end

%--------------------------------------------------------------------------