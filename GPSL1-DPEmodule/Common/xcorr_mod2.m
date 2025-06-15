function [c,lags] = xcorr_mod2(x,varargin)
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
        [c,lags] = xcorr_mod2(x,varargin{:});
    else
        c = xcorr_mod2(x,varargin{:});
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
    y = varargin{1}(:); % y was validated to be a vector, make it a column.
    
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
        c = crosscorr_mod(x,y,maxlag);
        
        % Scale the output.
        c = scaleOutput_mod(scale,c,x,y);
    end
else
    % Perform all auto- and cross-correlations of the columns of x.
    narginchk(1,3);
    
    maxlagDefault = size(x,1) - 1;
    if isempty(maxlagInput)
        maxlag = maxlagDefault;
    else
        maxlag = maxlagInput;
    end
    
    if isempty(x)
        [~,n] = size(x);
        if isa(x,'single')
            c = zeros(2*maxlag+1,n^2,'single');
        else
            c = zeros(2*maxlag+1,n^2);
        end
        lags = -maxlag:maxlag;
        return;
    else
        % Peform the auto- and cross-correlations.
        c = autocorr(x,maxlag);
        
        % Scale the output.
        c = scaleOutput_mod(scale,c,x);
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
