function X = forceFFTSymmetry(X)
% forceFFTSymmetry  A function to force conjugate symmetry on an FFT such that when an
% IFFT is performed the result is a real signal.

% The function has been written to replace MATLAB's ifft(X,'symmetric'), as this function
% is not compatible with MATLAB Coder.

% Licensed under Creative Commons Zero (CC0) so use freely.

XStartFlipped = fliplr(X(2:floor(end/2)));
X(ceil(end/2)+2:end) = real(XStartFlipped) - sqrt(complex(-1))*imag(XStartFlipped);

% Or
% X(ceil(end/2)+2:end) = conj(XStartFlipped);

end