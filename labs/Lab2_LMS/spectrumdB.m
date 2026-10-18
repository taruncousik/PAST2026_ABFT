function [f, PdB] = spectrumdB(y, fs)
% SPECTRUMDB  Power spectrum of y in dB, scaled so that a tone's peak
% equals the tone's power (same scale as measureTone).
%   [f, PdB] = spectrumdB(y, fs)     f in Hz, from -fs/2 to fs/2
y = y(:);
K = length(y);
win = 0.5 - 0.5*cos(2*pi*(0:K-1)'/K);          % Hann window
Y = fft(y .* win) / sum(win);
P = abs(Y).^2;
PdB = 10*log10(fftshift(P) + 1e-20);
f = ((0:K-1)' - floor(K/2)) * fs / K;
end
