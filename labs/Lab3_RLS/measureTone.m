function P = measureTone(y, fs, fTone)
% MEASURETONE  Power of the tone at frequency fTone in the signal y.
%
%   P = measureTone(y, fs, fTone)
%
%   Correlates y with a unit complex exponential at fTone:
%     a = (1/K) * sum_k  y(k) * exp(-j 2 pi fTone k / fs),   P = |a|^2
%   Only that tone contributes; other tones and noise average out.
K = length(y);
a = 0;
for k = 1:K
    a = a + y(k) * exp(-1j*2*pi*fTone*(k-1)/fs);
end
a = a / K;
P = abs(a)^2;
end
