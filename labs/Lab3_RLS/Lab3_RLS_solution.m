%% SOLUTION: LAB 3: RLS on real array data
% Same reference-driven setup as Lab 2: the SOI's transmitted IQ, corrected
% for the frequency offset, is the reference d(k). RLS solves the
% least-squares problem exactly at every snapshot instead of taking one
% small LMS step.
%
% Part A (fill in): the four RLS lines, RLS vs LMS on the same data, and the
%                   forgetting factor on a long capture.
% Part B (harder):  tracking an interferer that moves, and what the speed costs.
%
% Lines marked TODO are yours. Each TODO line already runs (with a
% placeholder), so you can run the script section by section.

clear; close all;
[cap, geom] = loadLabData(3);
Ny = geom.Ny;  Nz = geom.Nz;  N = geom.N;  dOverLambda = geom.dOverLambda;
for i = 1:numel(cap)
    fprintf('%d: %-5s  %s\n', i, cap(i).name, cap(i).description);
end
c = cap(3);                                  % s1_4: SOI + interferer (Labs 1 and 2)
X = c.X;  K = size(X, 2);
d = makeReference(c);                        % corrected reference, as in Lab 2 A1

%% A1. RLS, four lines per snapshot (slides, Module 3)
%   G(k) = P x / ( gamma + x^H P x )          gain vector
%   e(k) = d(k) - w^H x                       error against the reference
%   w    = w + G(k) e*(k)                     weight update
%   P    = ( P - G(k) x^H P ) / gamma         running inverse covariance
% Start: w = 0, P = I / delta with a small delta.
gamma = 0.999;
trR = 0;
for k = 1:K
    trR = trR + real(X(:, k)' * X(:, k)) / K;
end
delta = 0.01 * trR / N;                      % small compared with one element's power
P = eye(N) / delta;
w = zeros(N, 1);
eR = zeros(1, K);
for k = 1:K
    x = X(:, k);
    G = (P * x) / (gamma + x' * P * x);
    eR(k) = d(k) - w' * x;
    w = w + G * conj(eR(k));
    P = (P - G * (x' * P)) / gamma;
end
% checkAnswer('w_rls', w);
y = zeros(1, K);
for k = 1:K
    y(k) = w' * X(:, k);
end
fprintf('\nRLS:  SOI %.1f dB, INT %.1f dB, SIR %.1f dB\n', 10*log10(measureTone(y, c.fs, c.fTone(1))), ...
        10*log10(measureTone(y, c.fs, c.fTone(2))), 10*log10(measureTone(y, c.fs, c.fTone(1)) / measureTone(y, c.fs, c.fTone(2))));

%% A2. RLS vs LMS on the same data
eL = runLMS(X, d, 0.2 / trR);                % Lab 2's LMS, same step size
figure('Name', 'A2: RLS vs LMS');
semilogy(1:K, movingAverage(abs(eL).^2, 10), 1:K, movingAverage(abs(eR).^2, 10)); grid on;
xlim([1 300]); ylim([1e-3 10]); hold on; plot([2*N 2*N], [1e-3 10], 'k:'); hold off;
xlabel('snapshot k'); ylabel('|e(k)|^2, averaged over 10 snapshots');
legend('LMS, \mu = 0.2/tr(R)', 'RLS, \gamma = 0.999', 'k = 2N', 'Location', 'northeast');
title('First 300 snapshots');
% Questions: about how many snapshots does each need? Compare RLS with the
% rule of thumb 2N = 32.

%% A3. The forgetting factor on a long capture (Setup 3, 8184 snapshots)
% Run RLS with gamma = 1, 0.999 and 0.99 on capture 1 and compare.
% Note: this capture is stitched from 16384-sample hardware blocks with
% short gaps, so the SOI's phase jumps every ~1000 snapshots (Lab 2, B3).
c3 = cap(1);  X3 = c3.X;  K3 = size(X3, 2);  d3 = makeReference(c3);
gammas = [1 0.999 0.99];
figure('Name', 'A3: forgetting factor');
for g = 1:numel(gammas)
    e3 = runRLS(X3, d3, gammas(g));
    semilogy(1:K3, movingAverage(abs(e3).^2, 50)); hold on;
end
hold off; grid on; ylim([1e-3 10]); xlabel('snapshot k'); ylabel('|e(k)|^2, averaged over 50 snapshots');
legend('\gamma = 1 (remember everything)', '\gamma = 0.999 (memory ~1000)', '\gamma = 0.99 (memory ~100)', 'Location', 'northwest');
title('RLS with three forgetting factors');
% Questions: with gamma = 1, why does the error creep up and recover so slowly
% after each seam? Memory is about 1/(1 - gamma) snapshots. Which gamma
% would you pick here, and what does a shorter memory cost?

%% ===================== PART B: harder =====================
%% B1. Tracking an interferer that moves
% Captures 2, 3, 4 have the interferer at logged az = +50, +30, +10.
% Play them back to back, as if the interferer jumped twice, and keep one
% RLS running across all three. Compare gamma = 1 with gamma = 0.99.
% (These are separate captures, so the SOI's phase also jumps at each seam.)
seg = [2 3 4];
figure('Name', 'B1: tracking');
for g = [1 0.99]
    eAll = [];
    P = eye(N) / delta;  w = zeros(N, 1);
    for sIdx = seg
        cs = cap(sIdx);  Xs = cs.X;  ds = makeReference(cs);
        for k = 1:size(Xs, 2)
            x = Xs(:, k);
            G = (P * x) / (g + x' * P * x);
            ek = ds(k) - w' * x;
            w = w + G * conj(ek);
            P = (P - G * (x' * P)) / g;
            eAll(end+1) = ek;
        end
    end
    semilogy(1:numel(eAll), movingAverage(abs(eAll).^2, 25)); hold on;
end
Kseg = size(cap(2).X, 2);
plot([Kseg Kseg], [1e-3 10], 'k:', [2*Kseg 2*Kseg], [1e-3 10], 'k:'); hold off;
grid on; ylim([1e-3 10]); xlabel('snapshot k'); ylabel('|e(k)|^2, averaged over 25');
legend('\gamma = 1', '\gamma = 0.99', 'interferer jumps', 'Location', 'northwest');
title('Interferer at +50, then +30, then +10 deg');

%% B2. What the speed costs
% Time one pass of LMS and of RLS over the same capture. Per snapshot LMS
% needs about 2N multiplies, RLS about N^2 (slides, Module 3).
tic; runLMS(X3, d3, 0.2 / trR); tLMS = toc;
tic; runRLS(X3, d3, 0.999); tRLS = toc;
fprintf('\nTime per snapshot:  LMS %.2f us,  RLS %.2f us  (ratio %.1f)\n', 1e6*tLMS/K3, 1e6*tRLS/K3, tRLS/tLMS);
% Question: what would this ratio become for a 64-element array?

%% Local functions
function d = makeReference(c)
% Lab 2, A1: the transmitted SOI, scaled to unit amplitude and shifted to
% the frequency at which it was received.
dTx = c.dTx / mean(abs(c.dTx));
df = c.fTone(1) - c.fTxNominal(1);
d = zeros(1, numel(dTx));
for k = 1:numel(dTx)
    d(k) = dTx(k) * exp(1j*2*pi*df*(k-1)/c.fs);
end
end

function e = runLMS(X, d, mu)
% Lab 2's LMS loop; returns the error at every snapshot.
[N, K] = size(X);
w = zeros(N, 1);  e = zeros(1, K);
for k = 1:K
    e(k) = d(k) - w' * X(:, k);
    w = w + mu * X(:, k) * conj(e(k));
end
end

function e = runRLS(X, d, gamma)
% The A1 loop as a function; returns the error at every snapshot.
[N, K] = size(X);
trR = 0;
for k = 1:K
    trR = trR + real(X(:, k)' * X(:, k)) / K;
end
P = eye(N) / (0.01 * trR / N);  w = zeros(N, 1);  e = zeros(1, K);
for k = 1:K
    x = X(:, k);
    G = (P * x) / (gamma + x' * P * x);
    e(k) = d(k) - w' * x;
    w = w + G * conj(e(k));
    P = (P - G * (x' * P)) / gamma;
end
end

function m = movingAverage(x, L)
m = zeros(size(x));
for k = 1:numel(x)
    m(k) = mean(x(max(1, k-L+1):k));
end
end
