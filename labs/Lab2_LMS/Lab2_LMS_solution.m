%% SOLUTION: LAB 2: LMS on real array data
% Same array as Lab 1. The SOI (Tx1) is a tone, and we know exactly what
% was transmitted: the transmitter's IQ samples are stored with each capture.
% LMS uses that known waveform as its reference d(k) and learns weights that
% make the array output look like it. Nobody tells it where the interferer is.
%
% Part A (fill in): build the reference from the transmitted IQ, run LMS,
%                   watch the beam learn, measure what it achieved.
% Part B (harder):  skip the frequency correction, choose the step size,
%                   and three interferers at once.
%
% Lines marked TODO are yours. Each TODO line already runs (with a
% placeholder), so you can run the script section by section.

clear; close all;
[cap, geom] = loadLabData(2);
Ny = geom.Ny;  Nz = geom.Nz;  N = geom.N;  dOverLambda = geom.dOverLambda;
for i = 1:numel(cap)
    fprintf('%d: %-5s  %s\n', i, cap(i).name, cap(i).description);
end
c = cap(1);                                  % s1_4: SOI + interferer (same scene as Lab 1)
X = c.X;  K = size(X, 2);

%% A1. The reference: what the transmitter sent
% c.dTx is the SOI's transmitted IQ, shifted and decimated exactly like
% the received data. Compare its spectrum with what element 1 received.
dTx = c.dTx / mean(abs(c.dTx));              % scale to unit amplitude
figure('Name', 'A1: reference vs received');
[f, Pr] = spectrumdB(X(1, :), c.fs);  [~, Pt] = spectrumdB(dTx, c.fs);
plot(f/1e6, Pr - max(Pr), f/1e6, Pt - max(Pt)); grid on; xlim([-4 4]); ylim([-60 5]);
xlabel('frequency (MHz)'); ylabel('dB (each normalized to its peak)'); labelTones(c);
legend('received, element 1', 'transmitted SOI (reference)', 'Location', 'southwest');
fprintf('SOI: transmitted at %.1f kHz, received at %.1f kHz -> offset %.1f kHz\n', ...
        c.fTxNominal(1)/1e3, c.fTone(1)/1e3, (c.fTone(1) - c.fTxNominal(1))/1e3);
% The transmitter and receiver run from separate oscillators, so the SOI
% arrives slightly off its nominal frequency. Correct the reference:
%   d(k) = dTx(k) * exp( j 2 pi df (k-1) / fs ),   df = fReceived - fTransmitted
df = c.fTone(1) - c.fTxNominal(1);
d = zeros(1, K);
for k = 1:K
    d(k) = dTx(k) * exp(1j*2*pi*df*(k-1)/c.fs);
end
% checkAnswer('d_first10', d(1:10));

%% A2. LMS, three lines per snapshot (slides, Module 2)
%   y(k) = w^H x(k)
%   e(k) = d(k) - y(k)
%   w    = w + mu x(k) e*(k)
% Step size rule:  mu < 2 / (total received power) = 2 / trace(R).
% We use mu = 0.2 / trace(R). (trace(R) = sum of the elements' powers.)
trR = 0;
for k = 1:K
    trR = trR + real(X(:, k)' * X(:, k)) / K;
end
mu = 0.2 / trR;
w = zeros(N, 1);
y = zeros(1, K);  e = zeros(1, K);
s_soi = steeringVector(c.txAz(1), c.txEl(1), Ny, Nz, dOverLambda);
azScan = -90:1:90;
figure('Name', 'A2: LMS learning (live)');
subplot(2, 1, 1); hP = plot(azScan, zeros(size(azScan))); grid on; ylim([-40 10]);
hold on; plot([c.txAz(1) c.txAz(1)], [-40 10], 'g--', [c.txAz(2) c.txAz(2)], [-40 10], 'r--'); hold off;
xlabel('azimuth (deg)'); ylabel('|AF| (dB)'); title('Beam pattern of the current weights (cut at the SOI elevation)');
subplot(2, 1, 2); hE = semilogy(1:K, nan(1, K)); grid on; xlim([1 K]); ylim([1e-3 10]);
xlabel('snapshot k'); ylabel('|e(k)|^2'); title('Squared error');
for k = 1:K
    y(k) = w' * X(:, k);
    e(k) = d(k) - y(k);
    w = w + mu * X(:, k) * conj(e(k));
    if mod(k, 25) == 0 || k == K             % redraw every 25 snapshots
        AF = zeros(size(azScan));
        for i = 1:numel(azScan)
            AF(i) = w' * steeringVector(azScan(i), c.txEl(1), Ny, Nz, dOverLambda);
        end
        set(hP, 'YData', 20*log10(abs(AF) / max(abs(w' * s_soi), 1e-9) + 1e-6));
        set(hE, 'YData', [abs(e(1:k)).^2, nan(1, K - k)]);
        drawnow; pause(0.06);                    % ~5x slower so the beam is watchable; raise to slow more
    end
end
% checkAnswer('w_final', w);
% The pattern is drawn relative to its gain toward the logged SOI direction.
% Green dashed: SOI (logged). Red dashed: interferer (logged).

%% A3. What did LMS achieve?
% Apply the final weights to the whole capture and compare with element 1
% and with Lab 1's conventional beam.
yFinal = zeros(1, K);  yConv = zeros(1, K);
for k = 1:K
    yFinal(k) = w' * X(:, k);
    yConv(k) = (s_soi / N)' * X(:, k);
end
fprintf('\n%-14s %9s %9s %9s\n', '', 'SOI (dB)', 'INT (dB)', 'SIR (dB)');
printRow('element 1', X(1, :), c);
printRow('conventional', yConv, c);
printRow('LMS', yFinal, c);
figure('Name', 'A3: spectra');
[f, P1] = spectrumdB(X(1, :), c.fs);  [~, Pc] = spectrumdB(yConv, c.fs);  [~, Pl] = spectrumdB(yFinal, c.fs);
subplot(3, 1, 1); plot(f/1e6, P1); grid on; ylim([-80 40]); labelTones(c); title('element 1'); ylabel('dB');
subplot(3, 1, 2); plot(f/1e6, Pc); grid on; ylim([-80 40]); labelTones(c); title('conventional beam'); ylabel('dB');
subplot(3, 1, 3); plot(f/1e6, Pl); grid on; ylim([-80 40]); labelTones(c); title('LMS'); ylabel('dB'); xlabel('frequency (MHz)');
% Questions: how much did the interferer drop compared with the conventional
% beam? LMS was never told where the interferer is, or even that it exists.
% (The LMS output is scaled to match the reference, amplitude 1, so its
% powers sit near 0 dB rather than on the ADC scale.)

%% ===================== PART B: harder =====================
%% B1. Curve ball: skip the frequency correction
% Run LMS with the transmitted IQ as it is (d = dTx), no correction.
w = zeros(N, 1);  e2 = zeros(1, K);
for k = 1:K
    yk = w' * X(:, k);
    e2(k) = dTx(k) - yk;
    w = w + mu * X(:, k) * conj(e2(k));
end
figure('Name', 'B1: unsynchronized reference');
semilogy(1:K, movingAverage(abs(e).^2, 25), 1:K, movingAverage(abs(e2).^2, 25)); grid on; ylim([1e-3 10]);
xlabel('snapshot k'); ylabel('|e(k)|^2, averaged over 25 snapshots');
legend('reference corrected by df', 'reference as transmitted', 'Location', 'southwest');
title('LMS needs a reference that is synchronized with the received SOI');
% Why does a ~120 kHz offset ruin it? (How far does the reference's phase
% drift relative to the SOI over 25 snapshots at 25 MS/s?)

%% B2. Choosing the step size
% Try mu = [0.02 0.2 1 2.5] / trace(R) and compare the learning curves.
% What is the price of a small mu? Of a large one? What happens above 2?
muList = [0.02 0.2 1 2.5] / trR;
figure('Name', 'B2: step size');
for m = 1:numel(muList)
    w = zeros(N, 1);  em = zeros(1, K);
    for k = 1:K
        yk = w' * X(:, k);
        em(k) = d(k) - yk;
        w = w + muList(m) * X(:, k) * conj(em(k));
    end
    semilogy(1:K, movingAverage(abs(em).^2, 25)); hold on;
end
hold off; grid on; ylim([1e-3 1e3]); xlabel('snapshot k'); ylabel('|e(k)|^2, averaged over 25 snapshots');
legend('\mu = 0.02/tr(R)', '\mu = 0.2/tr(R)', '\mu = 1/tr(R)', '\mu = 2.5/tr(R)', 'Location', 'northeast');
title('Learning curves for four step sizes');

%% B3. Three interferers at once (Setup 3): scenario -> live LMS -> result
% Capture 2 (s3_3) has the SOI plus three interferers. The same three LMS
% lines as A2 handle all of them at once.
c3 = cap(2);  X3 = c3.X;  K3 = size(X3, 2);
dTx3 = c3.dTx / mean(abs(c3.dTx));  df3 = c3.fTone(1) - c3.fTxNominal(1);
d3 = dTx3 .* exp(1j*2*pi*df3*(0:K3-1)/c3.fs);          % corrected reference (A1)
trR3 = 0;
for k = 1:K3, trR3 = trR3 + real(X3(:, k)' * X3(:, k)) / K3; end
mu3 = 0.2 / trR3;

% (1) the new scenario, in the spectrum
figure('Name', 'B3.1: three-interferer scenario');
[fS, PS] = spectrumdB(X3(1, :), c3.fs);
plot(fS/1e6, PS); grid on; ylim([-80 40]); labelTones(c3);
xlabel('frequency (MHz)'); ylabel('dB'); title('Element 1: SOI + three interferers');

% (2) LMS learning, live
s_soi3 = steeringVector(c3.txAz(1), c3.txEl(1), Ny, Nz, dOverLambda);
azScan = -90:1:90;
w = zeros(N, 1);  e3 = zeros(1, K3);
figure('Name', 'B3.2: LMS learning (live)');
subplot(2, 1, 1); hP = plot(azScan, zeros(size(azScan))); grid on; ylim([-40 10]); hold on;
for t = find(c3.txActive), if t == 1, col = 'g--'; else, col = 'r--'; end; plot([c3.txAz(t) c3.txAz(t)], [-40 10], col); end
hold off; xlabel('azimuth (deg)'); ylabel('|AF| (dB)'); title('Beam pattern (cut at the SOI elevation)');
subplot(2, 1, 2); hE = semilogy(1:K3, nan(1, K3)); grid on; xlim([1 K3]); ylim([1e-3 10]);
xlabel('snapshot k'); ylabel('|e(k)|^2'); title('Squared error');
for k = 1:K3
    y = w' * X3(:, k);
    e3(k) = d3(k) - y;
    w = w + mu3 * X3(:, k) * conj(e3(k));
    if mod(k, 25) == 0 || k == K3
        AF = zeros(size(azScan));
        for i = 1:numel(azScan), AF(i) = w' * steeringVector(azScan(i), c3.txEl(1), Ny, Nz, dOverLambda); end
        set(hP, 'YData', 20*log10(abs(AF) / max(abs(w' * s_soi3), 1e-9) + 1e-6));
        set(hE, 'YData', [abs(e3(1:k)).^2, nan(1, K3 - k)]);
        drawnow; pause(0.06);
    end
end

% (3) the result: before vs after
yLMS = zeros(1, K3);
for k = 1:K3, yLMS(k) = w' * X3(:, k); end
figure('Name', 'B3.3: result');
subplot(2, 1, 1);
[fB, Pin] = spectrumdB(X3(1, :), c3.fs);  [~, Pout] = spectrumdB(yLMS, c3.fs);
plot(fB/1e6, Pin - max(Pin), fB/1e6, Pout - max(Pout)); grid on; ylim([-80 5]); labelTones(c3);
xlabel('frequency (MHz)'); ylabel('dB (each normalized to its peak)');
legend('element 1 (before)', 'LMS output (after)', 'Location', 'southwest'); title('Before vs after');
subplot(2, 1, 2);
toneBar(c3, {X3(1, :), yLMS}, {'element 1', 'LMS'});
title('Interference cancellation: tone power before vs after');
% Questions: LMS got no angles at all. Did it null all three interferers?
% How many dB did each drop, and did the SOI survive?

%% Local functions
function m = movingAverage(x, L)
% Average of the last L values at each point (shorter at the start).
m = zeros(size(x));
for k = 1:numel(x)
    m(k) = mean(x(max(1, k-L+1):k));
end
end

function printRow(name, y, c)
pS = 10*log10(measureTone(y, c.fs, c.fTone(1)));
pI = 10*log10(measureTone(y, c.fs, c.fTone(2)));
fprintf('%-14s %9.1f %9.1f %9.1f\n', name, pS, pI, pS - pI);
end

function labelTones(c)
nm = {'SOI', 'INT', 'Tx3', 'Tx4'};
if sum(c.txActive) > 2, nm = {'SOI', 'Tx2', 'Tx3', 'Tx4'}; end
yl = ylim;
for t = find(c.txActive)
    text(c.fTone(t)/1e6, yl(2) - 0.08*diff(yl), nm{t}, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
end
end

function toneBar(c, Ys, labels)
% Grouped bar of tone power (dB) at each transmitter frequency, one bar per
% output in Ys. Interferers drop while the SOI stays -> cancellation at a glance.
tones = find(c.txActive);
nm = {'SOI', 'INT', 'Tx3', 'Tx4'};
if sum(c.txActive) > 2, nm = {'SOI', 'Tx2', 'Tx3', 'Tx4'}; end
P = zeros(numel(Ys), numel(tones));
for b = 1:numel(Ys)
    for t = 1:numel(tones)
        P(b, t) = 10*log10(measureTone(Ys{b}, c.fs, c.fTone(tones(t))));
    end
end
bar(P.'); grid on; ylabel('tone power (dB)');
set(gca, 'XTickLabel', nm(tones)); legend(labels, 'Location', 'northeast');
end

function s = steeringVector(az, el, Ny, Nz, dOverLambda)
% Steering vector from Lab 1 (angles in degrees)
u = cosd(el) * sind(az);
v = sind(el);
s = zeros(Ny*Nz, 1);
for ny = 1:Ny
    for nz = 1:Nz
        s((ny-1)*Nz + nz) = exp(-1j*2*pi*dOverLambda*((ny-1)*u + (nz-1)*v));
    end
end
end
