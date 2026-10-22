%% LAB 4: MVDR on real array data
% Same array and transmitters as Lab 1: Tx1 is the signal of interest (SOI),
% Tx2 is an interferer, both are tones at different frequencies.
%
% Part A (fill in): MVDR from the data itself (SOI + interferer), then MVDR
%                   trained on the interferer alone, side by side; then
%                   where the array thinks everything is.
% Part B (harder):  keeping the SOI with diagonal loading, and how SOI
%                   power changes the damage.
%
% Lines marked TODO are yours. Each TODO line already runs (with a
% placeholder), so you can run the script section by section.

clear; close all;
[cap, geom] = loadLabData(4);
Ny = geom.Ny;  Nz = geom.Nz;  N = geom.N;  dOverLambda = geom.dOverLambda;
for i = 1:numel(cap)
    fprintf('%d: %-5s  %s\n', i, cap(i).name, cap(i).description);
end
iTrain = 1;  iTest = 2;                      % s1_3: interferer only;  s1_4: SOI + interferer

%% A1. Covariance of the data we want to clean
% Capture s1_4 has the SOI and the interferer on together. Estimate its
% spatial covariance from the K snapshots:
%   R = (1/K) sum_k x(k) x(k)^H         (N x N)
c = cap(iTest);
X = c.X;  K = size(X, 2);
azSOI = c.txAz(1);  elSOI = c.txEl(1);       % logged SOI position
azINT = c.txAz(2);  elINT = c.txEl(2);       % logged interferer position
Rboth = zeros(N, N);
for k = 1:K
    Rboth = Rboth + 0;                      % TODO: add x(k) x(k)^H / K
end
checkAnswer('Rboth', Rboth);
lambdaBoth = sort(real(eig(Rboth)), 'descend');

%% A2. MVDR from that covariance
%   w = R^{-1} s(theta_0) / ( s^H(theta_0) R^{-1} s(theta_0) )
% Use R\s rather than inv(R)*s: it solves R u = s directly.
s_soi = steeringVector(azSOI, elSOI, Ny, Nz, dOverLambda);
u = zeros(N, 1);                       % TODO: R^{-1} s, using backslash
w_both = zeros(N, 1);                  % TODO: normalize so that w^H s = 1
checkAnswer('w_both', w_both);
w_conv = s_soi / N;                          % conventional beam (Lab 1), for comparison
y_conv = zeros(1, K);  y_both = zeros(1, K);
for k = 1:K
    y_conv(k) = w_conv' * X(:, k);
    y_both(k) = 0;                         % TODO: y(k) = w^H x(k)
end
fprintf('\n%-26s %9s %9s %9s\n', '', 'SOI (dB)', 'INT (dB)', 'SIR (dB)');
printRow('conventional', y_conv, c);
printRow('MVDR, R from SOI + INT', y_both, c);
% Question: the interferer went down, but what happened to the SOI?

%% A3. Train on the interference alone, then apply
% Capture s1_3 was taken with only the interferer on. Estimate its
% covariance, build MVDR weights from it, and apply them to s1_4.
Xtr = cap(iTrain).X;  Ktr = size(Xtr, 2);
Rint = zeros(N, N);
for k = 1:Ktr
    Rint = Rint + 0;                    % TODO: add x(k) x(k)^H / K
end
checkAnswer('Rint', Rint);
u = zeros(N, 1);                       % TODO: R^{-1} s
w_int = zeros(N, 1);                   % TODO: normalize so that w^H s = 1
checkAnswer('w_int', w_int);
y_int = zeros(1, K);
for k = 1:K
    y_int(k) = 0;                          % TODO: y(k) = w^H x(k)
end
printRow('MVDR, R from INT only', y_int, c);

% Eigenvalues of the two covariances
lambdaInt = sort(real(eig(Rint)), 'descend');
figure('Name', 'A3: eigenvalues');
stem(1:N, 10*log10(lambdaBoth / lambdaBoth(end)), 'filled'); hold on;
stem(1:N, 10*log10(lambdaInt / lambdaInt(end)), 'filled'); hold off; grid on;
xlabel('eigenvalue index i'); ylabel('10 log_{10}(\lambda_i / \lambda_{min})  (dB)');
legend('R from SOI + INT', 'R from INT only'); title('How many strong directions does each covariance see?');

% The two MVDR beams overlaid, and the tone powers of all three beams
figure('Name', 'A3: the value of training');
[f, Pc] = spectrumdB(y_conv, c.fs);  [~, Pb] = spectrumdB(y_both, c.fs);  [~, Pi] = spectrumdB(y_int, c.fs);
subplot(2, 1, 1);
plot(f/1e6, Pb, f/1e6, Pi); grid on; ylim([-70 40]);
xlabel('frequency (MHz)'); ylabel('dB re 1 ADC unit^2'); labelTones(c);
legend('MVDR, R from SOI + INT', 'MVDR, R from INT only', 'Location', 'southwest');
title('Same capture, both MVDR beams aimed at the SOI');
Pbar = zeros(3, 2);                          % rows: beams;  columns: SOI, INT
Y = {y_conv, y_both, y_int};
for b = 1:3
    Pbar(b, 1) = 10*log10(measureTone(Y{b}, c.fs, c.fTone(1)));
    Pbar(b, 2) = 10*log10(measureTone(Y{b}, c.fs, c.fTone(2)));
end
subplot(2, 1, 2);
bar(Pbar); grid on; ylabel('tone power (dB)'); legend('SOI', 'INT', 'Location', 'northeast');
set(gca, 'XTickLabel', {'conventional', 'MVDR: SOI+INT', 'MVDR: INT only'});
title('Tone power at each beam output');
% Questions:
%  - Which beam keeps the SOI, and which removes the interferer best?
%  - Why does the SOI vanish when it is in the covariance? (The logged SOI
%    direction is a few degrees off. MVDR only protects the exact direction
%    s(theta_0); a strong signal anywhere else is "interference" to it.)
%  - The training capture was recorded minutes before the test capture.
%    Why do its weights still work?

%% A4. Where are the sources? The Capon (MVDR) power map
% Steer MVDR at every direction and record the output power it would have:
%   P(az, el) = 1 / ( s^H R^{-1} s )
% R is the covariance with both sources on (Rboth, from A1). Compare with Lab 1's
% conventional map: same sources, much sharper peaks.
azGrid = -60:1:60;  elGrid = -30:1:40;
Pcap = zeros(numel(elGrid), numel(azGrid));
for ia = 1:numel(azGrid)
    for ie = 1:numel(elGrid)
        sv = steeringVector(azGrid(ia), elGrid(ie), Ny, Nz, dOverLambda);
        Pcap(ie, ia) = 0;          % TODO: 10 log10 of 1 / (s^H R^{-1} s)
    end
end
checkAnswer('Pcap', Pcap);
[pkAz, pkEl] = twoPeaks(Pcap, azGrid, elGrid, 20);
figure('Name', 'A4: Capon map');
imagesc(azGrid, elGrid, Pcap); axis xy; colorbar; hold on;
xlabel('azimuth (deg)'); ylabel('elevation (deg)');
names = cell(1, 2);
for p = 1:2
    w = steeringVector(pkAz(p), pkEl(p), Ny, Nz, dOverLambda) / N;
    y = zeros(1, K);
    for k = 1:K
        y(k) = w' * X(:, k);
    end
    if measureTone(y, c.fs, c.fTone(1)) > measureTone(y, c.fs, c.fTone(2)), names{p} = 'SOI'; else, names{p} = 'INT'; end
    plot(pkAz(p), pkEl(p), 'wo', 'MarkerSize', 12, 'LineWidth', 2);
    text(pkAz(p) - 3, pkEl(p) + 3, sprintf('%s found (%d, %d)', names{p}, pkAz(p), pkEl(p)), 'Color', 'w', 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
end
hold off; title('Capon power map (dB):  o = found');
iINTfound = find(strcmp(names, 'INT'), 1);
if isempty(iINTfound), iINTfound = 2; end      % (only happens before the TODOs are filled in)
azINTfound = pkAz(iINTfound);  elINTfound = pkEl(iINTfound);

%% A5. Where did MVDR put its null?
% Beam pattern of the trained MVDR weights from A3:  |AF(az, el)| = |w^H s(az, el)|
G = zeros(numel(elGrid), numel(azGrid));
for ia = 1:numel(azGrid)
    for ie = 1:numel(elGrid)
        G(ie, ia) = 0;     % TODO: 20 log10 |w^H s(az, el)|
    end
end
figure('Name', 'A5: MVDR pattern');
imagesc(azGrid, elGrid, max(G, -50)); axis xy; colorbar; hold on;
plot(azSOI, elSOI, 'wx', azINT, elINT, 'wx', azINTfound, elINTfound, 'wo', 'MarkerSize', 12, 'LineWidth', 2);
text(azSOI + 3, elSOI - 3, 'SOI logged', 'Color', 'w');  text(azINT + 3, elINT - 3, 'INT logged', 'Color', 'w');
text(azINTfound - 3, elINTfound + 3, 'INT found', 'Color', 'w', 'HorizontalAlignment', 'right');
hold off; xlabel('azimuth (deg)'); ylabel('elevation (deg)');
title('MVDR pattern |w^H s| (dB), clipped at -50 dB');
% Questions: where is the null? Is it on the logged position, the found
% position, or somewhere else? MVDR nulls the interferer's real wavefront
% as the array measured it, not an angle. When the array's real response
% differs a little from the ideal model s(az, el), the null drawn with the
% model can sit a few degrees away from either marker.

%% ===================== PART B: harder =====================
%% B1. Keep the SOI when you cannot switch it off: diagonal loading
% In practice you often cannot turn the SOI off to train (A2's situation).
% Fix: diagonal loading.  R  ->  R + beta I.  Sweep beta from tiny to large
% (in units of the average element power) and plot SOI power and SIR.
sigma2 = real(trace(Rboth)) / N;             % average power per element
betas = [0 0.001 0.01 0.1 1 10];
Psoi = zeros(size(betas));  SIRb = zeros(size(betas));
for b = 1:numel(betas)
    % ---- YOUR CODE HERE ----------------------------------------
    % Load the covariance, build MVDR weights, apply them, and store Psoi(b)
    % and SIRb(b) in dB.
    Psoi(b) = NaN;  SIRb(b) = NaN;   % placeholder: replace with your code
    % ------------------------------------------------------------
end
figure('Name', 'B1: diagonal loading');
semilogx(max(betas, 1e-4), Psoi, 'o-', max(betas, 1e-4), SIRb, 's-'); grid on; hold on;
Pconv = 10*log10(measureTone(y_conv, c.fs, c.fTone(1)));
semilogx([1e-4 10], [Pconv Pconv], 'k--'); hold off;
xlabel('\beta  (x average element power;  0 plotted at 10^{-4})'); ylabel('dB');
legend('SOI power at the output', 'SIR', 'SOI power, conventional beam', 'Location', 'southeast');
title('Trained with the SOI present: loading trades nulling for robustness');

%% B2. How strong is too strong? Setup 2: SOI power from -5 to -45 dBm
% Captures 3-7: the interferer stays at -5 dBm, the SOI gets 10 dB weaker
% each time. Train on each capture itself (SOI present), steer at the
% logged SOI, and compare the SOI power out of MVDR with the conventional
% beam. Repeat with loading beta = 0.1.
iSweep = 3:7;
txP = zeros(size(iSweep));  lossNoLoad = zeros(size(iSweep));  lossLoad = zeros(size(iSweep));
for j = 1:numel(iSweep)
    cj = cap(iSweep(j));
    txP(j) = cj.txPowdBm(1);
    % ---- YOUR CODE HERE ----------------------------------------
    % For capture cj: estimate R, build conventional, MVDR and loaded-MVDR
    % (beta = 0.1) weights, and store the SOI power of each MVDR relative to
    % the conventional beam (dB).
    lossNoLoad(j) = NaN;  lossLoad(j) = NaN;   % placeholder: replace with your code
    % ------------------------------------------------------------
end
figure('Name', 'B2: SOI power sweep');
plot(txP, lossNoLoad, 'o-', txP, lossLoad, 's-'); grid on;
xlabel('SOI transmit power (dBm)'); ylabel('SOI at MVDR output minus conventional (dB)');
legend('MVDR, no loading', 'MVDR, \beta = 0.1', 'Location', 'best');
title('Signal cancellation grows with SOI power');
% Questions: when is training with the SOI present harmless, and when does
% it hurt? (At -45 dBm the SOI is barely above the noise, so its numbers are noisy.)

%% Local functions
function printRow(name, y, c)
% One row of the SOI / INT / SIR table.
pS = 10*log10(measureTone(y, c.fs, c.fTone(1)));
pI = 10*log10(measureTone(y, c.fs, c.fTone(2)));
fprintf('%-26s %9.1f %9.1f %9.1f\n', name, pS, pI, pS - pI);
end

function labelTones(c)
% Mark the SOI and interferer tone frequencies on the current spectrum plot.
nm = {'SOI', 'INT', 'Tx3', 'Tx4'};
yl = ylim;
for t = find(c.txActive)
    text(c.fTone(t)/1e6, yl(2) - 0.08*diff(yl), nm{t}, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
end
end

function [pkAz, pkEl] = twoPeaks(P, azGrid, elGrid, minSep)
% Brightest cell, then the brightest cell at least minSep degrees away.
[~, iMax] = max(P(:));
[ie1, ia1] = ind2sub(size(P), iMax);
Pm = P;
for ia = 1:numel(azGrid)
    for ie = 1:numel(elGrid)
        if abs(azGrid(ia) - azGrid(ia1)) < minSep && abs(elGrid(ie) - elGrid(ie1)) < minSep
            Pm(ie, ia) = -inf;
        end
    end
end
[~, iMax] = max(Pm(:));
[ie2, ia2] = ind2sub(size(P), iMax);
pkAz = [azGrid(ia1), azGrid(ia2)];  pkEl = [elGrid(ie1), elGrid(ie2)];
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
