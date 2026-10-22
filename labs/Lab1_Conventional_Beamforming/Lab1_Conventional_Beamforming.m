%% LAB 1: Conventional beamforming on real array data
% ADI ADXBAND16EBZ: 16 receive elements in a 4 x 4 grid, 10.4 GHz.
% Tx1 is the signal of interest (SOI); Tx2 is an interferer.
% Both are tones at slightly different frequencies, so we can tell them
% apart in a spectrum, but the array sees both at the same time.
%
% Part A (fill in): look at the data, steer beams, draw patterns, scan the room.
% Part B (harder):  tapering, an interferer near the main beam, and
%                   losses off broadside.
%
% Lines marked TODO are yours. Each TODO line already runs (with a
% placeholder), so you can run the script section by section.

clear; close all;
[cap, geom] = loadLabData(1);
Ny = geom.Ny;  Nz = geom.Nz;  N = geom.N;  dOverLambda = geom.dOverLambda;
for i = 1:numel(cap)
    fprintf('%d: %-5s  %s\n', i, cap(i).name, cap(i).description);
end

%% A1. Signal discovery: what does one element receive?
% Three captures: the SOI alone, the interferer alone, and both together.
% We look at element 1 only: a time trace (animated) and its spectrum.
iSOI = 1;  iINT = 3;  iBOTH = 4;           % capture indices (see list above)
show = [iSOI, iINT, iBOTH];
figure('Name', 'A1: one element');
for p = 1:3
    c = cap(show(p));
    x1 = c.X(1, :);                          % element 1, all snapshots
    subplot(3, 2, 2*p - 1);
    hT = plot(1:200, real(x1(1:200)), 1:200, imag(x1(1:200)));
    ylim(1.2*max(abs(x1))*[-1 1]); xlabel('snapshot k'); ylabel('x_1(k)');
    title(sprintf('%s: %s', c.name, c.description), 'Interpreter', 'none'); legend('I', 'Q');
    subplot(3, 2, 2*p);
    [f, PdB] = spectrumdB(x1, c.fs);
    plot(f/1e6, PdB); xlabel('frequency (MHz)'); ylabel('dB re 1 ADC unit^2'); grid on;
    title('spectrum of element 1');
    labelTones(c);
    for start = 1:10:600                     % "play" the capture
        set(hT(1), 'YData', real(x1(start:start+199)));
        set(hT(2), 'YData', imag(x1(start:start+199)));
        drawnow; pause(0.02);
    end
end
% Received power of each tone at element 1, in dB:
for p = 1:3
    c = cap(show(p));
    for t = find(c.txActive)
        fprintf('%s  Tx%d: %.1f dB at element 1\n', c.name, t, 10*log10(measureTone(c.X(1,:), c.fs, c.fTone(t))));
    end
end
% Observe: how strong is each tone at one element? Does it change when the
% other transmitter is switched on?

%% A2. Steering vectors and conventional weights
% Steering vector for azimuth az and elevation el (slides, Module 1):
%   u = cos(el) sin(az),   v = sin(el)
%   s_n = exp(-j 2 pi (d/lambda) [ (ny-1) u + (nz-1) v ]),   n = (ny-1) Nz + nz
% Conventional (delay-and-sum) weights, two scalings of the same beam:
%   w = s / sqrt(N)   lossless N-way combiner (like a Wilkinson): the SOI
%                     power goes UP by N (12 dB), the noise power stays the same
%   w = s / N         average of the elements: the SOI power stays the same,
%                     the noise power goes DOWN by N (12 dB)
% Either way the SNR improves by N = 16 (12 dB). We use w = s/N from A3 on,
% because it puts the peak of the beam pattern at 0 dB.
% Beamformer output:   y(k) = w^H x(k)
c = cap(iBOTH);                              % SOI and interferer both on
X = c.X;  K = size(X, 2);
azSOI = c.txAz(1);  elSOI = c.txEl(1);       % logged positions
azINT = c.txAz(2);  elINT = c.txEl(2);

s_soi = steeringVector(azSOI, elSOI, Ny, Nz, dOverLambda);
checkAnswer('s_soi', s_soi);
w_wilk = zeros(N, 1);                     % TODO: SOI beam, lossless combiner scaling
w_soi = zeros(N, 1);                      % TODO: SOI beam, averaging scaling
w_int = zeros(N, 1);                      % TODO: same, for the interferer
w_bs  = zeros(N, 1);                      % TODO: same, for broadside (az = 0, el = 0)
checkAnswer('w_soi', w_soi);

% Beamformer output for each beam, one snapshot at a time
y_wilk = zeros(1, K);  y_soi = zeros(1, K);  y_int = zeros(1, K);  y_bs = zeros(1, K);
for k = 1:K
    y_wilk(k) = 0;                        % TODO: y(k) = w^H x(k)
    y_soi(k) = 0;                         % TODO: y(k) = w^H x(k)
    y_int(k) = 0;                         % TODO
    y_bs(k)  = 0;                         % TODO
end
checkAnswer('y_soi_first10', y_soi(1:10));

% Power of each tone in each beam, against element 1
% (all powers in dB relative to 1 ADC unit squared: the raw digital scale,
%  not dBm; the receiver chain is not calibrated in absolute power)
fprintf('\n%-12s %10s %10s\n', '', 'SOI (dB)', 'INT (dB)');
rows = {'element 1', X(1, :); 'SOI s/sqrtN', y_wilk; 'SOI s/N', y_soi; 'INT s/N', y_int; 'bside s/N', y_bs};
for r = 1:size(rows, 1)
    fprintf('%-12s %10.1f %10.1f\n', rows{r, 1}, 10*log10(measureTone(rows{r, 2}, c.fs, c.fTone(1))), ...
            10*log10(measureTone(rows{r, 2}, c.fs, c.fTone(2))));
end
% Spectrum: element 1 vs the SOI beam, with both scalings
figure('Name', 'A2: element vs beam');
[f, P1] = spectrumdB(X(1, :), c.fs);
[~, Pw] = spectrumdB(y_wilk, c.fs);
[~, Pa] = spectrumdB(y_soi, c.fs);
subplot(2, 1, 1);
plot(f/1e6, P1, f/1e6, Pw); grid on; ylim([-50 50]); ylabel('dB re 1 ADC unit^2');
labelTones(c); legend('element 1', 'SOI beam, w = s/sqrt(N)', 'Location', 'northwest');
title('Lossless combiner: the SOI rises about 12 dB, the noise floor stays put');
subplot(2, 1, 2);
plot(f/1e6, P1, f/1e6, Pa); grid on; ylim([-50 50]); ylabel('dB re 1 ADC unit^2');
labelTones(c); legend('element 1', 'SOI beam, w = s/N', 'Location', 'northwest'); xlabel('frequency (MHz)');
title('Average: the SOI stays put, the noise floor drops about 12 dB');
% Questions:
%  - In both plots, how far apart are the SOI peak and the noise floor,
%    compared with element 1? Compare with 10 log10(N).
%  - The SOI rises a little less than 12 dB in the top plot. Element 1 is
%    one of the stronger elements, and the logged direction is a few degrees off.
%  - How much did the interferer drop in the SOI beam? And the SOI in the interferer beam?

%% A3. Beam patterns: AF(az) = w^H s(az, el)
azScan = -90:0.5:90;
AF = zeros(3, numel(azScan));
W = [w_soi, w_int, w_bs];
elCut = [elSOI, elINT, 0];                   % cut each beam through its own peak
for b = 1:3
    for i = 1:numel(azScan)
        sScan = steeringVector(azScan(i), elCut(b), Ny, Nz, dOverLambda);
        AF(b, i) = 0;                         % TODO: array factor
    end
end
checkAnswer('AF_soi', AF(1, :));
figure('Name', 'A3: patterns');
plot(azScan, 20*log10(abs(AF) + 1e-6)); grid on; ylim([-40 2]);
xlabel('azimuth (deg)'); ylabel('|AF| (dB)'); legend('SOI beam', 'INT beam', 'broadside');
title('Azimuth cut of each beam, through its own elevation');
% Questions: which beam is widest, and why? Where does the INT beam's
% pattern sit at the SOI's azimuth (and vice versa)?

%% A4. Scan the room: 1-degree sweep in azimuth and elevation (live)
azGrid = -60:1:60;  elGrid = -30:1:40;
Pmap = nan(numel(elGrid), numel(azGrid));
figure('Name', 'A4: power map');
hImg = imagesc(azGrid, elGrid, Pmap); axis xy; colorbar;
xlabel('azimuth (deg)'); ylabel('elevation (deg)'); title('Beam output power (dB re 1 ADC unit^2), filling in...');
for ia = 1:numel(azGrid)
    for ie = 1:numel(elGrid)
        w = steeringVector(azGrid(ia), elGrid(ie), Ny, Nz, dOverLambda) / N;
        P = 0;
        for k = 1:K
            P = P + 0;                             % TODO: average output power, (1/K) sum |w^H x(k)|^2
        end
        Pmap(ie, ia) = 10*log10(P);
    end
    if mod(ia, 5) == 0 || ia == numel(azGrid)
        set(hImg, 'CData', Pmap); drawnow;
    end
end
checkAnswer('Pmap', Pmap);

% Find the two strongest spots: the brightest cell, then the brightest cell
% at least 20 deg away from it. Steer at each and ask which tone dominates.
[~, iMax] = max(Pmap(:));
[ie1, ia1] = ind2sub(size(Pmap), iMax);
Pmasked = Pmap;
for ia = 1:numel(azGrid)
    for ie = 1:numel(elGrid)
        if abs(azGrid(ia) - azGrid(ia1)) < 20 && abs(elGrid(ie) - elGrid(ie1)) < 20
            Pmasked(ie, ia) = -inf;
        end
    end
end
[~, iMax] = max(Pmasked(:));
[ie2, ia2] = ind2sub(size(Pmap), iMax);
peaksAz = [azGrid(ia1), azGrid(ia2)];  peaksEl = [elGrid(ie1), elGrid(ie2)];
hold on;
for p = 1:2
    w = steeringVector(peaksAz(p), peaksEl(p), Ny, Nz, dOverLambda) / N;
    y = zeros(1, K);
    for k = 1:K
        y(k) = w' * X(:, k);
    end
    if measureTone(y, c.fs, c.fTone(1)) > measureTone(y, c.fs, c.fTone(2))
        name = 'SOI';
    else
        name = 'INT';
    end
    plot(peaksAz(p), peaksEl(p), 'wo', 'MarkerSize', 12, 'LineWidth', 2);
    if peaksAz(p) > 20, side = 'right'; dx = -3; else, side = 'left'; dx = 3; end
    text(peaksAz(p) + dx, peaksEl(p) + 3, sprintf('%s measured (%d, %d)', name, peaksAz(p), peaksEl(p)), ...
         'Color', 'w', 'FontWeight', 'bold', 'HorizontalAlignment', side);
end
hold off; title('Beam output power (dB):  o = measured peak');
% Questions: how many sources do you see? Do the measured peaks sit on the
% logged positions? The tripods were placed by hand and the logged angles
% were estimates. In practice the geometry is never exactly what we think,
% and the methods we use have to tolerate that.

%% ===================== PART B: harder =====================
%% B1. Tapering: push the interferer down further
% A taper weights the elements unequally:  w_n = a_n s_n / sum(a).
% On the 4 x 4 grid we apply the same 4-point window t along the rows and
% along the columns:  a(ny, nz) = t(ny) * t(nz).
% Three windows:
%   rectangular  t = [1 1 1 1]                    (no taper)
%   Chebyshev    t = [0.43 1 1 0.43]              (30 dB sidelobes, chebwin(4,30))
%   Hamming      t(m) = 0.54 - 0.46 cos(2 pi (m-1)/3)  =  [0.08 0.77 0.77 0.08]
% For each: draw the beam (left) and the output spectrum (right), and read
% off the SOI-to-interferer ratio (SIR).
winNames = {'rectangular', 'Chebyshev 30 dB', 'Hamming'};
winT = zeros(3, 4);
winT(1, :) = [1 1 1 1];
winT(2, :) = [0.4290 1 1 0.4290];
for m = 1:4
    winT(3, m) = 0.54 - 0.46*cos(2*pi*(m-1)/3);
end
figure('Name', 'B1: tapers');
for e = 1:3
    t = winT(e, :);
    % ---- YOUR CODE HERE ----------------------------------------
    % Build the 2-D taper a (N x 1) from t, the tapered weights w, the output
    % y, and the SIR in dB.
    w = zeros(N, 1);  y = zeros(1, K);  SIR = NaN;   % placeholder: replace with your code
    % ------------------------------------------------------------
    AFt = zeros(1, numel(azScan));
    for i = 1:numel(azScan)
        AFt(i) = w' * steeringVector(azScan(i), elSOI, Ny, Nz, dOverLambda);
    end
    subplot(3, 2, 2*e - 1);
    plot(azScan, 20*log10(abs(AFt) + 1e-6)); grid on; ylim([-40 2]); hold on;
    plot([azSOI azSOI], [-40 2], 'g--', [azINT azINT], [-40 2], 'r--', [peaksAz(2) peaksAz(2)], [-40 2], 'r:'); hold off;
    title(sprintf('%s beam (el = %g deg)', winNames{e}, elSOI)); ylabel('|AF| (dB)');
    subplot(3, 2, 2*e);
    [f, Py] = spectrumdB(y, c.fs);
    plot(f/1e6, Py); grid on; ylim([-60 40]); labelTones(c);
    title(sprintf('%s:  SIR = %.1f dB', winNames{e}, SIR)); ylabel('dB');
end
subplot(3, 2, 5); xlabel('azimuth (deg)');  subplot(3, 2, 6); xlabel('frequency (MHz)');
% Green dashed: SOI. Red dashed: interferer where we logged it; red dotted:
% where the power map found it. Curve ball: is the strongest taper the best?

%% B2. Interferer near the main beam
% Repeat B1 on capture 5 (interferer logged at az = +10 deg, inside the
% SOI's main beam). Before you run it: will a taper help?
c5 = cap(5);
figure('Name', 'B2: interferer near the main beam');
for e = 1:3
    t = winT(e, :);
    % ---- YOUR CODE HERE ----------------------------------------
    % As in B1, for capture c5: build the taper a, the weights w (steered at
    % c5's SOI), the output y and the SIR.
    w = zeros(N, 1);  y = zeros(1, size(c5.X, 2));  SIR = NaN;   % placeholder: replace with your code
    % ------------------------------------------------------------
    AFt = zeros(1, numel(azScan));
    for i = 1:numel(azScan)
        AFt(i) = w' * steeringVector(azScan(i), c5.txEl(1), Ny, Nz, dOverLambda);
    end
    subplot(3, 2, 2*e - 1);
    plot(azScan, 20*log10(abs(AFt) + 1e-6)); grid on; ylim([-40 2]); hold on;
    plot([c5.txAz(1) c5.txAz(1)], [-40 2], 'g--', [c5.txAz(2) c5.txAz(2)], [-40 2], 'r--'); hold off;
    title(sprintf('%s beam (el = %g deg)', winNames{e}, c5.txEl(1))); ylabel('|AF| (dB)');
    subplot(3, 2, 2*e);
    [f, Py] = spectrumdB(y, c5.fs);
    plot(f/1e6, Py); grid on; ylim([-60 40]); labelTones(c5);
    title(sprintf('%s:  SIR = %.1f dB', winNames{e}, SIR)); ylabel('dB');
end
subplot(3, 2, 5); xlabel('azimuth (deg)');  subplot(3, 2, 6); xlabel('frequency (MHz)');
% Green dashed: SOI. Red dashed: interferer (logged at +10 deg).
% Fixed weights cannot separate two sources inside one main beam.
% That is the job of the adaptive methods in the next modules.

%% B3. Off broadside, the beam receives less
% The interferer (Tx2) was moved through az = 50, 30, 10, 0, -10 deg.
% Steer at Tx2 in each capture and plot its received power against azimuth.
% ---- YOUR CODE HERE ----------------------------------------
% Collect azTx2 (logged azimuth) and PTx2 (Tx2 power in dB in a beam
% steered at it) for every capture with Tx2 on.
azTx2 = nan;  PTx2 = nan;   % placeholder: replace with your code
% ------------------------------------------------------------
figure('Name', 'B3: off broadside');
plot(azTx2, PTx2, 'o-'); grid on; xlabel('Tx2 logged azimuth (deg)'); ylabel('Tx2 power in its beam (dB)');
title('Received power vs scan angle');
% The array factor's peak is the same at every scan angle. What else changes
% off broadside? (Think about the element pattern and the projected aperture.)

%% Local functions
function labelTones(c)
% Mark the SOI and interferer tone frequencies on the current spectrum plot.
names = {'SOI', 'INT', 'Tx3', 'Tx4'};
yl = ylim;
for t = find(c.txActive)
    text(c.fTone(t)/1e6, yl(2) - 0.08*diff(yl), names{t}, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
end
end

% Steering vector (yours to complete)
function s = steeringVector(az, el, Ny, Nz, dOverLambda)
% s = steeringVector(az, el, Ny, Nz, dOverLambda), angles in degrees
u = cosd(el) * sind(az);
v = sind(el);
s = zeros(Ny*Nz, 1);
for ny = 1:Ny
    for nz = 1:Nz
        n = (ny-1)*Nz + nz;
        s(n) = 1;   % TODO: s_n = exp(-j 2 pi (d/lambda) [ (ny-1) u + (nz-1) v ])
    end
end
end
