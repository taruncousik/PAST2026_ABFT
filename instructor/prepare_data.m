%% prepare_data.m: build labN_data.mat from the raw ADXBAND16EBZ captures
% Instructors only. Participants get the prepared labN_data.mat files.
%
% Raw captures: Runs_1015/<setup folder>/<name>_<date>.mat, each with a
% struct "capture" (400 MS/s, 16 channels). For every capture we
%   1. apply the phase of ADI's calibration weights (phase only, so the
%      noise stays equal on all channels),
%   2. shift the band centre fc to 0 Hz,
%   3. low-pass filter and keep every DEC-th sample,
%   4. estimate each tone's exact frequency (the transmit and receive
%      boards are not locked, so tones sit ~120-170 kHz off nominal),
%   5. (Labs 2, 3) pass the SOI's transmitted IQ through steps 2-3 as well,
%      to serve as the LMS / RLS reference.
%
% NOTE: the lab files shipped so far were made with a Python version of
% these exact steps; this script has not been run in MATLAB yet.

rawDir = 'Runs_1015';          % folder with the unzipped raw captures
DEC = 16;                      % 400 MS/s -> 25 MS/s
b = fir1(8*DEC, 0.8/DEC);      % low-pass, passband about +-10 MHz
Ny = 4;  Nz = 4;  dOverLambda = 13.6e-3 / (3e8/10.425e9);

% {name, description, band centre fc}
labs = {
  1, {'s1_1', 'SOI only', -7.5e6; 's1_2', 'SOI + interferer at az +50', -7.5e6; ...
      's1_3', 'interferer only, az +30', -7.5e6; 's1_4', 'SOI + interferer at az +30', -7.5e6; ...
      's1_5', 'SOI + interferer at az +10', -7.5e6; 's1_6', 'SOI + interferer at az 0', -7.5e6; ...
      's1_7', 'SOI + interferer at az -10 (same place as the SOI)', -7.5e6}, false;
  2, {'s1_4', 'SOI + interferer at logged az +30', -7.5e6; 's3_3', 'SOI + three interferers (Setup 3)', -12.5e6}, true;
  3, {'s3_3', 'SOI + three interferers (Setup 3)', -12.5e6; 's1_2', 'SOI + interferer at logged az +50', -7.5e6; ...
      's1_4', 'SOI + interferer at logged az +30', -7.5e6; 's1_5', 'SOI + interferer at logged az +10', -7.5e6}, true;
  4, {'s1_3', 'interferer only, logged az +30 (training data)', -7.5e6; 's1_4', 'SOI + interferer at logged az +30', -7.5e6; ...
      's2_1', 'SOI -5 dBm + interferer -5 dBm', -7.5e6; 's2_2', 'SOI -15 dBm + interferer -5 dBm', -7.5e6; ...
      's2_3', 'SOI -25 dBm + interferer -5 dBm', -7.5e6; 's2_4', 'SOI -35 dBm + interferer -5 dBm', -7.5e6; ...
      's2_5', 'SOI -45 dBm + interferer -5 dBm', -7.5e6; 's2_6', 'interferer only, -5 dBm', -7.5e6}, false};

for L = 1:size(labs, 1)
    list = labs{L, 2};  withRef = labs{L, 3};
    clear cap
    for i = 1:size(list, 1)
        cap(i) = prepCapture(rawDir, list{i, 1}, list{i, 2}, list{i, 3}, DEC, b, withRef); %#ok<SAGROW>
    end
    if labs{L, 1} == 4
        % s2_5: the SOI is at -45 dBm, too weak to find its frequency reliably;
        % borrow it from s2_4 (same transmitter, minutes earlier)
        i4 = find(strcmp({cap.name}, 's2_4'));  i5 = find(strcmp({cap.name}, 's2_5'));
        cap(i5).fTone(1) = cap(i4).fTone(1);
    end
    save(sprintf('lab%d_data.mat', labs{L, 1}), 'cap', 'Ny', 'Nz', 'dOverLambda');
end

function c = prepCapture(rawDir, name, descr, fc, DEC, b, withRef)
f = dir(fullfile(rawDir, '*', [name '_2*.mat']));
R = load(fullfile(f(1).folder, f(1).name));
cp = R.capture;
Xraw = cp.rxRaw;  if size(Xraw, 2) ~= 16, Xraw = Xraw.'; end      % K x 16
fs = cp.fsRxIQ;  K = size(Xraw, 1);  n = (0:K-1).';
Xcal = Xraw .* exp(1j*angle(cp.calweights(:).'));                   % 1. phase-only calibration
Y = filter(b, 1, Xcal .* exp(-1j*2*pi*fc/fs*n));                    % 2-3. shift, low-pass ...
Y = Y(8*DEC+1:DEC:end, :);                                          %      ... and decimate
bb = cp.txIQ.bbHz(:).';
fTone = nan(1, 4);
for t = 1:4                                                         % 4. exact tone frequencies
    if cp.txActive(t)
        fTone(t) = findTone(Xcal, fs, bb(t) + 150e3) - fc;
    end
end
c.name = name;  c.description = descr;
c.X = Y.';  c.fs = fs/DEC;
c.txActive = double(cp.txActive(:).');  c.txAz = cp.txAz(:).';
c.txEl = cp.txEl(:).';  c.txPowdBm = cp.txPowdBm(:).';  c.fTone = fTone;
if withRef                                                          % 5. SOI reference
    iq = cp.txIQ.iq;  if size(iq, 1) < size(iq, 2), iq = iq.'; end  % Ktx x rows
    m = (0:size(iq, 1)-1).';
    % the stored rows are not always in transmitter order: pick the one at Tx1's frequency
    score = zeros(1, size(iq, 2));
    for j = 1:size(iq, 2)
        score(j) = abs(mean(iq(:, j) .* exp(-1j*2*pi*bb(1)/fs*m)));
    end
    [~, r] = max(score);
    dTx = filter(b, 1, iq(:, r) .* exp(-1j*2*pi*fc/fs*m));
    dTx = dTx(8*DEC+1:DEC:end);
    c.dTx = dTx(1:size(Y, 1)).';
    c.fTxNominal = bb - fc;
end
end

function fBest = findTone(X, fs, f0)
% Search +-400 kHz around f0 for the frequency with the most power, summed over channels.
K = size(X, 1);  k = (0:K-1).';
fGrid = f0 + (-400e3:fs/K/8:400e3);
P = zeros(size(fGrid));
for g = 1:numel(fGrid)
    a = sum(X .* exp(-1j*2*pi*fGrid(g)/fs*k), 1) / K;
    P(g) = sum(abs(a).^2);
end
[~, iBest] = max(P);
fBest = fGrid(iBest);
end
