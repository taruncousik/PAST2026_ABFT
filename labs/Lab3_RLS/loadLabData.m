function [cap, geom] = loadLabData(labNumber)
% LOADLABDATA  Load the recorded captures for one lab.
%
%   [cap, geom] = loadLabData(1)
%
%   cap(i) is one capture (one row of the experiment spreadsheet):
%     cap(i).name         e.g. 's1_4'
%     cap(i).description  what was switched on, in words
%     cap(i).X            N x K complex snapshots: X(n,k) = element n, snapshot k
%     cap(i).fs           sample rate of X (Hz)
%     cap(i).txActive     1x4, which transmitters were on (Tx1 = SOI, Tx2 = interferer, ...)
%     cap(i).txAz         1x4, logged azimuth of each transmitter (deg)
%     cap(i).txEl         1x4, logged elevation of each transmitter (deg)
%     cap(i).txPowdBm     1x4, transmit power setting (dBm)
%     cap(i).fTone        1x4, frequency of each transmitter's tone in X (Hz)
%
%   geom describes the array:
%     geom.Ny, geom.Nz    elements along y (columns) and z (rows): 4 x 4
%     geom.N              Ny*Nz = 16
%     geom.dOverLambda    element spacing in wavelengths (13.6 mm at 10.42 GHz)
%     Element n sits at column ny and row nz with n = (ny-1)*Nz + nz,
%     the same ordering as the slides.
%
% -------------------------------------------------------------------------
% How the data was prepared (you do not need this for the lab):
%   The ADXBAND16EBZ samples at 400 MS/s. Each tone sits in one narrow
%   slice of that 400 MHz band, so we kept only that slice:
%     1. per-channel phase calibration from ADI's calibration file
%     2. shift the band centre (-7.5 MHz) down to 0 Hz
%     3. low-pass filter (fir1, 129 taps) and keep every 16th sample
%   Result: 25 MS/s, about 1000 snapshots per capture, and roughly 12 dB
%   more signal-to-noise than the full-band data. See prep/prepare_data.m.
% -------------------------------------------------------------------------

S = load(sprintf('lab%d_data.mat', labNumber));
cap = S.cap;
geom.Ny = double(S.Ny);
geom.Nz = double(S.Nz);
geom.N = geom.Ny * geom.Nz;
geom.dOverLambda = double(S.dOverLambda);
end
