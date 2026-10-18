# Hands-On Adaptive Beamforming

**IEEE PAST 2026 Tutorial**
Mark Yeary (University of Oklahoma) · Jon Kraft (Analog Devices) · Tarun Cousik (Virginia Tech)

From delay-and-sum to LMS, RLS and MVDR, run on real 16-channel X-band array data recorded with the Analog Devices ADXBAND16EBZ.

> The materials are updated continuously until the workshop and remain available to participants afterwards.

## Before the workshop

- Bring a laptop and its charger.
- All code runs in **MATLAB Online** in your web browser. MathWorks has generously provided licenses for all participants.
- Nothing to install: open the repository in MATLAB Online, open a lab folder, and run the script.

[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=taruncousik/PAST2026_ABFT)

## Labs

| Lab | Topic | What you do | Data |
|---|---|---|---|
| [Lab 1](labs/Lab1_Conventional_Beamforming) | Conventional beamforming | Explore the received signals, build steering vectors, steer beams, scan the room, apply tapers | Setup 1 |
| [Lab 2](labs/Lab2_LMS) | LMS | Build a reference from the transmitted IQ, run LMS, choose the step size, null three interferers | Setups 1, 3 |
| [Lab 3](labs/Lab3_RLS) | RLS | Write the RLS recursion, compare with LMS, choose the forgetting factor, track a moving interferer | Setups 1, 3 |
| [Lab 4](labs/Lab4_MVDR) | MVDR | Compare training with and without the signal of interest, map the sources (Capon), use diagonal loading | Setups 1, 2 |

Each lab folder is self-contained:

| File | Purpose |
|---|---|
| `LabN_*.m` | Participant version. Lines marked `TODO` are yours; each already runs with a placeholder, so you can go section by section. |
| `LabN_*_solution.m` | Complete version. |
| `labN_data.mat` | The prepared captures for this lab. |
| `loadLabData.m` | Loads the data; its header documents every field. |
| `measureTone.m`, `spectrumdB.m` | Small measurement helpers. |
| `checkAnswer.m`, `labN_ref.mat` | Self-checks: prints `[OK]` when your result matches the reference. |

Part A of each lab is fill-in-the-blanks, following the MATLAB snippets on the slides. Part B is harder: you write more of the code and meet some of the surprises real arrays bring.

To view a script as a Live Script, right-click it and choose **Open as Live Script**.

## The data

- **Array:** ADI ADXBAND16EBZ, 16 receive elements in a 4 × 4 grid, 13.6 mm spacing, 10.42 GHz.
- **Transmitters:** up to four tones on tripods. Tx1 is the signal of interest (SOI); the others are interferers.
  - **Setup 1:** one interferer moved toward the SOI.
  - **Setup 2:** the SOI's power stepped down against a fixed interferer.
  - **Setup 3:** the SOI with three interferers.
- **Preparation:** each capture was phase-calibrated, shifted to baseband around the tones, and decimated from 400 MS/s to 25 MS/s. The SOI's transmitted IQ went through the same steps, to serve as the reference for LMS and RLS. Details are in `loadLabData.m` and `instructor/prepare_data.m`.

## Repository layout

```
labs/
  Lab1_Conventional_Beamforming/
  Lab2_LMS/
  Lab3_RLS/
  Lab4_MVDR/
instructor/
  prepare_data.m   rebuilds labN_data.mat from the raw captures (not in this repository)
```
