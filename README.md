# Warper

`PV_PolyFreqWarp` is a SuperCollider UGen extension that performs arbitrary polynomial frequency warping on spectral data in the phase-vocoder domain. Peak tracking and region-based phase synchronization are used to maintain phase coherence and avoid phase-cancellation artifacts.

## Reference Paper

The phase-locking algorithm and peak-based spectral region processing implemented in this UGen are based on the seminal phase-vocoder techniques established by **Jean Laroche** and **Mark Dolson**:

> **Jean Laroche and Mark Dolson.**  
> *"Improved Phase Vocoder Time-Scale Modification of Audio."*  
> *IEEE Transactions on Speech and Audio Processing*, Vol. 7, Issue 3, pp. 323–332, May 1999.  
> DOI: [10.1109/89.759041](https://doi.org/10.1109/89.759041)

## Overview & Features

- **Polynomial Frequency Warping**: Apply 5th-order polynomial warping using coefficients `curve0` through `curve3`, alongside global `scale` (transposition) and `shift`.
- **Phase-Locked Region Processing**: Detects spectral peaks, segments the spectrum into regions bounded by local minima between peaks, and updates frame-to-frame phases in sync for each peak region.
- **Multichannel Phase Alignment (`detector`)**: An optional secondary FFT detector chain can be passed to drive peak detection. For stereo or multichannel signals, feeding a shared detector chain (e.g. summed mix) to all channel instances ensures identical region boundaries and preserves stereo imaging and inter-channel phase coherence.
- **Boundary Handling (`reflect`)**: Supports spectral reflection with phase conjugation at DC and Nyquist boundaries.

## Building and Installing

### Prerequisites
- [SuperCollider](https://supercollider.github.io/) source tree or SDK headers
- [CMake](https://cmake.org/) (3.16 or later)
- C++17 compiler (Clang / GCC / MSVC)

### Build Instructions

```sh
mkdir build && cd build
cmake .. -DSUPERCOLLIDER_SOURCE_DIR=/path/to/supercollider
cmake --build . --target install
```

On macOS, this installs the compiled plugin (`PolyFreqWarp.scx`), the SuperCollider class (`PolyFreqWarp.sc`), and the documentation (`PV_PolyFreqWarp.schelp`) into:
`~/Library/Application Support/SuperCollider/Extensions/Warper/`

## Quick Example (SuperCollider)

```supercollider
(
s.waitForBoot {
    b = Buffer.read(s, Platform.resourceDir +/+ "sounds/a11wlk01.wav");
    s.sync;

    x = {
        var in, chain, sig;
        in = PlayBuf.ar(1, b, BufRateScale.kr(b), loop: 1);
        chain = FFT(LocalBuf(4096), in, hop: 0.25, wintype: 1);
        chain = PV_PolyFreqWarp(
            buffer: chain,
            scale: MouseX.kr(0.5, 2.0, \exponential), // Pitch transposition
            shift: MouseY.kr(-0.2, 0.2)               // Frequency shift
        );
        sig = IFFT(chain);
        sig ! 2 * 0.2;
    }.play;
};
)
```