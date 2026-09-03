#include "SC_PlugIn.h"
#include "FFT_UGens.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstring>

InterfaceTable* ft;

struct PV_PolyFreqWarp : PV_Unit
{
    int mNumBins;
    std::complex<float>* mPhase;
};

// Define the reflection modes for handling negative frequencies and those above the Nyquist frequency. 

enum class ReflectMode { None, Negative, Reflect };

static void PolyFreqWarpFree(PV_PolyFreqWarp* unit)
{
    RTFree(unit->mWorld, unit->mPhase);
}

// Ensure that the phase arrays are allocated and initialized for the given number of bins

static bool PolyFreqWarpEnsureState(PV_PolyFreqWarp* unit, int numBins)
{
    // If the number of bins is already correct, no need to reallocate

    if (unit->mNumBins == numBins)
        return true;

    // Free any existing phase arrays before allocating new ones

    PolyFreqWarpFree(unit);
    unit->mPhase = static_cast<std::complex<float>*>(RTAlloc(unit->mWorld, numBins * sizeof(std::complex<float>)));

    // If allocation failed, free any allocated memory and return false

    if (!unit->mPhase)
    {
        PolyFreqWarpFree(unit);
        unit->mNumBins = 0;
        return false;
    }

    // Initialize the phase arrays to represent no rotation (real=1, imag=0)

    for (int i = 0; i < numBins; ++i)
        unit->mPhase[i] = std::complex<float>(1.f, 0.f);
    unit->mNumBins = numBins;
    
    return true;
}

// Find the next peak in the magnitudes array starting from the given index.
// A peak is defined as a value that is greater than its two immediate neighbors on both sides.

static int PolyFreqWarpPeak(const float* magnitudes, int numBins, int start)
{
    for (int i = start; i < numBins; ++i)
    {
        if (magnitudes[i] > magnitudes[i - 1] && magnitudes[i] > magnitudes[i - 2] &&
            magnitudes[i] > magnitudes[i + 1] && magnitudes[i] > magnitudes[i + 2])
            return i;
    }
    return -1;
}

// Interpolate a peak position using a parabolic fit to the logarithm of the magnitudes

static double PolyFreqWarpPeakPosition(const float* magnitudes, int peak)
{
    const double left = std::log(std::max(static_cast<double>(magnitudes[peak - 1]), 1.0e-20));
    const double centre = std::log(std::max(static_cast<double>(magnitudes[peak]), 1.0e-20));
    const double right = std::log(std::max(static_cast<double>(magnitudes[peak + 1]), 1.0e-20));
    const double divisor = left + right - 2.0 * centre;
   
    if (divisor == 0.0)
        return static_cast<double>(peak);
    return static_cast<double>(peak) + 0.5 * (left - right) / divisor;
}

// Paste a single bin into the output array at a fractional position,
// Distribute the value between the two nearest bins based on the fractional part of the position.

enum class PasteMode { Add, Clip, DC, Nyquist };

static inline void PolyFreqWarpAddNorm(SCComplex* output, int idx, float fraction, std::complex<float> value)
{
    output[idx].real += value.real() * (1.f - fraction);
    output[idx].imag += value.imag() * (1.f - fraction);
    output[idx + 1].real += value.real() * fraction;
    output[idx + 1].imag += value.imag() * fraction;
}

static inline void PolyFreqWarpAddClip(SCComplex* output, int idx, float fraction, std::complex<float> value)
{
    output[idx].real += value.real() * (1.f - fraction);
    output[idx].imag += value.imag() * (1.f - fraction);
}

static inline void PolyFreqWarpAddDC(SCComplex* output, int /*idx*/, float fraction, std::complex<float> value)
{
    output[0].real += std::abs(value * (1.f - fraction));
    output[1].real += value.real() * fraction;
    output[1].imag += value.imag() * fraction;
}

static inline void PolyFreqWarpAddNyquist(SCComplex* output, int idx, float fraction, std::complex<float> value)
{
    output[idx].real += value.real() * (1.f - fraction);
    output[idx].imag += value.imag() * (1.f - fraction);
    output[idx + 1].real += std::abs(value * fraction);
}

// Shift and paste a single bin

template <PasteMode mode>
static inline void PolyFreqWarpDoBin(const SCComplex* input, SCComplex* output, int bin, int idx, float fraction,
                                std::complex<float> rotate, std::complex<float>* phase, bool conjugate)
{
    // Rotate the input bin by the given complex value

    const std::complex<float> source = std::complex<float>(input[bin].real, input[bin].imag);
    std::complex<float> value = source * rotate;
    value = conjugate ? std::conj(value) : value;

    // Add to the output (use the correct mode)
    
    switch (mode)
    {
        case PasteMode::Add:        PolyFreqWarpAddNorm(output, idx, fraction, value);          break;
        case PasteMode::Clip:       PolyFreqWarpAddClip(output, idx, fraction, value);          break;
        case PasteMode::DC:         PolyFreqWarpAddDC(output, idx, fraction, value);            break;
        case PasteMode::Nyquist:    PolyFreqWarpAddNyquist(output, idx, fraction, value);       break;
    } 
               
    // Update the phase arrays for the next iteration

    phase[bin] = rotate;
}

// Calculate the warped frequency position for a given peak position 

static inline double PolyFreqWarpPolynomial(double peak, const float* p, double binShift, int numBins)
{
    const double pNorm = peak / static_cast<double>(numBins - 1);
    const double polynomial = p[4] + 10000.0 * pNorm * (p[3] + pNorm * (p[2] + pNorm * (p[1] + pNorm * p[0])));
    
    return peak * polynomial + binShift;
}

// Process the input and write to the output

static void PolyFreqWarpProcess(const SCComplex* input, SCComplex* output, std::complex<float>* phase, const float* magnitudes,
                                int numBins, const float* params, ReflectMode reflect, float overlap)
{
    const double phaseConst = 2.0 * M_PI / overlap;
    const double binShift = params[5] * (numBins - 1);

    // Zero the output

    for (int i = 0; i < numBins; ++i)
    {
        output[i].real = 0.f;
        output[i].imag = 0.f;
    }
    
    // Find the first peak

    int regionStart = 0;
    int peak = PolyFreqWarpPeak(magnitudes, numBins, 0);
    if (peak == -1)
        peak = 0;

    // Process one region at a time, defined by the peaks in the magnitudes array

    while (regionStart < numBins)
    {
        // Find the next peak after the current peak to aid in finding the region end

        int nextPeak = PolyFreqWarpPeak(magnitudes, numBins, peak + 1);
        int regionEnd;

        // Find the minimum magnitude between current peak and the next to define the end of the current region

        if (nextPeak != -1)
        {
            auto minIt = std::min_element(magnitudes + peak, magnitudes + nextPeak + 1);
            regionEnd = static_cast<int>(std::distance(magnitudes, minIt)) + 1;
        }
        else
            regionEnd = numBins;

        // Calculate the shift for this region based on the peak position and the polynomial

        const double peakPosition = PolyFreqWarpPeakPosition(magnitudes, peak);
        const double peakPolynomial = PolyFreqWarpPolynomial(peakPosition, params, binShift, numBins);
        const double shift = peakPolynomial - peakPosition;

        // Calculate the rotation for this region based on the shift and the previous phase values

        const double phaseAngle = shift * phaseConst;
        const float real = static_cast<float>(std::cos(phaseAngle));
        const float imag = static_cast<float>(std::sin(phaseAngle));
        const std::complex<float> rotate = std::complex<float>(real, imag) * phase[peak];

        // Calculate the first bin destination
        
        bool conjugate = false;
        float destination = static_cast<float>(regionStart) + shift;

        // Handle negative frequencies

        if (regionStart < -shift)
        {
            if (reflect != ReflectMode::None)
            {
                destination = -destination;
                conjugate = !conjugate;
            }
            else
            {
                // N.B. we allow an extra bin so that we paste the DC amount

                regionStart = static_cast<int>(std::ceil(-shift)) + 1;
                destination = static_cast<float>(regionStart) + shift;
            }
        }

        // Handle frequencies beyond the Nyquist frequency

        if (destination > static_cast<float>(numBins - 1))
        {
           if (reflect != ReflectMode::None)
            {
                float twoN = 2.f * static_cast<float>(numBins - 1);
                float wrapped = std::fmod(destination, twoN);
                if (wrapped > static_cast<float>(numBins - 1))
                {
                    wrapped = twoN - wrapped;
                    conjugate = !conjugate;
                }
                destination = wrapped;
            }
            else
            {
                // Nothing to shift, so we can skip this region entirely

                continue;
            }
        }

        int idx = static_cast<int>(std::floor(destination));
        float fraction = destination - static_cast<float>(idx);

        // Shift this region and paste the the output in one or two passes depending on the mode

        for (int i = regionStart; i < regionEnd; )
        {
            if (!conjugate)
            {
                // Set the loop limit for pasting bins in the reverse direction, taking into account the reflection mode

                int loop = std::min(regionEnd, regionStart + (numBins - (reflect != ReflectMode::Reflect ? 1 : 2)) - idx);

                // Paste bins in the forward direction, handling edge cases for DC and Nyquist bins

                if (reflect == ReflectMode::None && idx == -1)
                    PolyFreqWarpDoBin<PasteMode::Clip>(input, output, i++, idx++, fraction, rotate, phase, conjugate);
                
                if (idx == 0 && i < regionEnd)
                    PolyFreqWarpDoBin<PasteMode::DC>(input, output, i++, idx++, fraction, rotate, phase, conjugate);
            
                for ( ; i < loop; )
                    PolyFreqWarpDoBin<PasteMode::Add>(input, output, i++, idx++, fraction, rotate, phase, conjugate);

                if (idx == numBins - 2 && i < regionEnd)
                    PolyFreqWarpDoBin<PasteMode::Nyquist>(input, output, i++, idx++, fraction, rotate, phase, conjugate);

                if (reflect != ReflectMode::Reflect && idx == numBins - 1 && i < regionEnd)
                {
                    PolyFreqWarpDoBin<PasteMode::Clip>(input, output, i++, idx++, fraction, rotate, phase, conjugate);
                    break;
                }
            }
            else
            {
                // Set the loop limit for pasting bins in the reverse direction, taking into account the reflection mode

                int loop = std::min(regionEnd, regionStart + idx + (reflect == ReflectMode::None ? 1 : 2));
                
                // Paste bins in the reverse direction, handling edge cases for DC and Nyquist bins

                if (reflect != ReflectMode::Reflect && idx == numBins - 1)
                    PolyFreqWarpDoBin<PasteMode::Clip>(input, output, i++, idx--, fraction, rotate, phase, conjugate);
                
                if (idx == numBins - 2 && i < regionEnd)
                    PolyFreqWarpDoBin<PasteMode::Nyquist>(input, output, i++, idx--, fraction, rotate, phase, conjugate);
                
                for ( ; i < loop; )
                    PolyFreqWarpDoBin<PasteMode::Add>(input, output, i++, idx--, fraction, rotate, phase, conjugate);

                if (idx == 0 && i < regionEnd)
                    PolyFreqWarpDoBin<PasteMode::DC>(input, output, i++, idx--, fraction, rotate, phase, conjugate);

                if (reflect == ReflectMode::None && idx == -1 && i < regionEnd)
                {
                    PolyFreqWarpDoBin<PasteMode::Clip>(input, output, i++, idx--, fraction, rotate, phase, conjugate);
                    break;
                }
            }

            conjugate = !conjugate;
            fraction = 1.f - fraction;
        }

        // Update the region start and peak for the next iteration

        regionStart = regionEnd;
        peak = nextPeak;
    }
}

// Retrieve the SndBuf pointer from the unit and buffer number

static inline SndBuf* GetSndBuf(Unit* unit, float fbufnum)
{
    // If the buffer number is negative, return nullptr

    if (fbufnum < 0.f)
        return nullptr;

    uint32 ibufnum = (uint32)fbufnum;
    World* world = unit->mWorld;
    SndBuf* buf;

    // If the buffer number is greater than or equal to the number of sound buffers in the world,
    // check if it is a local buffer in the parent graph. If so, return the local buffer. 
    // Otherwise, return the main sound buffers.

    if (ibufnum >= world->mNumSndBufs)
    {
        int localBufNum = ibufnum - world->mNumSndBufs;
        Graph* parent = unit->mParent;
        if (localBufNum <= parent->localBufNum)
            buf = parent->mLocalSndBufs + localBufNum;
        else
            buf = world->mSndBufs;
    }
    else
        buf = world->mSndBufs + ibufnum;
    
    return buf;
}

static void PolyFreqWarpNext(PV_PolyFreqWarp* unit, int)
{
    // Get the main input buffer from the unit. If it is null, return early.

    const float fmain = ZIN0(0);                                                                                                         \

    SndBuf* buf = GetSndBuf(unit, fmain);

    if (!buf)
    {
        ZOUT0(0) = -1.f;
        return;
    }

    // Set the output to be the input buffer
 
    ZOUT0(0) = fmain; 

    // Calculate the number of bins based on the input buffer

    const int numBins = (buf->samples >> 1) + 1;

    if (numBins <= 0 || !PolyFreqWarpEnsureState(unit, numBins))
    {
        ZOUT0(0) = -1.f;
        return;
    }

    // Read parameters from the input buffer. The first six parameters are polynomial coefficients,
    // followed by a boolean for reflection, a detector buffer index, and an overlap value.

    float params[6];
    params[0] = ZIN0(1);
    params[1] = ZIN0(2);
    params[2] = ZIN0(3);
    params[3] = ZIN0(4);
    params[4] = ZIN0(5);
    params[5] = ZIN0(6);

    const float freflect = ZIN0(7);
    const float fdetector = ZIN0(8);
    const float foverlap = ZIN0(9);
    const float overlap = foverlap > 0.f ? foverlap : 4.f;

    const ReflectMode reflect = freflect == 0.f ? ReflectMode::None : freflect == 2.f ? ReflectMode::Reflect : ReflectMode::Negative;

    SndBuf* detBuf = buf;

    // If a detector buffer is specified, use it for magnitude analysis. Otherwise, use the main buffer.

    if (fdetector >= 0.f)
    {
        SndBuf* b = GetSndBuf(unit, fdetector);
        if (b && b->data && b->samples == buf->samples)
            detBuf = b;
    }

    LOCK_SNDBUF2(buf, detBuf);

    // Ensure that the buffers are in complex format. If not, convert them to complex format.

    SCComplexBuf* mainComplex = ToComplexApx(buf);
    SCComplexBuf* detComplex = ToComplexApx(detBuf);

    // Check that the buffers match in size.

    if (buf->samples != detBuf->samples)
        return;

    // Prepare the magnitudes array for peak detection. We pad the magnitudes to simplify edge handling during peak detection.
    // This gives accurate results when the peaks are near the edges of the frequency bins.

    const int padding = 4;
    float* paddedMagnitudes = static_cast<float*>(alloca((numBins + 2 * padding) * sizeof(float)));
    float* magnitudes = paddedMagnitudes + padding;

    // Calculate the magnitudes of the complex bins for peak detection.

    magnitudes[0] = detComplex->dc * detComplex->dc;
    for (int i = 0; i < numBins - 2; ++i)
    {
        const float real = detComplex->bin[i].real;
        const float imag = detComplex->bin[i].imag;
        magnitudes[i + 1] = real * real + imag * imag;
    }
    magnitudes[numBins - 1] = detComplex->nyq * detComplex->nyq;

    // Fold the magnitudes at the edges to handle boundary conditions for peak detection

    for (int j = 1; j <= padding; ++j)
    {
        magnitudes[-j] = magnitudes[j];
        magnitudes[numBins - 1 + j] = magnitudes[numBins - 1 - j];
    }

    // Make temporary arrays for the inputs/outputs

    SCComplex* input = static_cast<SCComplex*>(alloca(numBins * sizeof(SCComplex)));
    SCComplex* output = static_cast<SCComplex*>(alloca(numBins * sizeof(SCComplex)));

    // Copy input from the main complex buffer to the temporary input array

    input[0] = Complex(mainComplex->dc, 0.f);
    std::copy_n(mainComplex->bin, numBins - 2, input + 1);
    input[numBins - 1] = Complex(mainComplex->nyq, 0.f);

    // Process
    
    PolyFreqWarpProcess(input, output, unit->mPhase, magnitudes, numBins, params, reflect, overlap);

    // Copy output back to the main complex buffer

    mainComplex->dc = output[0].real;
    std::copy_n(output + 1, numBins - 2, mainComplex->bin); 
    mainComplex->nyq = output[numBins - 1].real;
}

static void PV_PolyFreqWarp_Ctor(PV_PolyFreqWarp* unit)
{
    SETCALC(PolyFreqWarpNext);
    ZOUT0(0) = ZIN0(0);
    unit->mNumBins = 0;
    unit->mPhase = nullptr;
}

static void PV_PolyFreqWarp_Dtor(PV_PolyFreqWarp* unit)
{
    PolyFreqWarpFree(unit);
}

PluginLoad(PolyFreqWarp)
{
    ft = inTable;
    DefineDtorUnit(PV_PolyFreqWarp);
}
