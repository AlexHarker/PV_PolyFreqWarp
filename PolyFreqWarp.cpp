#include "SC_PlugIn.h"
#include "FFT_UGens.h"

#include <algorithm>
#include <cmath>
#include <cstring>

InterfaceTable* ft;

struct PV_PolyFreqWarp : PV_Unit
{
    int mNumBins;
    float* mPrevious;
    float* mPhaseReal;
    float* mPhaseImag;
};

static void PolyFreqWarpFree(PV_PolyFreqWarp* unit)
{
    RTFree(unit->mWorld, unit->mPrevious);
    RTFree(unit->mWorld, unit->mPhaseReal);
    RTFree(unit->mWorld, unit->mPhaseImag);
}

static bool PolyFreqWarpEnsureState(PV_PolyFreqWarp* unit, int numBins)
{
    if (unit->mNumBins == numBins)
        return true;

    PolyFreqWarpFree(unit);
    unit->mPrevious = static_cast<float*>(RTAlloc(unit->mWorld, numBins * sizeof(float)));
    unit->mPhaseReal = static_cast<float*>(RTAlloc(unit->mWorld, numBins * sizeof(float)));
    unit->mPhaseImag = static_cast<float*>(RTAlloc(unit->mWorld, numBins * sizeof(float)));

    if (!unit->mPrevious || !unit->mPhaseReal || !unit->mPhaseImag)
    {
        PolyFreqWarpFree(unit);
        unit->mNumBins = 0;
        return false;
    }

    std::memset(unit->mPrevious, 0, numBins * sizeof(float));
    for (int i = 0; i < numBins; ++i)
    {
        unit->mPhaseReal[i] = 1.f;
        unit->mPhaseImag[i] = 0.f;
    }
    unit->mNumBins = numBins;
    return true;
}

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

static float PolyFreqWarpPeakPosition(const float* magnitudes, int peak)
{
    const float left = std::log(std::max(magnitudes[peak - 1], 1.0e-20f));
    const float centre = std::log(std::max(magnitudes[peak], 1.0e-20f));
    const float right = std::log(std::max(magnitudes[peak + 1], 1.0e-20f));
    const float divisor = left + right - 2.f * centre;
    if (divisor == 0.f)
        return static_cast<float>(peak);
    return static_cast<float>(peak) + 0.5f * (left - right) / divisor;
}

static void PolyFreqWarpAdd(SCComplex* output, int numBins, float position, float real, float imag)
{
    if (position < 0.f || position > static_cast<float>(numBins - 1))
        return;
    const int lower = static_cast<int>(std::floor(position));
    const float fraction = position - static_cast<float>(lower);
    output[lower].real += real * (1.f - fraction);
    output[lower].imag += imag * (1.f - fraction);
    if (lower + 1 < numBins)
    {
        output[lower + 1].real += real * fraction;
        output[lower + 1].imag += imag * fraction;
    }
}

static void PolyFreqWarpProcessChannel(const SCComplex* input, SCComplex* output, float* previous,
                                       float* phaseReal, float* phaseImag, const float* magnitudes,
                                       int numBins, const float* params, float normShift,
                                       bool reflect, float overlap)
{
    const float constMultVal = 2.f * static_cast<float>(M_PI) / overlap;

    for (int i = 0; i < numBins; ++i)
    {
        output[i].real = 0.f;
        output[i].imag = 0.f;
    }

    int regionStart = 0;
    int peak = PolyFreqWarpPeak(magnitudes, numBins, 0);
    if (peak == -1)
        peak = 0;

    while (regionStart < numBins)
    {
        int nextPeak = (peak < numBins - 1) ? PolyFreqWarpPeak(magnitudes, numBins, peak + 1) : -1;
        int regionEnd;

        if (nextPeak != -1)
        {
            auto minIt = std::min_element(magnitudes + peak, magnitudes + nextPeak + 1);
            regionEnd = static_cast<int>(std::distance(magnitudes, minIt));
            if (regionEnd <= regionStart)
                regionEnd = regionStart + 1;
        }
        else
        {
            regionEnd = numBins;
        }

        const float peakPosition = PolyFreqWarpPeakPosition(magnitudes, peak);
        const float normalizedPeak = peakPosition / static_cast<float>(numBins - 1);
        const float peakPolynomial = params[4] + 10000.f * normalizedPeak *
            (params[3] + normalizedPeak *
            (params[2] + normalizedPeak * (params[1] + normalizedPeak * params[0])));
        const float shift = peakPosition * peakPolynomial + normShift - peakPosition;

        const float phaseAngle = shift * constMultVal;
        const float tempReal = std::cos(phaseAngle);
        const float tempImag = std::sin(phaseAngle);

        const float rotateReal = tempReal * phaseReal[peak] - tempImag * phaseImag[peak];
        const float rotateImag = tempReal * phaseImag[peak] + tempImag * phaseReal[peak];

        for (int i = regionStart; i < regionEnd; ++i)
        {
            const float sourceReal = input[i].real;
            const float sourceImag = input[i].imag;
            const float real = sourceReal * rotateReal - sourceImag * rotateImag;
            const float imag = sourceReal * rotateImag + sourceImag * rotateReal;
            float destination = static_cast<float>(i) + shift;
            bool conjugate = false;

            if (reflect)
            {
                if (destination < 0.f)
                {
                    destination = -destination;
                    conjugate = !conjugate;
                }
                if (destination > static_cast<float>(numBins - 1))
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
            }
            float outReal = real;
            float outImag = conjugate ? -imag : imag;
            PolyFreqWarpAdd(output, numBins, destination, outReal, outImag);
            phaseReal[i] = rotateReal;
            phaseImag[i] = rotateImag;
            previous[i] = sourceReal * sourceReal + sourceImag * sourceImag;
        }
        regionStart = regionEnd;
        peak = nextPeak;
    }
}

static inline SndBuf* GetSndBuf(Unit* unit, float fbufnum)
{
    if (fbufnum < 0.f)
        return nullptr;

    uint32 ibufnum = (uint32)fbufnum;
    World* world = unit->mWorld;
    SndBuf* buf;

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
    PV_GET_BUF

    if (!buf)
        return;

    // `numbins` from the SC macro covers only bins other than DC/Nyquist. 
    // The DC and Nyquist bins are handled separately, so we add 2 to get the total number of bins.

    const int numBins = numbins + 2;

    if (numBins <= 0 || !PolyFreqWarpEnsureState(unit, numBins))
        return;

    float params[6];
    params[0] = ZIN0(1);
    params[1] = ZIN0(2);
    params[2] = ZIN0(3);
    params[3] = ZIN0(4);
    params[4] = ZIN0(5);
    params[5] = ZIN0(6);

    const bool reflect = ZIN0(7) != 0.f;
    const float fdetector = ZIN0(8);
    const float foverlap = ZIN0(9);
    const float overlap = foverlap > 0.f ? foverlap : 4.f;

    SndBuf* detBuf = buf;
    if (fdetector >= 0.f)
    {
        SndBuf* b = GetSndBuf(unit, fdetector);
        if (b && b->data && b->samples == buf->samples)
        {
            detBuf = b;
        }
    }

    const float normShift = params[5] * numBins;

    SCComplexBuf* mainComplex = ToComplexApx(buf);
    SCComplexBuf* detComplex = ToComplexApx(detBuf);

    const int padding = 4;
    float* paddedMagnitudes = static_cast<float*>(alloca((numBins + 2 * padding) * sizeof(float)));
    float* magnitudes = paddedMagnitudes + padding;

    magnitudes[0] = detComplex->dc * detComplex->dc;
    for (int i = 0; i < numBins - 2; ++i)
    {
        const float r = detComplex->bin[i].real;
        const float im = detComplex->bin[i].imag;
        magnitudes[i + 1] = r * r + im * im;
    }
    magnitudes[numBins - 1] = detComplex->nyq * detComplex->nyq;

    for (int j = 1; j <= padding; ++j)
    {
        magnitudes[-j] = magnitudes[j];
        magnitudes[numBins - 1 + j] = magnitudes[numBins - 1 - j];
    }

    // Make temporary arrays for the inputs/outputs

    SCComplex* input = static_cast<SCComplex*>(alloca(numBins * sizeof(SCComplex)));
    SCComplex* output = static_cast<SCComplex*>(alloca(numBins * sizeof(SCComplex)));

    // Copy input from the main complex buffer to the temporary input array

    input[0].real = mainComplex->dc;
    input[0].imag = 0.f;

    for (int i = 0; i < numBins - 2; ++i)
        input[i + 1] = mainComplex->bin[i];

    input[numBins - 1].real = mainComplex->nyq;
    input[numBins - 1].imag = 0.f;

    // Process the channel
    
    PolyFreqWarpProcessChannel(input, output, unit->mPrevious, unit->mPhaseReal,
                               unit->mPhaseImag, magnitudes, numBins, params, normShift,
                               reflect, overlap);

    // Copy output back to the main complex buffer

    mainComplex->dc = output[0].real;
    for (int i = 0; i < numBins - 2; ++i)
        mainComplex->bin[i] = output[i + 1];
    mainComplex->nyq = output[numBins - 1].real;
}

static void PV_PolyFreqWarp_Ctor(PV_PolyFreqWarp* unit)
{
    SETCALC(PolyFreqWarpNext);
    ZOUT0(0) = ZIN0(0);
    unit->mNumBins = 0;
    unit->mPrevious = nullptr;
    unit->mPhaseReal = nullptr;
    unit->mPhaseImag = nullptr;
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
