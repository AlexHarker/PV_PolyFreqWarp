#include "SC_PlugIn.h"
#include "FFT_UGens.h"

#include <cmath>
#include <cstring>

InterfaceTable* ft;

struct PV_PolyFreqWarp : PV_Unit
{
    int mNumBins;
    float* mPreviousLeft;
    float* mPreviousRight;
    float* mPhaseRealLeft;
    float* mPhaseImagLeft;
    float* mPhaseRealRight;
    float* mPhaseImagRight;
};

static void PolyFreqWarpFree(PV_PolyFreqWarp* unit)
{
    RTFree(unit->mWorld, unit->mPreviousLeft);
    RTFree(unit->mWorld, unit->mPreviousRight);
    RTFree(unit->mWorld, unit->mPhaseRealLeft);
    RTFree(unit->mWorld, unit->mPhaseImagLeft);
    RTFree(unit->mWorld, unit->mPhaseRealRight);
    RTFree(unit->mWorld, unit->mPhaseImagRight);
}

static bool PolyFreqWarpEnsureState(PV_PolyFreqWarp* unit, int numBins)
{
    if (unit->mNumBins == numBins)
        return true;

    PolyFreqWarpFree(unit);
    unit->mPreviousLeft = static_cast<float*>(RTAlloc(unit->mWorld, numBins * sizeof(float)));
    unit->mPreviousRight = static_cast<float*>(RTAlloc(unit->mWorld, numBins * sizeof(float)));
    unit->mPhaseRealLeft = static_cast<float*>(RTAlloc(unit->mWorld, numBins * sizeof(float)));
    unit->mPhaseImagLeft = static_cast<float*>(RTAlloc(unit->mWorld, numBins * sizeof(float)));
    unit->mPhaseRealRight = static_cast<float*>(RTAlloc(unit->mWorld, numBins * sizeof(float)));
    unit->mPhaseImagRight = static_cast<float*>(RTAlloc(unit->mWorld, numBins * sizeof(float)));

    if (!unit->mPreviousLeft || !unit->mPreviousRight || !unit->mPhaseRealLeft ||
        !unit->mPhaseImagLeft || !unit->mPhaseRealRight || !unit->mPhaseImagRight)
    {
        PolyFreqWarpFree(unit);
        unit->mNumBins = 0;
        return false;
    }

    std::memset(unit->mPreviousLeft, 0, numBins * sizeof(float));
    std::memset(unit->mPreviousRight, 0, numBins * sizeof(float));
    for (int i = 0; i < numBins; ++i)
    {
        unit->mPhaseRealLeft[i] = 1.f;
        unit->mPhaseRealRight[i] = 1.f;
        unit->mPhaseImagLeft[i] = 0.f;
        unit->mPhaseImagRight[i] = 0.f;
    }
    unit->mNumBins = numBins;
    return true;
}

static int PolyFreqWarpPeak(const float* magnitudes, int numBins, int start)
{
    for (int i = start; i < numBins; ++i)
    {
        const float left = i ? magnitudes[i - 1] : magnitudes[i];
        const float right = i + 1 < numBins ? magnitudes[i + 1] : magnitudes[i];
        if (magnitudes[i] >= left && magnitudes[i] >= right && magnitudes[i] > 0.f)
            return i;
    }
    return numBins - 1;
}

static float PolyFreqWarpPeakPosition(const float* magnitudes, int numBins, int peak)
{
    if (peak <= 0 || peak + 1 >= numBins)
        return static_cast<float>(peak);

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
                                       bool reflect)
{
    const float binRadians = 2.f * static_cast<float>(M_PI) / (2.f * (numBins - 1));

    for (int i = 0; i < numBins; ++i)
    {
        output[i].real = 0.f;
        output[i].imag = 0.f;
    }

    int regionStart = 0;
    while (regionStart < numBins)
    {
        const int peak = PolyFreqWarpPeak(magnitudes, numBins, regionStart);
        int regionEnd = peak + 1;
        while (regionEnd < numBins && magnitudes[regionEnd] <= magnitudes[regionEnd - 1])
            ++regionEnd;
        if (regionEnd <= regionStart)
            regionEnd = regionStart + 1;

        const float peakPosition = PolyFreqWarpPeakPosition(magnitudes, numBins, peak);
        const float normalizedPeak = peakPosition / static_cast<float>(numBins - 1);
        const float peakPolynomial = params[4] + 10000.f * normalizedPeak *
            (params[3] + normalizedPeak *
            (params[2] + normalizedPeak * (params[1] + normalizedPeak * params[0])));
        const float shift = peakPosition * peakPolynomial + normShift - peakPosition;
        const float mappedPeak = peakPosition + shift;
        const float phase = mappedPeak * binRadians;
        const float rotateReal = std::cos(phase) * phaseReal[peak] - std::sin(phase) * phaseImag[peak];
        const float rotateImag = std::sin(phase) * phaseReal[peak] + std::cos(phase) * phaseImag[peak];

        for (int i = regionStart; i < regionEnd; ++i)
        {
            const float sourceReal = input[i].real;
            const float sourceImag = input[i].imag;
            const float real = sourceReal * rotateReal - sourceImag * rotateImag;
            const float imag = sourceReal * rotateImag + sourceImag * rotateReal;
            float destination = static_cast<float>(i) + shift;

            if (reflect)
            {
                float wrapped = std::fmod(std::fabs(destination), 2.f * (numBins - 1));
                if (wrapped > numBins - 1)
                    wrapped = 2.f * (numBins - 1) - wrapped;
                destination = wrapped;
            }
            PolyFreqWarpAdd(output, numBins, destination, real, imag);
            phaseReal[i] = rotateReal;
            phaseImag[i] = rotateImag;
            previous[i] = sourceReal * sourceReal + sourceImag * sourceImag;
        }
        regionStart = regionEnd;
    }
}

static void PolyFreqWarpNext(PV_PolyFreqWarp* unit, int inNumSamples)
{
    PV_GET_BUF2
    (void)inNumSamples;
    if (!buf1 || !buf2 || buf1->samples != buf2->samples)
        return;

    // `numbins` from the SC macro covers only bins other than DC/Nyquist. 
    // The DC and Nyquist bins are handled separately, so we add 2 to get the total number of bins.

    const int numBins = numbins + 2;

    if (numBins <= 0 || !PolyFreqWarpEnsureState(unit, numBins))
        return;

    float params[6];
    params[0] = ZIN0(2);
    params[1] = ZIN0(3);
    params[2] = ZIN0(4);
    params[3] = ZIN0(5);
    params[4] = ZIN0(6);
    params[5] = ZIN0(7);

    const float link = ZIN0(8);
    const bool reflect = ZIN0(9) != 0.f;
    const float normShift = params[5] * numBins;

    SCComplexBuf* left = ToComplexApx(buf1);
    SCComplexBuf* right = ToComplexApx(buf2);

    float* leftMagnitudes = static_cast<float*>(alloca(numBins * sizeof(float)));
    float* rightMagnitudes = static_cast<float*>(alloca(numBins * sizeof(float)));
    
    leftMagnitudes[0] = link ? left->dc : left->dc + right->dc;
    rightMagnitudes[0] = link ? right->dc : left->dc + right->dc;

    for (int i = 0; i < numBins - 2; ++i)
    {
        const float leftReal = left->bin[i].real;
        const float rightReal = right->bin[i].real;
        const float leftImag = left->bin[i].imag;
        const float rightImag = right->bin[i].imag;
        if (link)
        {
            leftMagnitudes[i + 1] = rightMagnitudes[i + 1] =
                (leftReal + rightReal) * (leftReal + rightReal) +
                (leftImag + rightImag) * (leftImag + rightImag);
        }
        else
        {
            leftMagnitudes[i + 1] = leftReal * leftReal + leftImag * leftImag;
            rightMagnitudes[i + 1] = rightReal * rightReal + rightImag * rightImag;
        }
    }

    leftMagnitudes[numBins - 1] = link ? left->nyq : left->nyq + right->nyq;
    rightMagnitudes[numBins - 1] = link ? right->nyq : left->nyq + right->nyq;

    // Make temporary arrays for the inputs/outputs

    SCComplex* leftInput = static_cast<SCComplex*>(alloca(numBins * sizeof(SCComplex)));
    SCComplex* rightInput = static_cast<SCComplex*>(alloca(numBins * sizeof(SCComplex)));
    SCComplex* leftOutput = static_cast<SCComplex*>(alloca(numBins * sizeof(SCComplex)));
    SCComplex* rightOutput = static_cast<SCComplex*>(alloca(numBins * sizeof(SCComplex)));

    // Copy the input bins to the leftInput and rightInput arrays for processing.

    leftInput[0].real = left->dc;
    leftInput[0].imag = 0.f;
    rightInput[0].real = right->dc;
    rightInput[0].imag = 0.f;

    for (int i = 0; i < numBins - 2; ++i)
    {
        leftInput[i + 1] = left->bin[i];
        rightInput[i + 1] = right->bin[i];
    }

    leftInput[numBins - 1].real = left->nyq;
    leftInput[numBins - 1].imag = 0.f;
    rightInput[numBins - 1].real = right->nyq;
    rightInput[numBins - 1].imag = 0.f;

    PolyFreqWarpProcessChannel(leftInput, leftOutput, unit->mPreviousLeft, unit->mPhaseRealLeft,
                               unit->mPhaseImagLeft, leftMagnitudes, numBins, params, normShift,
                               reflect);
    PolyFreqWarpProcessChannel(rightInput, rightOutput, unit->mPreviousRight, unit->mPhaseRealRight,
                               unit->mPhaseImagRight, rightMagnitudes, numBins, params, normShift,
                               reflect);

    left->dc = leftOutput[0].real;
    right->dc = rightOutput[0].real;

    for (int i = 0; i < numBins - 2; ++i)
    {
        left->bin[i] = leftOutput[i + 1];
        right->bin[i] = rightOutput[i + 1];
    }

    left->nyq = leftOutput[numBins - 1].real;
    right->nyq = rightOutput[numBins - 1].real;
}

static void PV_PolyFreqWarp_Ctor(PV_PolyFreqWarp* unit)
{
    SETCALC(PolyFreqWarpNext);
    ZOUT0(0) = ZIN0(0);
    ZOUT0(1) = ZIN0(1);
    unit->mNumBins = 0;
    unit->mPreviousLeft = nullptr;
    unit->mPreviousRight = nullptr;
    unit->mPhaseRealLeft = nullptr;
    unit->mPhaseImagLeft = nullptr;
    unit->mPhaseRealRight = nullptr;
    unit->mPhaseImagRight = nullptr;
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
