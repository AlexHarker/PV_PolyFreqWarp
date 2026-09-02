
// NB - July 2008 - Really shouldn't be correcting zero amplitudes to 1 for log estimate parabolic interpolation - should be using a v. small value instead.....

// NB. By using Macros for Vec operations we can test the platform on compile and therefore compile the same code for altivec and sse wherever possible 

#if (TARGET_RT_LITTLE_ENDIAN)
#define VEC_ROT_OP _mm_shuffle_ps
#define VEC_ROT_VAL 0x4E
#else
#define VEC_ROT_OP vec_sld
#define VEC_ROT_VAL 8
#endif


#define FFTRESET FreqWarp_Reset 
#define PROCESSPERFORM FreqWarp_perform

#include <ext.h>
#include <z_dsp.h>
#include <stdlib.h>
#include <math.h>


#include "FFTDataStereo.h"

/* 
FreqWarp~ is an object to apply an abitary frequency warping function to incoming fft data.
*/

void *this_class;

typedef struct _FreqWarp
{
    t_pxobject x_obj;
	FFTDataSt;
	
	double ConstMultVal, Params[6];
	float *Data[9];
	char ChannelLink, Reflect, NoMemory;
	
	void *proxies[10];
	long InletNumber;
	
} t_FreqWarp;

/* 
	
	N.B. - Data Pointers are:
	
	[0] - In Real 1 (stored)
	[1] - In Complex 1 (stored)
	[2] - In Real 2 (stored)
	[3] - In Complex 2 (stored)
	[4] - In Amps Squared
	[5] - Phase Accum Real 1
	[6] - Phase Accum Complex 1
	[7] - Phase Accum Real 2
	[8] - Phase Accum Complex 2
	
*/


void *FreqWarp_new();
void FreqWarp_Reset (t_FreqWarp *x);
void FreqWarp_flt (t_FreqWarp *x, double NewVal);
void FreqWarp_int (t_FreqWarp *x, long Val);
void FreqWarp_Reflect(t_FreqWarp *x, long Reflect);
void FreqWarp_perform(t_FreqWarp *x, COMPLEX_SPLIT *FFTProcess, int vectsize);
void FreqWarp_dsp(t_FreqWarp *x, t_signal **sp, short *count);
void FreqWarp_free(t_FreqWarp *x);
void FreqWarp_assist (t_FreqWarp *x, void *b, long m, long a, char *s);

#define T_EXTNAME t_FreqWarp
#include "FFTStereo.h"		

void main(void)
{
    setup((t_messlist **)&this_class, (method) FreqWarp_new, (method)FreqWarp_free, (short)sizeof(t_FreqWarp), 0L, 0);
    addmess ((method)FreqWarp_assist, "assist", A_CANT, 0);
	addfloat((method)FreqWarp_flt);
	addint((method)FreqWarp_int);
	addmess((method)FreqWarp_Reflect, "Reflect", 	A_LONG, 0);
	FFTMain();
    dsp_initclass();
}


void *FreqWarp_new()
{
    t_FreqWarp *x = (t_FreqWarp *)newobject(this_class);
	
	long i;
    for (i = 10; i > 0; i--)												// Make Proxies
        x->proxies[i - 1] = proxy_new(x, i ,&x->InletNumber);
		
    dsp_setup((t_pxobject *)x, 2);
    outlet_new((t_object *)x,"signal");
    outlet_new((t_object *)x,"signal");
	
	for (i = 0; i < 6; i++)
		x->Params[i] = 0;
		
	x->Params[4] = 1;
	
	x->Reflect = 1;
	x->ChannelLink = 0;
	
	x->NoMemory = 0;
	
	x->Data[0] = (float *) malloc((5 * MaxFFTOVER2 + 16) * sizeof (float));	// N.B. Changed to + 16 rather than 8 due to new peak finding routine for linked stereo
	x->Data[1] = x->Data[0] + MaxFFTOVER2;
	x->Data[2] = x->Data[1] + MaxFFTOVER2;
	x->Data[3] = x->Data[2] + MaxFFTOVER2;
	x->Data[4] = x->Data[3] + MaxFFTOVER2;
	x->Data[5] = (float *) malloc(MaxFFTOVER2 * 4 * sizeof (float));
	x->Data[6] = x->Data[5] + MaxFFTOVER2;
	x->Data[7] = x->Data[6] + MaxFFTOVER2;
	x->Data[8] = x->Data[7] + MaxFFTOVER2;
	
	FFTInitSt(11, 2);
	
	if (x->Data[0] == 0 || x->Data[5] == 0)
	{
		post ("Could Not Allocate Memory for Freq Warper - Bailing Out....");
		x->FFTOff = 1;
		x->NoMemory = 1;
		goto getout;
	}
	
	getout:
    return (x);
}


void FreqWarp_Reset (t_FreqWarp *x)
{
	if (!x->NoMemory)
	{
		int i;
		float *PhaseAccumR1 = x->Data[5];
		float *PhaseAccumI1 = x->Data[6];
		float *PhaseAccumR2 = x->Data[7];
		float *PhaseAccumI2 = x->Data[8];
		
		for (i = 0; i < x->FFTParam[1]; i++) // x->FFTParam[1]
		{
			PhaseAccumR1[i] = PhaseAccumR2[i] = 1;
			PhaseAccumI1[i] = PhaseAccumI2[i] = 0;
		}
		
		x->ConstMultVal = (2 * PI / ((double)x->FFTParam[3] * ZeroPad[1]));
	}
}

void FreqWarp_free(t_FreqWarp *x)
{
	dsp_free(&x->x_obj);
	
	if (x->Data[0]) free (x->Data[0]);
	if (x->Data[5]) free (x->Data[5]);
	
	int i;
    for (i = 0; i < 10; i++)
        freeobject((t_object *) x->proxies[i]);
        
    FFTFree();
}


void FreqWarp_flt (t_FreqWarp *x, double NewVal)
{
	switch (proxy_getinlet((t_object *)x))
	{
		case 5:
			x->Params[0] = NewVal;
			break;
		case 6:
			x->Params[1] = NewVal;
			break;
		case 7:
			x->Params[2] = NewVal;
			break;
		case 8:
			x->Params[3] = NewVal;
			break;
		case 9:
			x->Params[4] = NewVal;
			break;
		case 10:
			x->Params[5] = NewVal;
			break;
	}
}


void FreqWarp_int (t_FreqWarp *x, long Val)
{
	if (!x->NoMemory)
	{
		switch (proxy_getinlet((t_object *)x))
		{
			case 1:
				FFTNewSizeSt(Val);
				break;
			case 2:
				FFTNewOverlapSt(Val);
				break;
			case 3:
				FFTNewWindow(Val);
				break;
			case 4:
				switch (Val)
				{
					case 0:
						x->ChannelLink = 0;
						break;
					case 1:
						x->ChannelLink = 1;
						break;
				}
				break;
		}
	}
}


void FreqWarp_Reflect(t_FreqWarp *x, long Reflect)
{
	if (Reflect)
		x->Reflect = 1;
	else
		x->Reflect = 0;
}


void FreqWarp_perform(t_FreqWarp *x, COMPLEX_SPLIT *FFTProcess, int vectsize)
{	
	vFloat *in1 = (vFloat *) FFTProcess[0].realp;
	vFloat *in2 = (vFloat *) FFTProcess[0].imagp;
	vFloat *in3 = (vFloat *) FFTProcess[1].realp;
	vFloat *in4 = (vFloat *) FFTProcess[1].imagp;
	
	/* 
	
	N.B. - Float Pointers are:					// NOTE - Think this is up to date - However, VecPointers are different!
	
	Left 
	
	[0] - Out Real
	[1] - Out Imaginary
	[2] - In Real (stored)
	[3] - In Imaginary (stored)
	[4]	- Phase Acuum R
	[5] - Phase Accum I
	
	Right 
	
	[6] - Out Real
	[7] - Out Imaginary
	[8] - In Real (stored)
	[9] - In Imaginary (stored)
	[10] - Phase Acuum R
	[11] - Phase Accum I
	
	*/
	
	float **Data = x->Data;
	
	float *FltPointers[12];
	FltPointers[0] = (float *) FFTProcess[2].realp;
	FltPointers[1] = (float *) FFTProcess[2].imagp;
	FltPointers[2] = Data[0];
	FltPointers[3] = Data[1];
	FltPointers[4] = Data[5];
	FltPointers[5] = Data[6];
	FltPointers[6] = (float *) FFTProcess[3].realp;
	FltPointers[7] = (float *) FFTProcess[3].imagp;
	FltPointers[8] = Data[2];
	FltPointers[9] = Data[3];
	FltPointers[10] = Data[7];
	FltPointers[11] = Data[8];
		
	vFloat *VecPointers[8];
	VecPointers[0] = (vFloat *) FltPointers[0];
	VecPointers[1] = (vFloat *) FltPointers[1];
	VecPointers[2] = (vFloat *) FltPointers[2];
	VecPointers[3] = (vFloat *) FltPointers[3];
	VecPointers[4] = (vFloat *) FltPointers[6];
	VecPointers[5] = (vFloat *) FltPointers[7];
	VecPointers[6] = (vFloat *) FltPointers[8];
	VecPointers[7] = (vFloat *) FltPointers[9];

	float *AmpsSq = Data[4];
	vFloat *VecAmpsSq = (vFloat *) AmpsSq;
	
	double *Params = x->Params;
	double ConstMultVal = x->ConstMultVal;
	
	VecSplat BinShiftFract, LogEstimate, PhaseRotate[3];
	double ActShift, Wrap, Temp[2]; 
	float PeakAmp, NextPeakAmp, MinAmp, PeakEstDivisor, PeakCorrection, AcPeakLoc;
	int BinShiftInt, NewBin, PeakLoc, NextPeakLoc, StartOfRegion, NextStartOfRegion, i;
	
	float VectSizeRecip = 1 / (double) vectsize;
	float NormShiftParam = Params[5] * (float) vectsize;
	
	vFloat TempVec[2];
	vFloat Zero = (vFloat) {0.0,0.0,0.0,0.0};

	int FromUp, ToUp, FromDown, ToDown, Conj, NewBin2;
	char DoTopDown, DoTopUp; 
	char Reflect = x->Reflect;
	int TwoVect = 2 * vectsize;
	double TwoVectRecip = (double) 1 / (double) TwoVect; 
	double WrapDivisor, NormPeakLoc;
	
	// Load in data and zero the output array
	
	for (i = 0; i < vectsize >> 2; i++)
	{
		VecPointers[2][i] = *in1++;
		VecPointers[3][i] = *in2++;
		VecPointers[6][i] = *in3++;
		VecPointers[7][i] = *in4++;
		VecPointers[0][i] = Zero;
		VecPointers[1][i] = Zero;
		VecPointers[4][i] = Zero;
		VecPointers[5][i] = Zero;
	}
			
	if (x->ChannelLink)
	{
		///////////////////////////////////////// Channels Linked////////////////////////////////////
		
		// Load in amplitudes 
		
		for (i = 0; i < vectsize >> 2; i++)
		{
			TempVec[0] = VEC_ADD_OP(VecPointers[2][i], VecPointers[6][i]);																// vec_add(VecPointers[2][i], VecPointers[6][i]);
			TempVec[1] = VEC_ADD_OP(VecPointers[3][i], VecPointers[7][i]);																// vec_add(VecPointers[3][i], VecPointers[7][i]);
			VecAmpsSq[i + 1] = VEC_ADD_OP(VEC_MUL_OP(TempVec[1],TempVec[1] ZEROARG), VEC_MUL_OP(TempVec[0], TempVec[0] ZEROARG));		// vec_madd(TempVec[0], TempVec[0], vec_madd (TempVec[1],TempVec[1], Zero));
		}		
		
		// Search for first peak
		
		for (i = 0; i < vectsize; i++)							
		{
			if (AmpsSq[i + 4] > AmpsSq[i + 3] && AmpsSq[i + 4] > AmpsSq[i + 2] && AmpsSq[i + 4] > AmpsSq[i + 5] && AmpsSq[i + 4] > AmpsSq[i + 6])
			{
				PeakAmp = AmpsSq[i + 4];
				PeakLoc = i;
				break;
			}
		}
		
		// Just in case
		
		if (i == vectsize)										
			PeakLoc = vectsize - 1;
		
		StartOfRegion = 0;
		while (StartOfRegion < vectsize)
		{
			MinAmp = PeakAmp;
			for (i = PeakLoc + 1; i < vectsize; i++)	// Search for next region
			{
				if (AmpsSq[i + 4] < MinAmp)				// Find minimums
				{
					MinAmp = AmpsSq[i + 4];
					NextStartOfRegion = i;
				}
				else									// Find peaks
				{
					if (AmpsSq[i + 4] > AmpsSq[i + 3] && AmpsSq[i + 4] > AmpsSq[i + 2] && AmpsSq[i + 4] > AmpsSq[i + 5] && AmpsSq[i + 4] > AmpsSq[i + 6])
					{
						NextPeakAmp = AmpsSq[i + 4];
						NextPeakLoc = i;
						break;
					}
				}
			}
			
			if (i >= vectsize)
				NextStartOfRegion = vectsize;
			
			// More accurate peak location (parabolic interpolation)
			
			if (AmpsSq[PeakLoc + 3]) LogEstimate.flt[0] = (AmpsSq[PeakLoc + 3]); else LogEstimate.flt[0] = 1;
			if (AmpsSq[PeakLoc + 4]) LogEstimate.flt[1] = (AmpsSq[PeakLoc + 4]); else LogEstimate.flt[1] = 1;
			if (AmpsSq[PeakLoc + 5]) LogEstimate.flt[2] = (AmpsSq[PeakLoc + 5]); else LogEstimate.flt[2] = 1;
			LogEstimate.vec = vlogf(LogEstimate.vec);
			
			PeakEstDivisor = LogEstimate.flt[0] + LogEstimate.flt[2] - (2.0 * LogEstimate.flt[1]);
			if (PeakEstDivisor == 0) PeakCorrection = 0;
			else PeakCorrection = (0.5 * (LogEstimate.flt[0] - LogEstimate.flt[2])) / PeakEstDivisor;
		
			AcPeakLoc = (float) PeakLoc + PeakCorrection;
			
			// Calculate shift 
		
			NormPeakLoc = (double) (AcPeakLoc * VectSizeRecip);
			ActShift = ((double) AcPeakLoc * (Params[4] + 10000 * (NormPeakLoc * (Params[3] + (NormPeakLoc * (Params[2] + (NormPeakLoc * (Params[1] + (NormPeakLoc * Params[0]))))))))) + NormShiftParam - AcPeakLoc;			
			
			if (ActShift == 0)
			{
				// Do Two Channels At Once
				
				for (i = StartOfRegion; i < NextStartOfRegion; i++)
				{
					FltPointers[0][i] += FltPointers[2][i];
					FltPointers[1][i] += FltPointers[3][i];
					FltPointers[4][i] = 1;
					FltPointers[5][i] = 0;
					FltPointers[6][i] += FltPointers[8][i];
					FltPointers[7][i] += FltPointers[9][i];
					FltPointers[10][i] = 1;
					FltPointers[11][i] = 0;
				}

			}
			else
			{
				if (Reflect)
				{
					// Mirror Output Around DC and Nyquist
					
					Wrap = (double) StartOfRegion + ActShift;
					Conj = 0;

					if (Wrap < 0)
					{
						Wrap = -Wrap;
						Conj = 1;
					}
					if (Wrap > (double) TwoVect)
					{
						WrapDivisor = Wrap * TwoVectRecip;
						Conj = ((int) floor (WrapDivisor) + Conj) % 2;
						Wrap -= (floor(WrapDivisor) * (double) TwoVect);
					}
					if (Wrap > (double) vectsize)
					{
						Wrap = (double) TwoVect - Wrap;
						Conj = 1 - Conj;
					}
					
					NewBin = (int) floor(Wrap);
					
					// Work Out Which Bits Go Which Way
					
					if (Conj)
					{
						FromDown = StartOfRegion;
						DoTopUp = DoTopDown = FromUp = ToUp = NewBin2 = 0;
						ToDown = StartOfRegion + NewBin + 1;
						if (ToDown < NextStartOfRegion)
						{
							FromUp = ToDown;
							ToUp = NextStartOfRegion;
						}
						else
							ToDown = NextStartOfRegion;
					}
					else
					{
						FromUp = StartOfRegion;
						DoTopUp = DoTopDown = FromDown = ToDown = 0;
						NewBin2 = NewBin;
						NewBin = vectsize - 2;
						ToUp = StartOfRegion + vectsize - 1 - NewBin2;
						if (NextStartOfRegion > ToUp)
						{
							DoTopUp = 1;
							if (NextStartOfRegion > ToUp + 1)
							{
								DoTopDown = 1;
								FromDown = ToUp + 2;
								ToDown = NextStartOfRegion;
							}
						}
						else
							ToUp = NextStartOfRegion;
					}
				}
				else
				{
					DoTopUp = DoTopDown = FromDown = ToDown = 0;
					FromUp = StartOfRegion;
					ToUp = NextStartOfRegion;
					
					BinShiftInt = (int) floor (ActShift);
					if (BinShiftInt < 0)
					{
						if (StartOfRegion + BinShiftInt < 0)
							FromUp = - BinShiftInt;
						if (BinShiftInt < -vectsize)
							FromUp = vectsize;
					}
					else
					{
						if (NextStartOfRegion + BinShiftInt > vectsize - 1)
						{
							ToUp = vectsize - BinShiftInt - 1;
							if (NextStartOfRegion + BinShiftInt > vectsize)
								DoTopUp = 1;
							if (ToUp < 0)
								ToUp = 0;
						}
					}
					NewBin2 = FromUp + BinShiftInt;
				}
				
				// Calculate shift vals and phase rotation
				
				BinShiftFract.flt[0] = BinShiftFract.flt[1] = (float) (ActShift - (double) floor(ActShift));
				BinShiftFract.flt[2] = BinShiftFract.flt[3] = 1 - BinShiftFract.flt[0];
				Temp[0] = cos((double) ActShift * ConstMultVal);
				Temp[1] = sin((double) ActShift * ConstMultVal);
				PhaseRotate[0].flt[0] = PhaseRotate[0].flt[3] = (Temp[0] * FltPointers[4][PeakLoc]) - (Temp[1] * FltPointers[5][PeakLoc]);
				PhaseRotate[0].flt[1] = (Temp[0] * FltPointers[5][PeakLoc]) + (Temp[1] * FltPointers[4][PeakLoc]);
				PhaseRotate[0].flt[2] = - PhaseRotate[0].flt[1];
				
				// Do the Shifts (Two Channels at Once)
								
				// Down
				
				for (i = FromDown; i < ToDown; i++)
				{
					PhaseRotate[1].flt[0] = PhaseRotate[1].flt[1] = FltPointers[2][i];
					PhaseRotate[1].flt[2] = PhaseRotate[1].flt[3] = FltPointers[3][i];
					PhaseRotate[1].vec = VEC_MUL_OP(PhaseRotate[1].vec, PhaseRotate[0].vec ZEROARG);											// vec_madd(PhaseRotate[1].vec, PhaseRotate[0].vec, Zero);
					PhaseRotate[1].vec = VEC_ADD_OP(PhaseRotate[1].vec, VEC_ROT_OP(PhaseRotate[1].vec, PhaseRotate[1].vec, VEC_ROT_VAL));		// vec_add(PhaseRotate[1].vec, vec_sld(PhaseRotate[1].vec, PhaseRotate[1].vec, 8));
					PhaseRotate[1].vec = VEC_MUL_OP(PhaseRotate[1].vec, BinShiftFract.vec ZEROARG);												// vec_madd(PhaseRotate[1].vec, BinShiftFract.vec, Zero);
					
					PhaseRotate[2].flt[0] = PhaseRotate[2].flt[1] = FltPointers[8][i];
					PhaseRotate[2].flt[2] = PhaseRotate[2].flt[3] = FltPointers[9][i];
					PhaseRotate[2].vec = VEC_MUL_OP(PhaseRotate[2].vec, PhaseRotate[0].vec ZEROARG);											// vec_madd(PhaseRotate[2].vec, PhaseRotate[0].vec, Zero);
					PhaseRotate[2].vec = VEC_ADD_OP(PhaseRotate[2].vec, VEC_ROT_OP(PhaseRotate[2].vec, PhaseRotate[2].vec, VEC_ROT_VAL));		// vec_add(PhaseRotate[2].vec, vec_sld(PhaseRotate[2].vec, PhaseRotate[2].vec, 8));
					PhaseRotate[2].vec = VEC_MUL_OP(PhaseRotate[2].vec, BinShiftFract.vec ZEROARG);												// vec_madd(PhaseRotate[2].vec, BinShiftFract.vec, Zero);
					
					// N.B. Conjugation + Reversal of Interpolation!
					
					FltPointers[0][NewBin + 1] += PhaseRotate[1].flt[2];
					FltPointers[1][NewBin + 1] -= PhaseRotate[1].flt[3];
					FltPointers[0][NewBin] += PhaseRotate[1].flt[0];
					FltPointers[1][NewBin] -= PhaseRotate[1].flt[1];
					FltPointers[6][NewBin + 1] += PhaseRotate[2].flt[2];
					FltPointers[7][NewBin + 1] -= PhaseRotate[2].flt[3];
					FltPointers[6][NewBin] += PhaseRotate[2].flt[0];
					FltPointers[7][NewBin] -= PhaseRotate[2].flt[1];
					
					FltPointers[4][i] = FltPointers[10][i] = PhaseRotate[0].flt[0];
					FltPointers[5][i] = FltPointers[11][i] = PhaseRotate[0].flt[1];
					
					NewBin--; 
				}
				
				// Up
				
				for (i = FromUp; i < ToUp; i++)
				{
					PhaseRotate[1].flt[0] = PhaseRotate[1].flt[1] = FltPointers[2][i];
					PhaseRotate[1].flt[2] = PhaseRotate[1].flt[3] = FltPointers[3][i];
					PhaseRotate[1].vec = VEC_MUL_OP(PhaseRotate[1].vec, PhaseRotate[0].vec ZEROARG);											// vec_madd(PhaseRotate[1].vec, PhaseRotate[0].vec, Zero);
					PhaseRotate[1].vec = VEC_ADD_OP(PhaseRotate[1].vec, VEC_ROT_OP(PhaseRotate[1].vec, PhaseRotate[1].vec, VEC_ROT_VAL));		// vec_add(PhaseRotate[1].vec, vec_sld(PhaseRotate[1].vec, PhaseRotate[1].vec, 8));
					PhaseRotate[1].vec = VEC_MUL_OP(PhaseRotate[1].vec, BinShiftFract.vec ZEROARG);												// vec_madd(PhaseRotate[1].vec, BinShiftFract.vec, Zero);
					
					PhaseRotate[2].flt[0] = PhaseRotate[2].flt[1] = FltPointers[8][i];
					PhaseRotate[2].flt[2] = PhaseRotate[2].flt[3] = FltPointers[9][i];
					PhaseRotate[2].vec = VEC_MUL_OP(PhaseRotate[2].vec, PhaseRotate[0].vec ZEROARG);											// vec_madd(PhaseRotate[2].vec, PhaseRotate[0].vec, Zero);
					PhaseRotate[2].vec = VEC_ADD_OP(PhaseRotate[2].vec, VEC_ROT_OP(PhaseRotate[2].vec, PhaseRotate[2].vec, VEC_ROT_VAL));		// vec_add(PhaseRotate[2].vec, vec_sld(PhaseRotate[2].vec, PhaseRotate[2].vec, 8));
					PhaseRotate[2].vec = VEC_MUL_OP(PhaseRotate[2].vec, BinShiftFract.vec ZEROARG);												// vec_madd(PhaseRotate[2].vec, BinShiftFract.vec, Zero);
					
					FltPointers[0][NewBin2] += PhaseRotate[1].flt[2];
					FltPointers[1][NewBin2] += PhaseRotate[1].flt[3];
					FltPointers[0][NewBin2 + 1] += PhaseRotate[1].flt[0];
					FltPointers[1][NewBin2 + 1] += PhaseRotate[1].flt[1];
					FltPointers[6][NewBin2] += PhaseRotate[2].flt[2];
					FltPointers[7][NewBin2] += PhaseRotate[2].flt[3];
					FltPointers[6][NewBin2 + 1] += PhaseRotate[2].flt[0];
					FltPointers[7][NewBin2 + 1] += PhaseRotate[2].flt[1];
					
					FltPointers[4][i] = FltPointers[10][i] = PhaseRotate[0].flt[0];
					FltPointers[5][i] = FltPointers[11][i] = PhaseRotate[0].flt[1];
					
					NewBin2++; 
				}
				
				if (DoTopUp)
				{
					PhaseRotate[1].flt[0] = PhaseRotate[1].flt[1] = FltPointers[2][ToUp];
					PhaseRotate[1].flt[2] = PhaseRotate[1].flt[3] = FltPointers[3][ToUp];
					PhaseRotate[1].vec = VEC_MUL_OP(PhaseRotate[1].vec, PhaseRotate[0].vec ZEROARG);											// vec_madd(PhaseRotate[1].vec, PhaseRotate[0].vec, Zero);
					PhaseRotate[1].vec = VEC_ADD_OP(PhaseRotate[1].vec, VEC_ROT_OP(PhaseRotate[1].vec, PhaseRotate[1].vec, VEC_ROT_VAL));		// vec_add(PhaseRotate[1].vec, vec_sld(PhaseRotate[1].vec, PhaseRotate[1].vec, 8));
					PhaseRotate[1].vec = VEC_MUL_OP(PhaseRotate[1].vec, BinShiftFract.vec ZEROARG);												// vec_madd(PhaseRotate[1].vec, BinShiftFract.vec, Zero);
					
					PhaseRotate[2].flt[0] = PhaseRotate[2].flt[1] = FltPointers[8][i];
					PhaseRotate[2].flt[2] = PhaseRotate[2].flt[3] = FltPointers[9][i];
					PhaseRotate[2].vec = VEC_MUL_OP(PhaseRotate[2].vec, PhaseRotate[0].vec ZEROARG);											// vec_madd(PhaseRotate[2].vec, PhaseRotate[0].vec, Zero);
					PhaseRotate[2].vec = VEC_ADD_OP(PhaseRotate[2].vec, VEC_ROT_OP(PhaseRotate[2].vec, PhaseRotate[2].vec, VEC_ROT_VAL));		// vec_add(PhaseRotate[2].vec, vec_sld(PhaseRotate[2].vec, PhaseRotate[2].vec, 8));
					PhaseRotate[2].vec = VEC_MUL_OP(PhaseRotate[2].vec, BinShiftFract.vec ZEROARG);												// vec_madd(PhaseRotate[2].vec, BinShiftFract.vec, Zero);
					
					FltPointers[0][vectsize - 1] += PhaseRotate[1].flt[2];
					FltPointers[1][vectsize - 1] += PhaseRotate[1].flt[3];
					FltPointers[6][vectsize - 1] += PhaseRotate[2].flt[2];
					FltPointers[7][vectsize - 1] += PhaseRotate[2].flt[3];
					
					FltPointers[4][ToUp] = FltPointers[10][ToUp] = PhaseRotate[0].flt[0];
					FltPointers[5][ToUp] = FltPointers[11][ToUp] = PhaseRotate[0].flt[1];
				}
				
				if (DoTopDown)
				{
					PhaseRotate[1].flt[0] = PhaseRotate[1].flt[1] = FltPointers[2][ToUp + 1];
					PhaseRotate[1].flt[2] = PhaseRotate[1].flt[3] = FltPointers[3][ToUp + 1];
					PhaseRotate[1].vec = VEC_MUL_OP(PhaseRotate[1].vec, PhaseRotate[0].vec ZEROARG);											// vec_madd(PhaseRotate[1].vec, PhaseRotate[0].vec, Zero);
					PhaseRotate[1].vec = VEC_ADD_OP(PhaseRotate[1].vec, VEC_ROT_OP(PhaseRotate[1].vec, PhaseRotate[1].vec, VEC_ROT_VAL));		// vec_add(PhaseRotate[1].vec, vec_sld(PhaseRotate[1].vec, PhaseRotate[1].vec, 8));
					PhaseRotate[1].vec = VEC_MUL_OP(PhaseRotate[1].vec, BinShiftFract.vec ZEROARG);												// vec_madd(PhaseRotate[1].vec, BinShiftFract.vec, Zero);
					
					PhaseRotate[2].flt[0] = PhaseRotate[2].flt[1] = FltPointers[8][i];
					PhaseRotate[2].flt[2] = PhaseRotate[2].flt[3] = FltPointers[9][i];
					PhaseRotate[2].vec = VEC_MUL_OP(PhaseRotate[2].vec, PhaseRotate[0].vec ZEROARG);											// vec_madd(PhaseRotate[2].vec, PhaseRotate[0].vec, Zero);
					PhaseRotate[2].vec = VEC_ADD_OP(PhaseRotate[2].vec, VEC_ROT_OP(PhaseRotate[2].vec, PhaseRotate[2].vec, VEC_ROT_VAL));		// vec_add(PhaseRotate[2].vec, vec_sld(PhaseRotate[2].vec, PhaseRotate[2].vec, 8));
					PhaseRotate[2].vec = VEC_MUL_OP(PhaseRotate[2].vec, BinShiftFract.vec ZEROARG);												// vec_madd(PhaseRotate[2].vec, BinShiftFract.vec, Zero);
					
					// N.B. Conjugation + Reversal of Interpolation!
					
					FltPointers[0][vectsize - 1] += PhaseRotate[1].flt[0];
					FltPointers[1][vectsize - 1] -= PhaseRotate[1].flt[1];
					FltPointers[6][vectsize - 1] += PhaseRotate[2].flt[0];
					FltPointers[7][vectsize - 1] -= PhaseRotate[2].flt[1];
					
					FltPointers[4][ToUp + 1] = FltPointers[10][ToUp + 1] = PhaseRotate[0].flt[0];
					FltPointers[5][ToUp + 1] = FltPointers[11][ToUp + 1] = PhaseRotate[0].flt[1];
				}
			}
						
			StartOfRegion = NextStartOfRegion;
			PeakAmp = NextPeakAmp;
			PeakLoc = NextPeakLoc;
		}
	}
	else
	{
	
		///////////////////////////////////////// Channels Separately////////////////////////////////////		
		/////////////////////////////////////////// Left Channel ////////////////////////////////////////
		
		// Load in amplitudes
		 
		for (i = 0; i < vectsize >> 2; i++)
			VecAmpsSq[i + 1] = VEC_ADD_OP(VEC_MUL_OP(VecPointers[2][i],VecPointers[2][i] ZEROARG), VEC_MUL_OP(VecPointers[3][i], VecPointers[3][i] ZEROARG)); // vec_madd(VecPointers[2][i], VecPointers[2][i], vec_madd (VecPointers[3][i], VecPointers[3][i], Zero));
		
		// Search for first peak
		
		for (i = 0; i < vectsize; i++)							
		{
			if (AmpsSq[i + 4] > AmpsSq[i + 3] && AmpsSq[i + 4] > AmpsSq[i + 2] && AmpsSq[i + 4] > AmpsSq[i + 5] && AmpsSq[i + 4] > AmpsSq[i + 6])
			{
				PeakAmp = AmpsSq[i + 4];
				PeakLoc = i;
				break;
			}
		}
		
		// Just in case
		
		if (i == vectsize)										
			PeakLoc = vectsize - 1;
		
		StartOfRegion = 0;
		while (StartOfRegion < vectsize)
		{
			MinAmp = PeakAmp;
			for (i = PeakLoc + 1; i < vectsize; i++)	// Search for next region
			{
				if (AmpsSq[i + 4] < MinAmp)				// Find minimums
				{
					MinAmp = AmpsSq[i + 4];
					NextStartOfRegion = i;
				}
				else									// Find peaks
				{
					if (AmpsSq[i + 4] > AmpsSq[i + 3] && AmpsSq[i + 4] > AmpsSq[i + 2] && AmpsSq[i + 4] > AmpsSq[i + 5] && AmpsSq[i + 4] > AmpsSq[i + 6])
					{
						NextPeakAmp = AmpsSq[i + 4];
						NextPeakLoc = i;
						break;
					}
				}
			}
			
			if (i >= vectsize)
				NextStartOfRegion = vectsize;
			
			// More accurate peak location (parabolic interpolation)
			
			if (AmpsSq[PeakLoc + 3]) LogEstimate.flt[0] = (AmpsSq[PeakLoc + 3]); else LogEstimate.flt[0] = 1;
			if (AmpsSq[PeakLoc + 4]) LogEstimate.flt[1] = (AmpsSq[PeakLoc + 4]); else LogEstimate.flt[1] = 1;
			if (AmpsSq[PeakLoc + 5]) LogEstimate.flt[2] = (AmpsSq[PeakLoc + 5]); else LogEstimate.flt[2] = 1;
			LogEstimate.vec = vlogf(LogEstimate.vec);
			
			PeakEstDivisor = LogEstimate.flt[0] + LogEstimate.flt[2] - (2.0 * LogEstimate.flt[1]);
			if (PeakEstDivisor == 0) PeakCorrection = 0;
			else PeakCorrection = (0.5 * (LogEstimate.flt[0] - LogEstimate.flt[2])) / PeakEstDivisor;
		
			AcPeakLoc = (float) PeakLoc + PeakCorrection;
			
			// Calculate shift 
		
			NormPeakLoc = (double) (AcPeakLoc * VectSizeRecip);
			ActShift = ((double) AcPeakLoc * (Params[4] + 10000 * (NormPeakLoc * (Params[3] + (NormPeakLoc * (Params[2] + (NormPeakLoc * (Params[1] + (NormPeakLoc * Params[0]))))))))) + NormShiftParam - AcPeakLoc;			
			
			if (ActShift == 0)
			{
				// Left Channel Only
				for (i = StartOfRegion; i < NextStartOfRegion; i++)
				{
					FltPointers[0][i] += FltPointers[2][i];
					FltPointers[1][i] += FltPointers[3][i];
					FltPointers[4][i] = 1;
					FltPointers[5][i] = 0;
				}
			}
			else
			{
				if (Reflect)
				{
					// Mirror Output Around DC and Nyquist
					
					Wrap = (double) StartOfRegion + ActShift;
					Conj = 0;

					if (Wrap < 0)
					{
						Wrap = -Wrap;
						Conj = 1;
					}
					if (Wrap > (double) TwoVect)
					{
						WrapDivisor = Wrap * TwoVectRecip;
						Conj = ((int) floor (WrapDivisor) + Conj) % 2;
						Wrap -= (floor(WrapDivisor) * (double) TwoVect);
					}
					if (Wrap > (double) vectsize)
					{
						Wrap = (double) TwoVect - Wrap;
						Conj = 1 - Conj;
					}
					
					NewBin = (int) floor(Wrap);
					
					// Work Out Which Bits Go Which Way
					
					if (Conj)
					{
						FromDown = StartOfRegion;
						DoTopUp = DoTopDown = FromUp = ToUp = NewBin2 = 0;
						ToDown = StartOfRegion + NewBin + 1;
						if (ToDown < NextStartOfRegion)
						{
							FromUp = ToDown;
							ToUp = NextStartOfRegion;
						}
						else
							ToDown = NextStartOfRegion;
					}
					else
					{
						FromUp = StartOfRegion;
						DoTopUp = DoTopDown = FromDown = ToDown = 0;
						NewBin2 = NewBin;
						NewBin = vectsize - 2;
						ToUp = StartOfRegion + vectsize - 1 - NewBin2;
						if (NextStartOfRegion > ToUp)
						{
							DoTopUp = 1;
							if (NextStartOfRegion > ToUp + 1)
							{
								DoTopDown = 1;
								FromDown = ToUp + 2;
								ToDown = NextStartOfRegion;
							}
						}
						else
							ToUp = NextStartOfRegion;
					}
				}
				else
				{
					DoTopUp = DoTopDown = FromDown = ToDown = 0;
					FromUp = StartOfRegion;
					ToUp = NextStartOfRegion;
					
					BinShiftInt = (int) floor (ActShift);
					if (BinShiftInt < 0)
					{
						if (StartOfRegion + BinShiftInt < 0)
							FromUp = - BinShiftInt;
						if (BinShiftInt < -vectsize)
							FromUp = vectsize;
					}
					else
					{
						if (NextStartOfRegion + BinShiftInt > vectsize - 1)
						{
							ToUp = vectsize - BinShiftInt - 1;
							if (NextStartOfRegion + BinShiftInt > vectsize)
								DoTopUp = 1;
							if (ToUp < 0)
								ToUp = 0;
						}
					}
					NewBin2 = FromUp + BinShiftInt;
				}
				
				// Calculate shift vals and phase rotation
				
				BinShiftFract.flt[0] = BinShiftFract.flt[1] = (float) (ActShift - (double) floor(ActShift));
				BinShiftFract.flt[2] = BinShiftFract.flt[3] = 1 - BinShiftFract.flt[0];
				Temp[0] = cos((double) ActShift * ConstMultVal);
				Temp[1] = sin((double) ActShift * ConstMultVal);
				PhaseRotate[0].flt[0] = PhaseRotate[0].flt[3] = (Temp[0] * FltPointers[4][PeakLoc]) - (Temp[1] * FltPointers[5][PeakLoc]);
				PhaseRotate[0].flt[1] = (Temp[0] * FltPointers[5][PeakLoc]) + (Temp[1] * FltPointers[4][PeakLoc]);
				PhaseRotate[0].flt[2] = - PhaseRotate[0].flt[1];
				
				// Do the Shifts - Left Channel
				
				// Down
				
				for (i = FromDown; i < ToDown; i++)
				{
					PhaseRotate[1].flt[0] = PhaseRotate[1].flt[1] = FltPointers[2][i];
					PhaseRotate[1].flt[2] = PhaseRotate[1].flt[3] = FltPointers[3][i];
					PhaseRotate[1].vec = VEC_MUL_OP(PhaseRotate[1].vec, PhaseRotate[0].vec ZEROARG);											// vec_madd(PhaseRotate[1].vec, PhaseRotate[0].vec, Zero);
					PhaseRotate[1].vec = VEC_ADD_OP(PhaseRotate[1].vec, VEC_ROT_OP(PhaseRotate[1].vec, PhaseRotate[1].vec, VEC_ROT_VAL));		// vec_add(PhaseRotate[1].vec, vec_sld(PhaseRotate[1].vec, PhaseRotate[1].vec, 8));
					PhaseRotate[1].vec = VEC_MUL_OP(PhaseRotate[1].vec, BinShiftFract.vec ZEROARG);												// vec_madd(PhaseRotate[1].vec, BinShiftFract.vec, Zero);
					
					// N.B. Conjugation + Reversal of Interpolation!
					
					FltPointers[0][NewBin + 1] += PhaseRotate[1].flt[2];
					FltPointers[1][NewBin + 1] -= PhaseRotate[1].flt[3];
					FltPointers[0][NewBin] += PhaseRotate[1].flt[0];
					FltPointers[1][NewBin] -= PhaseRotate[1].flt[1];
					
					FltPointers[4][i] = PhaseRotate[0].flt[0];
					FltPointers[5][i] = PhaseRotate[0].flt[1];
					
					NewBin--; 
				}
				
				// Up
				
				for (i = FromUp; i < ToUp; i++)
				{
					PhaseRotate[1].flt[0] = PhaseRotate[1].flt[1] = FltPointers[2][i];
					PhaseRotate[1].flt[2] = PhaseRotate[1].flt[3] = FltPointers[3][i];
					PhaseRotate[1].vec = VEC_MUL_OP(PhaseRotate[1].vec, PhaseRotate[0].vec ZEROARG);											// vec_madd(PhaseRotate[1].vec, PhaseRotate[0].vec, Zero);
					PhaseRotate[1].vec = VEC_ADD_OP(PhaseRotate[1].vec, VEC_ROT_OP(PhaseRotate[1].vec, PhaseRotate[1].vec, VEC_ROT_VAL));		// vec_add(PhaseRotate[1].vec, vec_sld(PhaseRotate[1].vec, PhaseRotate[1].vec, 8));
					PhaseRotate[1].vec = VEC_MUL_OP(PhaseRotate[1].vec, BinShiftFract.vec ZEROARG);												// vec_madd(PhaseRotate[1].vec, BinShiftFract.vec, Zero);
					
					FltPointers[0][NewBin2] += PhaseRotate[1].flt[2];
					FltPointers[1][NewBin2] += PhaseRotate[1].flt[3];
					FltPointers[0][NewBin2 + 1] += PhaseRotate[1].flt[0];
					FltPointers[1][NewBin2 + 1] += PhaseRotate[1].flt[1];
					
					FltPointers[4][i] = PhaseRotate[0].flt[0];
					FltPointers[5][i] = PhaseRotate[0].flt[1];
					
					NewBin2++; 
				}
				
				if (DoTopUp)
				{
					PhaseRotate[1].flt[0] = PhaseRotate[1].flt[1] = FltPointers[2][ToUp];
					PhaseRotate[1].flt[2] = PhaseRotate[1].flt[3] = FltPointers[3][ToUp];
					PhaseRotate[1].vec = VEC_MUL_OP(PhaseRotate[1].vec, PhaseRotate[0].vec ZEROARG);											// vec_madd(PhaseRotate[1].vec, PhaseRotate[0].vec, Zero);
					PhaseRotate[1].vec = VEC_ADD_OP(PhaseRotate[1].vec, VEC_ROT_OP(PhaseRotate[1].vec, PhaseRotate[1].vec, VEC_ROT_VAL));		// vec_add(PhaseRotate[1].vec, vec_sld(PhaseRotate[1].vec, PhaseRotate[1].vec, 8));
					PhaseRotate[1].vec = VEC_MUL_OP(PhaseRotate[1].vec, BinShiftFract.vec ZEROARG);												// vec_madd(PhaseRotate[1].vec, BinShiftFract.vec, Zero);
					
					FltPointers[0][vectsize - 1] += PhaseRotate[1].flt[2];
					FltPointers[1][vectsize - 1] += PhaseRotate[1].flt[3];
					
					FltPointers[4][ToUp] = PhaseRotate[0].flt[0];
					FltPointers[5][ToUp] = PhaseRotate[0].flt[1];
				}
				
				if (DoTopDown)
				{
					PhaseRotate[1].flt[0] = PhaseRotate[1].flt[1] = FltPointers[2][ToUp + 1];
					PhaseRotate[1].flt[2] = PhaseRotate[1].flt[3] = FltPointers[3][ToUp + 1];
					PhaseRotate[1].vec = VEC_MUL_OP(PhaseRotate[1].vec, PhaseRotate[0].vec ZEROARG);											// vec_madd(PhaseRotate[1].vec, PhaseRotate[0].vec, Zero);
					PhaseRotate[1].vec = VEC_ADD_OP(PhaseRotate[1].vec, VEC_ROT_OP(PhaseRotate[1].vec, PhaseRotate[1].vec, VEC_ROT_VAL));		// vec_add(PhaseRotate[1].vec, vec_sld(PhaseRotate[1].vec, PhaseRotate[1].vec, 8));
					PhaseRotate[1].vec = VEC_MUL_OP(PhaseRotate[1].vec, BinShiftFract.vec ZEROARG);												// vec_madd(PhaseRotate[1].vec, BinShiftFract.vec, Zero);

					
					// N.B. Conjugation + Reversal of Interpolation!
					
					FltPointers[0][vectsize - 1] += PhaseRotate[1].flt[0];
					FltPointers[1][vectsize - 1] -= PhaseRotate[1].flt[1];
					
					FltPointers[4][ToUp + 1] = PhaseRotate[0].flt[0];
					FltPointers[5][ToUp + 1] = PhaseRotate[0].flt[1];
				}
			}
			
			StartOfRegion = NextStartOfRegion;
			PeakAmp = NextPeakAmp;
			PeakLoc = NextPeakLoc;
		}
		
		///////////////////////////////////////// Channels Separately////////////////////////////////////
		/////////////////////////////////////////// Right Channel ///////////////////////////////////////
		
		// Load in amplitudes
		 
		for (i = 0; i < vectsize >> 2; i++)
			VecAmpsSq[i + 1] = VecAmpsSq[i + 1] = VEC_ADD_OP(VEC_MUL_OP(VecPointers[6][i],VecPointers[6][i] ZEROARG), VEC_MUL_OP(VecPointers[7][i], VecPointers[7][i] ZEROARG)); // vec_madd(VecPointers[6][i], VecPointers[6][i], vec_madd (VecPointers[7][i], VecPointers[7][i], Zero));
		
		// Search for first peak
		
		for (i = 0; i < vectsize; i++)							
		{
			if (AmpsSq[i + 4] > AmpsSq[i + 3] && AmpsSq[i + 4] > AmpsSq[i + 2] && AmpsSq[i + 4] > AmpsSq[i + 5] && AmpsSq[i + 4] > AmpsSq[i + 6])
			{
				PeakAmp = AmpsSq[i + 4];
				PeakLoc = i;
				break;
			}
		}
		
		// Just in case
		
		if (i == vectsize)										
			PeakLoc = vectsize - 1;
		
		StartOfRegion = 0;
		while (StartOfRegion < vectsize)
		{
			MinAmp = PeakAmp;
			for (i = PeakLoc + 1; i < vectsize; i++)	// Search for next region
			{
				if (AmpsSq[i + 4] < MinAmp)				// Find minimums
				{
					MinAmp = AmpsSq[i + 4];
					NextStartOfRegion = i;
				}
				else									// Find peaks
				{
					if (AmpsSq[i + 4] > AmpsSq[i + 3] && AmpsSq[i + 4] > AmpsSq[i + 2] && AmpsSq[i + 4] > AmpsSq[i + 5] && AmpsSq[i + 4] > AmpsSq[i + 6])
					{
						NextPeakAmp = AmpsSq[i + 4];
						NextPeakLoc = i;
						break;
					}
				}
			}
			
			if (i >= vectsize)
				NextStartOfRegion = vectsize;
			
			// More accurate peak location (parabolic interpolation)
			
			if (AmpsSq[PeakLoc + 3]) LogEstimate.flt[0] = (AmpsSq[PeakLoc + 3]); else LogEstimate.flt[0] = 1;
			if (AmpsSq[PeakLoc + 4]) LogEstimate.flt[1] = (AmpsSq[PeakLoc + 4]); else LogEstimate.flt[1] = 1;
			if (AmpsSq[PeakLoc + 5]) LogEstimate.flt[2] = (AmpsSq[PeakLoc + 5]); else LogEstimate.flt[2] = 1;
			LogEstimate.vec = vlogf(LogEstimate.vec);
			
			PeakEstDivisor = LogEstimate.flt[0] + LogEstimate.flt[2] - (2.0 * LogEstimate.flt[1]);
			if (PeakEstDivisor == 0) PeakCorrection = 0;
			else PeakCorrection = (0.5 * (LogEstimate.flt[0] - LogEstimate.flt[2])) / PeakEstDivisor;
		
			AcPeakLoc = (float) PeakLoc + PeakCorrection;
			
			// Calculate shift 
		
			NormPeakLoc = (double) (AcPeakLoc * VectSizeRecip);
			ActShift = ((double) AcPeakLoc * (Params[4] + 10000 * (NormPeakLoc * (Params[3] + (NormPeakLoc * (Params[2] + (NormPeakLoc * (Params[1] + (NormPeakLoc * Params[0]))))))))) + NormShiftParam - AcPeakLoc;			
			
			if (ActShift == 0)
			{
				// Right Channel Only
				for (i = StartOfRegion; i < NextStartOfRegion; i++)
				{
					FltPointers[6][i] += FltPointers[8][i];
					FltPointers[7][i] += FltPointers[9][i];
					FltPointers[10][i] = 1;
					FltPointers[11][i] = 0;
				}
			}
			else
			{
				if (Reflect)
				{
					// Mirror Output Around DC and Nyquist
					
					Wrap = (double) StartOfRegion + ActShift;
					Conj = 0;

					if (Wrap < 0)
					{
						Wrap = -Wrap;
						Conj = 1;
					}
					if (Wrap > (double) TwoVect)
					{
						WrapDivisor = Wrap * TwoVectRecip;
						Conj = ((int) floor (WrapDivisor) + Conj) % 2;
						Wrap -= (floor(WrapDivisor) * (double) TwoVect);
					}
					if (Wrap > (double) vectsize)
					{
						Wrap = (double) TwoVect - Wrap;
						Conj = 1 - Conj;
					}
					
					NewBin = (int) floor(Wrap);
					
					// Work Out Which Bits Go Which Way
					
					if (Conj)
					{
						FromDown = StartOfRegion;
						DoTopUp = DoTopDown = FromUp = ToUp = NewBin2 = 0;
						ToDown = StartOfRegion + NewBin + 1;
						if (ToDown < NextStartOfRegion)
						{
							FromUp = ToDown;
							ToUp = NextStartOfRegion;
						}
						else
							ToDown = NextStartOfRegion;
					}
					else
					{
						FromUp = StartOfRegion;
						DoTopUp = DoTopDown = FromDown = ToDown = 0;
						NewBin2 = NewBin;
						NewBin = vectsize - 2;
						ToUp = StartOfRegion + vectsize - 1 - NewBin2;
						if (NextStartOfRegion > ToUp)
						{
							DoTopUp = 1;
							if (NextStartOfRegion > ToUp + 1)
							{
								DoTopDown = 1;
								FromDown = ToUp + 2;
								ToDown = NextStartOfRegion;
							}
						}
						else
							ToUp = NextStartOfRegion;
					}
				}
				else
				{
					DoTopUp = DoTopDown = FromDown = ToDown = 0;
					FromUp = StartOfRegion;
					ToUp = NextStartOfRegion;
					
					BinShiftInt = (int) floor (ActShift);
					if (BinShiftInt < 0)
					{
						if (StartOfRegion + BinShiftInt < 0)
							FromUp = - BinShiftInt;
						if (BinShiftInt < -vectsize)
							FromUp = vectsize;
					}
					else
					{
						if (NextStartOfRegion + BinShiftInt > vectsize - 1)
						{
							ToUp = vectsize - BinShiftInt - 1;
							if (NextStartOfRegion + BinShiftInt > vectsize)
								DoTopUp = 1;
							if (ToUp < 0)
								ToUp = 0;
						}
					}
					NewBin2 = FromUp + BinShiftInt;
				}
				
				// Calculate shift vals and phase rotation
				
				BinShiftFract.flt[0] = BinShiftFract.flt[1] = (float) (ActShift - (double) floor(ActShift));
				BinShiftFract.flt[2] = BinShiftFract.flt[3] = 1 - BinShiftFract.flt[0];
				Temp[0] = cos((double) ActShift * ConstMultVal);
				Temp[1] = sin((double) ActShift * ConstMultVal);
				PhaseRotate[0].flt[0] = PhaseRotate[0].flt[3] = (Temp[0] * FltPointers[4][PeakLoc]) - (Temp[1] * FltPointers[5][PeakLoc]);
				PhaseRotate[0].flt[1] = (Temp[0] * FltPointers[5][PeakLoc]) + (Temp[1] * FltPointers[4][PeakLoc]);
				PhaseRotate[0].flt[2] = - PhaseRotate[0].flt[1];
				
				// Do the Shifts - Right Channel
				
				// Down
				
				for (i = FromDown; i < ToDown; i++)
				{
					PhaseRotate[1].flt[0] = PhaseRotate[1].flt[1] = FltPointers[8][i];
					PhaseRotate[1].flt[2] = PhaseRotate[1].flt[3] = FltPointers[9][i];
					PhaseRotate[1].vec = VEC_MUL_OP(PhaseRotate[1].vec, PhaseRotate[0].vec ZEROARG);											// vec_madd(PhaseRotate[1].vec, PhaseRotate[0].vec, Zero);
					PhaseRotate[1].vec = VEC_ADD_OP(PhaseRotate[1].vec, VEC_ROT_OP(PhaseRotate[1].vec, PhaseRotate[1].vec, VEC_ROT_VAL));		// vec_add(PhaseRotate[1].vec, vec_sld(PhaseRotate[1].vec, PhaseRotate[1].vec, 8));
					PhaseRotate[1].vec = VEC_MUL_OP(PhaseRotate[1].vec, BinShiftFract.vec ZEROARG);												// vec_madd(PhaseRotate[1].vec, BinShiftFract.vec, Zero);
					
					// N.B. Conjugation + Reversal of Interpolation!
					
					FltPointers[6][NewBin + 1] += PhaseRotate[1].flt[2];
					FltPointers[7][NewBin + 1] -= PhaseRotate[1].flt[3];
					FltPointers[6][NewBin] += PhaseRotate[1].flt[0];
					FltPointers[7][NewBin] -= PhaseRotate[1].flt[1];
					
					FltPointers[10][i] = PhaseRotate[0].flt[0];
					FltPointers[11][i] = PhaseRotate[0].flt[1];
					
					NewBin--; 
				}
				
				// Up
				
				for (i = FromUp; i < ToUp; i++)
				{
					PhaseRotate[1].flt[0] = PhaseRotate[1].flt[1] = FltPointers[8][i];
					PhaseRotate[1].flt[2] = PhaseRotate[1].flt[3] = FltPointers[9][i];
					PhaseRotate[1].vec = VEC_MUL_OP(PhaseRotate[1].vec, PhaseRotate[0].vec ZEROARG);											// vec_madd(PhaseRotate[1].vec, PhaseRotate[0].vec, Zero);
					PhaseRotate[1].vec = VEC_ADD_OP(PhaseRotate[1].vec, VEC_ROT_OP(PhaseRotate[1].vec, PhaseRotate[1].vec, VEC_ROT_VAL));		// vec_add(PhaseRotate[1].vec, vec_sld(PhaseRotate[1].vec, PhaseRotate[1].vec, 8));
					PhaseRotate[1].vec = VEC_MUL_OP(PhaseRotate[1].vec, BinShiftFract.vec ZEROARG);												// vec_madd(PhaseRotate[1].vec, BinShiftFract.vec, Zero);
					
					FltPointers[6][NewBin2] += PhaseRotate[1].flt[2];
					FltPointers[7][NewBin2] += PhaseRotate[1].flt[3];
					FltPointers[6][NewBin2 + 1] += PhaseRotate[1].flt[0];
					FltPointers[7][NewBin2 + 1] += PhaseRotate[1].flt[1];
					
					FltPointers[10][i] = PhaseRotate[0].flt[0];
					FltPointers[11][i] = PhaseRotate[0].flt[1];
					
					NewBin2++; 
				}
				
				if (DoTopUp)
				{
					PhaseRotate[1].flt[0] = PhaseRotate[1].flt[1] = FltPointers[8][ToUp];
					PhaseRotate[1].flt[2] = PhaseRotate[1].flt[3] = FltPointers[9][ToUp];
					PhaseRotate[1].vec = VEC_MUL_OP(PhaseRotate[1].vec, PhaseRotate[0].vec ZEROARG);											// vec_madd(PhaseRotate[1].vec, PhaseRotate[0].vec, Zero);
					PhaseRotate[1].vec = VEC_ADD_OP(PhaseRotate[1].vec, VEC_ROT_OP(PhaseRotate[1].vec, PhaseRotate[1].vec, VEC_ROT_VAL));		// vec_add(PhaseRotate[1].vec, vec_sld(PhaseRotate[1].vec, PhaseRotate[1].vec, 8));
					PhaseRotate[1].vec = VEC_MUL_OP(PhaseRotate[1].vec, BinShiftFract.vec ZEROARG);												// vec_madd(PhaseRotate[1].vec, BinShiftFract.vec, Zero);
					
					FltPointers[6][vectsize - 1] += PhaseRotate[1].flt[2];
					FltPointers[7][vectsize - 1] += PhaseRotate[1].flt[3];
					
					FltPointers[10][ToUp] = PhaseRotate[0].flt[0];
					FltPointers[11][ToUp] = PhaseRotate[0].flt[1];
				}
				
				if (DoTopDown)
				{
					PhaseRotate[1].flt[0] = PhaseRotate[1].flt[1] = FltPointers[8][ToUp + 1];
					PhaseRotate[1].flt[2] = PhaseRotate[1].flt[3] = FltPointers[9][ToUp + 1];
					PhaseRotate[1].vec = VEC_MUL_OP(PhaseRotate[1].vec, PhaseRotate[0].vec ZEROARG);											// vec_madd(PhaseRotate[1].vec, PhaseRotate[0].vec, Zero);
					PhaseRotate[1].vec = VEC_ADD_OP(PhaseRotate[1].vec, VEC_ROT_OP(PhaseRotate[1].vec, PhaseRotate[1].vec, VEC_ROT_VAL));		// vec_add(PhaseRotate[1].vec, vec_sld(PhaseRotate[1].vec, PhaseRotate[1].vec, 8));
					PhaseRotate[1].vec = VEC_MUL_OP(PhaseRotate[1].vec, BinShiftFract.vec ZEROARG);												// vec_madd(PhaseRotate[1].vec, BinShiftFract.vec, Zero);
					
					// N.B. Conjugation + Reversal of Interpolation!
					
					FltPointers[6][vectsize - 1] += PhaseRotate[1].flt[0];
					FltPointers[7][vectsize - 1] -= PhaseRotate[1].flt[1];
					
					FltPointers[10][ToUp + 1] = PhaseRotate[0].flt[0];
					FltPointers[11][ToUp + 1] = PhaseRotate[0].flt[1];
				}
			}
			StartOfRegion = NextStartOfRegion;
			PeakAmp = NextPeakAmp;
			PeakLoc = NextPeakLoc;
		}
	}
}




void FreqWarp_assist(t_FreqWarp *x, void *b, long m, long a, char *s)
{
	if (m == ASSIST_OUTLET) 
	{
		switch (a) 
		{	
		case 0:
			sprintf(s,"(signal) Output Left");
			break;
		case 1:
			sprintf(s,"(signal) Output Right");
			break;

		}
	}
	else 
	{
		switch (a) 
		{	
		case 0:
			sprintf(s,"(signal) Input Left");
			break;
		case 1:
			sprintf(s,"(signal) Input Right");
			break;
		case 2:
			sprintf(s,"FFT Size");
			break;
		case 3:
			sprintf(s,"FFT Overlap");
			break;
		case 4:
			sprintf(s,"FFT Window");
			break;
		case 5:
			sprintf(s,"Channel Link");
			break;
		case 6:
			sprintf(s,"Parameter A In");
			break;
		case 7:
			sprintf(s,"Parameter B In");
			break;
		case 8:
			sprintf(s,"Parameter C In");
			break;
		case 9:
			sprintf(s,"Parameter D In");
			break;
		case 10:
			sprintf(s,"Parameter E In");
			break;
		case 11:
			sprintf(s,"Parameter F In");
			break;
		}
	}
}




