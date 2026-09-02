
// NB. By using Macros for Vec operations we can test the platform on compile and therefore compile the same code for altivec and sse wherever possible 

#if (TARGET_RT_LITTLE_ENDIAN)
#define VEC_MUL_OP _mm_mul_ps
#define VEC_ADD_OP _mm_add_ps
#define ZEROARG
#define ZEROINIT
#else
#define VEC_MUL_OP vec_madd
#define VEC_ADD_OP vec_add
#define ZEROARG , Zero
#define ZEROINIT vFloat Zero = {0.,0.,0.,0.};
#endif

//////////////////////////////////////////// Define FFT Limits /////////////////////////////////////////////////

#define MaxFFTSizelog2 14						

// NB. Max FFT defined by 2^MaxFFTSizeLog2 - 14 gives 16384

static const int MAXFFTSIZE = 2 << (MaxFFTSizelog2 - 1);
static const int MaxFFTOVER2 = (2 << (MaxFFTSizelog2 - 1)) >> 1;
static const int MaxFFTOVER4 = (2 << (MaxFFTSizelog2 - 1)) >> 2;

//////////////////////////////////////// Macro for Main FFT Stereo /////////////////////////////////////////////


#define FFTMain()										\
addmess((method)FFTExtSt_dsp, "dsp", A_CANT, 0);		\
addmess((method)FFTExt_DC, "DC", 	A_LONG, 0);


//////////////////////////////////// Macro for Initialising FFT Stereo /////////////////////////////////////////


#define FFTInitSt(FFTSizeLog2, OverlapLog2)											\
x->FFTBuffers[0] = (vFloat *) malloc ((16384 * 16 * sizeof(float)));			    \
if (x->FFTBuffers[0] == 0)															\
{																					\
	post ("Could Not Allocate Memory for FFT - Bailing Out....");					\
		x->FFTMemAlloc = 0;															\
}																					\
else																				\
x->FFTMemAlloc = 1;																	\
x->FFTBuffers[1] = x->FFTBuffers[0] + 8192;											\
x->FFTBuffers[2] = x->FFTBuffers[1] + 4096;											\
x->FFTBuffers[3] = x->FFTBuffers[2] + 4096;											\
x->FFTBuffers[4] = x->FFTBuffers[3] + 8192;											\
x->FFTBuffers[5] = x->FFTBuffers[4] + 4096;											\
x->FFTBuffers[6] = x->FFTBuffers[5] + 4096;											\
x->FFTBuffers[7] = x->FFTBuffers[6] + 4096;											\
x->FFTWindows[0] = x->FFTBuffers[7] + 4096;											\
x->FFTWindows[1] = x->FFTWindows[0] + 4096;											\
x->FFTProcess[0].realp = (float *) (x->FFTWindows[1] + 4096);						\
x->FFTProcess[0].imagp = x->FFTProcess[0].realp + 8192;								\
x->FFTProcess[1].realp = x->FFTProcess[0].imagp + 8192;								\
x->FFTProcess[1].imagp = x->FFTProcess[1].realp + 8192;								\
x->FFTProcess[2].realp = x->FFTProcess[1].imagp + 8192;								\
x->FFTProcess[2].imagp = x->FFTProcess[2].realp + 8192;								\
x->FFTProcess[3].realp = x->FFTProcess[2].imagp + 8192;								\
x->FFTProcess[3].imagp = x->FFTProcess[3].realp + 8192;								\
x->FFTWindowType = 0;																\
x->FFTParam[13] = 0;																\
FFTExtSt_ChangeFFTParams (x, FFTSizeLog2, OverlapLog2);								\
x->FFTSetupReal = create_fftsetup (14, 0)


////////////////////////////////////////// Macro for Freeing FFT //////////////////////////////////////////////


#define FFTFree()												\
if (x->FFTMemAlloc)												\
free (x->FFTBuffers[0]);										\
destroy_fftsetup(x->FFTSetupReal)


////////////////////////////////////////// Macro for New FFT Size /////////////////////////////////////////////


#define FFTNewSizeSt(NewVal)									\
if (NewVal >= 3 && NewVal <= 14 && NewVal != x->FFTParam[2])	\
FFTExtSt_ChangeFFTParams (x, NewVal, x->FFTParam[4])


////////////////////////////////////////// Macro for New Overlap //////////////////////////////////////////////


#define FFTNewOverlapSt(NewVal)									\
if (NewVal != x->FFTParam[4])									\
FFTExtSt_ChangeFFTParams (x, x->FFTParam[2], NewVal)			\


////////////////////////////////////////// Macro for New Window //////////////////////////////////////////////


#define FFTNewWindow(NewVal)																\
if (NewVal >= 0 && NewVal <= 6 && NewVal != x->FFTWindowType)								\
{																							\
	x->FFTWindowType = NewVal;																\
		FFTExt_GenerateWindow(x, (float *) x->FFTWindows[0], x->FFTParam[0] * 4, NewVal);	\
}						


////////////////////////////////////////// FFT Stereo Functions //////////////////////////////////////////////


void FFTExtSt_ChangeFFTParams (T_EXTNAME *x, long FFTlog2, long OverlapLog2);
void FFTExt_GenerateWindow(T_EXTNAME *x, float *FFTWindow, int FFTSize, int FFTWindowSel);
void FFTExt_DC(T_EXTNAME *x, long DC);
void FFTExtSt_dsp(T_EXTNAME *x, t_signal **sp, short *count);
t_int *FFTExtSt_perform(t_int *w);


void FFTExtSt_ChangeFFTParams (T_EXTNAME *x, long FFTlog2, long OverlapLog2)
{
	int i, HopSize;
	
	vFloat *Buffer = x->FFTBuffers[0];
	vFloat *Buffer2 = x->FFTBuffers[3];
	vFloat Zero = (vFloat) {0.0,0.0,0.0,0.0};
	
	char FFTOffPrev = x->FFTOff;
	x->FFTOff = 1;
	// Initialise FFTPointers and fft data
	
	if (FFTlog2 == 3)
		OverlapLog2 = 1;
	if (FFTlog2 == 4 && OverlapLog2 == 3)
		OverlapLog2 = 2;
	
	if (OverlapLog2 < 1 || OverlapLog2 > 3)
		OverlapLog2 = 1;
	
	x->FFTParam[3] = 1 << OverlapLog2;
	x->FFTParam[4] = OverlapLog2;
	
	HopSize	= 1 << (FFTlog2 - 2 - OverlapLog2);
	
	x->FFTPointers[0] = HopSize;
	x->FFTPointers[1] = 1;
	x->FFTPointers[2] = 0;
	x->FFTPointers[3] = 1 << (FFTlog2 - 2);
	x->FFTPointers[4] = 0;
	
	x->FFTParam[0] = 1 << (FFTlog2 - 2);
	x->FFTParam[1] = 1 << (FFTlog2 - 1);
	x->FFTParam[2] = FFTlog2;
	
	for (i = 0; i < x->FFTParam[3]; i++)
		x->FFTParam[i + 5] = i * HopSize;
	
	FFTExt_GenerateWindow(x, (float *) x->FFTWindows[0], x->FFTParam[0] * 4, x->FFTWindowType);
	
	// Zero FFTBuffers
	
	if (x->FFTMemAlloc)
	{
		for (i = 0; i < 12288; i++)
		{
			Buffer[i] = Zero;
			Buffer2[i] = Zero;
		}
	}
	
	// Any Additional Reset Stuff for the specific object
	FFTRESET (x);
	
	x->FFTOff = FFTOffPrev;
}

#define FFTW_PI			 3.14159265358979323846
#define FFTW_TWOPI		 6.28318530717958647692
#define FFTW_FOURPI		12.56637061435817295384
#define FFTW_SIXPI		18.84955592153875943076

void FFTExt_GenerateWindow(T_EXTNAME *x, float *FFTWindow, int FFTSize, int FFTWindowSel)
{
	long i, j;
	float AlphaBesselRecip, NewTerm, XSq, BFunc, Temp;
	
	if (x->FFTMemAlloc)
	{
		switch (FFTWindowSel)
		{
			///////////////////////// von Hann Window //////////////////////////
			case 0:
				for (i = 0; i < FFTSize; i++)
					FFTWindow[i] = 0.5 - (0.5 * cos(FFTW_TWOPI * ((float) i / (float) FFTSize)));
				break;
				///////////////////////// Hamming Window //////////////////////////
			case 1:
				for (i = 0; i < FFTSize; i++)
					FFTWindow[i] = 0.54347826 - (0.45652174 * cos(FFTW_TWOPI * ((float) i / (float) FFTSize)));
				break;
				/////////////////// Kaiser Window (Alpha = 6.8) ///////////////////
			case 2:
				
				// First Find Bessel Function of Alpha
				
				XSq = 46.24;
				NewTerm = 0.25 * XSq;
				BFunc = 1.0;
				j = 2;
				
				while (NewTerm)	// Gives Maximum Accuracy
				{
					BFunc += NewTerm;
					AlphaBesselRecip = (1.0 / (4.0 * (float) j * (float) j));
					NewTerm = NewTerm * XSq * (1.0 / (4.0 * (float) j * (float) j));
					j++;
				}
					
					AlphaBesselRecip = 1 / BFunc;
				
				// Now Create Kaiser Window
				
				for (i = 0; i < FFTSize; i++)
				{
					Temp = ((2.0 * (float) i) - ((float) FFTSize - 1.0));
					Temp = Temp / FFTSize;
					Temp *= Temp;
					XSq = (1 - Temp) * 46.24;
					NewTerm = 0.25 * XSq;
					BFunc = 1;
					j = 2;
					
					while (NewTerm)	// Gives Maximum Accuracy
					{
						BFunc += NewTerm;
						NewTerm = NewTerm * XSq * (1.0 / (4.0 * (float) j * (float) j));
						j++;
					}
					FFTWindow[i] = BFunc * AlphaBesselRecip;
				}
					break;
				//////////////////////// Triangular Window ////////////////////////
			case 3:
				for (i = 0; i < (FFTSize >> 1); i++)
					FFTWindow[i] = (float) i / (float) (FFTSize >> 1);
				for (; i < FFTSize; i++)
					FFTWindow[i] = (float) (((float) FFTSize - 1) - (float) i) / (float) (FFTSize >> 1);
					break;
				//////////////////////// Blackman Window //////////////////////////
			case 4:
				for (i = 0; i < FFTSize; i++)
					FFTWindow[i] = 0.42659071 - (0.49656062 * cos(FFTW_TWOPI * ((float) i / (float) FFTSize))) + (0.07684867 * cos(FFTW_FOURPI * ((float) i / (float) FFTSize)));
				break;
				///////////////////// Blackman-Harris Window /////////////////////
			case 5:
				for (i = 0; i < FFTSize; i++)
					FFTWindow[i] = 0.35875 - (0.48829 * cos(FFTW_TWOPI * ((float) i / (float) FFTSize))) + (0.14128 * cos(FFTW_FOURPI * ((float) i / (float) FFTSize))) - (0.01168 * cos(FFTW_SIXPI * ((float) i / (float) FFTSize)));
				break;
				///////////////////////// Flat Top Window ////////////////////////
			case 6:
				for (i = 0; i < FFTSize; i++)
					FFTWindow[i] = 0.2810639 - (0.5208972 * cos(FFTW_TWOPI * ((float) i / (float) FFTSize))) + (0.1980399 * cos(FFTW_FOURPI * ((float) i / (float) FFTSize)));
				break;
		}
		
		// Sqrt Window if overlap is 2x (Unless Flat Top)
		
		if (x->FFTParam[4] == 1 && FFTWindowSel != 6)
		{
			for (i = 0; i < FFTSize; i++)
				FFTWindow[i] = sqrt(FFTWindow[i]);
		}
		
		// Calculate the gain of the window
		
		double GainOUT = 0;
		double FFTSizeRecip = (double) 1 / (double) FFTSize;
		double Overlap = (double) x->FFTParam[3];
		for (i = 0; i < FFTSize; i++)
			GainOUT += (double) FFTWindow[i] * (double) FFTWindow[i] * FFTSizeRecip;
		GainOUT *= Overlap * (double) FFTSize;
		
		// Do Scaling
		
		float *InWin = (float *) x->FFTWindows[0];
		float *OutWin = (float *) x->FFTWindows[1];
		
		float Scale = (float) 1 / ((float) GainOUT);
		float Scale2 = (float) 0.5;
		
		for (i = 0; i < FFTSize; i++)
		{
			OutWin[i] = InWin[i] * Scale;
			InWin[i] *= Scale2;
		}
	}
}


void FFTExt_DC(T_EXTNAME *x, long DC)
{
	if (DC)
		x->FFTParam[13] = 0;
	else
		x->FFTParam[13] = 1;
}


void FFTExtSt_dsp(T_EXTNAME *x, t_signal **sp, short *count)
{
	if (x->FFTMemAlloc)
		dsp_add(FFTExtSt_perform, 6, sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, sp[3]->s_vec, sp[0]->s_n, x);
}


t_int *FFTExtSt_perform(t_int *w)
{	
    int vectsize;
	T_EXTNAME *x;
    
    vFloat *in1 = (vFloat *)(w[1]);
	vFloat *in2 = (vFloat *)(w[2]);
    vFloat *out1 = (vFloat *)(w[3]);
	vFloat *out2 = (vFloat *)(w[4]);
    vectsize = w[5];
    x = (T_EXTNAME *)(w[6]);
	
    COMPLEX_SPLIT *FFTProcess = x->FFTProcess;
    FFTSetup FFTSetupReal = x->FFTSetupReal;
	
    vFloat **FFTBuffers = x->FFTBuffers;
    vFloat **FFTWindows = x->FFTWindows;
    int *FFTPointers = x->FFTPointers;
	int *FFTParam = x->FFTParam;
	
	ZEROINIT

	vFloat *BufPointer, *BufPointer2;
    int FFTTest, Loop, NextPoint, i;
	int HalfFFTOver4 = FFTParam[1] >> 2;
	int VectRemain = vectsize >> 2;
	
	if (x->FFTOff || x->x_obj.z_disabled)
		goto out;
	
	while (VectRemain)
	{
		// Check whether there will be another FFT this vector
		
		FFTTest = VectRemain - FFTPointers[0];
		if (FFTTest > 0)
		{
			Loop = FFTPointers[0];
			VectRemain -= FFTPointers[0];
			FFTPointers[0] = 0;
		}
		else
		{
			Loop = VectRemain;
			FFTPointers[0] -= VectRemain;
			VectRemain = 0;
		}
		
		// Load input into buffer (twice) and output from the output buffer
		
		// Left
		
		for (i = 0; i < Loop; i++)
		{
			// Left
			
			*(FFTBuffers[0] + FFTPointers[2]) = *in1; 
			*(FFTBuffers[0] + FFTPointers[3]) = *in1;
			in1++;
			
			// Right
			
			*(FFTBuffers[3] + FFTPointers[2]++) = *in2; 
			*(FFTBuffers[3] + FFTPointers[3]++) = *in2;
			in2++;
			
			// Left and Right 
			
			*out1++ = *(FFTBuffers[1] + FFTPointers[4]);	
			*out2++ = *(FFTBuffers[4] + FFTPointers[4]++);	
		}
		
		
		// Check that there is a new FFTs worth of buffer
		
		if (FFTPointers[0] == 0)
		{
			NextPoint = FFTPointers[1] + 5;
			
			// Window the data
			
			// Left
			
			BufPointer = FFTBuffers[0] + FFTParam[NextPoint];
			BufPointer2 = FFTBuffers[2] + HalfFFTOver4;
			for (i = 0; i < HalfFFTOver4; i++)
				*(BufPointer2 + i) = VEC_MUL_OP(*(BufPointer + i),  FFTWindows[0][i] ZEROARG); // *(BufPointer2 + i) = _mm_mul_ps(*(BufPointer + i),  FFTWindows[0][i]); // vec_madd(*(BufPointer + i),  FFTWindows[0][i], Zero);
			BufPointer2 = FFTBuffers[2] - HalfFFTOver4;
			for (i = HalfFFTOver4; i < FFTParam[0]; i++)
				*(BufPointer2 + i) = VEC_MUL_OP(*(BufPointer + i),  FFTWindows[0][i] ZEROARG); // *(BufPointer2 + i) = _mm_mul_ps(*(BufPointer + i),  FFTWindows[0][i]); // vec_madd(*(BufPointer + i),  FFTWindows[0][i], Zero);
			
			// Right
			
			BufPointer = FFTBuffers[3] + FFTParam[NextPoint];
			BufPointer2 = FFTBuffers[5] + HalfFFTOver4;
			for (i = 0; i < HalfFFTOver4; i++)
				*(BufPointer2 + i) = VEC_MUL_OP(*(BufPointer + i),  FFTWindows[0][i] ZEROARG); // *(BufPointer2 + i) = _mm_mul_ps(*(BufPointer + i),  FFTWindows[0][i]); // vec_madd(*(BufPointer + i),  FFTWindows[0][i], Zero);
			BufPointer2 = FFTBuffers[5] - HalfFFTOver4;
			for (i = HalfFFTOver4; i < FFTParam[0]; i++)
				*(BufPointer2 + i) = VEC_MUL_OP(*(BufPointer + i),  FFTWindows[0][i] ZEROARG); // *(BufPointer2 + i) = _mm_mul_ps(*(BufPointer + i),  FFTWindows[0][i]); // vec_madd(*(BufPointer + i),  FFTWindows[0][i], Zero);
			
			// Do the FFT, Process and Convert back with scaling (Left and Right)
			
			ctoz ((COMPLEX *) FFTBuffers[2], 2, &FFTProcess[0], 1, FFTParam[1]);
			fft_zrip (FFTSetupReal, &FFTProcess[0], 1, FFTParam[2], FFT_FORWARD);
			ctoz ((COMPLEX *) FFTBuffers[5], 2, &FFTProcess[1], 1, FFTParam[1]);
			fft_zrip (FFTSetupReal, &FFTProcess[1], 1, FFTParam[2], FFT_FORWARD);
			// Zero Nyquist + DC bins (In)
			
			FFTProcess[0].imagp[0] = (float) 0.;
			FFTProcess[1].imagp[0] = (float) 0.;
			if (FFTParam[13])
			{
				FFTProcess[0].realp[0] = (float) 0.;
				FFTProcess[1].realp[0] = (float) 0.;
			}
			
			// Processing Done Here
			PROCESSPERFORM (x, FFTProcess, FFTParam[1]);
			
			// Zero Nyquist + DC bins (Out)
			
			FFTProcess[2].imagp[0] = (float) 0.;
			FFTProcess[3].imagp[0] = (float) 0.;
			if (FFTParam[13])
			{
				FFTProcess[2].realp[0] = (float) 0.;
				FFTProcess[3].realp[0] = (float) 0.;
			}
			
			fft_zrip (FFTSetupReal, &FFTProcess[2], 1, FFTParam[2], FFT_INVERSE);
			ztoc (&FFTProcess[2], 1, (COMPLEX *) FFTBuffers[6], 2, FFTParam[1]);
			fft_zrip (FFTSetupReal, &FFTProcess[3], 1, FFTParam[2], FFT_INVERSE);
			ztoc (&FFTProcess[3], 1, (COMPLEX *) FFTBuffers[7], 2, FFTParam[1]);
			
			// Window the output
			
			// Left
			
			BufPointer2 = FFTBuffers[6] + HalfFFTOver4;
			for (i = 0; i < HalfFFTOver4; i++)
				*(FFTBuffers[2] + i) = VEC_MUL_OP(*(BufPointer2 + i),  FFTWindows[1][i] ZEROARG); // *(FFTBuffers[2] + i) = _mm_mul_ps(*(BufPointer2 + i),  FFTWindows[1][i]); // vec_madd(*(BufPointer2 + i),  FFTWindows[1][i], Zero);
			BufPointer2 = FFTBuffers[6] - HalfFFTOver4;
			for (i = HalfFFTOver4; i < FFTParam[0]; i++)
				*(FFTBuffers[2] + i) = VEC_MUL_OP(*(BufPointer2 + i),  FFTWindows[1][i] ZEROARG); // *(FFTBuffers[2] + i) = _mm_mul_ps(*(BufPointer2 + i),  FFTWindows[1][i]); // vec_madd(*(BufPointer2 + i),  FFTWindows[1][i], Zero);
			
			// Right
			
			BufPointer2 = FFTBuffers[7] + HalfFFTOver4;
			for (i = 0; i < HalfFFTOver4; i++)
				*(FFTBuffers[5] + i) = VEC_MUL_OP(*(BufPointer2 + i),  FFTWindows[1][i] ZEROARG); // *(FFTBuffers[5] + i) = _mm_mul_ps(*(BufPointer2 + i),  FFTWindows[1][i]); // vec_madd(*(BufPointer2 + i),  FFTWindows[1][i], Zero);
			BufPointer2 = FFTBuffers[7] - HalfFFTOver4;
			for (i = HalfFFTOver4; i < FFTParam[0]; i++)
				*(FFTBuffers[5] + i) = VEC_MUL_OP(*(BufPointer2 + i),  FFTWindows[1][i] ZEROARG); // *(FFTBuffers[5] + i) = _mm_mul_ps(*(BufPointer2 + i),  FFTWindows[1][i]); // vec_madd(*(BufPointer2 + i),  FFTWindows[1][i], Zero);	
			
			// Overlap Add the Result
			
			switch (FFTPointers[1])
			{
				case 0:
					Loop = FFTParam[0] - FFTParam[6];
					// Left
					for (i = 0; i < Loop; i++)
						*(FFTBuffers[1] + i) = VEC_ADD_OP(*(FFTBuffers[1] + i), *(FFTBuffers[2] + i)); // *(FFTBuffers[1] + i) = _mm_add_ps(*(FFTBuffers[1] + i), *(FFTBuffers[2] + i)); // vec_add(*(FFTBuffers[1] + i), *(FFTBuffers[2] + i));
					for (i = Loop; i < FFTParam[0]; i++)
						*(FFTBuffers[1] + i) = *(FFTBuffers[2] + i);
					//Right
					for (i = 0; i < Loop; i++)
						*(FFTBuffers[4] + i) = VEC_ADD_OP(*(FFTBuffers[4] + i), *(FFTBuffers[5] + i)); // *(FFTBuffers[4] + i) = _mm_add_ps(*(FFTBuffers[4] + i), *(FFTBuffers[5] + i));//vec_add(*(FFTBuffers[4] + i), *(FFTBuffers[5] + i));
					for (i = Loop; i < FFTParam[0]; i++)
						*(FFTBuffers[4] + i) = *(FFTBuffers[5] + i);
									
					// Reset I/O Buffer FFTPointers
								
					FFTPointers[2] = 0;
					FFTPointers[3] = FFTParam[0];
					FFTPointers[4] = 0;
					break;
					
				case 1:
					// Left
					BufPointer = FFTBuffers[2] - FFTParam[6];
					for (i = FFTParam[6]; i < FFTParam[0]; i++)
						*(FFTBuffers[1] + i) = VEC_ADD_OP(*(FFTBuffers[1] + i), *(BufPointer + i)); // *(FFTBuffers[1] + i) = _mm_add_ps(*(FFTBuffers[1] + i), *(BufPointer + i)); // vec_add(*(FFTBuffers[1] + i), *(BufPointer + i));
						BufPointer = FFTBuffers[2] + (FFTParam[0] - FFTParam[6]);
					for (i = 0; i < FFTParam[6]; i++)
						*(FFTBuffers[1] + i) = *(BufPointer + i);
						// Right
						BufPointer = FFTBuffers[5] - FFTParam[6];
					for (i = FFTParam[6]; i < FFTParam[0]; i++)
						*(FFTBuffers[4] + i) = VEC_ADD_OP(*(FFTBuffers[4] + i), *(BufPointer + i)); // *(FFTBuffers[4] + i) = _mm_add_ps(*(FFTBuffers[4] + i), *(BufPointer + i)); // vec_add(*(FFTBuffers[4] + i), *(BufPointer + i));
						BufPointer = FFTBuffers[5] + (FFTParam[0] - FFTParam[6]);
					for (i = 0; i < FFTParam[6]; i++)
						*(FFTBuffers[4] + i) = *(BufPointer + i);
						break;
					
				default:
					// Left
					Loop = FFTParam[NextPoint] - FFTParam[6];
					BufPointer = FFTBuffers[2] - FFTParam[NextPoint];
					for (i = FFTParam[NextPoint]; i < FFTParam[0]; i++)
						*(FFTBuffers[1] + i) =  VEC_ADD_OP(*(FFTBuffers[1] + i), *(BufPointer + i)); // *(FFTBuffers[1] + i) =  _mm_add_ps(*(FFTBuffers[1] + i), *(BufPointer + i)); // vec_add(*(FFTBuffers[1] + i), *(BufPointer + i));
						BufPointer = FFTBuffers[2] + (FFTParam[0] - FFTParam[NextPoint]);
					for (i = 0; i < Loop; i++)
						*(FFTBuffers[1] + i) =  VEC_ADD_OP(*(FFTBuffers[1] + i), *(BufPointer + i)); // *(FFTBuffers[1] + i) =  _mm_add_ps(*(FFTBuffers[1] + i), *(BufPointer + i)); // vec_add(*(FFTBuffers[1] + i), *(BufPointer + i));
						for (i = Loop; i < FFTParam[NextPoint]; i++)
							*(FFTBuffers[1] + i) = *(BufPointer + i);
							// Right
							Loop = FFTParam[NextPoint] - FFTParam[6];
					BufPointer = FFTBuffers[5] - FFTParam[NextPoint];
					for (i = FFTParam[NextPoint]; i < FFTParam[0]; i++)
						*(FFTBuffers[4] + i) = VEC_ADD_OP(*(FFTBuffers[4] + i), *(BufPointer + i)); // *(FFTBuffers[4] + i) = _mm_add_ps(*(FFTBuffers[4] + i), *(BufPointer + i)); // vec_add(*(FFTBuffers[4] + i), *(BufPointer + i));
						BufPointer = FFTBuffers[5] + (FFTParam[0] - FFTParam[NextPoint]);
					for (i = 0; i < Loop; i++)
						*(FFTBuffers[4] + i) = VEC_ADD_OP(*(FFTBuffers[4] + i), *(BufPointer + i)); // *(FFTBuffers[4] + i) = _mm_add_ps(*(FFTBuffers[4] + i), *(BufPointer + i)); // vec_add(*(FFTBuffers[4] + i), *(BufPointer + i));
						for (i = Loop; i < FFTParam[NextPoint]; i++)
							*(FFTBuffers[4] + i) = *(BufPointer + i);
							break;
			}
			
			// Set wait till next data to hop time and next offset
			
			
			FFTPointers[0] = FFTParam[6];	
			if (++FFTPointers[1] >= FFTParam[3])
				FFTPointers[1] = 0;
		}
	}
	
out:
		return w + 7;
}




