#include <ext.h>
#include <z_dsp.h>
#include <Accelerate.h>

typedef union 
{
	float	flt[4];
	vFloat	vec;
} VecSplat;


/////////////////////////////////////// Macro for FFT Data Stereo ////////////////////////////////////////////


#define FFTDataSt												\
vFloat *FFTBuffers[8], *FFTWindows[2];					        \
int FFTPointers[5], FFTParam[14], FFTWindowType;				\
char FFTOff, FFTMemAlloc;										\
COMPLEX_SPLIT FFTProcess[4];									\
FFTSetup FFTSetupReal

/* N.B.

-- FFTBuffers

[0] In Buffer			--Left
[1] Out Buffer			--Left
[2] Processing Buffer	--Left
[3] In Buffer			--Right
[4] Out Buffer			--Right
[5] Processing Buffer	--Right
[6] Processing Buffer 2	--Left
[7] Processing Buffer 2	--Right

-- FFTPointers

[0]	Vectors till next FFT
[1]	Next Offset
[2]	Write (In) Pointer 1
[3]	Write (In) Pointer 2
[4]	Read (Out) Pointer

-- FFTParam

[0]	FFT Size in vectors
[1]	Half FFT Size
[2] FFT Size (log 2)
[3]	Overlap
[4]	Overlap (log 2)
[5] - [12] Offsets
[13] DCOff

Notes:
[5] is always 0
[6] is always the Hop Size

-- FFTProcess

[0] Process In Left
[1] Process In Right
[2] Process Out Left
[3] Process Out Right


*/



