(
//var soundFile = Platform.resourceDir +/+ "sounds/a11wlk01.wav";
var soundFile = "/Users/alexharker/Desktop/- EVERYTHING/- Audio Test Files/05 She Treats Me Well.wav";
//var soundFile = "/Users/alexharker/Desktop/- EVERYTHING/- Audio Test Files/01 Old Pine.aif";
var buffer;

s.waitForBoot {
    buffer = Buffer.read(s, soundFile, action: {
        {
            SynthDef(\polyFreqWarpTest, {
            var source, leftFFT, rightFFT, warped, trigger;
            var curve0, curve1, curve2, curve3, scale, shift;

            source = PlayBuf.ar(2, \buffer.kr(0), BufRateScale.kr(\buffer.kr(0)), loop: 1);
            //source = source ! 2;

            leftFFT = FFT(LocalBuf(4096), source[0], hop: 0.25, wintype: 1);
            rightFFT = FFT(LocalBuf(4096), source[1], hop: 0.25, wintype: 1);

            trigger = Impulse.kr(2);
            curve0 = 0.0;//TRand.kr(-0.0000002, 0.0000002, trigger);
            curve1 = 0.0;//TRand.kr(-0.000002, 0.000002, trigger);
            curve2 = 0.0001;//TRand.kr(-0.00002, 0.00002, trigger);
            curve3 = 0.0001;//TRand.kr(-0.0002, 0.0002, trigger);
            scale = 1;//TRand.kr(0.998, 1.002, trigger);
            shift = 0.001;

            // link = 0 and reflect = 1 stay fixed for this test.
            warped = PV_PolyFreqWarp(
					[leftFFT, rightFFT],
                curve0,
                curve1,
                curve2,
                curve3,
                scale,
                shift,
                1,
				PV_Add(leftFFT, rightFFT)
            );

            Out.ar(0, IFFT(warped));
			}).add;

            s.sync;
            Synth(\polyFreqWarpTest, [\buffer, buffer]);
        }.fork;
    });

    CmdPeriod.doOnce({ buffer.free });
};
)
