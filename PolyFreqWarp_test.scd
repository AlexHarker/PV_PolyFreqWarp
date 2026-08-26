(
var soundFile = Platform.resourceDir +/+ "sounds/a11wlk01.wav";
var buffer;

s.waitForBoot {
    buffer = Buffer.read(s, soundFile, action: {
        {
            SynthDef(\polyFreqWarpTest, {
            var source, leftFFT, rightFFT, warped;
            var curve0, curve1, curve2, curve3, scale, shift;

            source = PlayBuf.ar(1, \buffer.kr(0), BufRateScale.kr(\buffer.kr(0)), loop: 1);
            source = source ! 2;

            leftFFT = FFT(LocalBuf(4096), source[0], hop: 0.25, wintype: 0);
            rightFFT = FFT(LocalBuf(4096), source[1], hop: 0.25, wintype: 0);

            curve0 = Rand(-0.0000002, 0.0000002);
            curve1 = Rand(-0.000002, 0.000002);
            curve2 = Rand(-0.00002, 0.00002);
            curve3 = Rand(-0.0002, 0.0002);
            scale = Rand(0.998, 1.002);
            shift = Rand(-0.00002, 0.00002);

            // link = 0 and reflect = 1 stay fixed for this test.
            warped = PV_PolyFreqWarp(
                leftFFT,
                rightFFT,
                curve0,
                curve1,
                curve2,
                curve3,
                scale,
                shift,
                0,
                1
            );

            Out.ar(0, [IFFT(warped[0]), IFFT(warped[1])] * 0.2);
            }).add;

            s.sync;
            Synth(\polyFreqWarpTest, [\buffer, buffer]);
        }.fork;
    });

    CmdPeriod.doOnce({ buffer.free });
};
)
