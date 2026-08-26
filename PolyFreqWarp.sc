PV_PolyFreqWarp : PV_ChainUGen {
    var <channels;

    *new { |bufferLeft, bufferRight, curve0 = 0.0, curve1 = 0.0, curve2 = 0.0,
        curve3 = 0.0, scale = 1.0, shift = 0.0, link = 0, reflect = 1|
        ^this.multiNew('control', bufferLeft, bufferRight, curve0, curve1, curve2,
            curve3, scale, shift, link, reflect)
    }

    init { |...theInputs|
        super.init(*theInputs);
        channels = Array.fill(2, { |index| OutputProxy(rate, this, index) });
        ^channels
    }
    
    numOutputs {
        ^2
    }

    writeOutputSpecs { |file|
        channels.do { |output| output.writeOutputSpec(file) };
    }
}
