PV_PolyFreqWarp : PV_ChainUGen {
    *new { |buffer, curve0 = 0.0, curve1 = 0.0, curve2 = 0.0, curve3 = 0.0,
        scale = 1.0, shift = 0.0, reflect = 1, detector, overlap = 4|
        ^this.multiNew('control', buffer, curve0, curve1, curve2, curve3,
            scale, shift, reflect, detector ? -1, overlap)
    }
}
