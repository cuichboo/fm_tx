#include <vector>

std::vector<double> optfir_lowpass(
        double gain, double fs, \
        double freq1, double freq2, \
        double passband_ripple_db, \
        double stopband_atten_db, \
        int nextra_taps=2);
