#include <iostream>
#include <vector>
#include <math.h>
#include <algorithm>
#include <gnuradio/filter/pm_remez.h>
#include <optfir.h>

using namespace std;

double lporder (double freq1, double freq2, double delta_p, double delta_s);
void remezord(const std::vector<double> fcuts, \
         const std::vector<double> mags, \
         const std::vector<double> devs,
         int &n,
         std::vector<double> &fo,
         std::vector<double> &ao,
         std::vector<double> &w,
         int fsamp = 2);

double stopband_atten_to_dev (double atten_db)
{
    //Convert a stopband attenuation in dB to an absolute value
    return pow(10, (-atten_db/20));
}

double passband_ripple_to_dev (double ripple_db)
{
    //"""Convert passband ripple spec expressed in dB to an absolute value"""
    return (pow(10, (ripple_db/20))-1.0f) / (pow(10, (ripple_db/20))+1.0f);
}

std::vector<double> optfir_lowpass(double gain, double fs, \
                                   double freq1, double freq2, \
                                   double passband_ripple_db, \
                                   double stopband_atten_db, \
                                   int nextra_taps)
{
    int n = 0;
    std::vector<double> fo, ao, w;
    std::vector<double> desired_ampls(2), band(2), band_dev(2);
    
    double passband_dev = passband_ripple_to_dev(passband_ripple_db);
    double stopband_dev = stopband_atten_to_dev(stopband_atten_db);

    cout << "passband: " << passband_dev << " stopband: " << stopband_dev << endl;
    
    desired_ampls[0] = gain;
    desired_ampls[1] = 0;

    band[0] = freq1;
    band[1] = freq2;

    band_dev[0] = passband_dev;
    band_dev[1] = stopband_dev;

    /*
    (n, fo, ao, w) = remezord ([freq1, freq2], desired_ampls,
                               [passband_dev, stopband_dev], fs);
                               */

    remezord(band, desired_ampls, band_dev, n, fo, ao, w, fs);

    return gr::filter::pm_remez (n + nextra_taps, fo, ao, w, "bandpass");
}

void remezord(const std::vector<double> fcuts, \
         const std::vector<double> mags, \
         const std::vector<double> devs,
         int &n,
         std::vector<double> &fo,
         std::vector<double> &ao,
         std::vector<double> &w,
         int fsamp)
{
    std::vector<double> l_fcuts(fcuts), l_mags(mags), l_devs(devs);
    unsigned int nbands = 0;
    unsigned int i = 0;
    double l, max_dev, min_delta;
    std::vector<double> f1, f2;
    
    for (std::vector<double>::iterator it = l_fcuts.begin();
            it != l_fcuts.end(); it++) {

        *it = (float)(*it) / fsamp;
    }

    if (l_mags.size() != l_devs.size()) {
        std::cout << "Length of mags and devs must be equal" << std::endl;
    }

    nbands = l_mags.size();
    if (l_fcuts.size() != 2*(nbands - 1)) {
        std::cout << "Length of f must be 2 * len (mags) - 2" << std::endl;
    }

    for (i=0; i<nbands; i++) {
        if (l_mags[i] != 0)
            l_devs[i] = l_devs[i] / l_mags[i];
    }

    for (i=0; i<l_fcuts.size(); i+=2) {
        f1.push_back(l_fcuts[i]);
        f2.push_back(l_fcuts[i+1]);
    }

#ifdef DEBUG
    cout << "f1 "; 
    for (std::vector<double>::iterator it = f1.begin(); it != f1.end(); it++)
        std::cout << *it << ' ';
    cout << endl;

    cout << "f2 "; 
    for (std::vector<double>::iterator it = f2.begin(); it != f2.end(); it++)
        std::cout << *it << ' ';
    cout << endl;
#endif

    n = 0;
    min_delta = 2;

    for (i=0; i<f1.size(); i++) {
        if (f2[i] - f1[i] < min_delta) {
            n = i;
            min_delta = f2[i] - f1[i];
        }
    }

    if (nbands == 2) {
        l = lporder (f1[n], f2[n], l_devs[0], l_devs[1]);
    } else {
        l = 0;
        for (i=1; i<nbands - 1; i++) {
            double l1, l2;
            l1 = lporder(f1[i-1], f2[i-1], l_devs[i], l_devs[i-1]);
            l2 = lporder(f1[i], f2[i], devs[i], devs[i+1]);
            std::initializer_list<double> list = {l, l1, l2};
            l = std::max(list);
        }
    }

    n = (int)ceil(l) - 1;

    //std::vector<double> ff;
    
    fo.clear();
    fo.push_back(0);
    for (i=0; i<l_fcuts.size(); i++) {
        fo.push_back(l_fcuts[i]);
    }
    fo.push_back(1);

    for (i=0; i<fo.size() - 1; i++) {
        fo[i] *= 2;
    }

    for (std::vector<double>::iterator it=l_mags.begin(); it != l_mags.end(); it++) {
        ao.push_back(*it);
        ao.push_back(*it);
    }

    max_dev = *std::max_element(l_devs.begin(), l_devs.end());
    std::vector<double> wts(l_devs.size(), 1);

    for (i=0; i<wts.size(); i++) {
        wts[i] = max_dev / l_devs[i];
    }

    w = wts;
}

double lporder (double freq1, double freq2, double delta_p, double delta_s)
{
    /*
    FIR lowpass filter length estimator.  freq1 and freq2 are
    normalized to the sampling frequency.  delta_p is the passband
    deviation (ripple), delta_s is the stopband deviation (ripple).

    Note, this works for high pass filters too (freq1 > freq2), but
    doesnt work well if the transition is near f == 0 or f == fs/2

    From Herrmann et al (1973), Practical design rules for optimum
    finite impulse response filters.  Bell System Technical J., 52, 769-99
    */
    
    double df = fabs(freq2 - freq1);
    double ddp = log10 (delta_p);
    double dds = log10 (delta_s);

    double a1, a2, a3, a4, a5, a6, b1, b2, t1, t2, t3, t4;
    double dinf, ff, n;

    a1 = 5.309e-3;
    a2 = 7.114e-2;
    a3 = -4.761e-1;
    a4 = -2.66e-3;
    a5 = -5.941e-1;
    a6 = -4.278e-1;

    b1 = 11.01217;
    b2 = 0.5124401;

    t1 = a1 * ddp * ddp;
    t2 = a2 * ddp;
    t3 = a4 * ddp * ddp;
    t4 = a5 * ddp;

    dinf=((t1 + t2 + a3) * dds) + (t3 + t4 + a6);
    ff = b1 + b2 * (ddp - dds);
    n = dinf / df - ff * df + 1;
    return n;
}
