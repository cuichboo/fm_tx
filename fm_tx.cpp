#include <iostream>
#include <gnuradio/blocks/multiply_const_ff.h>
#include <gnuradio/blocks/wavfile_source.h>
#include <gnuradio/analog/frequency_modulator_fc.h>
#include <gnuradio/filter/iir_filter_ffd.h>
#include <gnuradio/filter/interp_fir_filter_fff.h>
#include <gnuradio/prefs.h>
#include <gnuradio/top_block.h>
#include <osmosdr/source.h>
#include <osmosdr/ranges.h>
#include <gnuradio/gr_complex.h>
#include <math.h>
#include "optfir.h"
#include "fm_tx.h"

fm_tx::fm_tx(const std::string filename,
             const std::string output_device,
             unsigned int decimation,
             double audio_rate,
             double quad_rate,
             double out_rate
             )
    :d_output_rate(out_rate),
    d_quad_rate(quad_rate),
    d_audio_rate(audio_rate),
    d_decim(decimation),
    d_rf_freq(433000000.0),
    d_max_dev(75e3)
{
    bool do_interp = false;
    float k, fh, fs, w_ch, w_cl, w_cla, w_cha, tau;
    float kl, kh, z1, p1, b0;
    float w_0dB, g;
    std::vector<double> btaps(2), ataps(2);
    bool input_is_file = true;

    tb = gr::make_top_block("wbfm_tx");
    
    if (!filename.empty()) {
        wav_src = gr::blocks::wavfile_source::make(filename.c_str(), true);
        input_is_file = true;
    } else {
        audio_src = make_pa_source("", d_audio_rate, 1, "wbfm_tx", "record");
        input_is_file = false;
    }

    output = osmosdr::sink::make(output_device);
    resampler = make_resampler_cc((float)d_output_rate/(float)d_quad_rate);

    do_interp = d_audio_rate != d_quad_rate;
    if (do_interp) {
        float interp_factor = d_quad_rate / d_audio_rate;
        std::vector<double> interp_taps = optfir_lowpass(interp_factor,   // gain
                                             d_quad_rate,       // Fs
                                             16000,           // passband cutoff
                                             18000,           // stopband cutoff
                                             0.1,             // passband ripple dB
                                             40);             // stopband atten dB
        std::vector<float> interp(interp_taps.size());
        std::cout << "n_taps: " << interp_taps.size() << std::endl;
        
        for (std::vector<double>::iterator it=interp_taps.begin(); it!=interp_taps.end(); it++) {
            interp.push_back(*it);
        }
        interpolator = gr::filter::interp_fir_filter_fff::make(interp_factor, interp);
    }

    fh = -1.0;
    tau = 300e-6;
    fs = d_quad_rate;
    k = 2 * M_PI * d_max_dev / d_quad_rate;
    modulator = gr::analog::frequency_modulator_fc::make(k);

    if ((fh <= 0.0f) || (fh >= fs/2.0f)) {
        fh = 0.925 * fs/2.0;
    }

	/* Digital corner frequencies */
	w_cl = 1.0 / tau;
	w_ch = 2.0 * M_PI * fh;

	/* Prewarped analog corner frequencies */
	w_cla = 2.0 * fs * tanf(w_cl / (2.0 * fs));
	w_cha = 2.0 * fs * tanf(w_ch / (2.0 * fs));

	/* Resulting digital pole, zero, and gain term from the bilinear
     * transformation of H(s) = (s + w_cla) / (s + w_cha) to
     * H(z) = b0 (1 - z1 z^-1)/(1 - p1 z^-1)
     * */
	kl = -w_cla / (2.0 * fs);
	kh = -w_cha / (2.0 * fs);
	z1 = (1.0 + kl) / (1.0 - kl);
	p1 = (1.0 + kh) / (1.0 - kh);
	b0 = (1.0 - kl) / (1.0 - kh);

	/* # Since H(s = infinity) = 1.0, then H(z = -1) = 1.0 and
     * # this filter  has 0 dB gain at fs/2.0.
     * # That isn't what users are going to expect, so adjust with a
     * # gain, g, so that H(z = 1) = 1.0 for 0 dB gain at DC.
     * */

	w_0dB = 2.0 * M_PI * 0.0;
#if 0
	g =        fabs(1.0 - p1 * cmath.rect(1.0, -w_0dB))  \
	   / (b0 * fabs(1.0 - z1 * cmath.rect(1.0, -w_0dB)));
#endif

    gr_complex tmp1 = gr_complex(1.0, 0.0) - p1 * gr_complex(cos(-w_0dB), sin(-w_0dB));
    gr_complex tmp2 = gr_complex(1.0, 0.0) - z1 * gr_complex(cos(-w_0dB), sin(-w_0dB));
    g = sqrt(pow(tmp1.real(), 2.0) + pow(tmp1.imag(), 2)) \
        / sqrt(pow(tmp2.real(), 2.0) + pow(tmp2.imag(), 2));

    btaps[0] = g * b0 * 1.0;
    btaps[1] = g * b0 * -z1;
    ataps[0] = 1.0;
    ataps[1] = -p1;

    preemph = gr::filter::iir_filter_ffd::make(btaps, ataps, false);

    if (input_is_file) {
        tb->connect(wav_src, 0, interpolator, 0);
    } else {
        tb->connect(audio_src, 0, interpolator, 0);
    }
    tb->connect(interpolator, 0, preemph, 0);
    tb->connect(preemph, 0, modulator, 0);


#if (0 == 1)
    tb->connect(modulator, 0, output, 0);
#elif (1 == 1)
    std::cout << "use rational_resampler" << std::endl;
    tb->connect(modulator, 0, rational_resampler, 0);
    tb->connect(rational_resampler, 0, output, 0);
#else
    tb->connect(modulator, 0, resampler, 0);
    tb->connect(resampler, 0, output, 0);
#endif
    
    output->set_bb_gain(20);
    output->set_if_gain(20);
    output->set_sample_rate(d_output_rate);
    output->set_center_freq(d_rf_freq);
}

fm_tx::~fm_tx(void)
{
    tb->stop();
}

void fm_tx::start(void)
{
    return tb->start();
}

void fm_tx::stop(void)
{
    return tb->stop();
}

fm_tx::status fm_tx::set_rf_freq(double freq_hz)
{
    d_rf_freq = freq_hz;
    output->set_center_freq(freq_hz);

    return fm_tx::STATUS_OK;
}
