#ifndef _FM_TX_H
#define _FM_TX_H

#include <gnuradio/blocks/multiply_const_ff.h>
#include <gnuradio/blocks/multiply_cc.h>
#include <gnuradio/blocks/wavfile_source.h>
#include <gnuradio/blocks/null_sink.h>
#include <gnuradio/analog/sig_source_c.h>
#include <gnuradio/analog/frequency_modulator_fc.h>
#include <gnuradio/filter/iir_filter_ffd.h>
#include <gnuradio/filter/interp_fir_filter_fff.h>
#include <gnuradio/prefs.h>
#include <gnuradio/top_block.h>
#include <osmosdr/sink.h>
#include <osmosdr/ranges.h>
#include <string>
#include "resampler_xx.h"

class fm_tx
{
    public:
        enum status {
            STATUS_OK    = 0, /*!< Operation was successful. */
            STATUS_ERROR = 1  /*!< There was an error. */
        };
        
        enum filter_shape {
            FILTER_SHAPE_SOFT = 0,   /*!< Soft: Transition band is TBD of width. */
            FILTER_SHAPE_NORMAL = 1, /*!< Normal: Transition band is TBD of width. */
            FILTER_SHAPE_SHARP = 2   /*!< Sharp: Transition band is TBD of width. */
        };

        fm_tx(const std::string filename,
             const std::string output_device,
             unsigned int decimation);
        void start(void);
        void stop(void);
        void set_output_device(const std::string device);
        void set_antenna(const std::string &antenna);
        double set_input_rate(double rate);
        double get_input_rate(void) const;

        unsigned int set_input_decim(unsigned int decim);
        unsigned int get_input_decim(void) const { return d_decim; }

        double set_analog_bandwidth(double bw);
        double get_analog_bandwidth(void) const;

        status      set_rf_freq(double freq_hz);
        double      get_rf_freq(void);
        status      get_rf_range(double *start, double *stop, double *step);

        status      set_auto_gain(bool automatic);
        status      set_gain(std::string name, double value);
        double      get_gain(std::string name) const;

        status      set_filter(double low, double high, filter_shape shape);
        status      set_freq_corr(double ppm);
        
    private:
        bool        d_running;          /*!< Whether receiver is running or not. */
        double      d_output_rate;       /*!< Input sample rate. */
        double      d_quad_rate;        /*!< Quadrature rate (input_rate / decim) */
        double      d_audio_rate;       /*!< Audio output rate. */
        unsigned int    d_decim;        /*!< input decimation. */
        double      d_rf_freq;          /*!< Current RF frequency. */
        double      d_sample_rate;
        double      d_filter_offset;    /*!< Current filter offset */
        double      d_max_dev;

        std::string output_devstr; /*!< Current output device string. */
        gr::top_block_sptr        tb;        /*!< The GNU Radio top block. */

        osmosdr::sink::sptr       output;       /*!< Real time I/Q sink. */

        gr::analog::sig_source_c::sptr      lo;  /*!< oscillator used for tuning. */
        gr::blocks::multiply_cc::sptr       mixer;

        gr::blocks::multiply_const_ff::sptr audio_gain0; /*!< Audio gain block. */
        gr::blocks::multiply_const_ff::sptr audio_gain1; /*!< Audio gain block. */

        gr::blocks::wavfile_source::sptr    wav_src;    /*!< WAV file source for playback. */
        gr::blocks::null_sink::sptr         audio_null_sink0; /*!< Audio null sink used during playback. */
        gr::blocks::null_sink::sptr         audio_null_sink1; /*!< Audio null sink used during playback. */

        resampler_cc_sptr resample; /*!< Sniffer resampler. */
        gr::analog::frequency_modulator_fc::sptr modulator;
        gr::filter::iir_filter_ffd::sptr    preemph;
        gr::filter::interp_fir_filter_fff::sptr interpolator;

#if 0
#ifdef WITH_PULSEAUDIO
        pa_sink_sptr              audio_snk;  /*!< Pulse audio sink. */
#else
        gr::audio::sink::sptr     audio_snk;  /*!< gr audio sink */
#endif
#endif
};

#endif
