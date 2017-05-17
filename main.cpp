#include "fm_tx.h"
#include <signal.h>
#include <getopt.h>
#include <stdlib.h>

using namespace std;

fm_tx *fm = NULL;

void signal_handler(int signum)
{
        cout << "signal caught: " << signum << endl;
        delete fm;
        exit(0);
}

int main(int argc, char *argv[])
{
    std::string file = "";
    std::string out_device = "hackrf=228acf";
    int opt = 0;
    float freq = 433e6;
    float quad_rate = 320000.0f;
    float audio_rate = 32000.0f;
    float output_rate = 320000.0f;

    while ((opt = getopt(argc, argv, "f:p:hq:a:o:")) > 0) {
        switch (opt) {
            case 'f':
                freq = strtof(optarg, NULL);
                break;
            case 'p':
                file = optarg;
                break;
            case 'q':
                quad_rate = strtof(optarg, NULL);
                break;
            case 'a':
                audio_rate = strtof(optarg, NULL);
            case 'o':
                output_rate = strtof(optarg, NULL);
            case 'h':
            default:
                cout << "usage: " << basename(argv[0]) << " [-f freq] -p wavfile(32k)" << endl;
                return 0;
            }
    }

    fm = new fm_tx(file, out_device, 1, audio_rate, quad_rate, output_rate);

    signal(SIGINT, signal_handler);
    signal(SIGTERM, signal_handler);
    signal(SIGQUIT, signal_handler);
    signal(SIGSEGV, signal_handler);
    
    fm->set_rf_freq(freq);
    fm->start();

    while (1) {
        usleep(100000);
    }

    return 0;
}
