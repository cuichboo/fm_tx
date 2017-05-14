#include "fm_tx.h"
#include <signal.h>

using namespace std;

fm_tx *fm = NULL;

void signal_handler(int signum)
{
        cout << "signal caught: " << signum << endl;
        fm->stop();
        exit(0);
}

int main(int argc, char **argv)
{
    std::string file = argv[1];
    std::string out_device = "hackrf=228acf";
    fm_tx *fm = new fm_tx(file, out_device, 1);
    double freq = 433e6;

    fm->start();

    while (1) {
        usleep(100000);
        /*
        fm->set_rf_freq(freq);
        freq+=10e3;
        cout << freq << endl;
        if (freq >= 434e6)
            freq = 433e6;
            */
    }

    return 0;
}
