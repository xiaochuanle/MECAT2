#include "overlap_show.hpp"


int main(int argc, char *argv[]) {
    OverlapShow os;

    if (os.ParseArgument(argc, argv)) {
        os.Run();
    }
    else {
        os.Usage();
    }
    return 0;

}
