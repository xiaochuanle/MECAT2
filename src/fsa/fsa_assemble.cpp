#include "assembly.hpp"

int main(int argc, char *argv[])
{
    Assembly ass;
    if (ass.ParseArgument(argc, argv)) {
        ass.Run();
    }
    else {
        ass.Usage();
    }

	return 0;

}