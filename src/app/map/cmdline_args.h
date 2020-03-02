#ifndef __CMDLINE_ARGS_H
#define __CMDLINE_ARGS_H

#include "hbn_options.h"

#ifdef __cplusplus
extern "C" {
#endif

void
ParseHbnProgramCmdLineArguments(int argc, char* argv[], HbnProgramOptions* opts);

char*
HbnProgramOptions2String(const HbnProgramOptions* opts);

#ifdef __cplusplus
}
#endif

#endif // __CMDLINE_ARGS_H