#ifdef __INTEL_COMPILER
#include <mathimf.h>
#else // other C compilers
#include <complex.h>
#include <math.h>
#endif // ?__INTEL_COMPILER
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "load_data.c"
#include "split_data.c"
#include "main.c"
