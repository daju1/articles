
#include "tzap.h"
#include <assert.h>
//#define USE_DEBUG
#include "dbg_info.h"
#include "backtrace.h"


int main()
{
	signal(SIGSEGV, catch_signal);
#ifndef _MSC_VER
	signal(SIGBUS, catch_signal);
#endif
#ifdef ALGORITHM_VERSION_0
	return test_v1();
	return test_v0();
#endif
#ifdef ALGORITHM_VERSION_1
	return test_v1();
#endif
}