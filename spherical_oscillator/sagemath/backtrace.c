#include <stdio.h>
#include <stdlib.h>
#include "backtrace.h"

#define ARR_SIZE 50
#define BACKTRACE_SIZ 64

void show_backtrace(void)
{
#if 0
#ifndef _MSC_VER
    void    *array[BACKTRACE_SIZ];
    size_t   size, i;
    char   **strings;

    size = backtrace(array, BACKTRACE_SIZ);
    strings = backtrace_symbols(array, size);

    for (i = 0; i < size; i++) {
        printf("%p : %s\n", array[i], strings[i]);
    }

    free(strings);  // malloced by backtrace_symbols
#endif
#endif
}

void catch_signal(int sig)
{
	printf("Show backtrace Segmentation fault\n");
	show_backtrace();
	exit(0);
}