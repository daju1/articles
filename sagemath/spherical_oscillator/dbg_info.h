
#ifdef USE_DEBUG
#define DBG_INFO(fmt, args, ...) \
    do \
    { \
        printf(fmt, ## args); \
        /*printf("\n");*/ \
    } \
    while(0)
#else
#define DBG_INFO(fmt, args, ...)
#endif

#ifdef USE_DEBUG
#define DBG_ERROR(fmt, args, ...) \
    do \
    { \
        printf(fmt, ## args); \
        /*printf("\n");*/ \
    } \
    while(0)
#else
#define DBG_ERROR(fmt, args, ...)
#endif