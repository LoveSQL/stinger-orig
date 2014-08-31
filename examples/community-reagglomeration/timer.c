#define _XOPEN_SOURCE 600
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#ifdef __MTA__

#include <sys/mta_task.h>
#include <machine/runtime.h>

void
init_timer (void)
{
    /* Empty. */
}

double
timer (void)
{
    return((double)mta_get_clock(0) / mta_clock_freq());
}

double
timer_getres (void)
{
    /* A guess. */
    return 1.0 / mta_clock_freq();
}

#else /* Elsewhere */

static clockid_t clockid;

#if defined(CLOCK_REALTIME_ID)
#define CLKID CLOCK_REALTIME_ID
#define CLKIDNAME "CLOCK_REALTIME_ID"
#elif defined(CLOCK_THREAD_CPUTIME_ID)
#define CLKID CLOCK_THREAD_CPUTIME_ID
#define CLKIDNAME "CLOCK_THREAD_CPUTIME_ID"
#elif defined(CLOCK_REALTIME_ID)
#warning "Falling back to realtime clock."
#define CLKID CLOCK_REALTIME_ID
#define CLKIDNAME "CLOCK_REALTIME_ID"
#else
#error "Cannot find a clock!"
#endif

void
init_timer (void)
{
    int err;
    err = clock_getcpuclockid (0, &clockid);
    if (err >= 0) return;
    fprintf (stderr, "Unable to find CPU clock, falling back to "
	     CLKIDNAME "\n");
    clockid = CLKID;
}

double
timer (void)
{
    struct timespec tp;
    clock_gettime(clockid, &tp);
    return (double)tp.tv_sec + 1.0e-9 * (double)tp.tv_nsec;
}

double
timer_getres (void)
{
    struct timespec tp;
    clock_getres(clockid, &tp);
    return (double)tp.tv_sec + 1.0e-9 * (double)tp.tv_nsec;
}

#endif

static double last_tic = -1.0;

void
tic (void)
{
    last_tic = timer ();
}

double
toc (void)
{
    const double t = timer ();
    const double out = t - last_tic;
    last_tic = t;
    return out;
}
