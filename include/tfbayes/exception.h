/* Copyright (C) 2010 Philipp Benner
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef EXCEPTION_H
# define EXCEPTION_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <stdarg.h>


#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

extern int verbose;

extern void std_warn(int mode, const char *msg, ...)
__attribute__ ((format(printf,2,3)));

extern void vstd_warn(int mode, const char *msg, va_list az)
__attribute__ ((format(printf,2,0)));

extern void std_err(int mode, const char *msg, ...)
__attribute__ ((format(printf,2,3), noreturn));

extern void vstd_err(int mode, const char *msg, va_list az)
__attribute__ ((format(printf,2,0), noreturn));

extern void s_warn(char *buf, size_t bufsize, int mode, const char *msg, ...)
__attribute__ ((format(printf,4,5)));

extern void vs_warn(char *buf, size_t bufsize, int mode, const char *msg, va_list az)
__attribute__ ((format(printf,4,0)));

extern void log_notice(int mode, const char *msg, ...)
__attribute__ ((format(printf,2,3)));

extern void vlog_notice(int mode, const char *msg, va_list az)
__attribute__ ((format(printf,2,0)));

extern void std_notice(int mode, const char *msg, ...)
__attribute__ ((format(printf,2,3)));

extern void vstd_notice(int mode, const char *msg, va_list az)
__attribute__ ((format(printf,2,0)));

extern void log_warn(int mode, const char *msg, ...)
__attribute__ ((format(printf,2,3)));

extern void vlog_warn(int mode, const char *msg, va_list az)
__attribute__ ((format(printf,2,0)));

extern void log_err(int mode, const char *msg, ...)
__attribute__ ((format(printf,2,3), noreturn));

extern void vlog_err(int mode, const char *msg, va_list az)
__attribute__ ((format(printf,2,0), noreturn));

#define DEBUG_COLOR(x)	"\033[0;31;1m"x"\033[0m"

#ifdef DEBUG
#  define print_debug(fmt, args...)                                             \
                fprintf(stderr, DEBUG_COLOR("DEBUG("__FILE__":%i): " fmt),      \
                        __LINE__, ##args)
#else
#  define print_debug(fmt, args...)
#endif

/* do exit(EXIT_CODE); for *_err functions*/

#define EXIT_CODE(a)	(a<<2)
#define EXIT		EXIT_CODE(1)

/* message type for all functions */

#define NONE		0	/* no extra error message */
#define PERR		1	/* sys error message */
#define HERR		2	/* network error */

/* measure execution time */
#include <time.h>

#define MET_INIT clock_t exception_met_start__;

#define MET(MSG, ...)                           \
        exception_met_start__ = clock();        \
        __VA_ARGS__;                            \
        fprintf(stderr, DEBUG_COLOR(""__FILE__":%i: " MSG " took %f seconds\n"), \
                __LINE__, ((double)clock()-exception_met_start__)/CLOCKS_PER_SEC);

__END_DECLS

#endif /* EXCEPTION_H */
