/* Copyright (C) 2013, 2014 Philipp Benner
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

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <iostream>

#include <cstdio>
#include <cstdlib>
#include <fcntl.h>
#include <getopt.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/types.h>
#include <cerrno>
#include <csignal>

#include <boost/format.hpp>

using boost::format;
using namespace std;

// Options
////////////////////////////////////////////////////////////////////////////////

typedef struct _options_t {
        size_t verbose;
        _options_t()
                : verbose(0)
                { }
} options_t;

static options_t options;

// Usage
////////////////////////////////////////////////////////////////////////////////

static
void print_usage(char *pname, FILE *fp)
{
        (void)fprintf(fp, "\nUsage: mkfifo FILE; COMMAND 2>&1 | %s [OPTION] FILE\n\n", pname);
        (void)fprintf(fp,
                      "Options:\n"
                      "             -v              - be verbose\n"
                      "   --help                    - print help and exit\n"
                      "   --version                 - print version information and exit\n\n");
}

void wrong_usage(char *pname, const char *msg)
{

        if(msg != NULL) {
                (void)fprintf(stderr, "%s\n", msg);
        }
        (void)fprintf(stderr,
                      "Try `%s --help' for more information.\n", pname);

        exit(EXIT_FAILURE);

}

static
void print_version(FILE *fp)
{
        (void)fprintf(fp,
                      "This is free software, and you are welcome to redistribute it\n"
                      "under certain conditions; see the source for copying conditions.\n"
                      "There is NO warranty; not even for MERCHANTABILITY or FITNESS\n"
                      "FOR A PARTICULAR PURPOSE.\n\n");
}

// Main
////////////////////////////////////////////////////////////////////////////////

int print_fifo(string filename)
{
        int readfd, writefd;
        struct stat status;
        char buffer[BUFSIZ];
        ssize_t bytes;

        signal(SIGPIPE, SIG_IGN);

        readfd = open(filename.c_str(), O_RDONLY | O_NONBLOCK);
        if(readfd == -1) {
                perror(str(format("open() %s failed") % filename).c_str());
                exit(EXIT_FAILURE);
        }
        if(fstat(readfd, &status) == -1) {
                perror("fstat() failed");
                close(readfd);
                exit(EXIT_FAILURE);
        }
        if(!S_ISFIFO(status.st_mode)) {
                cerr << format("%s is not a fifo!\n") % filename
                     << endl;
                close(readfd);
                exit(EXIT_FAILURE);
        }
        writefd = open(filename.c_str(), O_WRONLY | O_NONBLOCK);
        if(writefd == -1) {
                perror("open() failed");
                close(readfd);
                exit(EXIT_FAILURE);
        }
        close(readfd);

        while(true) {
                bytes = read(STDIN_FILENO, buffer, sizeof(buffer));
                if (bytes < 0 && errno == EINTR)
                        continue;
                if (bytes <= 0)
                        break;
                
                bytes = write(STDOUT_FILENO, buffer, bytes);
                if(bytes == -1) {
                        perror("writing to stdout failed");
                }
                bytes = write(writefd, buffer, bytes);
        }
        close(writefd); 

        return 0;
}

int main(int argc, char *argv[])
{
        for(;;) {
                int c, option_index = 0;
                static struct option long_options[] = {
                        { "help",              0, 0, 'h' },
                        { "version",           0, 0, 'q' }
                };

                c = getopt_long(argc, argv, "v",
                                long_options, &option_index);

                if(c == -1) {
                        break;
                }

                switch(c) {
                case 'v':
                        options.verbose = 1;
                        break;
                case 'h':
                        print_usage(argv[0], stdout);
                        exit(EXIT_SUCCESS);
                case 'q':
                        print_version(stdout);
                        exit(EXIT_SUCCESS);
                default:
                        wrong_usage(argv[0], NULL);
                        exit(EXIT_FAILURE);
                }
        }
        if(optind+1 != argc) {
                wrong_usage(argv[0], "Wrong number of arguments.");
                exit(EXIT_FAILURE);
        }

        return print_fifo(argv[optind]);
}
