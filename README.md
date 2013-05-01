## Configuration of local installations

	export PATH=$HOME/.usr/bin:$PATH
	export CPATH=$HOME/.usr/include
	export LD_LIBRARY_PATH=$HOME/.usr/lib:$LD_LIBRARY_PATH
	export LIBRARY_PATH=$LD_LIBRARY_PATH
	export MANPATH=$HOME/.usr/share/man:$MANPATH

## Install

First create all autoconf and automake files with

	autoreconf

For local installations use

	./configure --prefix=$HOME/.usr

to install tfbayes to *$HOME/.usr*. Now the source can be compiled and installed with

	make
	make install
