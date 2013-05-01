## Configuration of local installations

export PATH="$HOME/.usr/bin:$PATH"
export CPATH=$HOME/.usr/include
export LD_LIBRARY_PATH=$HOME/.usr/lib:$LD_LIBRARY_PATH
export LIBRARY_PATH=$LD_LIBRARY_PATH
export MANPATH=$HOME/.usr/share/man:$MANPATH

## Install

autoreconf
./configure --prefix=$HOME/.usr
make
make install
