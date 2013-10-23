#!/usr/bin/env python

# -------------------------------- WebLogo --------------------------------

#  Copyright (c) 2003-2004 The Regents of the University of California.
#  Copyright (c) 2005 Gavin E. Crooks
#  Copyright (c) 2006, The Regents of the University of California, through 
#  Lawrence Berkeley National Laboratory (subject to receipt of any required
#  approvals from the U.S. Dept. of Energy).  All rights reserved.

#  This software is distributed under the new BSD Open Source License.
#  <http://www.opensource.org/licenses/bsd-license.html>
#
#  Redistribution and use in source and binary forms, with or without 
#  modification, are permitted provided that the following conditions are met: 
#
#  (1) Redistributions of source code must retain the above copyright notice, 
#  this list of conditions and the following disclaimer. 
#
#  (2) Redistributions in binary form must reproduce the above copyright 
#  notice, this list of conditions and the following disclaimer in the 
#  documentation and or other materials provided with the distribution. 
#
#  (3) Neither the name of the University of California, Lawrence Berkeley 
#  National Laboratory, U.S. Dept. of Energy nor the names of its contributors 
#  may be used to endorse or promote products derived from this software 
#  without specific prior written permission. 
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
#  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
#  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
#  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
#  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
#  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
#  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
#  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
#  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
#  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
#  POSSIBILITY OF SUCH DAMAGE. 

# Replicates README.txt

"""
WebLogo (http://code.google.com/p/weblogo/) is a tool for creating sequence 
logos from biological sequence alignments.  It can be run on the command line,
as a standalone webserver, as a CGI webapp, or as a python library.

The main WebLogo webserver is located at http://bespoke.lbl.gov/weblogo/

Please consult the manual for installation instructions and more information:
(Also located in the weblogolib/htdocs subdirectory.)

    http://bespoke.lbl.gov/weblogo/manual.html

For help on the command line interface run
    ./weblogo --help

To build a simple logo run
    ./weblogo  < cap.fa > logo0.eps
    
To run as a standalone webserver at localhost:8080 
    ./weblogo --server

To create a logo in python code:
    >>> from weblogolib import *
    >>> fin = open('cap.fa')
    >>> seqs = read_seq_data(fin) 
    >>> data = LogoData.from_seqs(seqs)
    >>> options = LogoOptions()
    >>> options.title = "A Logo Title"
    >>> format = LogoFormat(data, options)
    >>> fout = open('cap.eps', 'w') 
    >>> eps_formatter( data, format, fout)


-- Distribution and Modification --
This package is distributed under the new BSD Open Source License. 
Please see the LICENSE.txt file for details on copyright and licensing.
The WebLogo source code can be downloaded from http://code.google.com/p/weblogo/
WebLogo requires Python 2.3, 2.4 or 2.5, the corebio python toolkit for computational 
biology (http://code.google.com/p/corebio), and the python array package 
'numpy' (http://www.scipy.org/Download)
"""

import sys
import stat
import copy
import os

from datetime import datetime
from StringIO import StringIO
from subprocess import Popen
from subprocess import PIPE

import random

from math import log, sqrt
from string import Template
from pkg_resources import resource_string

# Avoid 'from numpy import *' since numpy has lots of names defined
from numpy import array, asarray, float64, ones, zeros, int32, all, any, shape
import numpy as na

from color import *
from colorscheme import *

from ..corebio.seq import Alphabet, Seq, SeqList
from ..corebio.moremath import *
from ..corebio.data import amino_acid_composition
from ..corebio.seq import unambiguous_rna_alphabet, unambiguous_dna_alphabet, unambiguous_protein_alphabet

# ------ META DATA ------

__all__ = ['LogoSize', 
            'LogoOptions', 
            'description', 
            '__version__', 
            'LogoFormat',
            'LogoData',
            'Dirichlet',
            'GhostscriptAPI',
            'std_color_schemes',
            'default_color_schemes',
            'classic',
            'std_units',
            'std_sizes',
            'std_alphabets',
            'std_percentCG',
            'pdf_formatter',
            'jpeg_formatter',
            'png_formatter',
            'png_print_formatter',
            'txt_formatter',
            'eps_formatter',
            'formatters',
            'default_formatter',
            'base_distribution',
            'equiprobable_distribution',
            'which_alphabet',
            'color',
            'colorscheme',
            ]

description  = "Create sequence logos from biological sequence alignments." 

__version__ = "3.0"

# These keywords are subsituted by subversion.
# The date and revision will only  tell the truth after a branch or tag,
# since different files in trunk will have been changed at different times
release_date ="$Date: 2008-10-14 16:38:51 -0700 (Tue, 14 Oct 2008) $".split()[1]
release_build = "$Revision: 53 $".split()[1]
release_description = "WebLogo %s (%s)" % (__version__,  release_date)



def cgi(htdocs_directory) :
    import weblogolib._cgi
    weblogolib._cgi.main(htdocs_directory)
        
class GhostscriptAPI(object) :
    """Interface to the command line program Ghostscript ('gs')"""
    
    formats = ('png', 'pdf', 'jpeg')
    
    def __init__(self, path=None) :
        try:
            command = find_command('gs', path=path)
        except EnvironmentError:
            try:
                command = find_command('gswin32c.exe', path=path)
            except EnvironmentError:
                raise EnvironmentError("Could not find Ghostscript on path."
                " There should be either a gs executable or a gswin32c.exe on your system's path")
        
        self.command = command        
    
    def version(self) :
        args = [self.command, '--version']
        try :
            p = Popen(args, stdout=PIPE)
            (out,err) = p.communicate() 
        except OSError :
            raise RuntimeError("Cannot communicate with ghostscript.")  
        return out.strip()
       
    def convert(self, format, fin, fout,  width,  height, resolution=300) :
        device_map = { 'png':'png16m',  'pdf':'pdfwrite', 'jpeg':'jpeg'}
       
        try :
            device = device_map[format]
        except KeyError:
            raise ValueError("Unsupported format.")
        
        args = [self.command, 
            "-sDEVICE=%s" % device, 
            "-dPDFSETTINGS=/printer",
            #"-q",   # Quite: Do not dump messages to stdout.
            "-sstdout=%stderr", # Redirect messages and errors to stderr
            "-sOutputFile=-", # Stdout
            "-dDEVICEWIDTHPOINTS=%s" % str(width),  
            "-dDEVICEHEIGHTPOINTS=%s" % str(height),  
            "-dSAFER",  # For added security
            "-dNOPAUSE",]
            
        if device != 'pdf' :
            args.append("-r%s" % str(resolution) ) 
            if resolution < 300 : # Antialias if resolution is Less than 300 DPI
                args.append("-dGraphicsAlphaBits=4")
                args.append("-dTextAlphaBits=4")
                args.append("-dAlignToPixels=0")
        
        args.append("-")  # Read from stdin. Must be last argument.
        
        error_msg = "Unrecoverable error : Ghostscript conversion failed " \
                    "(Invalid postscript?). %s" % " ".join(args) 

        source = fin.read()
        try :
            p = Popen(args, stdin=PIPE, stdout = PIPE, stderr= PIPE) 
            (out,err) = p.communicate(source) 
        except OSError :
            raise RuntimeError(error_msg)
    
        if p.returncode != 0 : 
            error_msg += '\nReturn code: %i\n' % p.returncode 
            if err is not None : error_msg += err
            raise RuntimeError(error_msg)
        
        print >>fout, out
# end class Ghostscript


aa_composition = [ amino_acid_composition[_k] for _k in
                    unambiguous_protein_alphabet]



# ------  DATA ------

classic = ColorScheme([
    ColorGroup("G",  "orange" ),
    ColorGroup("TU", "red"),
    ColorGroup("C",  "blue"),
    ColorGroup("A",  "green")
    ] )
    
std_color_schemes = {"auto":            None,   # Depends on sequence type
                     "monochrome":      monochrome,
                     "base pairing":    base_pairing,
                     "classic":         classic,
                     "hydrophobicity" : hydrophobicity,  
                     "chemistry" :      chemistry,
                     "charge" :         charge,
                     }#

default_color_schemes = {
    unambiguous_protein_alphabet: hydrophobicity, 
    unambiguous_rna_alphabet: base_pairing,
    unambiguous_dna_alphabet: base_pairing 
}


std_units = {
    "bits"    : 1./log(2),
    "nats"    : 1.,
    "digits"  : 1./log(10),
    "kT"      : 1.,
    "kJ/mol"  : 8.314472 *298.15 /1000.,
    "kcal/mol": 1.987 *298.15  /1000.,
    "probability" : None,
}

class LogoSize(object) :
    def __init__(self, stack_width, stack_height) :
         self.stack_width = stack_width
         self.stack_height = stack_height
         
    def __repr__(self):
        return stdrepr(self)

# The base stack width is set equal to 9pt Courier. 
# (Courier has a width equal to 3/5 of the point size.)
# Check that can get 80 characters in journal page @small
# 40 chacaters in a journal column
std_sizes = {
    "small" : LogoSize( stack_width = 5.4, stack_height = 5.4*1*5), 
    "medium" : LogoSize( stack_width = 5.4*2, stack_height = 5.4*2*5),
    "large"  : LogoSize( stack_width = 5.4*3, stack_height = 5.4*3*5),    
}
            

            
std_alphabets = {
    'protein': unambiguous_protein_alphabet, 
    'rna': unambiguous_rna_alphabet,
    'dna': unambiguous_dna_alphabet}

std_percentCG = {
    'H. sapiens'    : 40.,
    'E. coli'       : 50.5,
    'S. cerevisiae' : 38.,
    'C. elegans'    : 36.,
    'D. melanogaster': 43.,
    'M. musculus'   :  42.,
    'T. thermophilus' : 69.4,
}

# Thermus thermophilus: Henne A, Bruggemann H, Raasch C, Wiezer A, Hartsch T,
# Liesegang H, Johann A, Lienard T, Gohl O, Martinez-Arias R, Jacobi C, 
# Starkuviene V, Schlenczeck S, Dencker S, Huber R, Klenk HP, Kramer W, 
# Merkl R, Gottschalk G, Fritz HJ: The genome sequence of the extreme 
# thermophile Thermus thermophilus.
# Nat Biotechnol 2004, 22:547-53

def update(obj, **entries):
    """Update an instance with new values. 

    >>> update({'a': 1}, a=10, b=20)
    {'a': 10, 'b': 20}
    """
    if hasattr(obj, 'update') :
        obj.update( entries) 
    else :
        for k, v in entries.iteritems() :
            setattr(obj, k, v)
    return obj

def stdrepr(obj) :
    attr = vars(obj).items()
    attr.sort()
    args = []
    for a in attr :
        if a[0][0]=='_' : continue
        args.append( '%s=%s' % ( a[0], repr(a[1])) )
    args = ',\n'.join(args).replace('\n', '\n    ')
    return '%s(\n    %s\n)' % (obj.__class__.__name__, args)

def _cull(potential, matches, verbose=0):
    """Cull inappropriate matches. Possible reasons:
        - a duplicate of a previous match
        - not a disk file
        - not executable (non-Windows)
    If 'potential' is approved it is returned and added to 'matches'.
    Otherwise, None is returned.
    """
    for match in matches:  # don't yield duplicates
        if _samefile(potential[0], match[0]):
            if verbose:
                sys.stderr.write("duplicate: %s (%s)\n" % potential)
            return None
    else:
        if not stat.S_ISREG(os.stat(potential[0]).st_mode):
            if verbose:
                sys.stderr.write("not a regular file: %s (%s)\n" % potential)
        elif not os.access(potential[0], os.X_OK):
            if verbose:
                sys.stderr.write("no executable access: %s (%s)\n"\
                                 % potential)
        else:
            matches.append(potential)
            return potential


def whichgen(command, path=None, verbose=0, exts=None):
    """Return a generator of full paths to the given command.
    
    "command" is a the name of the executable to search for.
    "path" is an optional alternate path list to search. The default it
        to use the PATH environment variable.
    "verbose", if true, will cause a 2-tuple to be returned for each
        match. The second element is a textual description of where the
        match was found.
    "exts" optionally allows one to specify a list of extensions to use
        instead of the standard list for this system. This can
        effectively be used as an optimization to, for example, avoid
        stat's of "foo.vbs" when searching for "foo" and you know it is
        not a VisualBasic script but ".vbs" is on PATHEXT. This option
        is only supported on Windows.

    This method returns a generator which yields either full paths to
    the given command or, if verbose, tuples of the form (<path to
    command>, <where path found>).
    """
    matches = []
    if path is None:
        usingGivenPath = 0
        path = os.environ.get("PATH", "").split(os.pathsep)
        if sys.platform.startswith("win"):
            path.insert(0, os.curdir)  # implied by Windows shell
    else:
        usingGivenPath = 1

    # Windows has the concept of a list of extensions (PATHEXT env var).
    if sys.platform.startswith("win"):
        if exts is None:
            exts = os.environ.get("PATHEXT", "").split(os.pathsep)
            # If '.exe' is not in exts then obviously this is Win9x and
            # or a bogus PATHEXT, then use a reasonable default.
            for ext in exts:
                if ext.lower() == ".exe":
                    break
            else:
                exts = ['.COM', '.EXE', '.BAT']
        elif not isinstance(exts, list):
            raise TypeError("'exts' argument must be a list or None")
    else:
        if exts is not None:
            raise WhichError("'exts' argument is not supported on "\
                             "platform '%s'" % sys.platform)
        exts = []

    # File name cannot have path separators because PATH lookup does not
    # work that way.
    if os.sep in command or os.altsep and os.altsep in command:
        pass
    else:
        for i in range(len(path)):
            dirName = path[i]
            # On windows the dirName *could* be quoted, drop the quotes
            if sys.platform.startswith("win") and len(dirName) >= 2\
               and dirName[0] == '"' and dirName[-1] == '"':
                dirName = dirName[1:-1]
            for ext in ['']+exts:
                absName = os.path.abspath(
                    os.path.normpath(os.path.join(dirName, command+ext)))
                if os.path.isfile(absName):
                    if usingGivenPath:
                        fromWhere = "from given path element %d" % i
                    elif not sys.platform.startswith("win"):
                        fromWhere = "from PATH element %d" % i
                    elif i == 0:
                        fromWhere = "from current directory"
                    else:
                        fromWhere = "from PATH element %d" % (i-1)
                    match = _cull((absName, fromWhere), matches, verbose)
                    if match:
                        if verbose:
                            yield match
                        else:
                            yield match[0]
        match = _getRegisteredExecutable(command)
        if match is not None:
            match = _cull(match, matches, verbose)
            if match:
                if verbose:
                    yield match
                else:
                    yield match[0]


def find_command(command, path=None):
    """Return the full path to the first match of the given command on
    the path.
    
    Arguments:
    - command -- is a the name of the executable to search for.
    - path -- is an optional alternate path list to search. The default is
        to use the COREBIOPATH environment variable, if it exists, else the 
        PATH environment variable.
        
    Raises:
    - EnvironmentError -- If no match is found for the command.
    
    By default the COREBIOPATH or PATH environment variable is searched (as well
    as, on Windows, the AppPaths key in the registry), but a specific 'path'
    list to search may be specified as well.  
        
    Author: Adapted from code by Trent Mick (TrentM@ActiveState.com)
    See: http://trentm.com/projects/which/
    """
    if path is None :
        path = os.environ.get("COREBIOPATH", "").split(os.pathsep)
        if path==['']: path = None   

    try :
        match = whichgen(command, path).next()
    except StopIteration, _which.WhichError:
        raise EnvironmentError("Could not find '%s' on the path." % command)
    return match


class LogoOptions(object) :
    """ A container for all logo formating options. Not all of these
    are directly accessible through the CLI or web interfaces. 
    
    To display LogoOption defaults:
    >>> from weblogolib import *
    >>> LogoOptions()
    
    
    Attributes:
        o alphabet
        o creator_text           -- Embedded as comment in figures.
        o logo_title             
        o logo_label
        o stacks_per_line
        o unit_name  
        o show_yaxis 
        o yaxis_label             -- Default depends on other settings.
        o yaxis_tic_interval 
        o yaxis_minor_tic_ratio 
        o yaxis_scale
        o show_xaxis 
        o xaxis_label
        o xaxis_tic_interval 
        o rotate_numbers
        o number_interval
        o show_ends 
        o show_fineprint
        o fineprint
        o show_boxes 
        o shrink_fraction
        o show_errorbars 
        o errorbar_fraction 
        o errorbar_width_fraction 
        o errorbar_gray 
        o resolution             -- Dots per inch
        o default_color 
        o color_scheme 
        o debug
        o logo_margin
        o stroke_width
        o tic_length
        o size
        o stack_margin
        o pad_right
        o small_fontsize
        o fontsize
        o title_fontsize
        o number_fontsize
        o text_font
        o logo_font
        o title_font
        o first_index
        o logo_start
        o logo_end
        o scale_width
    """       

    def __init__(self, **kwargs) :
        """ Create a new LogoOptions instance.
        
        >>> L = LogoOptions(logo_title = "Some Title String")
        >>> L.show_yaxis = False
        >>> repr(L)
        """

        self.creator_text = release_description,
        self.alphabet = None
        
        self.logo_title = ""
        self.logo_label = ""
        self.stacks_per_line = 40
        
        self.unit_name = "bits"
     
        self.show_yaxis = True
        # yaxis_lable default depends on other settings. See LogoFormat
        self.yaxis_label = None
        self.yaxis_tic_interval = 1.
        self.yaxis_minor_tic_ratio = 5
        self.yaxis_scale = None

        self.show_xaxis = True
        self.xaxis_label = ""
        self.xaxis_tic_interval =1
        self.rotate_numbers = False
        self.number_interval = 5            
        self.show_ends = False
          
        self.show_fineprint = True
        self.fineprint =  "WebLogo "+__version__ 
    
        self.show_boxes = False
        self.shrink_fraction = 0.5
  
        self.show_errorbars = True
        self.errorbar_fraction = 0.90
        self.errorbar_width_fraction = 0.25      
        self.errorbar_gray = 0.75
 
        self.resolution  = 96.     # Dots per inch
           
        self.default_color = Color.by_name("black")    
        self.color_scheme = None
        #self.show_color_key = False # NOT yet implemented
        
        self.debug = False
        
        self.logo_margin = 2
        self.stroke_width = 0.5
        self.tic_length = 5
        
        self.size = std_sizes["medium"]        
        
        self.stack_margin = 0.5
        self.pad_right = False  

        self.small_fontsize = 6
        self.fontsize = 10
        self.title_fontsize = 12
        self.number_fontsize = 8

        self.text_font    = "ArialMT"
        self.logo_font    = "Arial-BoldMT"
        self.title_font   = "ArialMT"

        self.first_index = 1        
        self.logo_start = None      
        self.logo_end=None          

        # Scale width of characters proportional to gaps
        self.scale_width = True

        #from corebio.utils import update
        update(self, **kwargs)

    def __repr__(self) :
        attr = vars(self).items()
        attr.sort()
        args = []
        for a in attr :
            if a[0][0]=='_' : continue
            args.append( '%s=%s' % ( a[0], repr(a[1])) )
        args = ',\n'.join(args).replace('\n', '\n    ')
        return '%s(\n    %s\n)' % (self.__class__.__name__, args)
# End class LogoOptions


        

class LogoFormat(LogoOptions) :     
    """ Specifies the format of the logo. Requires a LogoData and LogoOptions 
    objects.
    
    >>> data = LogoData.from_seqs(seqs )
    >>> options = LogoOptions()
    >>> options.title = "A Logo Title"
    >>> format = LogoFormat(data, options) 
    """
    # TODO: Raise ArgumentErrors instead of ValueError and document
    def __init__(self, data, options= None) :
        LogoOptions.__init__(self)
        
        if options is not None :
            self.__dict__.update(options.__dict__)
                 
        self.alphabet = data.alphabet
        self.seqlen = data.length
        
        self.show_title = False
        self.show_xaxis_label = False
        self.yaxis_minor_tic_interval = None
        self.lines_per_logo       = None
        self.char_width       = None
        self.line_margin_left = None
        self.line_margin_right    = None
        self.line_margin_bottom = None
        self.line_margin_top = None
        self.title_height = None
        self.xaxis_label_height = None
        self.line_height = None
        self.line_width = None
        self.logo_height = None
        self.logo_width = None
        self.creation_date = None
        self.end_type = None

        if self.stacks_per_line< 1 :
            raise ArgumentError("Stacks per line should be greater than zero.",
                                "stacks_per_line" )
        
        if self.size.stack_height<=0.0 : 
            raise ArgumentError( 
                "Stack height must be greater than zero.", "stack_height")
        
        if (self.small_fontsize <= 0 or self.fontsize <=0 or    
                self.title_fontsize<=0 ):
            raise ValueError("Font sizes must be positive.")
        
        if self.errorbar_fraction<0.0 or self.errorbar_fraction>1.0 :
            raise ValueError(
        "The visible fraction of the error bar must be between zero and one.")
        
        if self.yaxis_tic_interval<=0.0 :
            raise ArgumentError( "The yaxis tic interval cannot be negative.",
                'yaxis_tic_interval')
        
        if self.size.stack_width <= 0.0 :
            raise ValueError(
                "The width of a stack should be a positive number.")
        
        if self.yaxis_minor_tic_interval and \
                self.yaxis_minor_tic_interval<=0.0 : 
            raise ValueError("Distances cannot be negative.")
        
        if self.xaxis_tic_interval<=0 :
            raise ValueError("Tic interval must be greater than zero.")
        
        if self.number_interval<=0 :
            raise ValueError("Invalid interval between numbers.")
        
        if self.shrink_fraction<0.0 or self.shrink_fraction>1.0 :
            raise ValueError("Invalid shrink fraction.")
        
        if self.stack_margin<=0.0 : 
            raise ValueError("Invalid stack margin."  )
        
        if self.logo_margin<=0.0 : 
            raise ValueError("Invalid logo margin."  )
        
        if self.stroke_width<=0.0 : 
            raise ValueError("Invalid stroke width.")  
        
        if self.tic_length<=0.0 : 
            raise ValueError("Invalid tic length.") 

        # FIXME: More validation
         
        # Inclusive upper and lower bounds
        # FIXME: Validate here. Move from eps_formatter        
        if self.logo_start is  None: self.logo_start = self.first_index
        
        if self.logo_end is  None : 
            self.logo_end = self.seqlen + self.first_index -1 
        
        self.total_stacks = self.logo_end - self.logo_start +1

        if self.logo_start - self.first_index <0 :
            raise ArgumentError(
                "Logo range extends before start of available sequence.",
                'logo_range')
        
        if self.logo_end - self.first_index  >= self.seqlen  : 
            raise ArgumentError(
                "Logo range extends beyond end of available sequence.",
                'logo_range')
    
        if self.logo_title      : self.show_title = True
        if not self.fineprint   : self.show_fineprint = False
        if self.xaxis_label     : self.show_xaxis_label = True

        if self.yaxis_label is None : 
            self.yaxis_label = self.unit_name
        
        if self.yaxis_label : 
            self.show_yaxis_label = True
        else :
            self.show_yaxis_label = False
            self.show_ends = False
              
        if not self.yaxis_scale :
            conversion_factor = std_units[self.unit_name]
            if conversion_factor :
                self.yaxis_scale=log(len(self.alphabet))*conversion_factor
            else :
                self.yaxis_scale=1.0    # probability units

        if self.yaxis_scale<=0.0 : 
            raise ValueError(('yaxis_scale', "Invalid yaxis scale"))

        if self.yaxis_tic_interval >= self.yaxis_scale:
            self.yaxis_tic_interval /= 2. 

        self.yaxis_minor_tic_interval \
            = float(self.yaxis_tic_interval)/self.yaxis_minor_tic_ratio
                      
        if self.color_scheme is None :
            if self.alphabet in default_color_schemes :
                self.color_scheme = default_color_schemes[self.alphabet]
            else :
                self.color_scheme = monochrome

        self.lines_per_logo = 1+ ( (self.total_stacks-1) / self.stacks_per_line)
    
        if self.lines_per_logo==1 and not self.pad_right:
            self.stacks_per_line = min(self.stacks_per_line, self.total_stacks)

        self.char_width = self.size.stack_width - 2* self.stack_margin
    
    
        if self.show_yaxis :
            self.line_margin_left = self.fontsize * 3.0
        else :
            self.line_margin_left = 0

        if self.show_ends :
            self.line_margin_right = self.fontsize *1.5 
        else :
            self.line_margin_right = self.fontsize             

        if self.show_xaxis :
            if self.rotate_numbers :
                self.line_margin_bottom = self.number_fontsize *2.5
            else:
                self.line_margin_bottom = self.number_fontsize *1.5
        else :
            self.line_margin_bottom = 4

        self.line_margin_top = 4
    
        if self.show_title :
            self.title_height = self.title_fontsize 
        else :
            self.title_height = 0

        self.xaxis_label_height =0.
        if self.show_xaxis_label :
            self.xaxis_label_height += self.fontsize
        if self.show_fineprint :
            self.xaxis_label_height += self.small_fontsize

        self.line_height = (self.size.stack_height + self.line_margin_top +    
                            self.line_margin_bottom )
        self.line_width  = (self.size.stack_width*self.stacks_per_line + 
                            self.line_margin_left + self.line_margin_right )

        self.logo_height = int(2*self.logo_margin + self.title_height \
            + self.xaxis_label_height + self.line_height*self.lines_per_logo) 
        self.logo_width = int(2*self.logo_margin + self.line_width )


        self.creation_date = datetime.now().isoformat(' ') 

        end_type = '-'
        end_types = {
            unambiguous_protein_alphabet: 'p', 
            unambiguous_rna_alphabet: '-',
            unambiguous_dna_alphabet: 'd' 
        }
        if self.show_ends and self.alphabet in end_types:
            end_type = end_types[self.alphabet]
        self.end_type = end_type
    # End __init__
# End class LogoFormat



# ------ Logo Formaters ------
# Each formatter is a function f(LogoData, LogoFormat, output file).
# that draws a represntation of the logo into the given file.
# The main graphical formatter is eps_formatter. A mapping 'formatters'
# containing all available formatters is located after the formatter
# definitions. 

def pdf_formatter(data, format, fout) :
    """ Generate a logo in PDF format."""
    
    feps = StringIO()
    eps_formatter(data, format, feps)
    feps.seek(0)
    
    gs = GhostscriptAPI()    
    gs.convert('pdf', feps, fout, format.logo_width, format.logo_height)


def _bitmap_formatter(data, format, fout, device) :
    feps = StringIO()
    eps_formatter(data, format, feps)
    feps.seek(0)
    
    gs = GhostscriptAPI()    
    gs.convert(device, feps, fout, 
        format.logo_width, format.logo_height, format.resolution)


def jpeg_formatter(data, format, fout) : 
    """ Generate a logo in JPEG format."""
    _bitmap_formatter(data, format, fout, device="jpeg")


def png_formatter(data, format, fout) : 
    """ Generate a logo in PNG format."""
    _bitmap_formatter(data, format, fout, device="png")


def png_print_formatter(data, format, fout) : 
    """ Generate a logo in PNG format with print quality (600 DPI) resolution."""
    format.resolution = 600
    _bitmap_formatter(data, format, fout, device="png")


def txt_formatter( logodata, format, fout) :
    """ Create a text representation of the logo data. 
    """
    print >>fout, str(logodata)

    
def eps_formatter( logodata, format, fout) :
    """ Generate a logo in Encapsulated Postscript (EPS)"""
    
    subsitutions = {}
    from_format =[
        "creation_date",    "logo_width",           "logo_height",      
        "lines_per_logo",   "line_width",           "line_height",
        "line_margin_right","line_margin_left",     "line_margin_bottom",
        "line_margin_top",  "title_height",         "xaxis_label_height",
        "creator_text",     "logo_title",           "logo_margin",
        "stroke_width",     "tic_length",           
        "stacks_per_line",  "stack_margin",
        "yaxis_label",      "yaxis_tic_interval",   "yaxis_minor_tic_interval",
        "xaxis_label",      "xaxis_tic_interval",   "number_interval",
        "fineprint",        "shrink_fraction",      "errorbar_fraction",
        "errorbar_width_fraction",
        "errorbar_gray",    "small_fontsize",       "fontsize",
        "title_fontsize",   "number_fontsize",      "text_font",
        "logo_font",        "title_font",          
        "logo_label",       "yaxis_scale",          "end_type",
        "debug",            "show_title",           "show_xaxis",
        "show_xaxis_label", "show_yaxis",           "show_yaxis_label",
        "show_boxes",       "show_errorbars",       "show_fineprint",
        "rotate_numbers",   "show_ends",            
        ]
   
    for s in from_format :
        subsitutions[s] = getattr(format,s)


    from_format_size = ["stack_height", "stack_width"]
    for s in from_format_size :
        subsitutions[s] = getattr(format.size,s)

    subsitutions["shrink"] = str(format.show_boxes).lower()


    # --------- COLORS --------------
    def format_color(color):
        return  " ".join( ("[",str(color.red) , str(color.green), 
            str(color.blue), "]"))  

    subsitutions["default_color"] = format_color(format.default_color)

    colors = []  
    for group in format.color_scheme.groups :
        cf = format_color(group.color)
        for s in group.symbols :
            colors.append( "  ("+s+") " + cf )
    subsitutions["color_dict"] = "\n".join(colors)
        
    data = []
    
    # Unit conversion. 'None' for probability units
    conv_factor = std_units[format.unit_name]
    
    data.append("StartLine")
    

    seq_from = format.logo_start- format.first_index
    seq_to = format.logo_end - format.first_index +1

    # seq_index : zero based index into sequence data
    # logo_index : User visible coordinate, first_index based
    # stack_index : zero based index of visible stacks
    for seq_index in range(seq_from, seq_to) :
        logo_index = seq_index + format.first_index 
        stack_index = seq_index - seq_from
        
        if stack_index!=0 and (stack_index % format.stacks_per_line) ==0 :
            data.append("")
            data.append("EndLine")
            data.append("StartLine")
            data.append("")
        
        if logo_index % format.number_interval == 0 : 
            data.append("(%d) StartStack" % logo_index)
        else :            
            data.append("() StartStack" )

        if conv_factor: 
            stack_height = logodata.entropy[seq_index] * std_units[format.unit_name]
        else :
            stack_height = 1.0 # Probability

      #  if logodata.entropy_interval is not None and conv_factor:
            # Draw Error bars
       #     low, high = logodata.entropy_interval[seq_index]
       #     center = logodata.entropy[seq_index]


       #     down = (center - low) * conv_factor
       #     up   = (high - center) * conv_factor
       #     data.append(" %f %f %f DrawErrorbarFirst" % (down, up, stack_height) )
        
        s = zip(logodata.counts[seq_index], logodata.alphabet)
        def mycmp( c1, c2 ) :
            # Sort by frequency. If equal frequency then reverse alphabetic
            if c1[0] == c2[0] : return cmp(c2[1], c1[1])
            return cmp(c1[0], c2[0])
        
        s.sort(mycmp)

        C = float(sum(logodata.counts[seq_index])) 
        if C > 0.0 :
            fraction_width = 1.0
            if format.scale_width :
                fraction_width = logodata.weight[seq_index] 
            # print >>sys.stderr, fraction_width
            for c in s:
                data.append(" %f %f (%s) ShowSymbol" % (fraction_width, c[0]*stack_height/C, c[1]) )

        # Draw error bar on top of logo. Replaced by DrawErrorbarFirst above.
        if logodata.entropy_interval is not None and conv_factor:
            low, high = logodata.entropy_interval[seq_index]
            center = logodata.entropy[seq_index]

            down = (center - low) * conv_factor
            up   = (high - center) * conv_factor
            data.append(" %f %f DrawErrorbar" % (down, up) )
            
        data.append("EndStack")
        data.append("")
               
    data.append("EndLine")
    subsitutions["logo_data"] = "\n".join(data)  


    # Create and output logo
    template = resource_string( __name__, 'template.eps',)
    logo = Template(template).substitute(subsitutions)
    print >>fout, logo
 

# map between output format names and logo  
formatters = {
    'eps': eps_formatter, 
    'pdf': pdf_formatter,
    'png': png_formatter,
    'png_print' : png_print_formatter,
    'jpeg'  : jpeg_formatter,
    'txt' : txt_formatter,          
    }     
    
default_formatter = eps_formatter






    
def base_distribution(percentCG) :
    A = (1. - (percentCG/100.))/2.
    C = (percentCG/100.)/2.
    G = (percentCG/100.)/2.
    T = (1. - (percentCG/100))/2.
    return asarray((A,C,G,T), float64)   

def equiprobable_distribution( length) :
    return ones( (length), float64) /length   
  




#TODO: Move to seq_io? 
#   Would have to check that no symbol outside of full alphabet?
def which_alphabet(seqs) :
    """ Returns the most appropriate unambiguous protien, rna or dna alphabet
    for a Seq or SeqList.
    """
    alphabets = (unambiguous_protein_alphabet,
                unambiguous_rna_alphabet,
                unambiguous_dna_alphabet,
            )
    # Heuristic
    # Count occurances of each letter. Downweight longer alphabet.
    score = [1.0*asarray(seqs.tally(a)).sum()/sqrt(len(a)) for a in alphabets]
    #print score
    best = argmax(score) # Ties go to last in list.
    a = alphabets[best]
    return a

    
    
class LogoData(object) :
    """The data needed to generate a sequence logo.
       
    - alphabet 
    - length
    - counts  -- An array of character counts
    - entropy -- The relative entropy of each column
    - entropy_interval -- entropy confidence interval
     """
     
    def __init__(self, length=None, alphabet = None, counts =None, 
            entropy =None, entropy_interval = None, weight=None) :
        """Creates a new LogoData object"""
        self.length = length
        self.alphabet = alphabet
        self.counts = counts
        self.entropy = entropy
        self.entropy_interval = entropy_interval
        self.weight = weight
        
    
    #@classmethod    
    def from_counts(cls, alphabet, counts, prior= None):
        """Build a logodata object from counts."""
        seq_length, A = counts.shape
        
        if prior is not None: prior = array(prior, float64)
        
        if prior is None :
            R = log(A)
            ent = zeros(  seq_length, float64)
            entropy_interval = None    
            for i in range (0, seq_length) :
                C = sum(counts[i]) 
                #FIXME: fixup corebio.moremath.entropy()?
                if C == 0 :
                    ent[i] = 0.0
                else :
                    ent[i] = R - entropy(counts[i])
        else :
            ent = zeros(  seq_length, float64)
            entropy_interval = zeros( (seq_length,2) , float64)
        
            R = log(A)
            
            for i in range (0, seq_length) :
                alpha = array(counts[i] , float64)
                alpha += prior
                
                posterior = Dirichlet(alpha)
                ent[i] = posterior.mean_relative_entropy(prior/sum(prior)) 
                entropy_interval[i][0], entropy_interval[i][1] = \
                    posterior.interval_relative_entropy(prior/sum(prior), 0.95) 
 
        weight = array( na.sum(counts,axis=1) , float) 
        weight /= max(weight)
 
        return cls(seq_length, alphabet, counts, ent, entropy_interval, weight)
    from_counts = classmethod(from_counts)


    #@classmethod    
    def from_seqs(cls, seqs, prior= None):
        """Build a LogoData object form a SeqList, a list of sequences."""
        # --- VALIDATE DATA ---
        # check that at least one sequence of length at least 1 long
        if len(seqs)==0 or len(seqs[0]) ==0:
            raise ValueError("No sequence data found.")   
    
        # Check sequence lengths    
        seq_length = len(seqs[0])
        for i,s in enumerate(seqs) :
            #print i,s, len(s)
            if seq_length != len(s) :
                raise ArgumentError(
    "Sequence number %d differs in length from the previous sequences" % (i+1) ,
                        'sequences')

        # FIXME: Check seqs.alphabet?
        
        counts = asarray(seqs.tally())
        return cls.from_counts(seqs.alphabet, counts, prior)
    from_seqs = classmethod(from_seqs)

    def __str__(self) :
        out = StringIO()
        print >>out, '## LogoData'
        print >>out, '# First column is position number, couting from zero'
        print >>out, '# Subsequent columns are raw symbol counts'
        print >>out, '# Entropy is mean entropy measured in nats.' 
        print >>out, '# Low and High are the 95% confidence limits.'
        print >>out, '# Weight is the fraction of non-gap symbols in the column.'
        print >>out, '#\t'
        print >>out, '#\t',
        for a in self.alphabet :
            print >>out, a, '\t',
        print >>out, 'Entropy\tLow\tHigh\tWeight'
        
        for i in range(self.length) :
            print >>out, i, '\t',
            for c in self.counts[i] : print >>out, c, '\t',
            print >>out, self.entropy[i], '\t',
            if self.entropy_interval is not None:
                print >>out, self.entropy_interval[i][0], '\t', 
                print >>out, self.entropy_interval[i][1], '\t',
            else :
                print >>out, '\t','\t',
            if self.weight is not None :
                print >>out, self.weight[i],
            print >>out, ''
        print >>out, '# End LogoData'
        
        return out.getvalue()

##############################################################

class Dirichlet(object) :
    """The Dirichlet probability distribution. The Dirichlet is a continuous 
    multivariate probability distribution across non-negative unit length
    vectors. In other words, the Dirichlet is a probability distribution of 
    probability distributions. It is conjugate to the multinomial
    distribution and is widely used in Bayesian statistics.
    
    The Dirichlet probability distribution of order K-1 is 

     p(theta_1,...,theta_K) d theta_1 ... d theta_K = 
        (1/Z) prod_i=1,K theta_i^{alpha_i - 1} delta(1 -sum_i=1,K theta_i)

    The normalization factor Z can be expressed in terms of gamma functions:

      Z = {prod_i=1,K Gamma(alpha_i)} / {Gamma( sum_i=1,K alpha_i)}  

    The K constants, alpha_1,...,alpha_K, must be positive. The K parameters, 
    theta_1,...,theta_K are nonnegative and sum to 1.
    
    Status:
        Alpha
    """
    __slots__ = 'alpha', '_total', '_mean', 
    
    
    

    def __init__(self, alpha) :
        """
        Args:
            - alpha  -- The parameters of the Dirichlet prior distribution.
                        A vector of non-negative real numbers.  
        """
        # TODO: Check that alphas are positive
        #TODO : what if alpha's not one dimensional?
        self.alpha = asarray(alpha, float64)
        
        self._total = sum(alpha)
        self._mean = None
        
        
    def sample(self) :
        """Return a randomly generated probability vector.
        
        Random samples are generated by sampling K values from gamma
        distributions with parameters a=\alpha_i, b=1, and renormalizing. 
    
        Ref:
            A.M. Law, W.D. Kelton, Simulation Modeling and Analysis (1991).
        Authors:
            Gavin E. Crooks <gec@compbio.berkeley.edu> (2002)
        """
        alpha = self.alpha
        K = len(alpha)
        theta = zeros( (K,), float64)

        for k in range(K):
            theta[k] = random.gammavariate(alpha[k], 1.0) 
        theta /= sum(theta)

        return theta

    def mean(self) :
        if  self._mean ==None:
            self._mean = self.alpha / self._total
        return self._mean
    
    def covariance(self) :
        alpha = self.alpha
        A = sum(alpha)
        #A2 = A * A
        K = len(alpha)
        cv = zeros( (K,K), float64) 
        
        for i in range(K) :
            cv[i,i] = alpha[i] * (1. - alpha[i]/A) / (A * (A+1.) )
        
        for i in range(K) :
            for j in range(i+1,K) :
                v = - alpha[i] * alpha[j] / (A * A * (A+1.) )
                cv[i,j] = v
                cv[j,i] = v
        return cv
        
    def mean_x(self, x) :
        x = asarray(x, float64)
        if shape(x) != shape(self.alpha) :
            raise ValueError("Argument must be same dimension as Dirichlet")
        return sum( x * self.mean()) 

    def variance_x(self, x) :
        x = asarray(x, float64)
        if shape(x) != shape(self.alpha) :
            raise ValueError("Argument must be same dimension as Dirichlet")        
            
        cv = self.covariance()
        var = na.dot(na.dot(na.transpose( x), cv), x)
        return var


    def mean_entropy(self) :
        """Calculate the average entropy of probabilities sampled
        from this Dirichlet distribution. 
        
        Returns:
            The average entropy.
            
        Ref:
            Wolpert & Wolf, PRE 53:6841-6854 (1996) Theorem 7
            (Warning: this paper contains typos.)
        Status:
            Alpha
        Authors:
            GEC 2005
    
        """
        # TODO: Optimize
        alpha = self.alpha
        A = float(sum(alpha))
        ent = 0.0
        for a in alpha:
            if a>0 : ent += - 1.0 * a * digamma( 1.0+a) # FIXME: Check
        ent /= A
        ent += digamma(A+1.0)
        return ent    



    def variance_entropy(self):
        """Calculate the variance of the Dirichlet entropy. 

        Ref:
            Wolpert & Wolf, PRE 53:6841-6854 (1996) Theorem 8
            (Warning: this paper contains typos.)
        """
        alpha = self.alpha
        A = float(sum(alpha))
        A2 = A * (A+1)
        L = len(alpha)
        
        dg1 = zeros( (L) , float64)
        dg2 = zeros( (L) , float64)
        tg2 = zeros( (L) , float64)        
        
        for i in range(L) :
            dg1[i] = digamma(alpha[i] + 1.0)
            dg2[i] = digamma(alpha[i] + 2.0)
            tg2[i] = trigamma(alpha[i] + 2.0)

        dg_Ap2 = digamma( A+2. )
        tg_Ap2 = trigamma( A+2. )
        
        mean = self.mean_entropy()
        var = 0.0
    
        for i in range(L) :
            for j in range(L) :
                if i != j :
                    var += (
                        ( dg1[i] - dg_Ap2 ) * (dg1[j] - dg_Ap2 ) - tg_Ap2 
                        ) * (alpha[i] * alpha[j] ) / A2
                else : 
                    var += (
                        ( dg2[i] - dg_Ap2 ) **2 + ( tg2[i] - tg_Ap2 )
                        ) * ( alpha[i] * (alpha[i]+1.) ) / A2

        var -= mean**2
        return var
        
        
        
    def mean_relative_entropy(self, pvec) :
        ln_p = na.log(pvec)
        return - self.mean_x(ln_p) - self.mean_entropy() 
    
    
    def variance_relative_entropy(self, pvec) :
        ln_p = na.log(pvec)
        return self.variance_x(ln_p) + self.variance_entropy()
    
    
    def interval_relative_entropy(self, pvec, frac) :
        mean = self.mean_relative_entropy(pvec) 
        variance = self.variance_relative_entropy(pvec) 
        # If the variance is small, use the standard 95% 
        # confidence interval: mean +/- 1.96 * sd
        if variance< 0.1 :
            sd = sqrt(variance)
            return max(0.0, mean - sd*1.96), mean + sd*1.96
        
        g = Gamma.from_mean_variance(mean, variance)
        low_limit = g.inverse_cdf( (1.-frac)/2.)
        high_limit = g.inverse_cdf( 1. - (1.-frac)/2. )
        
        return low_limit, high_limit
