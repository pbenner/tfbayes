# Copyright (C) 2013 Philipp Benner
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

from Bio.AlignIO import *
from Bio.AlignIO import _FormatToIterator
from Bio.AlignIO import _FormatToWriter

# implant MafIO into Bio package
import MafIO

# register MafIO
_FormatToIterator["maf"] = MafIO.MafIterator
_FormatToWriter  ["maf"] = MafIO.MafWriter