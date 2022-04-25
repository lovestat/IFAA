# Note .XXXXXX means "import XXXXXX.py from current directory"

from .comb import IFAA
from .loadData import load_dataM, load_dataC

__version__ = '0.0.1'
__author__ = 'Shangchen Song'
__all__ = ['IFAA', 'load_dataM', 'load_dataC']