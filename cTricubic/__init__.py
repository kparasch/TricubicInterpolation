import sys
if sys.version_info[0] < 3:
    from .Tricubic2_c import *
else:
    from .Tricubic_c import *
