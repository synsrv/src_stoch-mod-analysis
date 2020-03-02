
import matlab, matlab.engine

import numpy as np

eng = matlab.engine.start_matlab()

eng.addpath('/home/hoffmann/lab/stoch-mod/analysis-dev/m/',nargout=0)

# eng.triarea(nargout=0)

alpha, bmin, L = eng.bplfit(matlab.double([900, 90, 9]),
                            matlab.double([1, 10, 100, 1000]),
                            'range',
                            matlab.double(list(np.arange(1.5,2.5,0.1))),
                            nargout=3)
                            
print(alpha, bmin, L)
                            
