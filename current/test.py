import glob
import sys
import time

import handeler as dh

start = time.time()
files = glob.glob(sys.argv[1] + "*.root")
data = dh.handeler(files, "new_test.root")
data.run()

end = time.time()
print("{} Sec".format(end - start))
