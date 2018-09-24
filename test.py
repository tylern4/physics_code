import time

import handeler as dh

start = time.time()
files = [
    "/Users/tylern/Data/e1d/h10_22904.root",
    "/Users/tylern/Data/e1d/h10_22905.root",
    "/Users/tylern/Data/e1d/h10_22906.root",
    "/Users/tylern/Data/e1d/h10_22907.root",
    "/Users/tylern/Data/e1d/h10_22908.root"
]

data = dh.handeler(files, "new_test.root")

data.run()

end = time.time()
print("{} Sec".format(end - start))
