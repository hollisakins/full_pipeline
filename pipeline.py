import sys

vers = '%s.%s' % (sys.version_info[0],sys.version_info[1])
if not vers=='2.7':
    raise Exception("Must be using Python 2.7")

import observatory

observatory.makeMasters()

f = observatory.Field()
f.initialize()

f.aperture_size = 100 ## pixel radius for aperture photometry, should be around 80-100 for bin1x1 but less for 2x2 or 3x3

for filename in f.list_of_files:
    f.openFits(filename,calibrated=False)
    f.Reduce()
    if f.calibrated:
        f.openFits(filename,calibrated=True)
        f.Extract()





