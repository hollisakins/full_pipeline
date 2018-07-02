import observatory 

observatory.makeMasters('ArchCal/')

f = observatory.Field()
f.initialize()

f.aperture_size = 80.0 ## pixel radius for aperture photometry, should be around 80-100 for bin1x1

for filename in f.list_of_files:
    f.openFits(filename,calibrated=False)
    f.Reduce()
    if f.calibrated:
        f.openFits(filename,calibrated=True)
        f.Extract()





