import observatory 

f = observatory.Field()
slow = False 
f.initialize()

for filename in f.list_of_files:
    f.openFits(filename,calibrated=False)
    f.Reduce()
    if f.calibrated():
        f.openFits(filename,calibrated=True)
        f.Extract()



