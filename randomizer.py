#mass, x, y, z(cartiensian coordinates), u, v, w (velocity vector)
import random
import sys

numBodies = int(sys.argv[1])

fileName = "nBodyInput" + str(numBodies) + ".txt"

#writing here so every time we run this script we
#have different random numbers
f = open(fileName, "w")

f.write(str(numBodies) + "\n")

#Our input file is in the format 
#number of bodies, followed by
#each body's mass, x, y, and z coordinates
#and u, v, and w velocity vector,
#all of which is seperated by a newline.
#This script generates random numbers so that
#we perform simulations on unique data.


f.write(str(1.))
f.write("\n")
for i in range (6):
    f.write(str(0.))
    f.write("\n")
    

for i in range(numBodies - 1):

    f.write(str(random.uniform(1., 1000000000.)))
    f.write("\n")
    for i in range(6):
        f.write(str(random.uniform(-1000000000., 1000000000.)))
        f.write("\n")

f.close()
