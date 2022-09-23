
from Vector import *

f = open("../animations/00_noanim.pv","w")

time_step = 1.0
nb_points = 4
nb_frames = 2

tous_les_points = []
tous_les_points += [Vector(-1,0,-1)]
tous_les_points += [Vector( 1,0,-1)]
tous_les_points += [Vector( 1,0, 1)]
tous_les_points += [Vector(-1,0, 1)]

f.write("#PV==\n" + str(time_step) + "\n" + str(nb_points) + "\n" + str(nb_frames) + "\n")
f.write("oLine 0:2\n")
f.write("Sphere 2:2 0.5\n")
f.write("Point 2:3\n")
f.write("====\n")

for fr in range(nb_frames):
    for n in range(nb_points):
        (x,y,z)=tous_les_points[n].get()
        f.write(str(x) + " " + str(y) + " " + str(z) + "\n")

f.close()
