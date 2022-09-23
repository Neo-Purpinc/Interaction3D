
import math
from Vector import *

f = open("../animations/01_battement.pv","w")

time_step = 0.1
nb_points = 4
nb_frames = 100

# placement de tous les points sur un carré
tous_les_points=[]
tous_les_points += [Vector(-1,0,-1)]
tous_les_points += [Vector( 1,0,-1)]
tous_les_points += [Vector( 1,0, 1)]
tous_les_points += [Vector(-1,0, 1)]

# ecriture de l'entete
f.write("#PV==\n" + str(time_step) + "\n" + str(nb_points) + "\n" + str(nb_frames) + "\n")
f.write("cLine 0:3\n")
f.write("====\n")

#-------------- animation
pulsation=0.5
amplitude=0.3

# parcours des frames (images)
for fr in range(nb_frames):
    # parcours des points
    for n in range(nb_points):
        (x,y,z)=tous_les_points[n].get()
        y=amplitude*math.sin(pulsation*fr - 0.5*(x-z))
        f.write(str(x) + " " + str(y) + " " + str(z) + "\n")

f.close()
