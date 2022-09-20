
import math
from Vector import *

f = open("animations/01_vagues.pv","w")

time_step = 0.1
nb_points = 100
nb_frames = 150

# place 10 lignes et 10 colonnes de points
tous_les_points=[]
for lgn in range(10):
    for col in range(10):
        tous_les_points += [Vector(col-5,0,lgn-5)]

# ecriture de l'entete
f.write("#PV==\n" + str(time_step) + "\n" + str(nb_points) + "\n" + str(nb_frames) + "\n")
for lgn in range(10):
    f.write("oLine "+str(10*lgn)+":"+str(10*(lgn+1)-1)+"\n")
f.write("====\n")


#-------------- animation
pulsation=0.15
amplitude=0.3

# parcours des frames (images)
for fr in range(nb_frames):
    # parcours des points
    for n in range(nb_points):
        (x,y,z)=tous_les_points[n].get()
        y=amplitude*math.sin(pulsation*fr - 0.5*(x-z))
        f.write(str(x) + " " + str(y) + " " + str(z) + "\n")

f.close()
