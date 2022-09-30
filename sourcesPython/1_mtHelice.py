
import math
from time import time
from Vector import *
import random
from Helice import *
f = open("animations/1_mtHelice.pv","w")

time_step = 0.1
nb_points = 1000
nb_frames = 600
R = 3 
omega = 1
pente = 5
tous_les_points = []
h = Helice(R,1,1)
h0 = h.courbe(0)
for i in range(nb_points):
    rayon = random.random()*R
    angle = random.random()*2*math.pi
    x=rayon*math.cos(angle)
    z=rayon*math.sin(angle)
    tous_les_points.append(Vector(x,0,z))


# ecriture de l'entete
f.write("#PV==\n" + str(time_step) + "\n" + str(nb_points) + "\n" + str(nb_frames) + "\n")
f.write("Point 0:"+str(nb_points-1)+"\n")
f.write("====\n")


#-------------- animation
# parcours des frames (images)
position_depart = Vector(0,0,0)
for frame in range(nb_frames):
    for i in range(nb_points):
        (x,y,z)=tous_les_points[i].get()
        f.write(str(x) + " " + str(y) + " " + str(z) + "\n")
        vitesse = h.tang_a_la_hauteur(y)
        distance_au_point = h.distance_avec_point(tous_les_points[i])
        if distance_au_point == R:
            tous_les_points[i] = tous_les_points[i] + 0*time_step
        if distance_au_point == 0:
            tous_les_points[i] = tous_les_points[i] + h.tang_a_la_hauteur(y)*time_step
        else:
            ## manque de temps
f.close()
    