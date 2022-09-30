
import math
from Vector import *
import random

f = open("animations/03_tourbillon.pv","w")

time_step = 0.1
nb_points = 600
nb_frames = 600
max_dist = 3
tous_les_points = []  

for i in range(nb_points):
    rayon = random.random()*max_dist
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
        # update position inversement proportionnelle au rayon
        rayon = tous_les_points[i].module()
        angle = math.atan2(z,x) # angle entre le point et l'axe des x
        angle += 0.05/rayon
        x=rayon*math.cos(angle)
        z=rayon*math.sin(angle)
        tous_les_points[i].set(x,y,z)



f.close()
    