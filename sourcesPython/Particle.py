# define class Particle
import Vector

class Particle:
    def __init__(self, m,p:Vector,v:Vector=Vector(0,0,0)):
        self.masse = m
        self.position = p
        self.vitesse = v
        self.forces = Vector(0, 0, 0)
    
    def addForce(self, f:Vector):
        self.forces += f
    
    def updateSpeed(self, dt):
        self.vitesse += self.forces * dt / self.masse
    
    def updatePosition(self, dt):
        self.position += self.vitesse * dt + self.forces * dt**2 / self.mass
    
