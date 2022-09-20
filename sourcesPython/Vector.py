import math


class Vector:
    def __init__(self,x=0,y=0,z=0):
        self.x=x
        self.y=y
        self.z=z

    def set(self,x,y,z):
        self.x=x
        self.y=y
        self.z=z

    def get(self):
        return (self.x,self.y,self.z)

    def show(self,label):
        print(label + " : (" + str(self.x) + "," + str(self.y) + "," + str(self.z) + ")")

    def __add__(self,autre):
        return Vector(self.x+autre.x,self.y+autre.y,self.z+autre.z)

    def __sub__(self,autre):
        return Vector(self.x-autre.x,self.y-autre.y,self.z-autre.z)

    def __mul__(self,nombre):
        return Vector(self.x*nombre,self.y*nombre,self.z*nombre)

    def __pow__(self,autre):
        return self.x*autre.x+self.y*autre.y+self.z*autre.z

    def module(self):
        return math.sqrt(self**self)

    def distanceAvec(self,autre):
        return (self-autre).module()

    def unit(self):
        norm = self.module()
        if norm!=0:
            return self*(1/norm)
        else:
            return Vector(1,0,0)

    def normalize(self):
        norm = self.module()
        if norm!=0:
            self = self * (1/norm)
        else:
            self = Vector(1,0,0)

if __name__ == "__main__":
    p = Vector(1,2,3)
    q = Vector(4,5,6)
    s = p - q
    s.show("diff")
