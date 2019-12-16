from .basic import GeometricWorld

class GasGiant(GeometricWorld):

    world_class = 'gasgiant'

    def __repr__(self):

        if 'name' in self.__dict__:
            if self.name is not None:
                return f'{self.name} {self.__class__} object at {hex(id(self))}'
        return f'{self.__class__} object at {hex(id(self))}'