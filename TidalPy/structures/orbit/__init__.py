from .base import OrbitBase
from .physics import PhysicsOrbit

# User will almost always want the Physics version of the orbit, so it is aliased as just `Orbit` here
Orbit = PhysicsOrbit
