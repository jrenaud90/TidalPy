import TidalPy

planet = TidalPy.build_planet('pluto_example', force_build=True, auto_save=False)
planet.paint(auto_show=True)