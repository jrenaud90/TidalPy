from TidalPy.Extending.burnman import burnman_installed, Mineral


class Pyrite(Mineral):

    def __init__(self):

        if not burnman_installed:
            raise ImportError('Burnman package not found.')

        """ Parameters from Thompson et al, 2016 in American Mineralogist Vol 101, Page 1046"""
        self.params = {
            'formula'          : {'Fe': 1., 'Si': 2.},
            'equation_of_state': 'bm3',
            'K_0'              : 140.2e9,
            'Kprime_0'         : 5.52,
            'grueneisen_0'     : 1.4,
            'Debye_0'          : 624.,
            'V_0'              : 2.393e-5,
            'q_0'              : 2.06,
            'molar_mass'       : 119.9750 / 1000.,
            'n'                : 3,
            # Shear information is from Whitaker & Wang 2010 Journal of Earth Science Vol 21 No 5 p. 792
            'G_0'              : 112.3e9,  # From Whitaker & Wang 2010 Journal of Earth Science Vol 21 No 5 p. 792
            'Gprime_0'         : 3.0,
            # Merkel et al (2002) in Physics and Chemistry of Minerals Vol 29 page 1) had slightly higher value (126e9)
            }

        Mineral.__init__(self)
