
def high_pressure_silicate_from_abundance(mg_number: float, mg_si_ratio: float):
    """ Determine the material composition for high pressure silicate based on Mg concentrations.

    Methods based on Sotin+2007

    Parameters
    ----------
    mg_number : float
        Magnesium to Iron number: Mg / (Mg + Fe)
    mg_si_ratio : float
        Magnesium to Silicon ratio: Mg / Si

    Returns
    -------

    """

    iron_number = 1. - mg_number
    #TODO