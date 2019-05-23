from send2trash import send2trash
from ..
def delete_planets(ask_prompt: bool = True):
    """ Will delete all planets in the dilled_planets folder """

    if ask_prompt:
        response = input("You are about to delete all planets from TidalPy's dilled planet folder!\n"
                         "    These files make loading planets much faster for TidalPy. Only delete as a last resort "
                         "troubleshooting process or if planets have become outdated.\n"
                         "Proceed with deletion? [y/n]: ")
        if response.lower().strip() in ['y', 'yes']:
            pass
        else:
            return False

    for file in

