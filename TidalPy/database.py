"""The pack database: the material, world, and system TOML files, read into memory when TidalPy is imported.

TidalPy's two packs (:mod:`TidalPy.Utilities.data_pack`) each keep their TOML files in memory: the MatPack's materials
(``<documents>/TidalPy/<version>/Materials``) and the WorldPack's worlds and systems (``.../Worlds``), the data
directory's copies preferred to the packaged files as on disk. :func:`load_database` installs both packs and reads
them; it runs when TidalPy is imported and again on ``TidalPy.reinit()``.

A material, world, or system named in a build then comes from memory, after a check that its file is unchanged on
disk (its modification time and size). A name the database does not hold, or whose file changed or disappeared since
it was read, is searched for on disk as before, and the database keeps what is found, so an edit made during a
session is picked up by the next build. A file added to a data directory during a session is found the first time it
is named.
"""

import os

import TidalPy
from TidalPy.paths import get_materials_dir, get_worlds_dir
from TidalPy.Utilities.data_pack import DataPack

# The packaged MatPack and WorldPack directories (read-only sources of the packs), relative to the package root.
PACKAGED_MATPACK_DIR = os.path.join(os.path.dirname(os.path.abspath(TidalPy.__file__)), "MatPack")
PACKAGED_WORLDPACK_DIR = os.path.join(os.path.dirname(os.path.abspath(TidalPy.__file__)), "WorldPack")

# File extensions installed into the user worlds directory: world and system TOMLs and their companion data files
# (e.g. PREM-like radial profiles). The database holds the TOMLs.
WORLDPACK_EXTENSIONS = (".toml", ".csv", ".txt", ".dat")

# The packs. TidalPy.Material.matpack and TidalPy.Structures.configs.worldpack point each pack's directory getter at
# their own get_materials_dir and get_worlds_dir, so a test that patches those redirects the pack.
MAT_PACK = DataPack(
    "MatPack",
    PACKAGED_MATPACK_DIR,
    get_materials_dir,
    (".toml",),
    "stale_matpack_copy",
    "TidalPy.Material.install_matpack(force=True)")

WORLD_PACK = DataPack(
    "WorldPack",
    PACKAGED_WORLDPACK_DIR,
    get_worlds_dir,
    WORLDPACK_EXTENSIONS,
    "stale_worldpack_copy",
    "TidalPy.Structures.install_worldpack(force=True)")

PACKS = (MAT_PACK, WORLD_PACK)


def load_database() -> None:
    """Install the MatPack and WorldPack and read every TOML file of both into memory.

    Called when TidalPy is imported and by ``TidalPy.reinit()``. A lookup notices a file edited, added, or deleted
    since, so calling it again is never needed to see a change. A pack that cannot be read is left for its first
    lookup to read, so importing TidalPy never fails here.
    """
    for pack in PACKS:
        try:
            pack.load()
        except OSError:
            pass
