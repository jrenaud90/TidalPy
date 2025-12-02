import os
from sphinx.ext.apidoc import main

package = os.path.join(os.pardir, os.pardir, "TidalPy")
output = "generated"

if not os.path.exists(output):
    os.makedirs(output)

main(["-f", "-o", output, package])
