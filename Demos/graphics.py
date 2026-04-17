import numpy as np
import matplotlib.pyplot as plt
import cmcrameri.cm as cmc
import cartopy.crs as ccrs

from TidalPy.utilities.graphics import projection_map

def tidal_heating_layers_plot(planet, arrays, radius_array, tidal_heating):
    colatitude, longitude, time, layer_indices = arrays

    for layer_i, layer in enumerate(planet.layers):

        # Create figure with 2 subplots
        fig = plt.figure(figsize=(10, 4))
        gs = fig.add_gridspec(1, 4, width_ratios=[55, 2, 15, (100-55-2-15)], wspace=0.2)
        ax_left = fig.add_subplot(gs[0, 0], projection=ccrs.Mollweide())
        cax = fig.add_subplot(gs[0, 1])
        # buffer_ax = fig.add_subplot(gs[0, 2])
        ax_right = fig.add_subplot(gs[0, 3])

        # Sum tidal heat in this layer's thickness
        layer_index = layer_indices[layer_i]
        tidal_heat_depth_sum = np.nansum(tidal_heating[layer_index, :, :], axis=0)

        # Sum heat across lat and long
        tidal_heat_latlong_sum = np.nansum(tidal_heating, axis=(1,2))

        projection_map(np.degrees(longitude), np.degrees(colatitude), tidal_heat_depth_sum,
                       cmap=cmc.lajolla,
                       projection="Mollweide",
                       zlabel='Tidal Heating [W]', zlog=False,
                       ax=ax_left, cbax=cax)
        # Plot other depths as black lines
        ax_right.plot(tidal_heat_latlong_sum[~layer_index], radius_array[~layer_index]/1000, c='k', linewidth=1.0)

        # Then plot this layers with a thicker line and different color
        ax_right.plot(tidal_heat_latlong_sum[layer_index], radius_array[layer_index]/1000, c='r', linewidth=3.0)
        
        ax_right.set_xlabel("Tidal Heating [W]")
        ax_right.set_xscale("log")
        ax_right.set_ylabel("Radius [km]")

        fig.suptitle(f"{planet.name} - {layer.name}", fontsize=16, color='k')

        plt.show()

