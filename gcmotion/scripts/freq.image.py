import numpy as np
import matplotlib.pyplot as plt

from skimage import measure
from matplotlib import transforms


from gcmotion.entities.profile import Profile
from gcmotion.plot._base._base_profile_contour import _base_profile_contour

from gcmotion.configuration.scripts_configuration import (
    FrequencyConfig as config,
)


def frequency(profile: Profile, **args):

    # Unpack parameters
    _thetalim = [-2 * np.pi, 2 * np.pi]
    args["psilim"] = args.get("psilim", config.psilim)
    args["flux_units"] = args.get("flux_units", config.flux_units)
    args["E_units"] = args.get("E_units", config.E_units)
    args["levels"] = args.get("levels", config.levels)
    args["grid_density"] = args.get("grid_density", config.grid_density)
    args["potential"] = args.get("potential", config.potential)

    thetalim = profile.Q(_thetalim, "radians")
    psilim = profile.Q(args["psilim"], "psi_wall").to(args["flux_units"])
    psigrid, thetagrid = np.meshgrid(
        np.linspace(psilim[0], psilim[1], args["grid_density"]),
        np.linspace(thetalim[0], thetalim[1], args["grid_density"]),
    )

    # Calculate Energy values
    Energy = profile.findEnergy(
        psigrid, thetagrid.m, args["E_units"], potential=args["potential"]
    )

    data = {
        "theta": thetagrid.m,
        "psi": psigrid.m,
        "Energy": Energy.m,
    }

    # Calculate Energy span
    Energyspan = [np.min(data["Energy"]), np.max(data["Energy"])]
    Energylevels = np.linspace(Energyspan[0], Energyspan[1], args["levels"])

    # Find all contours
    contour_groups = []
    for E in Energylevels:
        contour_groups.append(measure.find_contours(data["Energy"], E))

    # Display the image and plot all contours found
    fig, ax = plt.subplots()
    base = ax.transData
    rot = transforms.Affine2D().rotate_deg(90)
    trans = transforms.Affine2D().translate(args["grid_density"], 0)

    # ax.imshow(Energy.m, cmap="plasma", transform=rot + trans + base)
    _base_profile_contour(profile, ax)

    # Find contours at a constant value of 0.8
    for group in contour_groups:
        for contour in group:
            ax.plot(
                contour[:, 1],
                contour[:, 0],
                linewidth=2,
                transform=rot + trans + base,
            )
    ax.set_xlim([0, args["grid_density"]])
    ax.set_ylim([0, args["grid_density"]])
    plt.show()
    return contour_groups
