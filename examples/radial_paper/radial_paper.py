import gcmotion as gcm
import radial_paper_parameters as parameters

collection = gcm.Collection(parameters)
collection.run_all(orbit=True, terminal=5)

# gcm.tokamak_profile(collection.particles[0])

gcm.collection_energy_contour(
    collection,
    psi_lim="auto",
    plot_drift=True,
    contour_Phi=True,
    units="keV",
    levels=30,
    different_colors=False,
    adjust_cbar=False,
)

gcm.collection_poloidal_cut(collection)
