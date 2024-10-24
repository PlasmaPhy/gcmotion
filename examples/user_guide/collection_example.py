import gcmotion as gcm
import numpy as np
import parameters

collection = gcm.Collection(parameters)
collection.run_all(orbit=True, terminal=5)

# gcm.collection_drift(
#     collection, angle="theta", lim=[-np.pi, np.pi], plot_initial=True
# )

gcm.collection_energy_contour(
    collection,
    theta_lim=[-np.pi, np.pi],
    psi_lim="auto",
    plot_drift=True,
    contour_Phi=True,
    units="keV",
    levels=30,
    plot_initial=True,
    different_colors=False,
)

# gcm.collection_drifts(collection)
gcm.collection_poloidal_cut(collection)
gcm.collection_parabolas(collection, plot_label=False, autoscale=True)
