import gcmotion as gcm
import numpy as np
import guide_parameters as parameters

# --------------------Collection Setup--------------------
collection = gcm.Collection(parameters)
collection.run_all(orbit=True, terminal=8)

# -------------------------Plots-------------------------
# gcm.collection_drift(
#     collection, angle="theta", lim=[-np.pi, np.pi], plot_initial=True
# )

# gcm.collection_drifts(collection)

# gcm.collection_energy_contour(
#     collection,
#     psi_lim="auto",
#     plot_drift=True,
#     contour_Phi=True,
#     units="keV",
#     levels=30,
#     plot_initial=True,
#     different_colors=False,
#     adjust_cbar=False,
# )


gcm.collection_poloidal_cut(collection)
# gcm.collection_parabolas(collection, plot_label=False, autoscale=True)
