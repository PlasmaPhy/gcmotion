import gcmotion as gcm
import guide_parameters as parameters

# --------------------Collection Setup--------------------
collection = gcm.Collection(parameters)
collection.run_all(orbit=True)
# -------------------------Plots-------------------------

gcm.collection_energy_contour(
    collection,
    Ptheta_lim=[0, 0.4],
    plot_drift=True,
    contour_Phi=True,
    units="SI",
)
