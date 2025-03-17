import numpy as np
import pyvista as pv


def calculate_error_statistics(mesh: pv.UnstructuredGrid):
    """Calculate the error between analytical and coupled solution. The absolute
    error is stored as an additional data array on the mesh points.
    Additionally, the minimum, maximum, first quartile, third quartile, median,
    and mean L2 norm of the error is stored as field data on the mesh."""

    mesh["location_error"] = mesh.points - mesh["analytical_location"]
    location_error_norm = np.linalg.norm(mesh["location_error"], axis=1)
    mesh.field_data["location_error_L2_min"] = np.amin(location_error_norm)
    mesh.field_data["location_error_L2_max"] = np.amax(location_error_norm)
    mesh.field_data["location_error_L2_q1"] = np.percentile(location_error_norm, 25)
    mesh.field_data["location_error_L2_q3"] = np.percentile(location_error_norm, 75)
    mesh.field_data["location_error_L2_median"] = np.median(location_error_norm)
    mesh.field_data["location_error_L2_mean"] = np.mean(location_error_norm)

    mesh["velocity_error"] = mesh["velocity"] - mesh["analytical_velocity"]
    velocity_error_norm = np.linalg.norm(mesh["velocity_error"], axis=1)
    mesh.field_data["velocity_error_L2_min"] = np.amin(velocity_error_norm)
    mesh.field_data["velocity_error_L2_max"] = np.amax(velocity_error_norm)
    mesh.field_data["velocity_error_L2_q1"] = np.percentile(velocity_error_norm, 25)
    mesh.field_data["velocity_error_L2_q3"] = np.percentile(velocity_error_norm, 75)
    mesh.field_data["velocity_error_L2_median"] = np.median(velocity_error_norm)
    mesh.field_data["velocity_error_L2_mean"] = np.mean(velocity_error_norm)
