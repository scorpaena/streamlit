# import tempfile
# import os
# from functools import cache, cached_property
# from cadquery import Workplane, Edge, Wire, Vector, exporters, Sketch
# from math import radians, cos, tan, acos, degrees, sin
# import pyvista as pv
# from cadquery.occ_impl.shapes import TOLERANCE
# from cadquery.vis import show_object
# from stpyvista import stpyvista as stv
# import streamlit as st
# from splinecloud_scipy import load_spline, load_subset
#
# from OCP.GProp import GProp_GProps
# from OCP.BRepGProp import BRepGProp
# from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeEdge
# from OCP.gp import gp_Pnt
# from OCP.TColgp import TColgp_Array1OfPnt
# from OCP.TColStd import TColStd_Array1OfInteger, TColStd_Array1OfReal
# from OCP.Geom import Geom_BSplineCurve
#
# from scipy.interpolate import splprep, UnivariateSpline, interp1d
# import numpy as np
# import pandas as pd
#
#
# def build_occ_spline(knot_vector, control_points, degree):
#     ## assume there are no repeating knots inside
#     knot_mults = list(map(lambda k: degree + 1 if k in (0, 1) else 1, knot_vector))
#
#     poles = TColgp_Array1OfPnt(1, len(control_points))
#     for i, cp in enumerate(control_points):
#         poles.SetValue(i + 1, gp_Pnt(cp[0], cp[1], 0))
#
#     knots = TColStd_Array1OfReal(1, len(knot_vector))
#     for i, k in enumerate(knot_vector):
#         knots.SetValue(i + 1, k)
#
#     mults = TColStd_Array1OfInteger(1, len(knot_mults))
#     for i, m in enumerate(knot_mults):
#         mults.SetValue(i + 1, m)
#
#     print(">>>>>", poles, knots, mults, degree)
#     bspl = Geom_BSplineCurve(poles, knots, mults, degree)
#
#     return bspl
#
#
# class Airfoil:
#     def __init__(self, profile_data):
#         self.profile_data = profile_data
#         self.profile = load_subset(profile_data)
#
#
#
# class FanGeometry:
#     """Class to handle the generation of fan geometry."""
#
#     def __init__(self, airfoil, chord=1000, thicken_sharp_tail=True):
#         self.airfoil = airfoil
#         self.chord = chord
#         self.thicken_sharp_tail = thicken_sharp_tail
#
#         self._setup_profile(self.thicken_sharp_tail)
#         self.sketch = self.build_sketch()
#         self.cm = self.get_profile_center_of_mass()
#
#     @staticmethod
#     @cache
#     def get_spline(curve_id):
#         return load_spline(curve_id)
#
#     def get_profile_center_of_mass(self):
#         profile_face = Workplane().placeSketch(self.sketch.copy()).toOCC()
#         return profile_face.Center()
#
#     @staticmethod
#     def _has_sharp_tail(profile_points):
#         return np.array_equal(profile_points[0], profile_points[-1])
#
#     @staticmethod
#     def _patch_sharp_tail(profile_points):
#         """
#         Patch tail points to avoid malformed geometry construction
#         and enable shell creation in CadQuery
#         """
#         thicken_ratio = 0.3  # share of neighboring point offsets
#
#         ## height at the first point pair after the tail point
#         h1 = profile_points[1][1] - profile_points[-2][1]
#
#         h0 = min(h1 * thicken_ratio, 0.001)
#         profile_points[0][1] = profile_points[0][1] + h0 / 2
#         profile_points[-1][1] = profile_points[-1][1] - h0 / 2
#
#     def _setup_profile(self, thicken_sharp_tail):
#         # profile_points = self.airfoil.profile.to_numpy()
#         profile_points = pd.DataFrame.from_dict(self.airfoil.profile)
#         print("%%%%%%%", profile_points)
#
#         self.sharp_tail = self._has_sharp_tail(profile_points)
#         if self.sharp_tail and thicken_sharp_tail:
#             self._patch_sharp_tail(profile_points)
#             self.sharp_tail = False
#
#         t, c, k = self._bspl_profile_approx(profile_points)
#         self.profile_tck = t, c*self.chord, k
#
#         self.profile_points = [tuple(p) for p in profile_points*self.chord]
#         self._setup_profile_curves()
#
#     def _setup_profile_curves(self):
#         """
#         Inerpolate upper and butom profile countours with spline functions
#         to enable geometry analysis (these are not used in geometry consruction).
#         """
#         profile_points = np.array(self.profile_points)
#
#         nose_pnt_ind = None
#         min_x_val = np.min(profile_points[:, 0])
#         for i, p in enumerate(profile_points):
#             if p[0] == min_x_val:
#                 nose_pnt_ind = i
#                 break
#
#         upper_points = np.flip(profile_points[:nose_pnt_ind+1], axis=0)
#         x_top, y_top = upper_points[:, 0], upper_points[:, 1]
#         self.spl_top = UnivariateSpline(x_top, y_top, s=0.0)
#
#         bottom_points = np.array(profile_points[nose_pnt_ind:])
#         x_bot, y_bot = bottom_points[:, 0], bottom_points[:, 1]
#         self.spl_bottom = UnivariateSpline(x_bot, y_bot, s=0.0)
#
#     @staticmethod
#     def _bspl_profile_approx(profile_points, s=2e-5, k=3):
#         """
#         B-Spline approximation of the discrete profile data, represented in Selig format
#
#         Paameters:
#         'pp' - array of [x,y] coords
#         's' - smoothing factor
#         'k' - spline degree
#
#         Returns:
#         't', 'cx',  - B-Spline represntation, instances of the SplineCloud ParametricUnivariateSpline
#         """
#         xvals = [p[0] for p in profile_points]
#         yvals = [p[1] for p in profile_points]
#
#         ## set weights to critical profile points
#         # weight_mod = lambda x: 20 if x == 0.0 else (5 if x == 1.0 else 1)
#
#         # numeric_xvals = [x for x in xvals if isinstance(x, (int, float)) or (isinstance(x, str) and x.replace('.', '', 1).isdigit())]
#         # weight_mod = lambda x: 20 if np.isclose(float(x), 0.0) else (5 if np.isclose(float(x), 1.0) else 1)
#         # weights = list(map(weight_mod, xvals))
#         # weights = list(map(weight_mod, numeric_xvals))
#
#         print(">>>>>>>", xvals, yvals)
#
#         tck, u = splprep([xvals, yvals], s=s, k=k)#, w=weights)
#         t, c, k = tck
#         cx, cy = c
#
#         ## adjust spline start and end points to match corrsponding profile points
#         cx[0], cy[0] = profile_points[0]
#         cx[-1], cy[-1] = profile_points[-1]
#
#         print("number of profile data points:", len(profile_points))
#         print("number of brep control points:", len(cx))
#
#         return t[k:-k], np.array(list(zip(cx, cy))), k
#
#     def build_sketch(self, smooth_profile=True):
#         sketch = Sketch().spline(self.profile_points)
#         if smooth_profile:
#             bspl_occ = build_occ_spline(*self.profile_tck)
#             profile_edge = Edge(BRepBuilderAPI_MakeEdge(bspl_occ).Edge())
#             sketch = Sketch().edge(profile_edge)
#
#         if self.sharp_tail:
#             return sketch.assemble()
#         return sketch.close().assemble()
#
#
#
#
#
#
#     #
#     #
#     # def involute_curve(self, sign=1):
#     #     """Returns a function for generating involute curve points."""
#     #
#     #     def curve_function(radius):
#     #         alpha = sign * acos(self.base_radius / radius)
#     #         x = radius * cos(tan(alpha) - alpha)
#     #         y = radius * sin(tan(alpha) - alpha)
#     #         return x, y
#     #
#     #     return curve_function
#     #
#     # def create_tooth_profile(self):
#     #     """Creates a single tooth profile."""
#     #     # Angle calculations
#     #     angle_inv = tan(self.pressure_angle) - self.pressure_angle
#     #     angle_tip_inv = tan(acos(self.base_radius / self.top_radius)) - acos(self.base_radius / self.top_radius)
#     #     a = 90 / self.teeth + degrees(angle_inv)
#     #     a2 = 90 / self.teeth + degrees(angle_inv) - degrees(angle_tip_inv)
#     #     a3 = 360 / self.teeth - a
#     #
#     #     # Define tooth edges
#     #     right_edge = (
#     #         Workplane()
#     #         .transformed(rotate=(0, 0, a))
#     #         .parametricCurve(self.involute_curve(sign=-1), start=self.base_radius, stop=self.top_radius,
#     #                          N=self.points_per_curve, makeWire=False)
#     #         .val()
#     #     )
#     #
#     #     left_edge = (
#     #         Workplane()
#     #         .transformed(rotate=(0, 0, -a))
#     #         .parametricCurve(self.involute_curve(sign=1), start=self.base_radius, stop=self.top_radius,
#     #                          N=self.points_per_curve, makeWire=False)
#     #         .val()
#     #     )
#     #
#     #     # Define other edges
#     #     top_edge = Edge.makeCircle(self.top_radius, angle1=-a2, angle2=a2)
#     #     bottom_edge = Edge.makeCircle(self.dedendum_radius, angle1=-a3, angle2=-a)
#     #     side_edge = Edge.makeLine(Vector(self.dedendum_radius, 0), Vector(self.base_radius, 0))
#     #     side1 = side_edge.rotate(Vector(0, 0, 0), Vector(0, 0, 1), -a)
#     #     side2 = side_edge.rotate(Vector(0, 0, 0), Vector(0, 0, 1), -a3)
#     #
#     #     # Assemble and return tooth profile
#     #     tooth_profile = Wire.assembleEdges([left_edge, top_edge, right_edge, side1, bottom_edge, side2])
#     #     return tooth_profile
#     #
#     # def create_gear_profile(self):
#     #     """Creates the complete gear profile."""
#     #     tooth_profile = self.create_tooth_profile()
#     #     gear_profile = (
#     #         Workplane()
#     #         .polarArray(0, 0, 360, self.teeth)
#     #         .each(lambda loc: tooth_profile.located(loc))
#     #         .consolidateWires()
#     #     )
#     #     return gear_profile.val()
#
#
# # ----------------------- Gear Model Generation ----------------------- #
#
# # def generate_gear(airfoil):
# #     """Generates the complete gear model with specified parameters."""
# #     fan_geom = FanGeometry(airfoil)
# #     fan_profile = fan_geom.build_sketch()
# #
# #     # Create gear by duplicating tooth profile
# #     fan_blank = Workplane(obj=fan_profile).toPending().extrude(height)
# #
# #     # Cut center hole
# #     # gear_with_hole = gear_blank.faces(">Z").workplane().circle(center_hole_dia / 2).cutThruAll()
# #     #
# #     # key_slot = (
# #     #     Workplane("XY")
# #     #     .rect(center_hole_dia * 0.25, center_hole_dia * 0.25, centered=(True, True))  # Create a rectangular key slot profile
# #     #     .extrude(height)  # Extrude the slot through the gear height
# #     #     .translate((0, center_hole_dia / 2, 0))  # Position the slot on the edge of the center hole
# #     # )
# #
# #     # return gear_with_hole.cut(key_slot)
# #     return fan_blank
#
# def generate_gear(airfoil):
#     data = [
#         [1.00000, 0.00105], [0.95041, 0.0099], [0.90067, 0.01816], [0.80097, 0.03296], [0.70102, 0.04551],
#         [0.60085, 0.0558], [0.50049, 0.06356], [0.4, 0.06837], [0.29875, 0.06875], [0.24814, 0.06668],
#         [0.19761, 0.06276], [0.14722, 0.05665], [0.0971, 0.04766], [0.07217, 0.04169], [0.04742, 0.0342],
#         [0.02297, 0.02411], [0.01098, 0.01694], [0.00000, 0.00000], [0.01402, -0.01448], [0.02703, -0.01927],
#         [0.05258, -0.02482], [0.07783, -0.02809], [0.1029, -0.03016], [0.15278, -0.03227], [0.20239, -0.03276],
#         [0.25186, -0.0323], [0.30125, -0.03125], [0.40000, -0.02837], [0.49951, -0.02468], [0.59915, -0.02024],
#         [0.69898, -0.01551], [0.79903, -0.01074], [0.89933, -0.00594], [0.94959, -0.00352], [1.00000, -0.00105]
#     ]
#     points = np.array(data)  # Assuming your data is in [(x, y), ...]
#
#     x_vals, y_vals = points[:, 0], points[:, 1]
#
#     x_vals, unique_indices = np.unique(x_vals, return_index=True)
#     y_vals = y_vals[unique_indices]
#
#     # Create the spline interpolation (s=0 for exact fit, increase s for smoothing)
#     spline = UnivariateSpline(x_vals, y_vals, s=0)
#
#     # Generate finer points for a smooth curve
#     fine_x = np.linspace(min(x_vals), max(x_vals), 100)
#     fine_y = spline(fine_x)
#     cadquery_points = [(x, y) for x, y in zip(fine_x, fine_y)]
#
#     # print("@@@@@@@@@@", cadquery_points)
#
#     # points = [Vector(x, y, 0) for x, y in data]
#
#     curve = Workplane("XY").spline(points).close().extrude(2)
#     # show_object(curve)
#     return curve
#
# def refine_points(raw_data):
#     profile_points = np.array(raw_data)
#
#     nose_pnt_ind = np.argmin(profile_points[:, 0])
#
#     # Split into upper and lower surfaces
#     upper_points = np.flip(profile_points[:nose_pnt_ind + 1], axis=0)
#     bottom_points = profile_points[nose_pnt_ind:]
#
#     # Extract x and y values
#     x_top, y_top = upper_points[:, 0], upper_points[:, 1]
#     x_bot, y_bot = bottom_points[:, 0], bottom_points[:, 1]
#
#     # Ensure x-values are unique for spline interpolation
#     x_top, unique_top_idx = np.unique(x_top, return_index=True)
#     y_top = y_top[unique_top_idx]
#
#     x_bot, unique_bot_idx = np.unique(x_bot, return_index=True)
#     y_bot = y_bot[unique_bot_idx]
#
#     # Create spline interpolations
#     spl_top = UnivariateSpline(x_top, y_top, s=0.0)
#     spl_bottom = UnivariateSpline(x_bot, y_bot, s=0.0)
#
#     x_dense_top = np.linspace(x_top.min(), x_top.max(), 200)
#     x_dense_bot = np.linspace(x_bot.min(), x_bot.max(), 200)
#
#     y_dense_top = spl_top(x_dense_top)
#     y_dense_bot = spl_bottom(x_dense_bot)
#
#     points = [(x, y) for x, y in zip(x_dense_top, y_dense_top)] + [(x, y) for x, y in zip(x_dense_bot, y_dense_bot)]
#     return points
#
# def fix_knots(knot_vector, min_spacing=1e-6):
#     """ Ensure knots are strictly increasing with a minimum spacing. """
#     knots = np.array(knot_vector)
#     diffs = np.diff(knots)
#
#     if np.any(diffs < min_spacing):
#         print("Adjusting knots to prevent closeness issues...")
#         # Adjust knots by ensuring a minimum spacing
#         for i in range(1, len(knots)):
#             if knots[i] - knots[i - 1] < min_spacing:
#                 knots[i] = knots[i - 1] + min_spacing
#     return knots
#
#
# def get_curve(curve_id):
#     spline = load_spline(curve_id)
#     # t = fix_knots(spline.knots)
#     t = spline.knots
#     c = spline.load_data()[1]
#     k = spline.k
#     print (t, c, k)
#     return build_occ_spline(knot_vector=t, control_points=c, degree=k)
#
#
# def draw_turbofan_profile(points=None, curve=None):
#     # return Workplane("XY").spline(points).close()
#     profile_edge = Edge(BRepBuilderAPI_MakeEdge(curve).Edge())
#     sketch = Sketch().edge(profile_edge).close().assemble()
#     return Workplane().placeSketch(sketch).toOCC()
#
# def generate_turbofan_body(profile, height=2):
#     return profile.extrude(height)







import os
from enum import Enum
import tempfile
import json
from functools import cache, cached_property
from datetime import datetime

from OCP.gp import gp_Pnt
from OCP.TColgp import TColgp_Array1OfPnt
from OCP.TColStd import TColStd_Array1OfInteger, TColStd_Array1OfReal
from OCP.Geom import Geom_BSplineCurve
from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeEdge

import cadquery as cq
import streamlit as st
import pyvista as pv
import numpy as np
import splinecloud_scipy as sc
from stpyvista import stpyvista as stv

@cache
def get_refined_airfoils_collection(data = None, data_file = None):
    if not data:
        with open(data_file) as file:
            data = json.load(file)

    curves_dict = {}
    for category, airfoils in data.items():
        for airfoil_name, airfoil_data in airfoils.items():
            if "spline" in airfoil_data:
                curves_dict[airfoil_name] = airfoil_data["spline"]

    return curves_dict

TOLERANCE = 0.1
AIRFOILS_COLLECTION_FILE_LOCATION = f"{os.getcwd()}/airfoils/airfoils_collection.json"
REFINED_AIRFOILS_COLLECTION = get_refined_airfoils_collection(data_file=AIRFOILS_COLLECTION_FILE_LOCATION)


class Tessellation(Enum):
    COARSE = 1e-3
    MEDIUM = 1e-4
    FINE = 1e-5


class AirfoilSection:

    def __init__(self, profile_curve, chord, offset, twist_angle):
        self.profile_curve = profile_curve
        self.chord = chord
        self.offset = cq.Vector(0, 0, offset)
        self.twist_axis = cq.Vector(0, 0, 1)  # twist around Z axis
        self.twist_angle = twist_angle

    @cached_property
    def spline(self):
        if isinstance(self.profile_curve, sc.ParametricUnivariateSpline):
            return self.profile_curve
        return sc.load_spline(self.profile_curve)

    @staticmethod
    def _build_spline(knot_vector, control_points, degree):
        knot_mults = list(map(lambda k: degree + 1 if k in (0, 1) else 1, knot_vector))

        poles = TColgp_Array1OfPnt(1, len(control_points))
        for i, cp in enumerate(control_points):
            poles.SetValue(i + 1, gp_Pnt(cp[0], cp[1], 0))

        knots = TColStd_Array1OfReal(1, len(knot_vector))
        for i, k in enumerate(knot_vector):
            knots.SetValue(i + 1, k)

        mults = TColStd_Array1OfInteger(1, len(knot_mults))
        for i, m in enumerate(knot_mults):
            mults.SetValue(i + 1, m)

        return Geom_BSplineCurve(poles, knots, mults, degree)

    @cached_property
    def bspl(self):
        k = self.spline.k
        t = self.spline.knots[k:-k]
        cp = np.array([self.spline.coeffs_x, self.spline.coeffs_y]).T * self.chord
        return self._build_spline(t, cp, k)

    def build_sketch(self):
        edge = cq.Edge(BRepBuilderAPI_MakeEdge(self.bspl).Edge())
        return cq.Sketch().edge(edge).close().assemble()

    @cached_property
    def location(self):
        return cq.Location(self.offset, self.twist_axis, self.twist_angle)


class Turbofan:
    def __init__(self, sections: list, vanes_count: int, center_hole_diameter: int, hub_diameter: int):
        self.sections = sections
        self.vanes_count = vanes_count
        self.center_hole_diameter = center_hole_diameter
        self.hub_diameter = hub_diameter

    @cached_property
    def _vane_trailing_edge_max_coordinate(self):
        edges = self.sections[0].build_sketch().moved(self.sections[0].location).edges()
        vertices = [v for e in edges for v in e.Vertices()]
        vertex_coords = [v.toTuple() for v in vertices]
        trailing_edge_tip_coords = max(vertex_coords, key=lambda v: v[0])
        return max(trailing_edge_tip_coords)


    @cache
    def build_vane(self):
        if not self.sections:
            raise ValueError(f"The number of sections has to be > `0`, but you provided `{len(self.sections)}`")
        return cq.Workplane().placeSketch(
            *(section.build_sketch().moved(section.location) for section in self.sections)
        ).loft()

    @cached_property
    def hub_height(self):
        bbox = self.build_vane().val().BoundingBox()
        return bbox.ymax - bbox.ymin

    @cache
    def build_turbofan(self):
        vane_offset = ((self.hub_diameter / 2) ** 2 - self._vane_trailing_edge_max_coordinate ** 2) ** 0.5
        vanes = [
            self.build_vane().translate((0, 0, vane_offset)).rotate((0, 0, 0), (0, 1, 0), i * 360 / self.vanes_count)
            for i in range(self.vanes_count)]

        hub = cq.Workplane("XZ").circle(self.hub_diameter / 2).extrude(self.hub_height)
        hub = hub.faces("XZ").workplane().circle(self.hub_diameter / 2).extrude(-self.hub_height / 2)
        hub = hub.faces("+Y").edges().fillet(self.hub_height * 0.25)

        for vane in vanes:
            hub = hub.union(vane)

        return hub.faces("XZ").circle(self.center_hole_diameter / 2).cutThruAll()


# naca4421 = sc.load_spline('spl_QBhalXE844DC')
# naca4418 = sc.load_spline('spl_4OfePVC7ztg0')
# naca4412 = sc.load_spline('spl_oWHqlVtu67as')


def generate_temp_file(model, file_format, tessellation):
    """Generates a temporary file with the specified STL or STEP format."""
    if file_format.lower() == "stl":
        file_suffix = ".stl"
        export_type = cq.exporters.ExportTypes.STL

    elif file_format.lower() == "step":
        file_suffix = ".step"
        export_type = cq.exporters.ExportTypes.STEP
    else:
        raise ValueError(f"Unsupported file format: {file_format}")

    with tempfile.NamedTemporaryFile(delete=False, suffix=file_suffix) as tmpfile:
        cq.exporters.export(model, tmpfile.name, exportType=export_type, tolerance=getattr(Tessellation, tessellation.upper()).value)
        return tmpfile.name

@st.cache_data
def generate_and_export_turbofan_cached(
    _root_curve,
    _middle_curve,
    _tip_curve,
    root_chord_ratio,
    root_twist,
    middle_chord_ratio,
    middle_offset_distance,
    middle_twist,
    tip_chord_ratio,
    tip_offset_distance,
    tip_twist,
    file_format,
    tessellation,
    vanes_count,
    hub_diameter,
    center_hole_diameter
):
    sections = [
        AirfoilSection(_root_curve, root_chord_ratio, 0, root_twist),
        AirfoilSection(_middle_curve, middle_chord_ratio, middle_offset_distance, middle_twist),
        AirfoilSection(_tip_curve, tip_chord_ratio, middle_offset_distance + tip_offset_distance, tip_twist)
    ]
    turbofan = Turbofan(
        sections=sections,
        vanes_count=vanes_count,
        center_hole_diameter=center_hole_diameter,
        hub_diameter=hub_diameter
    ).build_turbofan()
    return generate_temp_file(turbofan, file_format, tessellation)


# ----------------------- Streamlit UI ----------------------- #

# Streamlit UI for parameter inputs
st.set_page_config(layout='wide')
col1, col2, col3 = st.columns([2, 2, 4])
with col1:
    st.header("Turbofan Model Generator")
    st.subheader("Root airfoil profile selection:")
    family_root = st.selectbox("Choose a profile family:", list(REFINED_AIRFOILS_COLLECTION.keys()), key="root_family")
    root_curve_id = REFINED_AIRFOILS_COLLECTION[family_root]
    st.divider()

    st.subheader("Middle airfoil profile selection:")
    family_middle = st.selectbox("Choose a profile family:", list(REFINED_AIRFOILS_COLLECTION.keys()), key="middle_family")
    middle_curve_id = REFINED_AIRFOILS_COLLECTION[family_middle]
    st.divider()

    st.subheader("Tip airfoil profile selection:")
    family_tip = st.selectbox("Choose a profile family:", list(REFINED_AIRFOILS_COLLECTION.keys()), key="tip_family")
    tip_curve_id = REFINED_AIRFOILS_COLLECTION[family_tip]

with col2:
    tessellation_value = st.select_slider('Surface quality', options=[t.name.lower() for t in Tessellation])
    blades_count = st.slider('Blades count (items)', min_value=2, max_value=10, value=2, step=1)
    hub_dia = st.slider('Turbofan hub diameter', min_value=2.0, max_value=5.0, value=2.0, step=0.25)
    center_hole_dia = st.slider('Turbofan center hole diameter', min_value=0.5, max_value=hub_dia-0.5, value=1.0, step=0.25)
    root_twist_angle = st.slider('Root section twist angle (degree)', min_value=-90, max_value=0, value=-20, step=1)
    middle_twist_angle = st.slider('Middle section twist angle (degree)', min_value=-90, max_value=0, value=-12, step=1)
    tip_twist_angle = st.slider('Tip section twist angle (degree)', min_value=-90, max_value=0, value=-3, step=1)

    root_chord = st.slider('Root section chord (ratio)', min_value=0.1, max_value=2.0, value=0.6, step=0.1)
    middle_chord = st.slider('Middle section chord (ratio)', min_value=0.1, max_value=2.0, value=1.0, step=0.1)
    tip_chord = st.slider('Tip section chord (ratio)', min_value=0.1, max_value=2.0, value=0.3, step=0.1)

    middle_offset = st.slider('Middle section offset', min_value=1, max_value=10, value=2, step=1)
    tip_offset = st.slider('Tip section offset', min_value=1, max_value=10, value=3, step=1)


# ----------------------- Visualization ----------------------- #
file_path_stl = generate_and_export_turbofan_cached(
    # _root_curve=naca4421,
    # _middle_curve=naca4418,
    # _tip_curve=naca4412,
    _root_curve=root_curve_id,
    _middle_curve=middle_curve_id,
    _tip_curve=tip_curve_id,
    root_chord_ratio=root_chord,
    root_twist=root_twist_angle,
    middle_chord_ratio=middle_chord,
    middle_offset_distance=middle_offset,
    middle_twist=middle_twist_angle,
    tip_chord_ratio=tip_chord,
    tip_offset_distance=tip_offset,
    tip_twist=tip_twist_angle,
    file_format="stl",
    tessellation=tessellation_value,
    vanes_count=blades_count,
    hub_diameter=hub_dia,
    center_hole_diameter=center_hole_dia,
)


file_path_step = generate_and_export_turbofan_cached(
    # _root_curve=naca4421,
    # _middle_curve=naca4418,
    # _tip_curve=naca4412,
    _root_curve=root_curve_id,
    _middle_curve=middle_curve_id,
    _tip_curve=tip_curve_id,
    root_chord_ratio=root_chord,
    root_twist=root_twist_angle,
    middle_chord_ratio=middle_chord,
    middle_offset_distance=middle_offset,
    middle_twist=middle_twist_angle,
    tip_chord_ratio=tip_chord,
    tip_offset_distance=tip_offset,
    tip_twist=tip_twist_angle,
    file_format="step",
    tessellation=tessellation_value,
    vanes_count=blades_count,
    hub_diameter=hub_dia,
    center_hole_diameter=center_hole_dia,
)

if os.getenv("OS_TYPE") != "windows":
    pv.start_xvfb()

mesh = pv.read(file_path_stl)

plotter = pv.Plotter(window_size=[500, 500])
plotter.add_mesh(mesh)
plotter.view_isometric()
plotter.set_background("#0e1117")

with col3:
    stv(plotter, key=f"turbofan_{datetime.now()}")
    with open(file_path_stl, 'rb') as stl_file:
        st.download_button(
            label="Download STL File",
            data=stl_file,
            file_name="gear.stl",
            mime="application/vnd.ms-pki.stl",
        )
    with open(file_path_step, 'rb') as step_file:
        st.download_button(
            label=f"Download STEP File",
            data=step_file,
            file_name=f"gear.step",
            mime="application/step",
        )
