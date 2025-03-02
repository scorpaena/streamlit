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
            if "profile_curve" in airfoil_data:
                if profile_curve := airfoil_data["profile_curve"]:
                    curves_dict[airfoil_name] = profile_curve[0]

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
