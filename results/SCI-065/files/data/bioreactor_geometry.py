#!/usr/bin/env python3
"""Generate parametric geometry and OpenFOAM-ready mesh metadata for a perfusion bioreactor."""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path
from typing import Any, Dict


def cylindrical_volume(radius: float, height: float) -> float:
    return math.pi * radius**2 * height


def make_wedge_vertices(
    inner_radius: float, outer_radius: float, height: float, wedge_angle_deg: float
) -> list[list[float]]:
    half_angle = math.radians(wedge_angle_deg / 2.0)
    cos_a = math.cos(half_angle)
    sin_a = math.sin(half_angle)
    return [
        [inner_radius * cos_a, -inner_radius * sin_a, 0.0],
        [outer_radius * cos_a, -outer_radius * sin_a, 0.0],
        [inner_radius * cos_a, inner_radius * sin_a, 0.0],
        [outer_radius * cos_a, outer_radius * sin_a, 0.0],
        [inner_radius * cos_a, -inner_radius * sin_a, height],
        [outer_radius * cos_a, -outer_radius * sin_a, height],
        [inner_radius * cos_a, inner_radius * sin_a, height],
        [outer_radius * cos_a, outer_radius * sin_a, height],
    ]


def build_geometry() -> Dict[str, Any]:
    vessel_diameter = 80e-3
    vessel_height = 120e-3
    port_diameter = 6e-3
    basket_diameter = 60e-3
    basket_height = 80e-3
    basket_z_min = 20e-3
    basket_z_max = basket_z_min + basket_height
    wedge_angle_deg = 5.0

    vessel_radius = vessel_diameter / 2.0
    basket_radius = basket_diameter / 2.0
    port_radius = port_diameter / 2.0
    axis_core_radius = 0.5e-3

    radial_cells = 48
    axial_cells = 96
    azimuthal_cells = 1

    vertices = make_wedge_vertices(axis_core_radius, vessel_radius, vessel_height, wedge_angle_deg)

    geometry = {
        "units": "SI",
        "vessel": {
            "diameter_m": vessel_diameter,
            "radius_m": vessel_radius,
            "height_m": vessel_height,
            "volume_m3": cylindrical_volume(vessel_radius, vessel_height),
        },
        "ports": {
            "count": 2,
            "diameter_m": port_diameter,
            "radius_m": port_radius,
            "cross_section_area_m2": math.pi * port_radius**2,
            "inlet_location": "bottom axial center",
            "outlet_location": "top axial center",
        },
        "basket": {
            "diameter_m": basket_diameter,
            "radius_m": basket_radius,
            "height_m": basket_height,
            "z_min_m": basket_z_min,
            "z_max_m": basket_z_max,
            "volume_m3": cylindrical_volume(basket_radius, basket_height),
        },
        "flow_rates_ml_min": [0.5, 1.0, 2.0, 5.0],
        "mesh": {
            "description": "Axisymmetric wedge mesh metadata for OpenFOAM blockMesh.",
            "convertToMeters": 1.0,
            "wedge_angle_deg": wedge_angle_deg,
            "axis_core_radius_m": axis_core_radius,
            "cells": {
                "radial": radial_cells,
                "axial": axial_cells,
                "azimuthal": azimuthal_cells,
            },
            "grading": {
                "radial": 1.0,
                "axial": 1.0,
                "azimuthal": 1.0,
            },
            "vertices": vertices,
            "blocks": [
                {
                    "type": "hex",
                    "vertices": [0, 1, 3, 2, 4, 5, 7, 6],
                    "cells": [radial_cells, azimuthal_cells, axial_cells],
                    "grading": [1.0, 1.0, 1.0],
                }
            ],
            "edges": [
                {"type": "arc", "vertices": [0, 2], "midpoint": [axis_core_radius, 0.0, 0.0]},
                {"type": "arc", "vertices": [1, 3], "midpoint": [vessel_radius, 0.0, 0.0]},
                {"type": "arc", "vertices": [4, 6], "midpoint": [axis_core_radius, 0.0, vessel_height]},
                {"type": "arc", "vertices": [5, 7], "midpoint": [vessel_radius, 0.0, vessel_height]},
            ],
            "boundaries": {
                "inlet": {"type": "patch", "faces": [[0, 2, 3, 1]]},
                "outlet": {"type": "patch", "faces": [[4, 5, 7, 6]]},
                "wall": {"type": "wall", "faces": [[1, 3, 7, 5]]},
                "core": {"type": "symmetryPlane", "faces": [[0, 4, 6, 2]]},
                "wedgeLow": {"type": "wedge", "faces": [[0, 1, 5, 4]]},
                "wedgeHigh": {"type": "wedge", "faces": [[2, 6, 7, 3]]},
            },
            "basket_zone": {
                "selection": "Cylindrical basket approximated by a wedge-aligned box for topoSet.",
                "box_min_m": [axis_core_radius * math.cos(math.radians(wedge_angle_deg / 2.0)), -basket_radius * math.sin(math.radians(wedge_angle_deg / 2.0)), basket_z_min],
                "box_max_m": [basket_radius * math.cos(math.radians(wedge_angle_deg / 2.0)), basket_radius * math.sin(math.radians(wedge_angle_deg / 2.0)), basket_z_max],
            },
        },
    }
    return geometry


def write_json(data: Dict[str, Any], destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", encoding="utf-8") as handle:
        json.dump(data, handle, indent=2)


def main() -> int:
    base_dir = Path(__file__).resolve().parent
    output_path = base_dir / "geometry_params.json"
    try:
        geometry = build_geometry()
        write_json(geometry, output_path)
        print(f"Geometry parameters written to {output_path}")
        return 0
    except Exception as exc:  # pragma: no cover - defensive runtime handling
        print(f"Failed to generate geometry parameters: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
