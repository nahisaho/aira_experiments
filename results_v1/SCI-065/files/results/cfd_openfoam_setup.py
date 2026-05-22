#!/usr/bin/env python3
"""Generate OpenFOAM-compatible case directories for perfusion bioreactor CFD studies."""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path
from typing import Dict, Iterable

DENSITY = 1007.0
VISCOSITY = 0.001
POROSITY = 0.4
PERMEABILITY = 1.0e-10
FORCHHEIMER = 1.0e5
FLOW_RATES_ML_MIN = [0.5, 1.0, 2.0, 5.0]


def load_geometry(geometry_path: Path) -> Dict:
    if not geometry_path.exists():
        raise FileNotFoundError(
            f"Missing geometry definition at {geometry_path}. Run data/bioreactor_geometry.py first."
        )
    with geometry_path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def write_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content.strip() + "\n", encoding="utf-8")


def inlet_velocity(flow_rate_ml_min: float, port_area: float) -> float:
    flow_rate_m3_s = flow_rate_ml_min * 1.0e-6 / 60.0
    return flow_rate_m3_s / port_area


def foam_header(class_name: str, object_name: str, location: str) -> str:
    return f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       {class_name};
    location    \"{location}\";
    object      {object_name};
}}
"""


def format_vector(values: Iterable[float]) -> str:
    return "(" + " ".join(f"{value:.8g}" for value in values) + ")"


def build_block_mesh_dict(geometry: Dict) -> str:
    mesh = geometry["mesh"]
    vertices = mesh["vertices"]
    block = mesh["blocks"][0]
    edges = mesh["edges"]
    boundary = mesh["boundaries"]

    vertices_text = "\n".join(f"    {format_vector(vertex)}" for vertex in vertices)
    block_text = (
        f"    {block['type']} ({' '.join(str(v) for v in block['vertices'])}) "
        f"({block['cells'][0]} {block['cells'][1]} {block['cells'][2]}) "
        f"simpleGrading {format_vector(block['grading'])}"
    )
    edges_text = "\n".join(
        f"    {edge['type']} {edge['vertices'][0]} {edge['vertices'][1]} {format_vector(edge['midpoint'])}"
        for edge in edges
    )

    boundary_entries = []
    for name, patch in boundary.items():
        faces_text = "\n".join(
            f"            ({' '.join(str(vertex) for vertex in face)})" for face in patch["faces"]
        )
        boundary_entries.append(
            f"    {name}\n    {{\n        type {patch['type']};\n        faces\n        (\n{faces_text}\n        );\n    }}"
        )
    boundary_text = "\n".join(boundary_entries)

    return (
        foam_header("dictionary", "blockMeshDict", "system")
        + f"convertToMeters {mesh['convertToMeters']};\n\n"
        + "vertices\n(\n"
        + vertices_text
        + "\n);\n\nblocks\n(\n"
        + block_text
        + "\n);\n\nedges\n(\n"
        + edges_text
        + "\n);\n\nboundary\n(\n"
        + boundary_text
        + "\n);\n\nmergePatchPairs\n(\n);"
    )


def build_topo_set_dict(geometry: Dict) -> str:
    box_min = format_vector(geometry["mesh"]["basket_zone"]["box_min_m"])
    box_max = format_vector(geometry["mesh"]["basket_zone"]["box_max_m"])
    return (
        foam_header("dictionary", "topoSetDict", "system")
        + f"actions\n(\n    {{\n        name basketCells;\n        type cellSet;\n        action new;\n        source boxToCell;\n        sourceInfo\n        {{\n            box {box_min} {box_max};\n        }}\n    }}\n    {{\n        name basketZone;\n        type cellZoneSet;\n        action new;\n        source setToCellZone;\n        sourceInfo\n        {{\n            set basketCells;\n        }}\n    }}\n);"
    )


def build_fv_options() -> str:
    darcy = VISCOSITY / PERMEABILITY
    return (
        foam_header("dictionary", "fvOptions", "constant")
        + f"basketResistance\n{{\n    type explicitPorositySource;\n    active yes;\n    explicitPorositySourceCoeffs\n    {{\n        selectionMode cellZone;\n        cellZone basketZone;\n        type DarcyForchheimer;\n        // porosity = {POROSITY}\n        DarcyForchheimerCoeffs\n        {{\n            d d [0 -2 0 0 0 0 0] ({darcy:.6e} {darcy:.6e} {darcy:.6e});\n            f f [0 -1 0 0 0 0 0] ({FORCHHEIMER:.6e} {FORCHHEIMER:.6e} {FORCHHEIMER:.6e});\n            coordinateSystem\n            {{\n                type cartesian;\n                origin (0 0 0);\n                coordinateRotation\n                {{\n                    type axesRotation;\n                    e1 (1 0 0);\n                    e2 (0 1 0);\n                }}\n            }}\n        }}\n    }}\n}}"
    )


def build_control_dict() -> str:
    return (
        foam_header("dictionary", "controlDict", "system")
        + "application     simpleFoam;\nstartFrom       startTime;\nstartTime       0;\nstopAt          endTime;\nendTime         2000;\ndeltaT          1;\nwriteControl    timeStep;\nwriteInterval   200;\npurgeWrite      0;\nwriteFormat     ascii;\nwritePrecision  8;\nwriteCompression off;\ntimeFormat      general;\ntimePrecision   6;\nrunTimeModifiable yes;"
    )


def build_fv_schemes() -> str:
    return (
        foam_header("dictionary", "fvSchemes", "system")
        + "ddtSchemes\n{\n    default steadyState;\n}\n\ngradSchemes\n{\n    default Gauss linear;\n}\n\ndivSchemes\n{\n    default none;\n    div(phi,U) Gauss linearUpwind grad(U);\n    div(phi,k) Gauss upwind;\n    div(phi,epsilon) Gauss upwind;\n    div((nuEff*dev2(T(grad(U))))) Gauss linear;\n}\n\nlaplacianSchemes\n{\n    default Gauss linear corrected;\n}\n\ninterpolationSchemes\n{\n    default linear;\n}\n\nsnGradSchemes\n{\n    default corrected;\n}"
    )


def build_fv_solution() -> str:
    return (
        foam_header("dictionary", "fvSolution", "system")
        + "solvers\n{\n    p\n    {\n        solver GAMG;\n        tolerance 1e-08;\n        relTol 0.05;\n        smoother GaussSeidel;\n    }\n\n    U\n    {\n        solver smoothSolver;\n        smoother symGaussSeidel;\n        tolerance 1e-09;\n        relTol 0.1;\n    }\n}\n\nSIMPLE\n{\n    nNonOrthogonalCorrectors 1;\n    consistent yes;\n    residualControl\n    {\n        p 1e-4;\n        U 1e-5;\n    }\n}\n\nrelaxationFactors\n{\n    fields\n    {\n        p 0.3;\n    }\n    equations\n    {\n        U 0.7;\n    }\n}"
    )


def build_transport_properties() -> str:
    return (
        foam_header("dictionary", "transportProperties", "constant")
        + f"transportModel  Newtonian;\nnu              nu [0 2 -1 0 0 0 0] {VISCOSITY / DENSITY:.10e};\nrho             rho [1 -3 0 0 0 0 0] {DENSITY:.6f};\ndynamicViscosity mu [1 -1 -1 0 0 0 0] {VISCOSITY:.6f};"
    )


def build_u_file(velocity: float) -> str:
    return (
        foam_header("volVectorField", "U", "0")
        + f"dimensions      [0 1 -1 0 0 0 0];\ninternalField   uniform (0 0 0);\nboundaryField\n{{\n    inlet\n    {{\n        type fixedValue;\n        value uniform (0 0 {velocity:.8e});\n    }}\n    outlet\n    {{\n        type zeroGradient;\n    }}\n    wall\n    {{\n        type noSlip;\n    }}\n    wedgeLow\n    {{\n        type wedge;\n    }}\n    wedgeHigh\n    {{\n        type wedge;\n    }}\n    core\n    {{\n        type symmetryPlane;\n    }}\n}}"
    )


def build_p_file() -> str:
    return (
        foam_header("volScalarField", "p", "0")
        + "dimensions      [0 2 -2 0 0 0 0];\ninternalField   uniform 0;\nboundaryField\n{\n    inlet\n    {\n        type zeroGradient;\n    }\n    outlet\n    {\n        type fixedValue;\n        value uniform 0;\n    }\n    wall\n    {\n        type zeroGradient;\n    }\n    wedgeLow\n    {\n        type wedge;\n    }\n    wedgeHigh\n    {\n        type wedge;\n    }\n    core\n    {\n        type symmetryPlane;\n    }\n}"
    )


def build_case(geometry: Dict, case_dir: Path, flow_rate_ml_min: float) -> None:
    port_area = geometry["ports"]["cross_section_area_m2"]
    velocity = inlet_velocity(flow_rate_ml_min, port_area)

    ensure_dir(case_dir / "0")
    ensure_dir(case_dir / "constant")
    ensure_dir(case_dir / "system")

    write_text(case_dir / "system" / "controlDict", build_control_dict())
    write_text(case_dir / "system" / "fvSchemes", build_fv_schemes())
    write_text(case_dir / "system" / "fvSolution", build_fv_solution())
    write_text(case_dir / "system" / "blockMeshDict", build_block_mesh_dict(geometry))
    write_text(case_dir / "system" / "topoSetDict", build_topo_set_dict(geometry))
    write_text(case_dir / "constant" / "fvOptions", build_fv_options())
    write_text(case_dir / "constant" / "transportProperties", build_transport_properties())
    write_text(case_dir / "0" / "U", build_u_file(velocity))
    write_text(case_dir / "0" / "p", build_p_file())
    readme = (
        f"Flow rate: {flow_rate_ml_min:.1f} mL/min\n"
        f"Inlet velocity: {velocity:.6e} m/s\n"
        f"Porosity: {POROSITY}\nPermeability: {PERMEABILITY:.2e} m^2\n"
        f"Forchheimer coefficient: {FORCHHEIMER:.2e} 1/m\n"
        "Run sequence:\n  blockMesh\n  topoSet\n  simpleFoam\n"
    )
    write_text(case_dir / "README.txt", readme)


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    geometry_path = root / "data" / "geometry_params.json"
    try:
        geometry = load_geometry(geometry_path)
        base_dir = root / "results" / "openfoam_cases"
        for flow_rate in FLOW_RATES_ML_MIN:
            case_name = f"flow_{str(flow_rate).replace('.', 'p')}_ml_min"
            build_case(geometry, base_dir / case_name, flow_rate)
        print(f"OpenFOAM cases created under {base_dir}")
        return 0
    except Exception as exc:  # pragma: no cover - runtime guard
        print(f"Failed to build OpenFOAM setup: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
