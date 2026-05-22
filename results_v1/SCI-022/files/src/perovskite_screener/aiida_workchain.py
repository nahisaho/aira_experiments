"""
AiiDA WorkChain for lead-free perovskite high-throughput screening.
Requires: aiida-core >= 2.0, aiida-vasp >= 3.0
"""
from aiida.engine import WorkChain, calcfunction, if_, while_, ToContext
from aiida.orm import Dict, StructureData, List, Float
from aiida.plugins import CalculationFactory, DataFactory


VaspCalculation   = CalculationFactory("vasp.vasp")
VaspNEBCalculation = CalculationFactory("vasp.neb")


class PerovskiteScreeningWorkChain(WorkChain):
    """
    AiiDA WorkChain for high-throughput Sn/Ge/Bi perovskite screening.

    Workflow outline:
      1. tolerance_filter  – discard non-perovskite compositions
      2. ml_prescreen      – ML band gap filter (0.9–2.5 eV window)
      3. dft_relax         – VASP structure relaxation (PBE+D3)
      4. dft_bands         – HSE06+SOC band structure
      5. defect_calc       – Defect formation energies (supercell)
      6. neb_calc          – CI-NEB ion migration barriers
      7. scaps_sim         – SCAPS-1D device simulation
      8. aggregate_rank    – Multi-criteria ranking
    """

    @classmethod
    def define(cls, spec):
        super().define(spec)
        spec.input("candidates", valid_type=List)
        spec.input("dft_parameters", valid_type=Dict)
        spec.input("screening_parameters", valid_type=Dict)
        spec.output("ranked_candidates", valid_type=List)
        spec.output("top_material", valid_type=Dict)
        spec.output("screening_report", valid_type=Dict)

        spec.outline(
            cls.tolerance_filter,
            if_(cls.has_candidates)(
                cls.ml_prescreen,
                if_(cls.has_ml_passed)(
                    cls.dft_relax,
                    cls.dft_bands,
                    cls.defect_calc,
                    cls.neb_calc,
                    cls.scaps_sim,
                )
            ),
            cls.aggregate_rank,
        )

        spec.exit_code(300, "NO_CANDIDATES_PASSED_FILTER",
                       "No candidates passed the tolerance factor filter")
        spec.exit_code(400, "DFT_FAILED",
                       "DFT calculation failed for all candidates")

    def tolerance_filter(self):
        from perovskite_screener.tolerance_factor import analyze_perovskite
        passed = []
        for c in self.inputs.candidates:
            res = analyze_perovskite(c["A"], c["B"], c["X"], c.get("B_ox", 2))
            if res.stability_class == "perovskite":
                passed.append(c)
        self.ctx.candidates = passed
        if not passed:
            return self.exit_codes.NO_CANDIDATES_PASSED_FILTER
        self.report(f"Tolerance filter: {len(passed)}/{len(self.inputs.candidates)} passed")

    def has_candidates(self):
        return len(self.ctx.candidates) > 0

    def ml_prescreen(self):
        from perovskite_screener.bandgap_ml import BandGapPredictor
        predictor = BandGapPredictor().fit(verbose=False)
        passed = []
        for c in self.ctx.candidates:
            res = predictor.predict(c["A"], c["B"], c["X"], c.get("B_ox", 2))
            if 0.9 <= res["Eg_predicted_eV"] <= 2.5:
                c["Eg_ml"] = res["Eg_predicted_eV"]
                passed.append(c)
        self.ctx.candidates = passed
        self.report(f"ML prescreening: {len(passed)} candidates in 0.9–2.5 eV window")

    def has_ml_passed(self):
        return len(self.ctx.candidates) > 0

    def dft_relax(self):
        calcs = {}
        for c in self.ctx.candidates:
            code     = self.inputs.dft_parameters["vasp_code"]
            builder  = VaspCalculation.get_builder()
            builder.structure = c["structure"]
            builder.parameters = Dict(dict=self.inputs.dft_parameters["relax_incar"])
            builder.kpoints = DataFactory("array.kpoints")()
            builder.kpoints.set_kpoints_mesh([6, 6, 6])
            future = self.submit(builder)
            calcs[f"relax_{c['formula']}"] = future
        return ToContext(**calcs)

    def dft_bands(self):
        # Similar to dft_relax but with HSE06+SOC settings
        pass

    def defect_calc(self):
        pass

    def neb_calc(self):
        pass

    def scaps_sim(self):
        from perovskite_screener.scaps_interface import run_scaps_simulation
        results = {}
        for c in self.ctx.candidates:
            res = run_scaps_simulation(c["A"], c["B"], c["X"], c.get("Eg_ml", 1.5))
            results[c["formula"]] = res
        self.ctx.scaps_results = results

    def aggregate_rank(self):
        from perovskite_screener.ranking import rank_candidates
        ranked = rank_candidates(
            self.ctx.candidates,
            scaps_results=getattr(self.ctx, "scaps_results", {}),
        )
        self.out("ranked_candidates", List(list=ranked))
        self.out("top_material", Dict(dict=ranked[0] if ranked else {}))
