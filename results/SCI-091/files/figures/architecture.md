```
┌─────────────────────────────────────────────────────────────────────────┐
│              Research Integrity Assessment System (RIAS)                │
│              NLP + Computer Vision 統合型 研究公正性評価                  │
└──────────────────────────────┬──────────────────────────────────────────┘
                               │
                        ┌──────▼──────┐
                        │ Paper Input │
                        │ (PDF/Text)  │
                        └──────┬──────┘
                               │
        ┌──────────┬───────────┼───────────┬──────────┬──────────┐
        │          │           │           │          │          │
   ┌────▼────┐ ┌──▼───┐ ┌────▼────┐ ┌───▼────┐ ┌──▼───┐ ┌───▼────┐
   │ Module1 │ │Mod 2 │ │ Module3 │ │Module4 │ │Mod 5 │ │Module6 │
   │ Image   │ │Stats │ │Plagiar- │ │P-hack/ │ │Repro │ │Extern  │
   │Forensic │ │Check │ │  ism    │ │HARKing │ │Score │ │Signal  │
   │ (0.20)  │ │(0.20)│ │ (0.15)  │ │ (0.20) │ │(0.15)│ │ (0.10) │
   └────┬────┘ └──┬───┘ └────┬────┘ └───┬────┘ └──┬───┘ └───┬────┘
        │         │          │          │         │          │
   ┌────┴────┐ ┌──┴───┐ ┌───┴────┐ ┌──┴────┐ ┌──┴───┐ ┌───┴────┐
   │• pHash  │ │•GRIM │ │•Winnow │ │•Pcurve│ │•Meth │ │•PubPeer│
   │• ELA    │ │•SPRIT│ │•Cite-  │ │•Zcurve│ │ Eval │ │•Retract│
   │• Block  │ │  E   │ │ Aware  │ │•Calip │ │•RIAS │ │ Watch  │
   │  Match  │ │•p-val│ │ Sim    │ │•HARK  │ │Score │ │•Ground │
   │• CNN    │ │ dist │ │•SimHash│ │ ing   │ │•Pred │ │ Truth  │
   └────┬────┘ └──┬───┘ └───┬────┘ └──┬────┘ └──┬───┘ └───┬────┘
        │         │          │         │         │          │
        └────┬────┴──────────┴────┬────┴─────────┴──────────┘
             │                    │
        ┌────▼────────────────────▼────┐
        │    Weighted Score Fusion     │
        │ S = Σ(wᵢ × moduleᵢ_risk)    │
        │ Integrity = 1 - S            │
        └──────────────┬───────────────┘
                       │
            ┌──────────▼──────────┐
            │  Integrity Report   │
            │ ├ Overall Score     │
            │ ├ Risk Level        │
            │ ├ Module Details    │
            │ ├ Concerns          │
            │ └ Recommendations   │
            └─────────────────────┘
```
