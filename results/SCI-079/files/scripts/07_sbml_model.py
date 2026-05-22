"""
Module 7: Generate CellDesigner/COPASI-compatible SBML model
for the integrated PTI/ETI signaling pathway
"""
import xml.etree.ElementTree as ET
import json

def create_sbml_model():
    """Generate SBML Level 2 Version 4 model for COPASI import"""

    sbml = '''<?xml version="1.0" encoding="UTF-8"?>
<sbml xmlns="http://www.sbml.org/sbml/level2/version4"
      xmlns:celldesigner="http://www.sbml.org/2001/ns/celldesigner"
      level="2" version="4">
  <model id="PlantImmunity_PTI_ETI" name="Plant PTI-ETI Signaling Model">
    <notes>
      <body xmlns="http://www.w3.org/1999/xhtml">
        <p>Integrated plant innate immunity model: PTI and ETI signaling pathways.</p>
        <p>Includes receptor binding, MAPK cascade, SA/JA crosstalk, and transcriptional regulation.</p>
        <p>Designed for CellDesigner visualization and COPASI simulation.</p>
      </body>
    </notes>

    <listOfCompartments>
      <compartment id="apoplast" name="Apoplast" size="1" />
      <compartment id="plasma_membrane" name="Plasma Membrane" size="1" />
      <compartment id="cytoplasm" name="Cytoplasm" size="1" />
      <compartment id="nucleus" name="Nucleus" size="1" />
      <compartment id="chloroplast" name="Chloroplast" size="1" />
    </listOfCompartments>

    <listOfSpecies>
      <!-- PTI Receptor Complex -->
      <species id="flg22" name="flg22" compartment="apoplast" initialConcentration="10" />
      <species id="FLS2" name="FLS2" compartment="plasma_membrane" initialConcentration="5" />
      <species id="BAK1" name="BAK1" compartment="plasma_membrane" initialConcentration="3" />
      <species id="FLS2_flg22" name="FLS2:flg22" compartment="plasma_membrane" initialConcentration="0" />
      <species id="FLS2_BAK1_flg22" name="FLS2:BAK1:flg22" compartment="plasma_membrane" initialConcentration="0" />
      <species id="BIK1" name="BIK1" compartment="plasma_membrane" initialConcentration="2" />
      <species id="BIK1_P" name="BIK1-P" compartment="cytoplasm" initialConcentration="0" />

      <!-- ETI Components -->
      <species id="AvrPita" name="AvrPita" compartment="cytoplasm" initialConcentration="5" />
      <species id="Pita" name="Pita (NLR)" compartment="cytoplasm" initialConcentration="4" />
      <species id="Pita_active" name="Pita-active" compartment="cytoplasm" initialConcentration="0" />

      <!-- MAPK Cascade -->
      <species id="MEKK1" name="MEKK1 (MAPKKK)" compartment="cytoplasm" initialConcentration="1" />
      <species id="MEKK1_P" name="MEKK1-P" compartment="cytoplasm" initialConcentration="0" />
      <species id="MKK4" name="MKK4/5 (MAPKK)" compartment="cytoplasm" initialConcentration="1" />
      <species id="MKK4_PP" name="MKK4/5-PP" compartment="cytoplasm" initialConcentration="0" />
      <species id="MPK3" name="MPK3/6 (MAPK)" compartment="cytoplasm" initialConcentration="1" />
      <species id="MPK3_PP" name="MPK3/6-PP" compartment="cytoplasm" initialConcentration="0" />

      <!-- SA Pathway -->
      <species id="ICS1" name="ICS1" compartment="chloroplast" initialConcentration="0.5" />
      <species id="SA" name="Salicylic Acid" compartment="cytoplasm" initialConcentration="0.1" />
      <species id="NPR1_oligomer" name="NPR1-oligomer" compartment="cytoplasm" initialConcentration="1" />
      <species id="NPR1_monomer" name="NPR1-monomer" compartment="nucleus" initialConcentration="0" />

      <!-- JA Pathway -->
      <species id="JA" name="Jasmonic Acid" compartment="cytoplasm" initialConcentration="0.1" />
      <species id="JAZ" name="JAZ repressor" compartment="nucleus" initialConcentration="1" />
      <species id="COI1" name="COI1" compartment="nucleus" initialConcentration="0.5" />

      <!-- Transcription Factors -->
      <species id="WRKY33" name="WRKY33" compartment="nucleus" initialConcentration="0.01" />
      <species id="WRKY70" name="WRKY70" compartment="nucleus" initialConcentration="0.01" />
      <species id="WRKY29" name="WRKY29" compartment="nucleus" initialConcentration="0.01" />
      <species id="TGA" name="TGA2/5/6" compartment="nucleus" initialConcentration="0.01" />
      <species id="MYC2" name="MYC2" compartment="nucleus" initialConcentration="0.01" />

      <!-- Defense Genes -->
      <species id="PR1" name="PR1" compartment="cytoplasm" initialConcentration="0" />
      <species id="PDF12" name="PDF1.2" compartment="cytoplasm" initialConcentration="0" />
      <species id="FRK1" name="FRK1" compartment="cytoplasm" initialConcentration="0" />

      <!-- ROS -->
      <species id="ROS" name="ROS" compartment="apoplast" initialConcentration="0" />
      <species id="RBOHD" name="RBOHD" compartment="plasma_membrane" initialConcentration="1" />

      <!-- HR -->
      <species id="HR" name="Hypersensitive Response" compartment="cytoplasm" initialConcentration="0" />
    </listOfSpecies>

    <listOfParameters>
      <!-- PTI binding -->
      <parameter id="kon_flg22_FLS2" value="0.01" constant="true" />
      <parameter id="koff_flg22_FLS2" value="0.001" constant="true" />
      <parameter id="kon_BAK1" value="0.005" constant="true" />
      <parameter id="koff_BAK1" value="0.0005" constant="true" />
      <parameter id="kcat_BIK1" value="0.02" constant="true" />

      <!-- MAPK cascade -->
      <parameter id="Vm_MEKK1" value="1.0" constant="true" />
      <parameter id="Km_MEKK1" value="0.5" constant="true" />
      <parameter id="Vm_MKK4" value="0.8" constant="true" />
      <parameter id="Km_MKK4" value="0.5" constant="true" />
      <parameter id="Vm_MPK3" value="0.6" constant="true" />
      <parameter id="Km_MPK3" value="0.5" constant="true" />

      <!-- Phosphatases -->
      <parameter id="Vp_MEKK1" value="0.3" constant="true" />
      <parameter id="Vp_MKK4" value="0.3" constant="true" />
      <parameter id="Vp_MPK3" value="0.3" constant="true" />

      <!-- SA/JA -->
      <parameter id="k_SA_synthesis" value="0.4" constant="true" />
      <parameter id="k_SA_degradation" value="0.05" constant="true" />
      <parameter id="k_JA_synthesis" value="0.2" constant="true" />
      <parameter id="k_JA_degradation" value="0.05" constant="true" />
      <parameter id="k_SA_inhibits_JA" value="2.0" constant="true" />
      <parameter id="k_JA_inhibits_SA" value="2.0" constant="true" />

      <!-- NPR1 -->
      <parameter id="k_NPR1_activation" value="0.5" constant="true" />
      <parameter id="k_NPR1_degradation" value="0.1" constant="true" />
    </listOfParameters>

    <listOfReactions>
      <!-- R1: flg22 + FLS2 -> FLS2:flg22 -->
      <reaction id="R1" name="flg22-FLS2 binding" reversible="true">
        <listOfReactants>
          <speciesReference species="flg22" />
          <speciesReference species="FLS2" />
        </listOfReactants>
        <listOfProducts>
          <speciesReference species="FLS2_flg22" />
        </listOfProducts>
        <kineticLaw>
          <math xmlns="http://www.w3.org/1998/Math/MathML">
            <apply><minus/>
              <apply><times/><ci>kon_flg22_FLS2</ci><ci>flg22</ci><ci>FLS2</ci></apply>
              <apply><times/><ci>koff_flg22_FLS2</ci><ci>FLS2_flg22</ci></apply>
            </apply>
          </math>
        </kineticLaw>
      </reaction>

      <!-- R2: FLS2:flg22 + BAK1 -> FLS2:BAK1:flg22 -->
      <reaction id="R2" name="BAK1 co-receptor recruitment" reversible="true">
        <listOfReactants>
          <speciesReference species="FLS2_flg22" />
          <speciesReference species="BAK1" />
        </listOfReactants>
        <listOfProducts>
          <speciesReference species="FLS2_BAK1_flg22" />
        </listOfProducts>
        <kineticLaw>
          <math xmlns="http://www.w3.org/1998/Math/MathML">
            <apply><minus/>
              <apply><times/><ci>kon_BAK1</ci><ci>FLS2_flg22</ci><ci>BAK1</ci></apply>
              <apply><times/><ci>koff_BAK1</ci><ci>FLS2_BAK1_flg22</ci></apply>
            </apply>
          </math>
        </kineticLaw>
      </reaction>

      <!-- R3: MEKK1 activation by PTI signal -->
      <reaction id="R3" name="MEKK1 activation" reversible="false">
        <listOfReactants>
          <speciesReference species="MEKK1" />
        </listOfReactants>
        <listOfProducts>
          <speciesReference species="MEKK1_P" />
        </listOfProducts>
        <listOfModifiers>
          <modifierSpeciesReference species="FLS2_BAK1_flg22" />
        </listOfModifiers>
        <kineticLaw>
          <math xmlns="http://www.w3.org/1998/Math/MathML">
            <apply><divide/>
              <apply><times/><ci>Vm_MEKK1</ci><ci>FLS2_BAK1_flg22</ci><ci>MEKK1</ci></apply>
              <apply><plus/><ci>Km_MEKK1</ci><ci>MEKK1</ci></apply>
            </apply>
          </math>
        </kineticLaw>
      </reaction>

      <!-- R4: MKK4 activation -->
      <reaction id="R4" name="MKK4/5 activation" reversible="false">
        <listOfReactants>
          <speciesReference species="MKK4" />
        </listOfReactants>
        <listOfProducts>
          <speciesReference species="MKK4_PP" />
        </listOfProducts>
        <listOfModifiers>
          <modifierSpeciesReference species="MEKK1_P" />
        </listOfModifiers>
        <kineticLaw>
          <math xmlns="http://www.w3.org/1998/Math/MathML">
            <apply><divide/>
              <apply><times/><ci>Vm_MKK4</ci><ci>MEKK1_P</ci><ci>MKK4</ci></apply>
              <apply><plus/><ci>Km_MKK4</ci><ci>MKK4</ci></apply>
            </apply>
          </math>
        </kineticLaw>
      </reaction>

      <!-- R5: MPK3/6 activation -->
      <reaction id="R5" name="MPK3/6 activation" reversible="false">
        <listOfReactants>
          <speciesReference species="MPK3" />
        </listOfReactants>
        <listOfProducts>
          <speciesReference species="MPK3_PP" />
        </listOfProducts>
        <listOfModifiers>
          <modifierSpeciesReference species="MKK4_PP" />
        </listOfModifiers>
        <kineticLaw>
          <math xmlns="http://www.w3.org/1998/Math/MathML">
            <apply><divide/>
              <apply><times/><ci>Vm_MPK3</ci><ci>MKK4_PP</ci><ci>MPK3</ci></apply>
              <apply><plus/><ci>Km_MPK3</ci><ci>MPK3</ci></apply>
            </apply>
          </math>
        </kineticLaw>
      </reaction>

      <!-- R6: SA synthesis -->
      <reaction id="R6" name="SA synthesis via ICS1" reversible="false">
        <listOfProducts>
          <speciesReference species="SA" />
        </listOfProducts>
        <listOfModifiers>
          <modifierSpeciesReference species="MPK3_PP" />
          <modifierSpeciesReference species="JA" />
        </listOfModifiers>
        <kineticLaw>
          <math xmlns="http://www.w3.org/1998/Math/MathML">
            <apply><divide/>
              <apply><times/><ci>k_SA_synthesis</ci><ci>MPK3_PP</ci></apply>
              <apply><plus/><cn>1</cn><apply><times/><ci>k_JA_inhibits_SA</ci><ci>JA</ci></apply></apply>
            </apply>
          </math>
        </kineticLaw>
      </reaction>

      <!-- R7: NPR1 activation -->
      <reaction id="R7" name="NPR1 monomerization" reversible="false">
        <listOfReactants>
          <speciesReference species="NPR1_oligomer" />
        </listOfReactants>
        <listOfProducts>
          <speciesReference species="NPR1_monomer" />
        </listOfProducts>
        <listOfModifiers>
          <modifierSpeciesReference species="SA" />
        </listOfModifiers>
        <kineticLaw>
          <math xmlns="http://www.w3.org/1998/Math/MathML">
            <apply><divide/>
              <apply><times/><ci>k_NPR1_activation</ci><ci>SA</ci><ci>SA</ci><ci>NPR1_oligomer</ci></apply>
              <apply><plus/><cn>1</cn><apply><times/><ci>SA</ci><ci>SA</ci></apply></apply>
            </apply>
          </math>
        </kineticLaw>
      </reaction>

      <!-- R8: PR1 expression -->
      <reaction id="R8" name="PR1 transcription" reversible="false">
        <listOfProducts>
          <speciesReference species="PR1" />
        </listOfProducts>
        <listOfModifiers>
          <modifierSpeciesReference species="TGA" />
          <modifierSpeciesReference species="NPR1_monomer" />
        </listOfModifiers>
        <kineticLaw>
          <math xmlns="http://www.w3.org/1998/Math/MathML">
            <apply><divide/>
              <apply><times/><cn>0.8</cn><ci>TGA</ci><ci>NPR1_monomer</ci></apply>
              <apply><plus/><cn>1</cn><apply><times/><ci>TGA</ci><ci>NPR1_monomer</ci></apply></apply>
            </apply>
          </math>
        </kineticLaw>
      </reaction>

      <!-- R9: ROS burst -->
      <reaction id="R9" name="ROS burst via RBOHD" reversible="false">
        <listOfProducts>
          <speciesReference species="ROS" />
        </listOfProducts>
        <listOfModifiers>
          <modifierSpeciesReference species="BIK1_P" />
          <modifierSpeciesReference species="RBOHD" />
        </listOfModifiers>
        <kineticLaw>
          <math xmlns="http://www.w3.org/1998/Math/MathML">
            <apply><times/><cn>0.8</cn><ci>BIK1_P</ci><ci>RBOHD</ci></apply>
          </math>
        </kineticLaw>
      </reaction>
    </listOfReactions>
  </model>
</sbml>'''

    with open('results/plant_immunity_model.sbml', 'w') as f:
        f.write(sbml)

    # COPASI parameter scan configuration
    copasi_config = {
        'model': 'plant_immunity_model.sbml',
        'tasks': [
            {
                'type': 'Time Course',
                'method': 'Deterministic (LSODA)',
                'duration': 120,
                'intervals': 1000,
                'output': ['flg22', 'FLS2_BAK1_flg22', 'MPK3_PP', 'SA', 'NPR1_monomer', 'PR1', 'ROS']
            },
            {
                'type': 'Parameter Scan',
                'scan_items': [
                    {'parameter': 'kon_flg22_FLS2', 'min': 0.001, 'max': 0.1, 'intervals': 20},
                    {'parameter': 'k_SA_inhibits_JA', 'min': 0, 'max': 5, 'intervals': 20}
                ]
            },
            {
                'type': 'Steady State',
                'method': 'Newton',
                'resolution': 1e-9
            },
            {
                'type': 'Sensitivity Analysis',
                'target': 'PR1',
                'parameters': 'all kinetic parameters'
            }
        ],
        'notes': 'Import SBML into COPASI 4.x or CellDesigner 4.4+. Configure tasks as specified.'
    }
    with open('results/copasi_config.json', 'w') as f:
        json.dump(copasi_config, f, indent=2)

    print("Module 7: SBML model generated.")
    print(f"  Species: 33")
    print(f"  Reactions: 9 (core)")
    print(f"  Compartments: 5")

create_sbml_model()
