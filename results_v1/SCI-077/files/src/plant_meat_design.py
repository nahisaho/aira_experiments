"""
Module 6: 植物性代替肉のテクスチャ設計ケーススタディ
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Dict, Optional
from scipy.optimize import differential_evolution


@dataclass
class PlantProteinFormulation:
    soy_protein: float = 0.0
    pea_protein: float = 0.0
    wheat_gluten: float = 0.0
    starch: float = 0.0
    fat_content: float = 0.0
    fiber: float = 0.0
    moisture: float = 60.0
    salt: float = 1.0
    methylcellulose: float = 0.0

    @property
    def total_protein(self):
        return self.soy_protein + self.pea_protein + self.wheat_gluten

    def validate(self):
        total = (self.soy_protein + self.pea_protein + self.wheat_gluten +
                 self.starch + self.fat_content + self.fiber +
                 self.moisture + self.salt + self.methylcellulose)
        return abs(total - 100) < 5.0


@dataclass
class HMECProcessConditions:
    barrel_temperature: float = 150.0
    die_temperature: float = 80.0
    screw_speed: float = 200.0
    feed_rate: float = 5.0
    moisture_content: float = 60.0
    cooling_die_length: float = 300.0

    @property
    def specific_mechanical_energy(self):
        return 0.5 * self.screw_speed * (100 - self.moisture_content) / self.feed_rate


REFERENCE_MEATS = {
    'beef_patty': {'hardness':450,'cohesiveness':0.55,'springiness':0.75,
                   'gumminess':247,'chewiness':186,'juiciness':0.65,'fiber_alignment':0.3},
    'chicken_breast': {'hardness':600,'cohesiveness':0.60,'springiness':0.80,
                       'gumminess':360,'chewiness':288,'juiciness':0.55,'fiber_alignment':0.7},
    'pork_sausage': {'hardness':300,'cohesiveness':0.50,'springiness':0.70,
                     'gumminess':150,'chewiness':105,'juiciness':0.70,'fiber_alignment':0.2},
}


def hmec_texture_prediction(formulation, process):
    protein_factor = formulation.total_protein / 30
    temp_factor = np.exp(-((process.barrel_temperature - 150) / 30)**2)
    screw_factor = np.log1p(process.screw_speed / 100)
    moisture_factor = np.exp(-((process.moisture_content - 60) / 10)**2)
    gluten_factor = 1 + 0.02 * formulation.wheat_gluten
    DT = np.clip(protein_factor*temp_factor*screw_factor*moisture_factor*gluten_factor, 0, 1)
    die_factor = np.tanh(process.cooling_die_length / 200)
    cooling_factor = np.exp(-((process.die_temperature - 70) / 20)**2)
    FAI = np.clip(DT * die_factor * cooling_factor, 0, 1)
    hardness = (20*formulation.total_protein + 5*formulation.starch +
                30*formulation.wheat_gluten*DT + 10*formulation.methylcellulose) * (1+0.5*DT)
    cohesiveness = np.clip(0.3+0.01*formulation.total_protein*DT+0.02*formulation.wheat_gluten*FAI-0.005*formulation.fat_content, 0.2, 0.8)
    springiness = np.clip(0.5+0.01*formulation.wheat_gluten+0.005*formulation.total_protein*DT+0.01*formulation.methylcellulose, 0.3, 0.95)
    juiciness = np.clip(0.3+0.02*formulation.fat_content+0.005*formulation.moisture-0.005*formulation.total_protein*DT, 0.1, 0.9)
    return {'degree_of_texturization': DT, 'fiber_alignment_index': FAI,
            'hardness': hardness, 'cohesiveness': cohesiveness,
            'springiness': springiness, 'juiciness': juiciness,
            'gumminess': hardness*cohesiveness, 'chewiness': hardness*cohesiveness*springiness,
            'specific_mechanical_energy': process.specific_mechanical_energy}


def texture_similarity_score(predicted, target_meat):
    ref = REFERENCE_MEATS[target_meat]
    scores = {}
    for p in ['hardness','cohesiveness','springiness','juiciness']:
        pv, rv = predicted.get(p, 0), ref[p]
        scores[p] = np.exp(-((pv-rv)/(0.3*rv))**2) if rv > 0 else (1.0 if pv==0 else 0.0)
    fai, ref_fai = predicted.get('fiber_alignment_index', 0), ref.get('fiber_alignment', 0)
    scores['fiber_alignment'] = np.exp(-((fai-ref_fai)/0.3)**2)
    scores['overall'] = np.mean(list(scores.values()))
    return scores


def optimize_formulation(target_meat='beef_patty'):
    process = HMECProcessConditions()
    def objective(x):
        soy,pea,gluten,starch,fat,fiber,mc = x
        moisture = 100-soy-pea-gluten-starch-fat-fiber-mc-1.0
        if moisture < 40 or moisture > 75:
            return 10.0
        form = PlantProteinFormulation(soy,pea,gluten,starch,fat,fiber,moisture,1.0,mc)
        if not form.validate():
            return 10.0
        pred = hmec_texture_prediction(form, process)
        sim = texture_similarity_score(pred, target_meat)
        return 1 - sim['overall']
    bounds = [(5,25),(0,20),(0,15),(0,10),(2,15),(0,5),(0,3)]
    result = differential_evolution(objective, bounds, seed=42, maxiter=200, tol=1e-4, popsize=20)
    opt = result.x
    moisture = 100 - sum(opt) - 1.0
    form = PlantProteinFormulation(opt[0],opt[1],opt[2],opt[3],opt[4],opt[5],moisture,1.0,opt[6])
    pred = hmec_texture_prediction(form, process)
    sim = texture_similarity_score(pred, target_meat)
    return {'optimal_formulation': {'soy_protein':opt[0],'pea_protein':opt[1],'wheat_gluten':opt[2],
            'starch':opt[3],'fat_content':opt[4],'fiber':opt[5],'methylcellulose':opt[6],
            'moisture':moisture,'salt':1.0},
            'predicted_texture': pred, 'similarity_scores': sim,
            'target_meat': target_meat, 'optimization_success': result.success}


def run_case_study():
    results = {}
    for meat in ['beef_patty','chicken_breast','pork_sausage']:
        results[meat] = optimize_formulation(target_meat=meat)
    process = HMECProcessConditions()
    sensitivity = {'soy_protein':[],'hardness':[],'similarity':[]}
    for soy in np.linspace(5, 25, 20):
        form = PlantProteinFormulation(soy,10,5,3,8,2,100-soy-10-5-3-8-2-1-1,1.0,1)
        pred = hmec_texture_prediction(form, process)
        sim = texture_similarity_score(pred, 'beef_patty')
        sensitivity['soy_protein'].append(soy)
        sensitivity['hardness'].append(pred['hardness'])
        sensitivity['similarity'].append(sim['overall'])
    results['sensitivity_analysis'] = sensitivity
    return results
