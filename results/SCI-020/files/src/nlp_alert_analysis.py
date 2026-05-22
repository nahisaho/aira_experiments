"""
NLP Alert Analysis Module.
Parses ProMED/WHO alerts using TF-IDF, NER, and DBSCAN clustering.
"""

import numpy as np
import pandas as pd
import re
from datetime import datetime, timedelta
from typing import Dict, List, Optional, Tuple
import warnings
warnings.filterwarnings("ignore")

try:
    from sklearn.feature_extraction.text import TfidfVectorizer
    from sklearn.cluster import DBSCAN
    from sklearn.preprocessing import normalize
    from sklearn.metrics.pairwise import cosine_similarity
    HAS_SKLEARN = True
except ImportError:
    HAS_SKLEARN = False


# ── Synthetic ProMED/WHO alert corpus ──────────────────────────────────────────
SAMPLE_ALERTS = [
    {
        "source": "ProMED",
        "date": "2025-05-01",
        "title": "UNDIAGNOSED RESPIRATORY ILLNESS - CHINA: (HUBEI) REQUEST FOR INFORMATION",
        "text": (
            "An undiagnosed respiratory illness has been reported in Wuhan, Hubei Province, China. "
            "As of 30 December, 27 cases of pneumonia of unknown etiology have been identified. "
            "Clinical presentation includes fever, dry cough, and dyspnea. No deaths have been "
            "reported so far. The source of infection and transmission route remain unknown. "
            "Local health authorities are conducting epidemiological investigations."
        ),
    },
    {
        "source": "WHO",
        "date": "2025-05-03",
        "title": "Novel coronavirus - China",
        "text": (
            "WHO has been informed of a cluster of cases of pneumonia of unknown etiology detected "
            "in Wuhan City, Hubei Province of China. As of 3 January 2025, 44 patients with "
            "pneumonia of unknown etiology have been officially reported. A novel coronavirus has "
            "been identified as the causative agent. WHO advises against the application of any "
            "travel or trade restrictions. Human-to-human transmission appears limited."
        ),
    },
    {
        "source": "ProMED",
        "date": "2025-05-05",
        "title": "INFLUENZA A (H5N1) - CAMBODIA: HUMAN CASE",
        "text": (
            "A fatal human case of influenza A (H5N1) has been reported in Cambodia. The patient, "
            "a 9-year-old female from Prey Veng Province, developed illness on 16 February 2025 "
            "with fever and cough. She was hospitalized on 21 February and died on 22 February 2025. "
            "Exposure to dead poultry was identified. No other cases have been identified through "
            "contact tracing. The risk of human-to-human transmission remains low."
        ),
    },
    {
        "source": "HealthMap",
        "date": "2025-05-07",
        "title": "Mpox outbreak update - Democratic Republic of Congo",
        "text": (
            "The Democratic Republic of Congo continues to report the highest burden of mpox globally. "
            "As of week 20 of 2025, 3,890 cases including 148 deaths (CFR 3.8%) have been reported "
            "from 22 of 26 provinces. Clade I remains predominant, with Clade Ib showing increased "
            "transmissibility. Cross-border spread to neighboring countries including Uganda, Rwanda, "
            "and Burundi has been confirmed. Vaccination campaigns with MVA-BN vaccine are ongoing."
        ),
    },
    {
        "source": "WHO",
        "date": "2025-05-09",
        "title": "Ebola virus disease - Uganda: update",
        "text": (
            "The Ministry of Health of Uganda declared an Ebola Virus Disease (EVD) outbreak on "
            "20 September 2022. As of 5 October 2022, a total of 63 confirmed and probable cases "
            "including 29 deaths (CFR 46%) have been reported from 6 districts. The outbreak is "
            "caused by the Sudan strain of Ebola virus. Contact tracing is ongoing for 1,658 contacts. "
            "No approved vaccine is available for Sudan Ebolavirus."
        ),
    },
    {
        "source": "ProMED",
        "date": "2025-05-11",
        "title": "CHOLERA - MULTIPLE COUNTRIES: INCREASED ACTIVITY",
        "text": (
            "Cholera outbreaks continue to be reported across multiple countries. Haiti reported 847 "
            "cases in the past week. Syria continues to experience widespread transmission with over "
            "90,000 cases since the current outbreak began. Malawi has reported more than 58,000 cases. "
            "The resurgence is attributed to disrupted water and sanitation infrastructure, displacement, "
            "and reduced vaccination coverage."
        ),
    },
    {
        "source": "WHO",
        "date": "2025-05-13",
        "title": "Marburg virus disease - Equatorial Guinea: update",
        "text": (
            "Equatorial Guinea is reporting a Marburg virus disease outbreak. As of 19 March 2023, "
            "9 confirmed cases and 11 probable cases, including 12 deaths, have been reported in "
            "Kie-Ntem Province. Marburg virus, a highly lethal hemorrhagic fever, has no approved "
            "vaccines or specific treatments. Response teams are conducting contact tracing of 4,000+ "
            "contacts. No spread to neighboring countries has been confirmed."
        ),
    },
    {
        "source": "HealthMap",
        "date": "2025-05-15",
        "title": "Novel variant SARS-CoV-2 JN.1 sublineage rapid spread",
        "text": (
            "A novel sublineage of JN.1, designated KP.2, is showing rapid spread across multiple "
            "countries. KP.2 carries additional mutations in the spike protein receptor binding domain "
            "at positions L455S and F456L, conferring increased immune escape from BA.2.86 antibodies. "
            "Current prevalence: Japan 34%, USA 22%, Germany 18%. Vaccine effectiveness against "
            "symptomatic disease estimated at 52-68%. Hospitalizations remain at pre-pandemic baseline."
        ),
    },
    {
        "source": "ProMED",
        "date": "2025-05-17",
        "title": "DENGUE FEVER - BRAZIL: (SAO PAULO) EMERGENCY DECLARED",
        "text": (
            "Sao Paulo state has declared a dengue fever public health emergency following a 200% "
            "increase in cases compared to the same period last year. Over 500,000 dengue cases have "
            "been reported nationwide in the first 5 weeks of 2025. The outbreak is driven primarily "
            "by DENV-3 serotype following years of low circulation and reduced population immunity. "
            "Vector control measures including elimination of Aedes aegypti breeding sites are being reinforced."
        ),
    },
    {
        "source": "WHO",
        "date": "2025-05-19",
        "title": "Antimicrobial resistance in Klebsiella pneumoniae: global update",
        "text": (
            "WHO has issued a global alert on the emergence of pan-drug-resistant Klebsiella pneumoniae "
            "strains carrying NDM-1, KPC-2, and OXA-48 carbapenemases simultaneously. These 'superbugs' "
            "have been identified in 47 countries. Mortality rates in ICU patients exceed 80%. "
            "The strains are plasmid-mediated and show horizontal gene transfer potential. WHO urges "
            "immediate strengthening of infection prevention and control measures globally."
        ),
    },
]


# ── NER patterns ────────────────────────────────────────────────────────────────
PATHOGEN_PATTERNS = [
    r'\b(SARS-CoV-2|COVID-19|coronavirus|influenza|H5N1|H1N1|ebola|marburg|mpox|monkeypox'
    r'|cholera|dengue|KP\.\d+|JN\.\d+|BA\.\d+|XBB|BQ\.\d+|Klebsiella|NDM-1|KPC)\b'
]
LOCATION_PATTERNS = [
    r'\b(China|Wuhan|Cambodia|Congo|DRC|Uganda|Haiti|Syria|Malawi|Brazil|Japan|USA|'
    r'Germany|Equatorial Guinea|Hubei|Sao Paulo|Africa|Asia|Europe|global|worldwide)\b'
]
CASE_COUNT_PATTERNS = [
    r'(\d[\d,]*)\s*(?:confirmed\s+)?cases',
    r'(\d[\d,]*)\s*deaths',
    r'CFR\s+([0-9.]+)%',
]


class ProMEDParser:
    """Parses ProMED and HealthMap alert text."""

    def extract_named_entities(self, text: str) -> Dict[str, List[str]]:
        """Simple regex-based NER for pathogens, locations, case counts."""
        pathogens = []
        for pattern in PATHOGEN_PATTERNS:
            pathogens.extend(re.findall(pattern, text, re.IGNORECASE))

        locations = []
        for pattern in LOCATION_PATTERNS:
            locations.extend(re.findall(pattern, text, re.IGNORECASE))

        case_counts = []
        for pattern in CASE_COUNT_PATTERNS:
            case_counts.extend(re.findall(pattern, text, re.IGNORECASE))

        return {
            "pathogens": list(set(pathogens)),
            "locations": list(set(locations)),
            "case_counts": case_counts,
        }

    def classify_alert_severity(self, text: str, title: str = "") -> str:
        """
        Rule-based severity classification.
        Returns: LOW / MEDIUM / HIGH / CRITICAL
        """
        combined = (title + " " + text).upper()

        # CRITICAL indicators
        critical_keywords = [
            "PANDEMIC", "GLOBAL EMERGENCY", "PHEIC", "NOVEL PATHOGEN",
            "UNKNOWN ETIOLOGY", "UNDIAGNOSED", "HIGHLY LETHAL", "PAN-DRUG-RESISTANT",
            "EBOLA", "MARBURG", "HEMORRHAGIC FEVER",
        ]
        if any(kw in combined for kw in critical_keywords):
            return "CRITICAL"

        # HIGH indicators
        high_keywords = [
            "OUTBREAK", "EMERGENCY DECLARED", "HUMAN CASE", "RAPID SPREAD",
            "IMMUNE ESCAPE", "CROSS-BORDER", "FATALITY", "DEATHS", "CFR",
            "NOVEL VARIANT", "NOVEL SUBLINEAGE",
        ]
        high_count = sum(1 for kw in high_keywords if kw in combined)
        if high_count >= 2:
            return "HIGH"

        # MEDIUM indicators
        medium_keywords = [
            "INCREASED ACTIVITY", "UPDATE", "SURVEILLANCE", "MONITORING",
            "CONTACT TRACING", "INVESTIGATION", "CLUSTER",
        ]
        medium_count = sum(1 for kw in medium_keywords if kw in combined)
        if medium_count >= 1 or high_count == 1:
            return "MEDIUM"

        return "LOW"

    def parse_alerts(self, alerts: List[Dict]) -> pd.DataFrame:
        """Parse a list of alert dicts into a structured DataFrame."""
        records = []
        for alert in alerts:
            entities = self.extract_named_entities(alert["text"])
            severity = self.classify_alert_severity(alert["text"], alert.get("title", ""))
            records.append({
                "source": alert.get("source", "Unknown"),
                "date": alert.get("date", ""),
                "title": alert.get("title", ""),
                "severity": severity,
                "pathogens": ", ".join(entities["pathogens"]) if entities["pathogens"] else "Unknown",
                "locations": ", ".join(entities["locations"]) if entities["locations"] else "Unknown",
                "n_pathogens": len(entities["pathogens"]),
                "n_locations": len(entities["locations"]),
                "text_length": len(alert.get("text", "")),
            })
        return pd.DataFrame(records)


class WHOAlertAnalyzer:
    """Analyzes WHO Disease Outbreak News."""

    def compute_novelty_score(self, alert_text: str,
                               corpus: List[str],
                               top_n: int = 20) -> float:
        """TF-IDF cosine distance from corpus mean = novelty score."""
        if not HAS_SKLEARN or len(corpus) < 2:
            return 0.5

        all_docs = corpus + [alert_text]
        try:
            vectorizer = TfidfVectorizer(max_features=500, stop_words="english",
                                          ngram_range=(1, 2))
            tfidf = vectorizer.fit_transform(all_docs)
            tfidf_norm = normalize(tfidf)
            corpus_centroid = tfidf_norm[:-1].mean(axis=0)
            query_vec = tfidf_norm[-1]
            sim = float(cosine_similarity(query_vec, corpus_centroid)[0, 0])
            return round(1.0 - sim, 4)  # Dissimilarity = novelty
        except Exception:
            return 0.5


class EpidemicSignalExtractor:
    """Clusters related alerts and extracts epidemic signals."""

    def cluster_alerts(self, alerts_df: pd.DataFrame) -> pd.DataFrame:
        """DBSCAN clustering of alerts by TF-IDF similarity."""
        if not HAS_SKLEARN or len(alerts_df) < 3:
            alerts_df = alerts_df.copy()
            alerts_df["cluster"] = 0
            return alerts_df

        texts = alerts_df["title"].tolist()
        try:
            vectorizer = TfidfVectorizer(max_features=200, stop_words="english")
            tfidf = vectorizer.fit_transform(texts)
            tfidf_norm = normalize(tfidf)
            dist_matrix = 1 - cosine_similarity(tfidf_norm)
            db = DBSCAN(eps=0.7, min_samples=2, metric="precomputed")
            labels = db.fit_predict(dist_matrix)
            alerts_df = alerts_df.copy()
            alerts_df["cluster"] = labels
        except Exception:
            alerts_df = alerts_df.copy()
            alerts_df["cluster"] = -1
        return alerts_df

    def generate_summary(self, alerts_df: pd.DataFrame) -> List[Dict]:
        """Generate extractive cluster summaries."""
        summaries = []
        for cluster_id in alerts_df["cluster"].unique():
            cluster = alerts_df[alerts_df["cluster"] == cluster_id]
            dominant_severity = cluster["severity"].value_counts().idxmax()
            all_pathogens = ", ".join(
                set(p.strip() for row in cluster["pathogens"]
                    for p in str(row).split(",") if p.strip() and p.strip() != "Unknown")
            ) or "Unknown"
            summaries.append({
                "cluster_id": cluster_id,
                "n_alerts": len(cluster),
                "dominant_severity": dominant_severity,
                "pathogens": all_pathogens,
                "sources": ", ".join(cluster["source"].unique()),
                "date_range": f"{cluster['date'].min()} to {cluster['date'].max()}",
            })
        return summaries

    def compute_signal_score(self, alerts_df: pd.DataFrame) -> float:
        """Aggregate NLP signal score from alert severities."""
        severity_weights = {"LOW": 10, "MEDIUM": 30, "HIGH": 60, "CRITICAL": 90}
        scores = alerts_df["severity"].map(severity_weights).fillna(0)
        return round(float(scores.mean()), 2)


def run_nlp_pipeline(alerts: Optional[List[Dict]] = None) -> Dict:
    """Run the full NLP alert analysis pipeline."""
    alerts = alerts or SAMPLE_ALERTS

    parser = ProMEDParser()
    who_analyzer = WHOAlertAnalyzer()
    extractor = EpidemicSignalExtractor()

    # Parse alerts
    alerts_df = parser.parse_alerts(alerts)

    # Compute novelty scores
    corpus = [a["text"] for a in alerts]
    novelty_scores = []
    for alert in alerts:
        score = who_analyzer.compute_novelty_score(alert["text"], corpus)
        novelty_scores.append(score)
    alerts_df["novelty_score"] = novelty_scores

    # Cluster alerts
    alerts_df = extractor.cluster_alerts(alerts_df)

    # Summarize clusters
    cluster_summaries = extractor.generate_summary(alerts_df)

    # Overall signal score
    signal_score = extractor.compute_signal_score(alerts_df)

    # Severity distribution
    severity_dist = alerts_df["severity"].value_counts().to_dict()

    return {
        "alerts_df": alerts_df,
        "cluster_summaries": cluster_summaries,
        "signal_score": signal_score,
        "severity_distribution": severity_dist,
        "n_alerts_processed": len(alerts_df),
        "n_critical": int(severity_dist.get("CRITICAL", 0)),
        "n_high": int(severity_dist.get("HIGH", 0)),
        "top_novelty_alert": alerts_df.nlargest(1, "novelty_score")["title"].iloc[0]
        if len(alerts_df) else "",
    }


if __name__ == "__main__":
    results = run_nlp_pipeline()
    print(f"Alerts processed: {results['n_alerts_processed']}")
    print(f"Signal score: {results['signal_score']}")
    print(f"Severity distribution: {results['severity_distribution']}")
    print(f"Top novel alert: {results['top_novelty_alert']}")
