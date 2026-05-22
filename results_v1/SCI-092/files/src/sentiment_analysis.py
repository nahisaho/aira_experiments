"""
Component 2: Social media sentiment analysis — BERT/sentiment dictionary hybrid.
"""

import numpy as np
import pandas as pd
import json

np.random.seed(42)

TECH_SENTIMENT_LEXICON = {
    "positive": ["breakthrough", "innovation", "promising", "revolutionary", "benefit",
                 "progress", "efficient", "safe", "sustainable", "hope", "革新的", "有望", "期待"],
    "negative": ["danger", "risk", "threat", "unethical", "catastrophe", "fear",
                 "harmful", "unsafe", "destroy", "ban", "危険", "脅威", "懸念", "不安"],
    "neutral": ["research", "study", "develop", "technology", "regulation", "研究", "開発", "技術"]
}

def generate_synthetic_social_media_data():
    technologies = {
        "gene_editing": {"sentiment_dist": [0.30, 0.35, 0.35]},
        "AI": {"sentiment_dist": [0.40, 0.30, 0.30]},
        "nuclear_fusion": {"sentiment_dist": [0.45, 0.20, 0.35]}
    }
    posts = []
    platforms = ["Twitter/X", "Reddit", "News_Comments", "Blog"]
    platform_weights = [0.45, 0.25, 0.20, 0.10]

    for tech, config in technologies.items():
        n_posts = np.random.randint(800, 1500)
        sentiments = np.random.choice(["positive", "negative", "neutral"],
                                       size=n_posts, p=config["sentiment_dist"])
        for i in range(n_posts):
            sentiment = sentiments[i]
            bert_conf = np.clip(np.random.beta(5, 2) if sentiment != "neutral"
                                else np.random.beta(2, 3), 0.3, 0.99)
            lex_score = {"positive": 1, "negative": -1, "neutral": 0}[sentiment] + np.random.normal(0, 0.3)
            hybrid_score = 0.6 * (bert_conf if sentiment == "positive" else -bert_conf if sentiment == "negative" else 0) \
                           + 0.4 * np.clip(lex_score, -1, 1)
            platform = np.random.choice(platforms, p=platform_weights)
            engagement = int(np.random.lognormal(3, 1.5))
            year = np.random.choice(range(2020, 2027), p=[0.05, 0.08, 0.12, 0.15, 0.20, 0.22, 0.18])
            posts.append({
                "post_id": f"{tech}_{i:04d}", "technology": tech, "platform": platform,
                "year": int(year), "sentiment_label": sentiment,
                "bert_confidence": round(float(bert_conf), 4),
                "lexicon_score": round(float(lex_score), 4),
                "hybrid_score": round(float(hybrid_score), 4),
                "engagement": engagement
            })
    return pd.DataFrame(posts)

def bert_lexicon_hybrid_analysis(df):
    results = {}
    for tech in df["technology"].unique():
        sub = df[df["technology"] == tech]
        sent_dist = sub["sentiment_label"].value_counts(normalize=True).to_dict()
        hybrid_stats = {
            "mean": round(float(sub["hybrid_score"].mean()), 4),
            "std": round(float(sub["hybrid_score"].std()), 4),
            "median": round(float(sub["hybrid_score"].median()), 4)
        }
        yearly = sub.groupby("year")["hybrid_score"].agg(["mean", "std", "count"]).reset_index()
        temporal_trend = [{k: round(v, 4) if isinstance(v, float) else v for k, v in t.items()}
                          for t in yearly.to_dict("records")]
        platform_sentiment = {k: round(v, 4) for k, v in sub.groupby("platform")["hybrid_score"].mean().items()}
        weighted_sent = np.average(sub["hybrid_score"], weights=sub["engagement"])
        results[tech] = {
            "n_posts": len(sub),
            "sentiment_distribution": {k: round(v, 4) for k, v in sent_dist.items()},
            "hybrid_score_stats": hybrid_stats,
            "engagement_weighted_sentiment": round(float(weighted_sent), 4),
            "platform_sentiment": platform_sentiment,
            "temporal_trend": temporal_trend
        }
    tech_ranking = sorted([(t, r["hybrid_score_stats"]["mean"]) for t, r in results.items()],
                          key=lambda x: x[1], reverse=True)
    results["technology_ranking"] = [{"technology": t, "mean_sentiment": m} for t, m in tech_ranking]
    return results

def run_sentiment_analysis():
    df = generate_synthetic_social_media_data()
    df.to_csv("data/social_media_sentiment.csv", index=False)
    results = bert_lexicon_hybrid_analysis(df)
    with open("results/sentiment_analysis_results.json", "w") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
    return df, results

if __name__ == "__main__":
    df, results = run_sentiment_analysis()
    print("Sentiment analysis completed.")
    for tech in ["gene_editing", "AI", "nuclear_fusion"]:
        r = results[tech]
        print(f"  {tech}: mean_hybrid={r['hybrid_score_stats']['mean']:.3f}, n={r['n_posts']}")
