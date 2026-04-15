import pandas as pd
import numpy as np
import os
from dot_chart import plot_dot_chart


def generate_sample_metabolite_data():
    """Generate sample metabolite correlation data for demonstration."""
    np.random.seed(42)
    
    pathways = [
        "Amino Acid Metabolism",
        "Lipid Metabolism", 
        "Carbohydrate Metabolism",
        "Nucleotide Metabolism",
        "Energy Metabolism"
    ]
    
    metabolites = [
        ("Leucine", "Leu", 0),
        ("Isoleucine", "Ile", 0),
        ("Valine", "Val", 0),
        ("Glutamate", "Glu", 0),
        ("Aspartate", "Asp", 0),
        ("Palmitic acid", "PA", 1),
        ("Stearic acid", "SA", 1),
        ("Oleic acid", "OA", 1),
        ("Linoleic acid", "LA", 1),
        ("Glucose", "Glc", 2),
        ("Fructose", "Fru", 2),
        ("Galactose", "Gal", 2),
        ("ATP", "ATP", 4),
        ("ADP", "ADP", 4),
        ("AMP", "AMP", 4),
        ("Adenine", "Ade", 3),
        ("Guanine", "Gua", 3),
        ("Cytosine", "Cyt", 3),
    ]
    
    n_metabolites = len(metabolites)
    
    r_values = np.array([
        0.65, 0.58, 0.72, 0.45, 0.52,
        0.68, 0.61, 0.55, 0.70,
        0.48, 0.53, 0.59,
        0.62, 0.57, 0.64,
        0.49, 0.56, 0.51
    ])
    
    p_values = np.array([
        0.001, 0.003, 0.0005, 0.02, 0.01,
        0.002, 0.004, 0.008, 0.001,
        0.015, 0.009, 0.006,
        0.003, 0.007, 0.002,
        0.012, 0.008, 0.011
    ])
    
    data = {
        'Metabolite_ID': [f"MET{i+1:03d}" for i in range(n_metabolites)],
        'Compound_Name': [name for name, _, _ in metabolites],
        'Compound_Abbr': [abbr for _, abbr, _ in metabolites],
        'SUPER_META_PATHWAY': [pathways[idx] for _, _, idx in metabolites],
        'SUB_META_PATHWAY': [pathways[idx] for _, _, idx in metabolites],
        'r': r_values,
        'p': p_values,
        'n_samples': [150] * n_metabolites,
    }
    
    df = pd.DataFrame(data)
    return df


def generate_sample_protein_data():
    """Generate sample protein correlation data for demonstration."""
    np.random.seed(42)
    
    groups = [
        "Extracellular Matrix",
        "Mitochondria",
        "Cytoplasm",
        "Nucleus",
        "Plasma Membrane"
    ]
    
    proteins = [
        ("Collagen I", "COL1A1", 0),
        ("Collagen III", "COL3A1", 0),
        ("Fibronectin", "FN1", 0),
        ("Laminin", "LAMA1", 0),
        ("NADH Dehydrogenase", "NDUFA1", 1),
        ("Cytochrome C", "CYCS", 1),
        ("ATP Synthase", "ATP5A", 1),
        ("SOD2", "SOD2", 1),
        ("GAPDH", "GAPDH", 2),
        ("ACTB", "ACTB", 2),
        ("TUBB", "TUBB", 2),
        ("HSP90", "HSP90AA1", 2),
        ("Histone H3", "H3F3A", 3),
        ("p53", "TP53", 3),
        ("NF-kB", "NFKB1", 3),
        ("Integrin b1", "ITGB1", 4),
        ("EGFR", "EGFR", 4),
        ("CD44", "CD44", 4),
    ]
    
    n_proteins = len(proteins)
    
    r_values = np.array([
        0.75, 0.68, 0.62, 0.59,
        0.71, 0.66, 0.73, 0.64,
        0.55, 0.58, 0.61, 0.52,
        0.69, 0.63, 0.67,
        0.56, 0.60, 0.54
    ])
    
    p_values = np.array([
        0.0001, 0.0003, 0.001, 0.002,
        0.0002, 0.0005, 0.0001, 0.0008,
        0.005, 0.003, 0.001, 0.009,
        0.0004, 0.0009, 0.0006,
        0.004, 0.002, 0.007
    ])
    
    data = {
        'Protein_ID': [f"PROT{i+1:03d}" for i in range(n_proteins)],
        'Gene_Name': [name for name, _, _ in proteins],
        'UniProt_ID': [uniprot for _, uniprot, _ in proteins],
        'Group': [groups[idx] for _, _, idx in proteins],
        'r': r_values,
        'p': p_values,
        'n_samples': [200] * n_proteins,
    }
    
    df = pd.DataFrame(data)
    return df


def example_metabolite_dot_chart():
    """Example: Plot metabolite dot chart."""
    print("[INFO] Generating sample metabolite data...")
    corr_df = generate_sample_metabolite_data()
    
    protein_name = "APOE"
    outdir = os.path.join("output", "metabolite_dotplot")
    
    print(f"[INFO] Plotting metabolite dot chart for {protein_name}...")
    plot_dot_chart(
        corr_df=corr_df,
        protein_name=protein_name,
        outdir=outdir,
        r_threshold=0.25,
        p_threshold=0.05,
        sex_tag="",
        font_size=8,
        fig_width=6.89,
        mode="metabolite"
    )
    print(f"[SUCCESS] Metabolite dot chart example completed")


def example_protein_dot_chart():
    """Example: Plot protein dot chart."""
    print("[INFO] Generating sample protein data...")
    corr_df = generate_sample_protein_data()
    
    protein_name = "TNF"
    outdir = os.path.join("output", "protein_dotplot")
    
    print(f"[INFO] Plotting protein dot chart for {protein_name}...")
    plot_dot_chart(
        corr_df=corr_df,
        protein_name=protein_name,
        outdir=outdir,
        r_threshold=0.25,
        p_threshold=0.05,
        sex_tag="",
        font_size=8,
        fig_width=6.89,
        mode="protein"
    )
    print(f"[SUCCESS] Protein dot chart example completed")


def example_with_sex_stratification():
    """Example: Plot sex-stratified dot charts."""
    print("[INFO] Generating sample data with sex stratification...")
    
    np.random.seed(42)
    pathways = ["Lipid Metabolism", "Amino Acid Metabolism", "Carbohydrate Metabolism"]
    
    metabolites_data = []
    for i, (name, abbr, pw_idx) in enumerate([
        ("Cholesterol", "Chol", 0),
        ("Triglyceride", "TG", 0),
        ("Phosphatidylcholine", "PC", 0),
        ("Alanine", "Ala", 1),
        ("Glycine", "Gly", 1),
        ("Serine", "Ser", 1),
        ("Pyruvate", "Pyr", 2),
        ("Lactate", "Lac", 2),
    ]):
        metabolites_data.append({
            'Metabolite_ID': f"MET{i+1:03d}",
            'Compound_Name': name,
            'Compound_Abbr': abbr,
            'SUPER_META_PATHWAY': pathways[pw_idx],
            'SUB_META_PATHWAY': pathways[pw_idx],
            'r': np.random.uniform(0.3, 0.8),
            'p': np.random.uniform(0.001, 0.04),
            'n_samples': 100,
        })
    
    corr_df = pd.DataFrame(metabolites_data)
    
    protein_name = "IL6"
    outdir = os.path.join("output", "sex_stratified_dotplot")
    
    print(f"[INFO] Plotting sex-stratified dot chart for {protein_name}...")
    
    for sex in ["", "Male", "Female"]:
        plot_dot_chart(
            corr_df=corr_df,
            protein_name=protein_name,
            outdir=outdir,
            r_threshold=0.25,
            p_threshold=0.05,
            sex_tag=sex,
            font_size=8,
            fig_width=6.89,
            mode="metabolite"
        )
    
    print(f"[SUCCESS] Sex stratification example completed")


if __name__ == "__main__":
    print("=" * 60)
    print("Dot Chart Examples")
    print("=" * 60)
    
    example_metabolite_dot_chart()
    print()
    
    example_protein_dot_chart()
    print()
    
    example_with_sex_stratification()
    
    print("\n" + "=" * 60)
    print("[SUCCESS] All examples completed. Check the 'output' directory for figures.")
    print("=" * 60)
