import sys
import subprocess
# התקנה אוטומטית עוקפת חסימות של adjustText
try:
    from adjustText import adjust_text
except ImportError:
    subprocess.run([sys.executable, "-m", "pip", "install", "adjustText", "--break-system-packages"])
    from adjustText import adjust_text
import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from adjustText import adjust_text

def get_project_paths():
    """מציאת נתיבי הפרויקט בצורה דינמית כדי למנוע נתיבים קשיחים"""
    base_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(base_dir, "data")
    output_dir = os.path.join(base_dir, "outputs")
    
    os.makedirs(output_dir, exist_ok=True)
    return data_dir, output_dir

def load_and_analyze_data():
    data_dir, output_dir = get_project_paths()
    
    # חיפוש דינמי של כל קובצי ה-CSV בתיקיית data
    csv_files = glob.glob(os.path.join(data_dir, "*.csv"))
    
    if not csv_files:
        raise FileNotFoundError("שגיאה: לא נמצאו קובצי נתונים (CSV) בתיקיית data.")
        
    all_dfs = []
    
    # טעינה דינמית של כל חלבון והגדרת השם שלו לפי שם הקובץ
    for file_path in csv_files:
        filename = os.path.basename(file_path)
        protein_name = os.path.splitext(filename)[0].replace("_", " ").title()
        
        df = pd.read_csv(file_path)
        df['protein'] = protein_name
        all_dfs.append(df)
        
    # איחוד כל החלבונים שמצאנו ל-DataFrame אחד
    combined_df = pd.concat(all_dfs, ignore_index=True)
    
    # 2. ניתוח סטטיסטי מורחב
    print(f"--- DNP Enhancement Summary Statistics ({len(csv_files)} Proteins) ---")
    summary = combined_df.groupby('protein')['enhancement'].describe()
    print(summary)
    print("\n-----------------------------------------")
    
    high_enhancement = combined_df[combined_df['enhancement'] > 2.5]
    print(f"Found {len(high_enhancement)} residues with exceptional DNP enhancement (> 2.5):")
    print(high_enhancement[['protein', 'residue', 'enhancement']].to_string(index=False))
    
    # -------------------------------------------------------------
    # 3. ויזואליזציה ושמירת תוצרים אוטומטית (Expanded Analysis Outputs)
    # -------------------------------------------------------------
    sns.set_theme(style="whitegrid")
    
    # --- גרף 1: מפת פיקים דו-ממדית (Scatter Plot) עם תוויות חומצות אמינו מעל 2.0 ---
    plt.figure(figsize=(12, 8))
    sns.scatterplot(
        data=combined_df, 
        x='h', 
        y='c', 
        hue='protein', 
        size='enhancement', 
        sizes=(40, 400), 
        alpha=0.7, 
        palette='Set1'
    )
    
    # הוספת שמות החומצות האמינו עבור פיקים עם Enhancement משמעותי
    for idx, row in combined_df.iterrows():
        if row['enhancement'] > 2.0:
            plt.text(
                x=row['h'] - 0.01,
                y=row['c'] + 0.15,
                s=row['residue'],
                fontdict={'size': 8, 'weight': 'bold'},
                alpha=0.8
            )
            
    plt.title("2D NMR Peak Map ($^1H$ vs $^{13}C$) with Significant Residue Labels")
    plt.xlabel("$^1H$ Chemical Shift (ppm)")
    plt.ylabel("$^{13}C$ Chemical Shift (ppm)")
    plt.gca().invert_xaxis()
    plt.gca().invert_yaxis()
    
    plot1_path = os.path.join(output_dir, "nmr_peak_map.png")
    plt.savefig(plot1_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    # --- גרף 2: ספקטרום NMR מקצועי (Contour Plot) עם כל החומצות אמינו בצורה חכמה ---
    def draw_nmr_contours(ax, df, color, label):
        """מייצרת Contours דמויי NMR עבור כל פיק דו-ממדי ב-DataFrame"""
        h_all = np.arange(-1.0, 3.0, 0.01)
        c_all = np.arange(5.0, 30.0, 0.1)
        H, C = np.meshgrid(h_all, c_all)
        spectrum = np.zeros_like(H)
        
        for _, row in df.iterrows():
            h_peak = row['h']
            c_peak = row['c']
            enh = row['enhancement']
            sigma_h = 0.04
            sigma_c = 0.5
            peak = enh * np.exp(-((H - h_peak)**2 / (2 * sigma_h**2) + (C - c_peak)**2 / (2 * sigma_c**2)))
            spectrum += peak
        
        levels = [df['enhancement'].min() * 0.8, df['enhancement'].mean(), df['enhancement'].max() * 1.5]
        ax.contour(H, C, spectrum, levels=levels, colors=[color], linewidths=1.0, alpha=0.9, linestyles='solid')
        
        from matplotlib.lines import Line2D
        return Line2D([0], [0], color=color, lw=1.5, label=label)

    fig, ax = plt.subplots(figsize=(14, 9))
    sns.set_theme(style="white") # רקע לבן נקי לספקטרום
    
    proxy_artists = []
    df1 = combined_df[combined_df['protein'] == 'Protein 1']
    proxy1 = draw_nmr_contours(ax, df1, color='#FF4500', label='Protein 1 (DNP Enhanced)')
    proxy_artists.append(proxy1)
    
    df2 = combined_df[combined_df['protein'] == 'Protein 2']
    proxy2 = draw_nmr_contours(ax, df2, color='#007ACC', label='Protein 2 (DNP Enhanced)')
    proxy_artists.append(proxy2)
    
    # הוספת התוויות של כל החומצות האמינו
    texts = []
    for idx, row in combined_df.iterrows():
        texts.append(plt.text(row['h'], row['c'], row['residue'], fontsize=9, fontweight='bold'))
    
    # סידור חכם למניעת חפיפה
    adjust_text(texts, arrowprops=dict(arrowstyle='-', color='black', alpha=0.3, lw=0.5))
    
    plt.title("2D $^1H$-$^{13}C$ HMQC Spectrum with DNP Enhancement Contours", fontsize=16)
    plt.xlabel("$^1H$ Chemical Shift (ppm)", fontsize=14)
    plt.ylabel("$^{13}C$ Chemical Shift (ppm)", fontsize=14)
    plt.gca().invert_xaxis()
    plt.gca().invert_yaxis()
    plt.legend(handles=proxy_artists, loc='upper right', frameon=True, shadow=True, title="Protein Sample", fontsize=11)
    
    plot2_path = os.path.join(output_dir, "nmr_spectrum_contours.png")
    plt.savefig(plot2_path, dpi=400, bbox_inches='tight')
    plt.close()
    
    # --- גרף 3: השוואת התפלגות ה-Enhancement (Boxplot) ---
    plt.figure(figsize=(8, 5))
    sns.set_theme(style="whitegrid")
    sns.boxplot(data=combined_df, x='protein', y='enhancement', palette='Pastel1')
    sns.stripplot(data=combined_df, x='protein', y='enhancement', color='black', alpha=0.5, jitter=0.2)
    plt.title("DNP Enhancement Distribution Comparison")
    plt.xlabel("Protein")
    plt.ylabel("Enhancement Factor")
    
    plot3_path = os.path.join(output_dir, "enhancement_comparison.png")
    plt.savefig(plot3_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    # שמירת דוח הסטטיסטיקה לקובץ CSV
    summary.to_csv(os.path.join(output_dir, "dnp_summary_statistics.csv"))
    
    print(f"\n[SUCCESS] Analysis complete! All 3 plots saved to: {output_dir}")

if __name__ == "__main__":
    load_and_analyze_data()