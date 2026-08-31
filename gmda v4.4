import tkinter as tk
from tkinter import filedialog, messagebox, ttk, scrolledtext
import openpyxl
from openpyxl.utils import get_column_letter
from openpyxl.styles import Font
import re
import os
import math
import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from scipy.cluster.hierarchy import cophenet
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import pandas as pd
from sklearn.metrics import pairwise_distances
from collections import defaultdict
import threading
import time
from packaging import version
from io import StringIO
import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
from sklearn.impute import SimpleImputer, KNNImputer

try:
    from statistics import mode, multimode
except ImportError:
    from statistics import mode
    def multimode(data):
        from collections import Counter
        count = Counter(data)
        max_count = max(count.values())
        return [x for x, c in count.items() if c == max_count]

from scipy import stats
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import string
from sklearn.manifold import MDS
import matplotlib.patches as mpatches
from scipy.spatial import procrustes
from scipy.stats import pearsonr
from scipy.stats import permutation_test
import warnings
from scipy.spatial.distance import cdist

MODERN_THEME = {
    "primary": "#2C3E50",
    "secondary": "#34495E",
    "accent": "#3498DB",
    "success": "#27AE60",
    "warning": "#F39C12",
    "danger": "#E74C3C",
    "light": "#ECF0F1",
    "dark": "#2C3E50",
    "text": "#2C3E50",
    "text_light": "#7F8C8D",
    "bg": "#F8F9FA",
    "card_bg": "#FFFFFF"
}

def mantel(matrix1, matrix2, method='pearson', permutations=999):
    import numpy as np
    from scipy.stats import pearsonr
    from scipy.spatial.distance import squareform

    if matrix1.shape != matrix2.shape:
        min_size = min(matrix1.shape[0], matrix2.shape[0])
        matrix1 = matrix1[:min_size, :min_size]
        matrix2 = matrix2[:min_size, :min_size]

    vec1 = squareform(matrix1, checks=False)
    vec2 = squareform(matrix2, checks=False)

    if len(vec1) != len(vec2):
        min_len = min(len(vec1), len(vec2))
        vec1 = vec1[:min_len]
        vec2 = vec2[:min_len]

    corr_original, _ = pearsonr(vec1, vec2)

    corr_permutations = []
    for _ in range(permutations):
        vec2_perm = vec2.copy()
        np.random.shuffle(vec2_perm)
        corr_perm, _ = pearsonr(vec1, vec2_perm)
        corr_permutations.append(corr_perm)

    p_value = np.mean(np.abs(corr_permutations) >= np.abs(corr_original))

    return corr_original, p_value, np.std(corr_permutations)

warnings.filterwarnings('ignore')


class ModernButton(tk.Button):
    def __init__(self, master=None, **kwargs):
        super().__init__(master, **kwargs)
        self.configure(
            relief="flat",
            bd=0,
            padx=15,
            pady=8,
            font=("Segoe UI", 9),
            cursor="hand2"
        )

class ModernEntry(tk.Entry):
    def __init__(self, master=None, **kwargs):
        super().__init__(master, **kwargs)
        self.configure(
            relief="flat",
            bd=1,
            highlightthickness=1,
            font=("Segoe UI", 9),
            highlightcolor=MODERN_THEME["accent"],
            highlightbackground="#CCCCCC"
        )

class ModernLabel(tk.Label):
    def __init__(self, master=None, **kwargs):
        super().__init__(master, **kwargs)
        self.configure(
            font=("Segoe UI", 9),
            bg=MODERN_THEME["bg"],
            fg=MODERN_THEME["text"]
        )

class ModernFrame(tk.Frame):
    def __init__(self, master=None, **kwargs):
        super().__init__(master, **kwargs)
        self.configure(bg=MODERN_THEME["bg"])

class CardFrame(tk.Frame):
    def __init__(self, master=None, **kwargs):
        super().__init__(master, **kwargs)
        self.configure(
            bg=MODERN_THEME["card_bg"],
            relief="raised",
            bd=1,
            padx=12,
            pady=12
        )


class GMDA:
    def __init__(self, root):
        self.root = root
        self.root.title("GMDA - Genetic and Morphometric Data Analysis")
        self.root.geometry("1400x900")
        self.root.configure(bg=MODERN_THEME["bg"])

        self.setup_styles()

        self.top3_results = []
        self.morpho_results = []
        self.anova_results = []
        self.tukey_results = []
        self.effect_size_results = []
        self.input_file_directory = ""
        self.replicates_var = tk.IntVar(value=3)
        self.morpho_means = []

        # Codominant
        self.genetic_linkage_matrix = None
        self.genetic_distance_matrix = None
        self.genetic_labels = None
        self.genetic_pcoa_coords = None
        self.genetic_cophenetic_correlation = None

        # Morphometric
        self.morpho_linkage_matrix = None
        self.morpho_distance_matrix = None
        self.morpho_labels = None
        self.morpho_pca_coords = None
        self.morpho_cophenetic_correlation = None

        # Dominant
        self.dominant_linkage_matrix = None
        self.dominant_distance_matrix = None
        self.dominant_labels = None
        self.dominant_pcoa_coords = None
        self.dominant_cophenetic_correlation = None

        # Haploid
        self.haploid_linkage_matrix = None
        self.haploid_distance_matrix = None
        self.haploid_labels = None
        self.haploid_pcoa_coords = None
        self.haploid_cophenetic_correlation = None

        # Query-only matrices
        self.genetic_query_distance_matrix = None
        self.genetic_query_labels = None
        self.genetic_query_pcoa_coords = None
        self.genetic_query_cophenetic_correlation = None

        self.dominant_query_distance_matrix = None
        self.dominant_query_labels = None
        self.dominant_query_pcoa_coords = None
        self.dominant_query_cophenetic_correlation = None

        self.morpho_query_distance_matrix = None
        self.morpho_query_labels = None
        self.morpho_query_pca_coords = None
        self.morpho_query_cophenetic_correlation = None

        self.haploid_query_distance_matrix = None
        self.haploid_query_labels = None
        self.haploid_query_pcoa_coords = None
        self.haploid_query_cophenetic_correlation = None

        self.dominant_top3_results = []
        self.haploid_top3_results = []

        self.setup_ui()
        self.setup_notebook()

    def setup_styles(self):
        style = ttk.Style()
        style.theme_use('clam')

        style.configure("TNotebook", background=MODERN_THEME["bg"])
        style.configure("TNotebook.Tab",
                       font=("Segoe UI", 9, "bold"),
                       padding=[15, 8],
                       background=MODERN_THEME["light"],
                       foreground=MODERN_THEME["text"])
        style.map("TNotebook.Tab",
                 background=[("selected", MODERN_THEME["accent"])],
                 foreground=[("selected", "white")])

        style.configure("TProgressbar",
                       thickness=20,
                       background=MODERN_THEME["success"])

    def setup_ui(self):
        header_frame = ModernFrame(self.root)
        header_frame.pack(fill="x", padx=20, pady=(15, 10))

        self.title_label = tk.Label(header_frame,
                                  text="GMDA 4.0",
                                  font=("Segoe UI", 24, "bold"),
                                  bg=MODERN_THEME["bg"],
                                  fg=MODERN_THEME["primary"])
        self.title_label.pack(pady=(0, 5))

        self.subtitle_label = tk.Label(header_frame,
                                     text="Genetic and Morphometric Data Analysis Platform",
                                     font=("Segoe UI", 11),
                                     bg=MODERN_THEME["bg"],
                                     fg=MODERN_THEME["text_light"])
        self.subtitle_label.pack(pady=(0, 10))

        self.status_frame = ModernFrame(self.root)
        self.status_frame.pack(fill="x", side="bottom", padx=20, pady=5)

        self.status_label = tk.Label(self.status_frame,
                                   text="Ready",
                                   font=("Segoe UI", 9),
                                   bg=MODERN_THEME["bg"],
                                   fg=MODERN_THEME["text_light"])
        self.status_label.pack(side="left")

        self.version_label = tk.Label(self.status_frame,
                                    text="v4.0",
                                    font=("Segoe UI", 9),
                                    bg=MODERN_THEME["bg"],
                                    fg=MODERN_THEME["text_light"])
        self.version_label.pack(side="right")

    def update_status(self, message):
        self.status_label.config(text=message)
        self.root.update_idletasks()

    def setup_notebook(self):
        self.notebook = ttk.Notebook(self.root)

        self.genetic_codominant_frame = ModernFrame(self.notebook)
        self.genetic_dominant_frame = ModernFrame(self.notebook)
        self.genetic_haploid_frame = ModernFrame(self.notebook)
        self.morpho_frame = ModernFrame(self.notebook)
        self.correlation_frame = ModernFrame(self.notebook)
        self.outlier_frame = ModernFrame(self.notebook)

        self.notebook.add(self.genetic_codominant_frame, text="🧬 Genetic Analysis - Codominant Markers")
        self.notebook.add(self.genetic_dominant_frame, text="🧬 Genetic Analysis - Dominant Markers")
        self.notebook.add(self.genetic_haploid_frame, text="🧬 Genetic Analysis - Haploid Markers")
        self.notebook.add(self.morpho_frame, text="📊 Morphometric Analysis")
        self.notebook.add(self.correlation_frame, text="🔗 Correlation Analysis")
        self.notebook.add(self.outlier_frame, text="📈 Outlier Detection")

        self.notebook.pack(expand=1, fill="both", padx=20, pady=10)

        self.setup_genetic_codominant_tab()
        self.setup_genetic_dominant_tab()
        self.setup_genetic_haploid_tab()
        self.setup_morpho_tab()
        self.setup_correlation_tab()
        self.setup_outlier_tab()

    # ==================== COPHENETIC CORRELATION FUNCTION ====================
    
    def calculate_cophenetic_correlation(self, distance_matrix_condensed, linkage_matrix):
        """
        Calculate the cophenetic correlation coefficient between original distances
        and cophenetic distances from the dendrogram.
        
        Args:
            distance_matrix_condensed: condensed distance matrix from original data
            linkage_matrix: linkage matrix from hierarchical clustering
            
        Returns:
            cophenetic_correlation: Pearson correlation between original and cophenetic distances
        """
        try:
            # Calculate cophenetic distances from the linkage matrix
            cophenetic_distances = cophenet(linkage_matrix, distance_matrix_condensed)[1]
            
            # Calculate Pearson correlation
            correlation = np.corrcoef(distance_matrix_condensed, cophenetic_distances)[0, 1]
            
            return correlation
        except Exception as e:
            print(f"Error calculating cophenetic correlation: {e}")
            return None

    # ==================== UTILITY FUNCTIONS ====================

    def read_genetic_file(self, filename):
        try:
            df = pd.read_excel(filename, engine='openpyxl', header=0)
            
            if len(df.columns) < 2:
                raise ValueError("The file must have at least 2 columns (individual + loci)")
            
            sample_col = df.columns[0]
            df = df.dropna(how='all')
            df = df[df[sample_col].notna()]
            df = df[df[sample_col].astype(str).str.strip() != '']
            
            if df.empty:
                raise ValueError("No valid individuals were found in the file")
            
            return df, sample_col
            
        except Exception as e:
            raise Exception(f"Erro ao ler arquivo {filename}: {str(e)}")

    def is_missing_codominant(self, value):
        try:
            if value is None or value == '':
                return True
            if isinstance(value, float) and np.isnan(value):
                return True
            if pd.isna(value):
                return True
            return False
        except (TypeError, ValueError):
            return False

    def is_missing_dominant(self, value):
        try:
            if value is None or value == '':
                return True
            if isinstance(value, float) and np.isnan(value):
                return True
            if pd.isna(value):
                return True
            return False
        except (TypeError, ValueError):
            return False

    def is_missing_haploid(self, value):
        try:
            if value is None or value == '':
                return True
            if isinstance(value, float) and np.isnan(value):
                return True
            if pd.isna(value):
                return True
            return False
        except (TypeError, ValueError):
            return False

    def convert_to_numeric(self, value):
        try:
            return float(value)
        except (ValueError, TypeError):
            return value

    def calcular_similaridade_codominante(self, query, reference):
        total_similarity = 0
        valid_loci_count = 0

        for q, r in zip(query, reference):
            q_valid = not self.is_missing_codominant(q)
            r_valid = not self.is_missing_codominant(r)

            if q_valid and r_valid:
                try:
                    q_val = float(q)
                    r_val = float(r)
                    if r_val != 0:
                        loco_similarity = max(0, 1 - abs(q_val - r_val) / r_val)
                    else:
                        loco_similarity = 0
                    total_similarity += loco_similarity
                    valid_loci_count += 1
                except (ValueError, TypeError):
                    if str(q) == str(r):
                        total_similarity += 1
                    else:
                        total_similarity += 0
                    valid_loci_count += 1

        if valid_loci_count == 0:
            return 0.0

        return total_similarity / valid_loci_count

    def calcular_similaridade_dominante(self, query, reference):
        a, b, c = 0, 0, 0

        for q, r in zip(query, reference):
            q_valid = not self.is_missing_dominant(q)
            r_valid = not self.is_missing_dominant(r)

            if q_valid and r_valid:
                try:
                    q_val = float(q)
                    r_val = float(r)
                except (ValueError, TypeError):
                    q_str = str(q).upper().strip()
                    r_str = str(r).upper().strip()
                    q_val = 1 if q_str in ['1', 'TRUE', 'YES', 'PRESENT', 'P', '+'] else 0
                    r_val = 1 if r_str in ['1', 'TRUE', 'YES', 'PRESENT', 'P', '+'] else 0

                if q_val == 1 and r_val == 1:
                    a += 1
                elif q_val == 1 and r_val == 0:
                    b += 1
                elif q_val == 0 and r_val == 1:
                    c += 1

        if (a + b + c) == 0:
            return 0.0

        return a / (a + b + c)

    def calcular_similaridade_haploid(self, query, reference):
        total_similarity = 0
        valid_loci_count = 0

        for q, r in zip(query, reference):
            q_valid = not self.is_missing_haploid(q)
            r_valid = not self.is_missing_haploid(r)

            if q_valid and r_valid:
                if str(q).strip() == str(r).strip():
                    total_similarity += 1
                else:
                    total_similarity += 0
                valid_loci_count += 1

        if valid_loci_count == 0:
            return 0.0

        return total_similarity / valid_loci_count

    def load_file(self, entry_widget):
        filename = filedialog.askopenfilename(
            title="Select File",
            filetypes=[("Excel files", "*.xlsx *.xls"), ("All files", "*.*")]
        )
        if filename:
            entry_widget.delete(0, tk.END)
            entry_widget.insert(0, filename)
            self.input_file_directory = os.path.dirname(filename)
            self.update_status(f"Loaded: {os.path.basename(filename)}")

    @staticmethod
    def isanumber(a):
        try:
            float(a)
            return not np.isnan(float(a))
        except (ValueError, TypeError):
            return False

    def get_cluster_members(self, linkage_matrix, n_samples):
        clusters = {}
        for i, (c1, c2, _, _) in enumerate(linkage_matrix):
            c1, c2 = int(c1), int(c2)
            members = set()
            if c1 < n_samples:
                members.add(c1)
            else:
                members.update(clusters.get(c1, set()))
            if c2 < n_samples:
                members.add(c2)
            else:
                members.update(clusters.get(c2, set()))
            clusters[i + n_samples] = members
        return clusters

    def create_output_directory(self, analysis_type):
        base_dir = self.input_file_directory
        if not base_dir:
            base_dir = os.getcwd()

        dir_name = f"GMDA_{analysis_type}_Outputs"
        output_dir = os.path.join(base_dir, dir_name)

        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        return output_dir

    # ==================== CODOMINANT TAB ====================

    def setup_genetic_codominant_tab(self):
        main_container = ModernFrame(self.genetic_codominant_frame)
        main_container.pack(fill="both", expand=True, padx=10, pady=10)

        input_card = CardFrame(main_container)
        input_card.pack(fill="x", pady=(0, 15))

        input_title = ModernLabel(input_card, text="Input Files (Format: Row1=locus names, Col1=sample names)", font=("Segoe UI", 12, "bold"))
        input_title.pack(anchor="w", pady=(0, 15))

        db_frame = ModernFrame(input_card)
        db_frame.pack(fill="x", pady=8)

        label_database = ModernLabel(db_frame, text="Reference Database:", font=("Segoe UI", 10, "bold"))
        label_database.grid(row=0, column=0, sticky="w", padx=(0, 10))

        self.entry_database = ModernEntry(db_frame, width=50)
        self.entry_database.grid(row=0, column=1, padx=5, sticky="ew")

        btn_browse_db = ModernButton(db_frame, text="Browse",
                                   command=lambda: self.load_file(self.entry_database),
                                   bg=MODERN_THEME["accent"], fg="white")
        btn_browse_db.grid(row=0, column=2, padx=5)

        query_frame = ModernFrame(input_card)
        query_frame.pack(fill="x", pady=8)

        label_query = ModernLabel(query_frame, text="Query Samples:", font=("Segoe UI", 10, "bold"))
        label_query.grid(row=0, column=0, sticky="w", padx=(0, 10))

        self.entry_query = ModernEntry(query_frame, width=50)
        self.entry_query.grid(row=0, column=1, padx=5, sticky="ew")

        btn_browse_query = ModernButton(query_frame, text="Browse",
                                      command=lambda: self.load_file(self.entry_query),
                                      bg=MODERN_THEME["accent"], fg="white")
        btn_browse_query.grid(row=0, column=2, padx=5)

        options_frame = ModernFrame(input_card)
        options_frame.pack(fill="x", pady=10)

        self.bootstrap_var_gen = tk.BooleanVar()
        check_bootstrap = tk.Checkbutton(options_frame, text="Use Bootstrap",
                                       variable=self.bootstrap_var_gen,
                                       bg=MODERN_THEME["card_bg"],
                                       font=("Segoe UI", 10),
                                       selectcolor=MODERN_THEME["light"])
        check_bootstrap.pack(anchor="w")

        db_frame.columnconfigure(1, weight=1)
        query_frame.columnconfigure(1, weight=1)

        action_card = CardFrame(main_container)
        action_card.pack(fill="x", pady=(0, 15))

        action_title = ModernLabel(action_card, text="Analysis Tools", font=("Segoe UI", 12, "bold"))
        action_title.pack(anchor="w", pady=(0, 15))

        btn_frame = ModernFrame(action_card)
        btn_frame.pack(fill="x")

        btn_top3 = ModernButton(btn_frame,
                              text="🧬 Generate Similarity Matrix & Dendrogram",
                              command=self.calculate_genetic,
                              bg=MODERN_THEME["success"], fg="white")
        btn_top3.grid(row=0, column=0, padx=5, pady=5, sticky="ew")

        btn_save = ModernButton(btn_frame,
                              text="💾 Save Top 3 Results",
                              command=self.save_genetic_results,
                              bg=MODERN_THEME["accent"], fg="white")
        btn_save.grid(row=0, column=1, padx=5, pady=5, sticky="ew")

        btn_query_analysis = ModernButton(btn_frame,
                                         text="📈 Query Data Analysis",
                                         command=self.calculate_genetic_indices,
                                         bg=MODERN_THEME["primary"], fg="white")
        btn_query_analysis.grid(row=1, column=0, columnspan=2, padx=5, pady=5, sticky="ew")

        btn_frame.columnconfigure(0, weight=1)
        btn_frame.columnconfigure(1, weight=1)

        progress_frame = ModernFrame(action_card)
        progress_frame.pack(fill="x", pady=(10, 0))

        self.progress_label_gen = ModernLabel(progress_frame, text="Progress:", font=("Segoe UI", 9, "bold"))
        self.progress_label_gen.pack(anchor="w")

        self.progress_bar_gen = ttk.Progressbar(progress_frame, mode='determinate', length=100)
        self.progress_bar_gen.pack(fill="x", pady=(5, 0))

        results_card = CardFrame(main_container)
        results_card.pack(fill="both", expand=True)

        results_title = ModernLabel(results_card, text="Results", font=("Segoe UI", 12, "bold"))
        results_title.pack(anchor="w", pady=(0, 10))

        text_frame = ModernFrame(results_card)
        text_frame.pack(fill="both", expand=True)

        self.results_text_gen = scrolledtext.ScrolledText(text_frame,
                                      bg="white",
                                      fg=MODERN_THEME["text"],
                                      font=("Consolas", 9),
                                      relief="flat",
                                      padx=10,
                                      pady=10)
        self.results_text_gen.pack(fill="both", expand=True)

    def calculate_genetic(self):
        self.update_status("Calculating genetic similarity for codominant markers...")
        self.progress_bar_gen['value'] = 0
        self.progress_bar_gen['maximum'] = 100
        self.root.update_idletasks()

        db_filename = self.entry_database.get().strip()
        query_filename = self.entry_query.get().strip()

        if not db_filename or not query_filename:
            messagebox.showwarning("Warning", "Please select both files.")
            self.update_status("Ready")
            return

        try:
            db_df, db_name_col = self.read_genetic_file(db_filename)
            query_df, query_name_col = self.read_genetic_file(query_filename)
        except Exception as e:
            messagebox.showerror("Error", f"Error reading files:\n{e}")
            self.update_status("Error reading files")
            return

        if db_df.empty or query_df.empty:
            messagebox.showerror("Error", "One or both files are empty.")
            self.update_status("Empty files")
            return

        self.progress_bar_gen['value'] = 10
        self.root.update_idletasks()

        output_dir = self.create_output_directory("Genetic_codominant")

        self.results_text_gen.delete(1.0, tk.END)
        self.results_text_gen.insert(tk.END, "🧬 GENETIC SIMILARITY - CODOMINANT MARKERS\n")
        self.results_text_gen.insert(tk.END, "=" * 60 + "\n\n")

        self.top3_results = []
        total_queries = len(query_df)
        queries_processed = 0

        for idx_query, query_row in query_df.iterrows():
            query_name = query_row[query_name_col]
            if pd.isna(query_name) or str(query_name).strip() == '':
                continue

            query_data = query_row.drop(query_name_col).values
            similarities = {}

            for idx_ref, ref_row in db_df.iterrows():
                ref_name = ref_row[db_name_col]
                if pd.isna(ref_name) or str(ref_name).strip() == '':
                    continue
                ref_data = ref_row.drop(db_name_col).values
                similarity = self.calcular_similaridade_codominante(query_data, ref_data)
                similarities[ref_name] = similarity

            if not similarities:
                continue

            top3 = sorted(similarities.items(), key=lambda x: x[1], reverse=True)[:3]
            self.results_text_gen.insert(tk.END, f"🔍 Query sample: {query_name}\n")
            for ref_name, score in top3:
                self.results_text_gen.insert(tk.END, f"   📌 {ref_name}: {score * 100:.2f}%\n")
                self.top3_results.append((query_name, ref_name, score))

            self.top3_results.append(("", "", ""))
            queries_processed += 1
            progress = 10 + (queries_processed / total_queries * 40)
            self.progress_bar_gen['value'] = progress
            self.root.update_idletasks()

        try:
            self.progress_bar_gen['value'] = 60
            self.root.update_idletasks()
            self.generate_genetic_matrix_and_dendrogram(
                query_df, db_df, query_name_col, db_name_col,
                self.bootstrap_var_gen.get(), output_dir)
        except Exception as e:
            self.results_text_gen.insert(tk.END, f"\n❌ Error generating matrix/dendrogram: {e}\n")

        self.progress_bar_gen['value'] = 100
        self.root.update_idletasks()
        self.update_status("Genetic analysis completed")

    def generate_genetic_matrix_and_dendrogram(self, query_df, db_df, query_name_col,
                                                db_name_col, bootstrap_enabled, output_dir):
        combined_df = pd.concat([query_df, db_df], ignore_index=True)
        name_col = query_name_col if query_name_col in combined_df.columns else db_name_col
        combined_df = combined_df.dropna(subset=[name_col])

        if combined_df.empty:
            self.results_text_gen.insert(tk.END, "❌ Error: No valid samples found.\n")
            return

        names = []
        for name in combined_df.iloc[:, 0]:
            if pd.notna(name):
                names.append(str(name).strip())
            else:
                names.append(None)

        data = combined_df.iloc[:, 1:]
        valid_indices = [i for i, name in enumerate(names) if name and name != '']
        names = [names[i] for i in valid_indices]
        data = data.iloc[valid_indices]

        if len(names) < 2:
            self.results_text_gen.insert(tk.END, "❌ Error: At least 2 samples required.\n")
            return

        missing_by_sample = {}
        for i, sample in enumerate(names):
            sample_data = data.iloc[i]
            n_missing = sum(1 for val in sample_data if self.is_missing_codominant(val))
            if n_missing > 0:
                missing_by_sample[sample] = n_missing

        if missing_by_sample:
            self.results_text_gen.insert(tk.END,
                f"\n⚠️ Missing data detected in {len(missing_by_sample)} samples\n")
            for sample, count in missing_by_sample.items():
                self.results_text_gen.insert(tk.END,
                    f"   • {sample}: {count} empty cells\n")

        matrix = pd.DataFrame(index=names, columns=names, dtype=float)
        np.fill_diagonal(matrix.values, 1.0)

        total_pairs = len(data) * (len(data) - 1) / 2
        pairs_done = 0

        for i in range(len(data)):
            for j in range(i + 1, len(data)):
                sim = self.calcular_similaridade_codominante(data.iloc[i].values, data.iloc[j].values)
                matrix.iat[i, j] = sim
                matrix.iat[j, i] = sim
                pairs_done += 1
                progress = 60 + (pairs_done / total_pairs * 20)
                self.progress_bar_gen['value'] = progress
                self.root.update_idletasks()

        matrix_filename = os.path.join(
            output_dir, "codominant_genetic_similarity_matrix.xlsx")
        try:
            with pd.ExcelWriter(matrix_filename, engine='openpyxl') as writer:
                title_df = pd.DataFrame({"": ["GMDA - Genetic Similarity Matrix"]})
                title_df.to_excel(writer, sheet_name='Similarity Matrix',
                                  index=False, header=False)
                matrix.to_excel(writer, sheet_name='Similarity Matrix', startrow=2)
            self.results_text_gen.insert(tk.END,
                f"\n💾 Similarity matrix saved as: {matrix_filename}\n")
        except Exception as e:
            self.results_text_gen.insert(tk.END, f"\n❌ Error saving matrix: {e}\n")

        distance_matrix = 1 - matrix.astype(float).values
        np.fill_diagonal(distance_matrix, 0)

        if np.any(np.isnan(distance_matrix)):
            self.results_text_gen.insert(tk.END,
                f"\n⚠️ Fixing distance matrix (replacing NaN with 1.0)...\n")
            distance_matrix = np.nan_to_num(distance_matrix, nan=1.0)

        try:
            condensed_dist = squareform(distance_matrix, checks=False)
            linkage_matrix = linkage(condensed_dist, method='average')

            self.genetic_linkage_matrix = linkage_matrix
            self.genetic_distance_matrix = condensed_dist
            self.genetic_labels = names.copy()
            
            # Calculate cophenetic correlation
            cophenetic_corr = self.calculate_cophenetic_correlation(condensed_dist, linkage_matrix)
            self.genetic_cophenetic_correlation = cophenetic_corr
            
            self.results_text_gen.insert(tk.END,
                f"\n📊 Cophenetic Correlation Coefficient: {cophenetic_corr:.4f}\n")
            self.results_text_gen.insert(tk.END,
                f"   (Measures how well the dendrogram represents the original distances)\n")
            
            if cophenetic_corr > 0.8:
                self.results_text_gen.insert(tk.END,
                    f"   ✅ Excellent representation (c > 0.8)\n")
            elif cophenetic_corr > 0.7:
                self.results_text_gen.insert(tk.END,
                    f"   ✓ Good representation (c > 0.7)\n")
            elif cophenetic_corr > 0.6:
                self.results_text_gen.insert(tk.END,
                    f"   ⚠️ Moderate representation (c > 0.6)\n")
            else:
                self.results_text_gen.insert(tk.END,
                    f"   ❌ Poor representation (c ≤ 0.6)\n")

            self.progress_bar_gen['value'] = 85
            self.root.update_idletasks()

            if bootstrap_enabled:
                self.results_text_gen.insert(tk.END,
                    "🔄 Calculating bootstrap supports...\n")
                self.root.update_idletasks()
                supports = self.calculate_bootstrap_support(data, linkage_matrix)
                self.plot_genetic_dendrogram(
                    linkage_matrix, names, supports, output_dir, missing_by_sample, cophenetic_corr)
            else:
                self.plot_genetic_dendrogram(
                    linkage_matrix, names, output_dir=output_dir,
                    missing_info=missing_by_sample, cophenetic_corr=cophenetic_corr)

            self.progress_bar_gen['value'] = 95
            self.root.update_idletasks()

        except Exception as e:
            self.results_text_gen.insert(tk.END, f"❌ Error computing dendrogram: {e}\n")

    def calculate_bootstrap_support(self, data, linkage_matrix, replicates=100):
        n = data.shape[0]
        clusters_suporte = defaultdict(int)
        original_clusters = self.get_cluster_members(linkage_matrix, n)

        for rep in range(replicates):
            try:
                sampled_cols = np.random.choice(
                    data.shape[1], size=data.shape[1], replace=True)
                data_boot = data.iloc[:, sampled_cols]

                matrix = pd.DataFrame(
                    index=range(len(data_boot)),
                    columns=range(len(data_boot)), dtype=float)
                np.fill_diagonal(matrix.values, 1.0)

                for i1 in range(len(data_boot)):
                    for i2 in range(i1 + 1, len(data_boot)):
                        sim = self.calcular_similaridade_codominante(
                            data_boot.iloc[i1].values, data_boot.iloc[i2].values)
                        matrix.iat[i1, i2] = sim
                        matrix.iat[i2, i1] = sim

                dist = 1 - matrix.values
                if np.any(np.isnan(dist)):
                    dist = np.nan_to_num(dist, nan=1.0)

                condensed_dist = squareform(dist, checks=False)
                linkage_boot = linkage(condensed_dist, method='average')
                bootstrap_clusters = self.get_cluster_members(linkage_boot, n)

                for cluster_id, members in original_clusters.items():
                    if any(members == b_members
                           for b_members in bootstrap_clusters.values()):
                        clusters_suporte[cluster_id] += 1

                progress = 85 + (rep / replicates * 10)
                self.progress_bar_gen['value'] = progress
                self.root.update_idletasks()

            except Exception:
                continue

        for c in clusters_suporte:
            clusters_suporte[c] = clusters_suporte[c] / replicates * 100

        return clusters_suporte

    def plot_genetic_dendrogram(self, linkage_matrix, labels, supports=None,
                                output_dir=None, missing_info=None, cophenetic_corr=None):
        if output_dir is None:
            output_dir = self.create_output_directory("Genetic_codominant")

        plt.figure(figsize=(12, 8))
        dendro = dendrogram(linkage_matrix, labels=labels,
                            leaf_rotation=90, leaf_font_size=10)

        if supports:
            icoord = np.array(dendro['icoord'])
            dcoord = np.array(dendro['dcoord'])
            for i, d in zip(icoord, dcoord):
                x = 0.5 * (i[1] + i[2])
                y = d[1]
                for cluster_idx, row in enumerate(linkage_matrix):
                    if abs(row[2] - y) < 1e-5:
                        support = supports.get(cluster_idx + len(labels), 0)
                        if support >= 50:
                            plt.text(x, y, f"{int(support)}%",
                                     ha='center', va='bottom',
                                     fontsize=8, color='red')
                        break

        title = "GMDA - Genetic UPGMA Dendrogram"
        if supports:
            title += " (with bootstrap)"
        if cophenetic_corr is not None:
            title += f"\nCophenetic Correlation: {cophenetic_corr:.4f}"
        plt.title(title, fontsize=14)
        plt.ylabel("Dissimilarity (1-S)", fontsize=14)
        plt.tight_layout()

        filename = "codominant_genetic_dendrogram"
        if supports:
            filename += "_with_bootstrap"
        plot_path = os.path.join(output_dir, f"{filename}.png")
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        plt.close()
        self.results_text_gen.insert(tk.END, f"📊 Dendrogram saved as: {plot_path}\n")

    def save_genetic_results(self):
        if not self.top3_results:
            messagebox.showinfo("No data", "No results available to save.")
            return

        output_dir = self.create_output_directory("Genetic_codominant")
        df_results = pd.DataFrame(
            self.top3_results, columns=["Query", "Reference", "Similarity"])
        df_results["Similarity (%)"] = pd.to_numeric(
            df_results["Similarity"], errors='coerce') * 100
        df_results.drop(columns=["Similarity"], inplace=True)

        save_path = os.path.join(output_dir, "Top3_codominant.tsv")
        try:
            with open(save_path, 'w', encoding='utf-8') as f:
                f.write("GMDA - Genetic Similarity Results\n\n")
                df_results.to_csv(f, sep="\t", index=False)
            self.results_text_gen.insert(tk.END,
                f"\n💾 Results saved as TSV: {save_path}\n")
        except Exception as e:
            self.results_text_gen.insert(tk.END, f"❌ Error saving TSV: {e}\n")

    def calculate_genetic_indices(self):
        self.update_status("Calculating genetic indices...")
        self.progress_bar_gen['value'] = 0
        self.progress_bar_gen['maximum'] = 100
        self.root.update_idletasks()

        query_filename = self.entry_query.get().strip()
        if not query_filename:
            messagebox.showwarning("Warning", "Please select the query samples file.")
            self.update_status("Ready")
            return

        try:
            query_df, query_name_col = self.read_genetic_file(query_filename)
        except Exception as e:
            messagebox.showerror("Error", f"Error reading file:\n{e}")
            self.update_status("Error reading file")
            return

        if query_df.empty:
            messagebox.showerror("Error", "The query file is empty.")
            self.update_status("Empty file")
            return

        self.progress_bar_gen['value'] = 5
        self.root.update_idletasks()

        output_dir = self.create_output_directory("Genetic_codominant")

        sample_results = []
        pop_A, pop_Ae, pop_He, pop_Ho, pop_F = [], [], [], [], []

        total_missing_loci = 0
        total_valid_loci = 0
        missing_by_sample = {}
        missing_by_locus = defaultdict(int)

        total_samples = sum(
            1 for _, row in query_df.iterrows()
            if not (pd.isna(row[query_name_col]) or
                    str(row[query_name_col]).strip() == ''))

        self.progress_bar_gen['value'] = 10
        self.root.update_idletasks()

        samples_processed = 0

        for idx, row in query_df.iterrows():
            sample = row[query_name_col]
            if pd.isna(sample) or str(sample).strip() == '':
                continue

            alleles = []
            heterozygotes = 0
            total_loci = 0
            sample_missing = 0

            data_cols = [c for c in query_df.columns if c != query_name_col]
            for locus_idx in range(0, len(data_cols) - 1, 2):
                col1 = data_cols[locus_idx]
                col2 = data_cols[locus_idx + 1] if locus_idx + 1 < len(data_cols) else None
                locus1 = row[col1]
                locus2 = row[col2] if col2 else None

                locus1_valid = not self.is_missing_codominant(locus1)
                locus2_valid = not self.is_missing_codominant(locus2) if col2 else False
                locus_num = (locus_idx // 2) + 1

                if locus1_valid:
                    try:
                        locus1_val = float(locus1)
                        alleles.append((locus_num, locus1_val))
                        total_valid_loci += 1
                    except (ValueError, TypeError):
                        alleles.append((locus_num, str(locus1)))
                        total_valid_loci += 1
                else:
                    sample_missing += 1
                    total_missing_loci += 1
                    missing_by_locus[f"Locus_{locus_num}"] += 1

                if locus2_valid and col2:
                    try:
                        locus2_val = float(locus2)
                        alleles.append((locus_num, locus2_val))
                        total_valid_loci += 1
                    except (ValueError, TypeError):
                        alleles.append((locus_num, str(locus2)))
                        total_valid_loci += 1
                elif col2:
                    sample_missing += 1
                    total_missing_loci += 1
                    missing_by_locus[f"Locus_{locus_num}"] += 1

                if locus1_valid and locus2_valid and col2:
                    try:
                        locus1_val = float(locus1)
                        locus2_val = float(locus2)
                        if locus1_val != locus2_val:
                            heterozygotes += 1
                    except (ValueError, TypeError):
                        if str(locus1) != str(locus2):
                            heterozygotes += 1
                    total_loci += 1

            if sample_missing > 0:
                missing_by_sample[sample] = sample_missing

            if not alleles:
                sample_results.append({
                    'Sample': sample, 'A': 0, 'Ho': 0.0,
                    'Missing_loci': sample_missing})
                continue

            unique_alleles = set(alleles)
            A = len(unique_alleles)
            Ho = heterozygotes / total_loci if total_loci > 0 else 0

            sample_results.append({
                'Sample': sample, 'A': A, 'Ho': round(Ho, 3),
                'Missing_loci': sample_missing})

            allele_counts = defaultdict(int)
            for allele in alleles:
                allele_counts[allele] += 1

            freqs = np.array([count / len(alleles) for count in allele_counts.values()])
            p_squared_sum = np.sum(freqs ** 2)

            Ae = 1 / p_squared_sum if p_squared_sum > 0 else 0
            He = 1 - p_squared_sum
            F = (He - Ho) / He if He > 0 else 0

            pop_A.append(A)
            pop_Ae.append(Ae)
            pop_He.append(He)
            pop_Ho.append(Ho)
            pop_F.append(F)

            samples_processed += 1
            progress = 10 + (samples_processed / total_samples * 30)
            self.progress_bar_gen['value'] = progress
            self.root.update_idletasks()

        if not sample_results:
            messagebox.showwarning("Warning", "No valid samples found.")
            self.update_status("No valid samples")
            return

        self.results_text_gen.insert(tk.END, "\n" + "=" * 60 + "\n")
        self.results_text_gen.insert(tk.END, "📊 MISSING DATA STATISTICS\n")
        self.results_text_gen.insert(tk.END, "=" * 60 + "\n")

        total_alleles = total_missing_loci + total_valid_loci
        missing_percentage = (total_missing_loci / total_alleles) * 100 \
            if total_alleles > 0 else 0
        self.results_text_gen.insert(tk.END,
            f"Total valid alleles: {total_valid_loci}\n")
        self.results_text_gen.insert(tk.END,
            f"Total empty cells: {total_missing_loci}\n")
        self.results_text_gen.insert(tk.END,
            f"Missing data percentage: {missing_percentage:.2f}%\n\n")

        if missing_by_sample:
            self.results_text_gen.insert(tk.END, "📊 Missing data by sample:\n")
            for sample, count in missing_by_sample.items():
                self.results_text_gen.insert(tk.END,
                    f"  • {sample}: {count} empty cells\n")

        if missing_by_locus:
            self.results_text_gen.insert(tk.END, "\n📊 Missing data by locus:\n")
            for locus, count in missing_by_locus.items():
                self.results_text_gen.insert(tk.END,
                    f"  • {locus}: {count} occurrences\n")

        self.results_text_gen.insert(tk.END, "\n" + "=" * 60 + "\n\n")

        mean_A = round(np.mean(pop_A), 3) if pop_A else 0
        mean_Ae = round(np.mean(pop_Ae), 3) if pop_Ae else 0
        mean_He = round(np.mean(pop_He), 3) if pop_He else 0
        mean_Ho = round(np.mean(pop_Ho), 3) if pop_Ho else 0
        mean_F = round(np.mean(pop_F), 3) if pop_F else 0

        df_samples = pd.DataFrame(sample_results)
        df_population = pd.DataFrame({
            'A': [mean_A], 'Ae': [mean_Ae], 'Ho': [mean_Ho],
            'He': [mean_He], 'F': [mean_F],
            'Missing_percentage': [round(missing_percentage, 2)]
        })

        tsv_path = os.path.join(output_dir, "codominant_genetic_indices.tsv")
        try:
            with open(tsv_path, 'w', encoding='utf-8') as f:
                f.write("GMDA - Genetic Indices (missing data excluded)\n\n")
                f.write("Samples\n")
                df_samples.to_csv(f, sep="\t", index=False)
                f.write("\nPopulation\n")
                df_population.to_csv(f, sep=" ", index=False)
        except Exception as e:
            messagebox.showerror("Error", f"Error saving TSV:\n{e}")
            self.update_status("Error saving results")
            return

        self.progress_bar_gen['value'] = 50
        self.root.update_idletasks()

        try:
            query_names = []
            for name in query_df.iloc[:, 0]:
                if pd.notna(name):
                    query_names.append(str(name).strip())
                else:
                    query_names.append(None)

            query_data = query_df.iloc[:, 1:]
            valid_indices = [i for i, name in enumerate(query_names)
                             if name and name != '']
            query_names = [query_names[i] for i in valid_indices]
            query_data = query_data.iloc[valid_indices]

            if len(query_names) < 2:
                self.results_text_gen.insert(tk.END,
                    "\n⚠️ At least 2 Query samples needed for dendrogram and PCoA\n")
                return

            matrix = pd.DataFrame(index=query_names, columns=query_names, dtype=float)
            np.fill_diagonal(matrix.values, 1.0)

            total_pairs = len(query_data) * (len(query_data) - 1) / 2
            pairs_done = 0

            for i in range(len(query_data)):
                for j in range(i + 1, len(query_data)):
                    sim = self.calcular_similaridade_codominante(
                        query_data.iloc[i].values, query_data.iloc[j].values)
                    matrix.iat[i, j] = sim
                    matrix.iat[j, i] = sim
                    pairs_done += 1
                    progress = 50 + (pairs_done / total_pairs * 20)
                    self.progress_bar_gen['value'] = progress
                    self.root.update_idletasks()

            distance_matrix = 1 - matrix.astype(float).values
            np.fill_diagonal(distance_matrix, 0)

            if np.any(np.isnan(distance_matrix)):
                self.results_text_gen.insert(tk.END,
                    f"\n⚠️ Fixing distance matrix (replacing NaN with 1.0)...\n")
                distance_matrix = np.nan_to_num(distance_matrix, nan=1.0)

            distance_matrix = (distance_matrix + distance_matrix.T) / 2
            condensed_dist = squareform(distance_matrix, checks=False)
            linkage_matrix = linkage(condensed_dist, method='average')

            self.genetic_query_distance_matrix = condensed_dist
            self.genetic_query_labels = query_names.copy()
            
            # Calculate cophenetic correlation for query dendrogram
            query_cophenetic = self.calculate_cophenetic_correlation(condensed_dist, linkage_matrix)
            self.genetic_query_cophenetic_correlation = query_cophenetic
            
            self.results_text_gen.insert(tk.END,
                f"\n📊 Query Dendrogram Cophenetic Correlation: {query_cophenetic:.4f}\n")
            if query_cophenetic > 0.8:
                self.results_text_gen.insert(tk.END,
                    f"   ✅ Excellent representation (c > 0.8)\n")
            elif query_cophenetic > 0.7:
                self.results_text_gen.insert(tk.END,
                    f"   ✓ Good representation (c > 0.7)\n")
            elif query_cophenetic > 0.6:
                self.results_text_gen.insert(tk.END,
                    f"   ⚠️ Moderate representation (c > 0.6)\n")
            else:
                self.results_text_gen.insert(tk.END,
                    f"   ❌ Poor representation (c ≤ 0.6)\n")

            self.progress_bar_gen['value'] = 75
            self.root.update_idletasks()

            if self.bootstrap_var_gen.get():
                self.results_text_gen.insert(tk.END,
                    "\n🔄 Calculating bootstrap supports for Query samples...\n")
                self.root.update_idletasks()
                supports = self.calculate_bootstrap_support(query_data, linkage_matrix)
                self.plot_genetic_query_dendrogram(
                    linkage_matrix, query_names, supports, output_dir, missing_by_sample, query_cophenetic)
            else:
                self.plot_genetic_query_dendrogram(
                    linkage_matrix, query_names, output_dir=output_dir,
                    missing_info=missing_by_sample, cophenetic_corr=query_cophenetic)

            self.progress_bar_gen['value'] = 85
            self.root.update_idletasks()

            self.genetic_query_pcoa_coords = self.generate_genetic_query_pcoa(
                distance_matrix, query_names, output_dir, missing_by_sample)

        except Exception as e:
            self.results_text_gen.insert(tk.END,
                f"\n❌ Error in Query-only analysis: {e}\n")
            import traceback
            traceback.print_exc()

        self.progress_bar_gen['value'] = 100
        self.root.update_idletasks()
        self.update_status("Genetic indices analysis completed")

    def plot_genetic_query_dendrogram(self, linkage_matrix, labels, supports=None,
                                      output_dir=None, missing_info=None, cophenetic_corr=None):
        if output_dir is None:
            output_dir = self.create_output_directory("Genetic_codominant")

        plt.figure(figsize=(12, 8))
        dendro = dendrogram(linkage_matrix, labels=labels,
                            leaf_rotation=90, leaf_font_size=14)

        if supports:
            icoord = np.array(dendro['icoord'])
            dcoord = np.array(dendro['dcoord'])
            for i, d in zip(icoord, dcoord):
                x = 0.5 * (i[1] + i[2])
                y = d[1]
                for cluster_idx, row in enumerate(linkage_matrix):
                    if abs(row[2] - y) < 1e-5:
                        support = supports.get(cluster_idx + len(labels), 0)
                        if support >= 50:
                            plt.text(x, y, f"{int(support)}%",
                                     ha='center', va='bottom',
                                     fontsize=8, color='red')
                        break

        title = "GMDA - Query Samples UPGMA Dendrogram"
        if supports:
            title += " (with bootstrap)"
        if cophenetic_corr is not None:
            title += f"\nCophenetic Correlation: {cophenetic_corr:.4f}"
        plt.title(title, fontsize=14)
        plt.ylabel("Dissimilarity (1-S)", fontsize=14)
        plt.tight_layout()

        filename = "codominant_genetic_query_dendrogram"
        if supports:
            filename += "_with_bootstrap"
        plot_path = os.path.join(output_dir, f"{filename}.png")
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        plt.close()
        self.results_text_gen.insert(tk.END,
            f"📊 Query dendrogram saved as: {plot_path}\n")

    def generate_genetic_query_pcoa(self, distance_matrix, labels, output_dir,
                                    missing_info=None):
        try:
            if np.any(np.isnan(distance_matrix)) or np.any(np.isinf(distance_matrix)):
                self.results_text_gen.insert(tk.END,
                    f"\n⚠️ Fixing distance matrix for PCoA...\n")
                distance_matrix = np.nan_to_num(
                    distance_matrix, nan=1.0, posinf=1.0, neginf=0.0)

            distance_matrix = (distance_matrix + distance_matrix.T) / 2
            np.fill_diagonal(distance_matrix, 0)
            distance_matrix += np.eye(len(distance_matrix)) * 1e-10

            if np.any(distance_matrix < 0):
                self.results_text_gen.insert(tk.END,
                    f"   ⚠️ Negative values detected, fixing...\n")
                distance_matrix = np.abs(distance_matrix)

            mds = MDS(n_components=2, dissimilarity='precomputed',
                      random_state=42, max_iter=300)
            coords = mds.fit_transform(distance_matrix)

            total_inertia = np.sum(distance_matrix ** 2) / (2 * len(labels))
            inertia_per_component = [
                np.sum((coords[:, i] - coords[:, i].mean()) ** 2)
                for i in range(2)
            ]
            var_explained = [
                inertia / total_inertia * 100
                for inertia in inertia_per_component
            ]

            plt.figure(figsize=(10, 8))
            plt.scatter(coords[:, 0], coords[:, 1],
                        c='steelblue', s=100, alpha=0.7)
            for i, label in enumerate(labels):
                plt.text(coords[i, 0], coords[i, 1], label,
                         fontsize=10, ha='center', va='bottom')
            plt.xlabel(f"PCo1 ({var_explained[0]:.1f}%)", fontsize=14)
            plt.ylabel(f"PCo2 ({var_explained[1]:.1f}%)", fontsize=14)
            plt.title("GMDA - Query Samples PCoA", fontsize=14)
            plt.grid(True, alpha=0.3, linestyle='--')
            plt.axhline(0, color='gray', linewidth=0.5)
            plt.axvline(0, color='gray', linewidth=0.5)
            plt.tight_layout()

            pcoa_path = os.path.join(output_dir, "codominant_genetic_query_PCoA.png")
            plt.savefig(pcoa_path, dpi=300, bbox_inches='tight')
            plt.close()

            self.results_text_gen.insert(tk.END,
                f"📈 Query PCoA plot saved as: {pcoa_path}\n")
            self.results_text_gen.insert(tk.END,
                f"📊 Variance explained: PCo1 = {var_explained[0]:.1f}%, "
                f"PCo2 = {var_explained[1]:.1f}%\n")

            has_missing_col = [
                label in missing_info for label in labels
            ] if missing_info else [False] * len(labels)
            pcoa_data = pd.DataFrame({
                'Sample': labels,
                'PCo1': coords[:, 0],
                'PCo2': coords[:, 1],
                'Has_Missing': has_missing_col
            })
            pcoa_excel_path = os.path.join(
                output_dir, "codominant_genetic_query_PCoA_coordinates.xlsx")
            pcoa_data.to_excel(pcoa_excel_path, index=False)
            self.results_text_gen.insert(tk.END,
                f"💾 PCoA coordinates saved as: {pcoa_excel_path}\n")

            return coords

        except Exception as e:
            self.results_text_gen.insert(tk.END,
                f"❌ Error generating PCoA: {e}\n")
            import traceback
            traceback.print_exc()
            return None

    # ==================== DOMINANT MARKER FUNCTIONS ====================

    def setup_genetic_dominant_tab(self):
        main_container = ModernFrame(self.genetic_dominant_frame)
        main_container.pack(fill="both", expand=True, padx=10, pady=10)

        input_card = CardFrame(main_container)
        input_card.pack(fill="x", pady=(0, 15))

        input_title = ModernLabel(input_card, text="Input Files (Dominant Markers - Presence/Absence)", font=("Segoe UI", 12, "bold"))
        input_title.pack(anchor="w", pady=(0, 15))

        info_label = ModernLabel(input_card,
                                 text="💡 Format: Row1=locus names, Col1=sample names. Use 1=Present, 0=Absent (or any text that converts to 1/0).",
                                 font=("Segoe UI", 9),
                                 fg=MODERN_THEME["text_light"])
        info_label.pack(anchor="w", pady=(0, 10))

        db_frame = ModernFrame(input_card)
        db_frame.pack(fill="x", pady=8)

        label_database = ModernLabel(db_frame, text="Reference Database:", font=("Segoe UI", 10, "bold"))
        label_database.grid(row=0, column=0, sticky="w", padx=(0, 10))

        self.entry_database_dominant = ModernEntry(db_frame, width=50)
        self.entry_database_dominant.grid(row=0, column=1, padx=5, sticky="ew")

        btn_browse_db = ModernButton(db_frame, text="Browse",
                                   command=lambda: self.load_file(self.entry_database_dominant),
                                   bg=MODERN_THEME["accent"], fg="white")
        btn_browse_db.grid(row=0, column=2, padx=5)

        query_frame = ModernFrame(input_card)
        query_frame.pack(fill="x", pady=8)

        label_query = ModernLabel(query_frame, text="Query Samples:", font=("Segoe UI", 10, "bold"))
        label_query.grid(row=0, column=0, sticky="w", padx=(0, 10))

        self.entry_query_dominant = ModernEntry(query_frame, width=50)
        self.entry_query_dominant.grid(row=0, column=1, padx=5, sticky="ew")

        btn_browse_query = ModernButton(query_frame, text="Browse",
                                      command=lambda: self.load_file(self.entry_query_dominant),
                                      bg=MODERN_THEME["accent"], fg="white")
        btn_browse_query.grid(row=0, column=2, padx=5)

        options_frame = ModernFrame(input_card)
        options_frame.pack(fill="x", pady=10)

        self.bootstrap_var_dominant = tk.BooleanVar()
        check_bootstrap = tk.Checkbutton(options_frame, text="Use Bootstrap",
                                       variable=self.bootstrap_var_dominant,
                                       bg=MODERN_THEME["card_bg"],
                                       font=("Segoe UI", 10),
                                       selectcolor=MODERN_THEME["light"])
        check_bootstrap.pack(anchor="w")

        db_frame.columnconfigure(1, weight=1)
        query_frame.columnconfigure(1, weight=1)

        action_card = CardFrame(main_container)
        action_card.pack(fill="x", pady=(0, 15))

        action_title = ModernLabel(action_card, text="Analysis Tools", font=("Segoe UI", 12, "bold"))
        action_title.pack(anchor="w", pady=(0, 15))

        btn_frame = ModernFrame(action_card)
        btn_frame.pack(fill="x")

        btn_top3 = ModernButton(btn_frame,
                              text="🧬 Generate Similarity Matrix & Dendrogram",
                              command=self.calculate_dominant_similarity,
                              bg=MODERN_THEME["success"], fg="white")
        btn_top3.grid(row=0, column=0, padx=5, pady=5, sticky="ew")

        btn_save = ModernButton(btn_frame,
                              text="💾 Save Top 3 Results",
                              command=self.save_dominant_results,
                              bg=MODERN_THEME["accent"], fg="white")
        btn_save.grid(row=0, column=1, padx=5, pady=5, sticky="ew")

        btn_query_analysis = ModernButton(btn_frame,
                                         text="📈 Query Data Analysis",
                                         command=self.calculate_dominant_query_analysis,
                                         bg=MODERN_THEME["primary"], fg="white")
        btn_query_analysis.grid(row=1, column=0, columnspan=2, padx=5, pady=5, sticky="ew")

        btn_frame.columnconfigure(0, weight=1)
        btn_frame.columnconfigure(1, weight=1)

        progress_frame = ModernFrame(action_card)
        progress_frame.pack(fill="x", pady=(10, 0))

        self.progress_label_dominant = ModernLabel(progress_frame, text="Progress:", font=("Segoe UI", 9, "bold"))
        self.progress_label_dominant.pack(anchor="w")

        self.progress_bar_dominant = ttk.Progressbar(progress_frame, mode='determinate', length=100)
        self.progress_bar_dominant.pack(fill="x", pady=(5, 0))

        results_card = CardFrame(main_container)
        results_card.pack(fill="both", expand=True)

        results_title = ModernLabel(results_card, text="Results", font=("Segoe UI", 12, "bold"))
        results_title.pack(anchor="w", pady=(0, 10))

        text_frame = ModernFrame(results_card)
        text_frame.pack(fill="both", expand=True)

        self.results_text_dominant = scrolledtext.ScrolledText(text_frame,
                                           bg="white",
                                           fg=MODERN_THEME["text"],
                                           font=("Consolas", 9),
                                           relief="flat",
                                           padx=10,
                                           pady=10)
        self.results_text_dominant.pack(fill="both", expand=True)

    def calculate_dominant_similarity(self):
        self.update_status("Calculating genetic similarity for dominant markers...")
        self.progress_bar_dominant['value'] = 0
        self.progress_bar_dominant['maximum'] = 100
        self.root.update_idletasks()

        db_filename = self.entry_database_dominant.get().strip()
        query_filename = self.entry_query_dominant.get().strip()

        if not db_filename or not query_filename:
            messagebox.showwarning("Warning", "Please select both files.")
            self.update_status("Ready")
            return

        try:
            db_df, db_name_col = self.read_genetic_file(db_filename)
            query_df, query_name_col = self.read_genetic_file(query_filename)
        except Exception as e:
            messagebox.showerror("Error", f"Error reading files:\n{e}")
            self.update_status("Error reading files")
            return

        if db_df.empty or query_df.empty:
            messagebox.showerror("Error", "One or both files are empty.")
            self.update_status("Empty files")
            return

        self.progress_bar_dominant['value'] = 10
        self.root.update_idletasks()

        output_dir = self.create_output_directory("Genetic_dominant")

        self.results_text_dominant.delete(1.0, tk.END)
        self.results_text_dominant.insert(tk.END,
            "🔗 GENETIC SIMILARITY - DOMINANT MARKERS (Jaccard Index)\n")
        self.results_text_dominant.insert(tk.END, "=" * 60 + "\n\n")

        self.dominant_top3_results = []
        total_queries = len(query_df)
        queries_processed = 0

        for idx_query, query_row in query_df.iterrows():
            query_name = query_row[query_name_col]
            if pd.isna(query_name) or str(query_name).strip() == '':
                continue

            query_data = query_row.drop(query_name_col).values
            similarities = {}

            for idx_ref, ref_row in db_df.iterrows():
                ref_name = ref_row[db_name_col]
                if pd.isna(ref_name) or str(ref_name).strip() == '':
                    continue
                ref_data = ref_row.drop(db_name_col).values
                similarity = self.calcular_similaridade_dominante(query_data, ref_data)
                similarities[ref_name] = similarity

            if not similarities:
                continue

            top3 = sorted(similarities.items(), key=lambda x: x[1], reverse=True)[:3]
            self.results_text_dominant.insert(tk.END, f"🔍 Query sample: {query_name}\n")
            for ref_name, score in top3:
                self.results_text_dominant.insert(tk.END,
                    f"   📌 {ref_name}: {score * 100:.2f}%\n")
                self.dominant_top3_results.append((query_name, ref_name, score))

            self.dominant_top3_results.append(("", "", ""))
            queries_processed += 1
            progress = 10 + (queries_processed / total_queries * 40)
            self.progress_bar_dominant['value'] = progress
            self.root.update_idletasks()

        try:
            self.progress_bar_dominant['value'] = 60
            self.root.update_idletasks()
            self.generate_dominant_matrix_and_dendrogram(
                query_df, db_df, query_name_col, db_name_col,
                self.bootstrap_var_dominant.get(), output_dir)
        except Exception as e:
            self.results_text_dominant.insert(tk.END,
                f"\n❌ Error generating matrix/dendrogram: {e}\n")

        self.progress_bar_dominant['value'] = 100
        self.root.update_idletasks()
        self.update_status("Genetic similarity analysis completed")

    def generate_dominant_matrix_and_dendrogram(self, query_df, db_df, query_name_col,
                                                 db_name_col, bootstrap_enabled, output_dir):
        combined_df = pd.concat([query_df, db_df], ignore_index=True)
        name_col = query_name_col if query_name_col in combined_df.columns else db_name_col
        combined_df = combined_df.dropna(subset=[name_col])

        if combined_df.empty:
            self.results_text_dominant.insert(tk.END, "❌ Error: No valid samples found.\n")
            return

        names = []
        for name in combined_df.iloc[:, 0]:
            if pd.notna(name):
                names.append(str(name).strip())
            else:
                names.append(None)

        data = combined_df.iloc[:, 1:]
        valid_indices = [i for i, name in enumerate(names) if name and name != '']
        names = [names[i] for i in valid_indices]
        data = data.iloc[valid_indices]

        if len(names) < 2:
            self.results_text_dominant.insert(tk.END,
                "❌ Error: At least 2 samples required.\n")
            return

        missing_by_sample = {}
        for i, sample in enumerate(names):
            sample_data = data.iloc[i]
            n_missing = sum(self.is_missing_dominant(val) for val in sample_data)
            if n_missing > 0:
                missing_by_sample[sample] = n_missing

        if missing_by_sample:
            self.results_text_dominant.insert(tk.END,
                f"\n⚠️ Missing data detected in {len(missing_by_sample)} samples\n")
            for sample, count in missing_by_sample.items():
                self.results_text_dominant.insert(tk.END,
                    f"   • {sample}: {count} empty cells\n")

        matrix = pd.DataFrame(index=names, columns=names, dtype=float)
        np.fill_diagonal(matrix.values, 1.0)

        total_pairs = len(data) * (len(data) - 1) / 2
        pairs_done = 0

        for i in range(len(data)):
            for j in range(i + 1, len(data)):
                sim = self.calcular_similaridade_dominante(
                    data.iloc[i].values, data.iloc[j].values)
                matrix.iat[i, j] = sim
                matrix.iat[j, i] = sim
                pairs_done += 1
                progress = 60 + (pairs_done / total_pairs * 20)
                self.progress_bar_dominant['value'] = progress
                self.root.update_idletasks()

        matrix_filename = os.path.join(
            output_dir, "dominant_genetic_similarity_matrix.xlsx")
        try:
            with pd.ExcelWriter(matrix_filename, engine='openpyxl') as writer:
                title_df = pd.DataFrame(
                    {"": ["GMDA - Dominant Genetic Similarity Matrix (Jaccard)"]})
                title_df.to_excel(writer, sheet_name='Similarity Matrix',
                                  index=False, header=False)
                matrix.to_excel(writer, sheet_name='Similarity Matrix', startrow=2)
            self.results_text_dominant.insert(tk.END,
                f"\n💾 Similarity matrix saved as: {matrix_filename}\n")
        except Exception as e:
            self.results_text_dominant.insert(tk.END,
                f"\n❌ Error saving matrix: {e}\n")

        distance_matrix = 1 - matrix.astype(float).values
        np.fill_diagonal(distance_matrix, 0)
        if np.any(np.isnan(distance_matrix)):
            distance_matrix = np.nan_to_num(distance_matrix, nan=1.0)

        try:
            condensed_dist = squareform(distance_matrix, checks=False)
            linkage_matrix = linkage(condensed_dist, method='average')

            self.dominant_linkage_matrix = linkage_matrix
            self.dominant_distance_matrix = condensed_dist
            self.dominant_labels = names.copy()
            
            # Calculate cophenetic correlation
            cophenetic_corr = self.calculate_cophenetic_correlation(condensed_dist, linkage_matrix)
            self.dominant_cophenetic_correlation = cophenetic_corr
            
            self.results_text_dominant.insert(tk.END,
                f"\n📊 Cophenetic Correlation Coefficient: {cophenetic_corr:.4f}\n")
            if cophenetic_corr > 0.8:
                self.results_text_dominant.insert(tk.END,
                    f"   ✅ Excellent representation (c > 0.8)\n")
            elif cophenetic_corr > 0.7:
                self.results_text_dominant.insert(tk.END,
                    f"   ✓ Good representation (c > 0.7)\n")
            elif cophenetic_corr > 0.6:
                self.results_text_dominant.insert(tk.END,
                    f"   ⚠️ Moderate representation (c > 0.6)\n")
            else:
                self.results_text_dominant.insert(tk.END,
                    f"   ❌ Poor representation (c ≤ 0.6)\n")

            self.progress_bar_dominant['value'] = 85
            self.root.update_idletasks()

            if bootstrap_enabled:
                self.results_text_dominant.insert(tk.END,
                    "🔄 Calculating bootstrap supports...\n")
                self.root.update_idletasks()
                supports = self.calculate_bootstrap_support_dominant(data, linkage_matrix)
                self.plot_dominant_dendrogram(
                    linkage_matrix, names, supports, output_dir, missing_by_sample, cophenetic_corr)
            else:
                self.plot_dominant_dendrogram(
                    linkage_matrix, names, output_dir=output_dir,
                    missing_info=missing_by_sample, cophenetic_corr=cophenetic_corr)

            self.progress_bar_dominant['value'] = 95
            self.root.update_idletasks()

        except Exception as e:
            self.results_text_dominant.insert(tk.END,
                f"❌ Error computing dendrogram: {e}\n")

    def calculate_bootstrap_support_dominant(self, data, linkage_matrix, replicates=100):
        n = data.shape[0]
        clusters_suporte = defaultdict(int)
        original_clusters = self.get_cluster_members(linkage_matrix, n)

        for rep in range(replicates):
            try:
                sampled_cols = np.random.choice(
                    data.shape[1], size=data.shape[1], replace=True)
                data_boot = data.iloc[:, sampled_cols]

                matrix_boot = pd.DataFrame(
                    index=range(len(data_boot)),
                    columns=range(len(data_boot)), dtype=float)
                np.fill_diagonal(matrix_boot.values, 1.0)

                for i1 in range(len(data_boot)):
                    for i2 in range(i1 + 1, len(data_boot)):
                        sim = self.calcular_similaridade_dominante(
                            data_boot.iloc[i1].values, data_boot.iloc[i2].values)
                        matrix_boot.iat[i1, i2] = sim
                        matrix_boot.iat[i2, i1] = sim

                dist = 1 - matrix_boot.values
                if np.any(np.isnan(dist)):
                    dist = np.nan_to_num(dist, nan=1.0)

                condensed_dist = squareform(dist, checks=False)
                linkage_boot = linkage(condensed_dist, method='average')
                bootstrap_clusters = self.get_cluster_members(linkage_boot, n)

                for cluster_id, members in original_clusters.items():
                    if any(members == b_members
                           for b_members in bootstrap_clusters.values()):
                        clusters_suporte[cluster_id] += 1

                progress = 85 + (rep / replicates * 10)
                self.progress_bar_dominant['value'] = progress
                self.root.update_idletasks()

            except Exception:
                continue

        for c in clusters_suporte:
            clusters_suporte[c] = clusters_suporte[c] / replicates * 100

        return clusters_suporte

    def plot_dominant_dendrogram(self, linkage_matrix, labels, supports=None,
                                 output_dir=None, missing_info=None, cophenetic_corr=None):
        if output_dir is None:
            output_dir = self.create_output_directory("Genetic_dominant")

        plt.figure(figsize=(12, 8))
        dendro = dendrogram(linkage_matrix, labels=labels,
                            leaf_rotation=90, leaf_font_size=10)

        if supports:
            icoord = np.array(dendro['icoord'])
            dcoord = np.array(dendro['dcoord'])
            for i, d in zip(icoord, dcoord):
                x = 0.5 * (i[1] + i[2])
                y = d[1]
                for cluster_idx, row in enumerate(linkage_matrix):
                    if abs(row[2] - y) < 1e-5:
                        support = supports.get(cluster_idx + len(labels), 0)
                        if support >= 50:
                            plt.text(x, y, f"{int(support)}%",
                                     ha='center', va='bottom',
                                     fontsize=8, color='red')
                        break

        title = "GMDA - Dominant Markers UPGMA Dendrogram (Jaccard)"
        if supports:
            title += " (with bootstrap)"
        if cophenetic_corr is not None:
            title += f"\nCophenetic Correlation: {cophenetic_corr:.4f}"
        plt.title(title, fontsize=14)
        plt.ylabel("Distance (1 - Jaccard Similarity)", fontsize=14)
        plt.tight_layout()

        filename = "dominant_markers_dendrogram"
        if supports:
            filename += "_with_bootstrap"
        plot_path = os.path.join(output_dir, f"{filename}.png")
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        plt.close()
        self.results_text_dominant.insert(tk.END, f"📊 Dendrogram saved as: {plot_path}\n")

    def save_dominant_results(self):
        if not self.dominant_top3_results:
            messagebox.showinfo("No data", "No results available to save.")
            return

        output_dir = self.create_output_directory("Genetic_dominant")
        df_results = pd.DataFrame(
            self.dominant_top3_results, columns=["Query", "Reference", "Similarity"])
        df_results["Similarity (%)"] = pd.to_numeric(
            df_results["Similarity"], errors='coerce') * 100
        df_results.drop(columns=["Similarity"], inplace=True)

        save_path = os.path.join(output_dir, "Top3_dominant.tsv")
        try:
            with open(save_path, 'w', encoding='utf-8') as f:
                f.write("GMDA - Dominant Genetic Similarity Results\n\n")
                df_results.to_csv(f, sep="\t", index=False)
            self.results_text_dominant.insert(tk.END,
                f"\n💾 Results saved as TSV: {save_path}\n")
        except Exception as e:
            self.results_text_dominant.insert(tk.END, f"❌ Error saving TSV: {e}\n")

    def calculate_dominant_query_analysis(self):
        self.update_status("Performing complete query analysis for dominant markers...")
        self.progress_bar_dominant['value'] = 0
        self.progress_bar_dominant['maximum'] = 100
        self.root.update_idletasks()

        filename = self.entry_query_dominant.get().strip()
        if not filename:
            messagebox.showwarning("Warning", "Please select a query data file.")
            self.update_status("Ready")
            return

        try:
            df, name_col = self.read_genetic_file(filename)
        except Exception as e:
            messagebox.showerror("Error", f"Error reading file:\n{e}")
            self.update_status("Error reading file")
            return

        if df.empty:
            messagebox.showerror("Error", "The file is empty.")
            self.update_status("Empty file")
            return

        self.progress_bar_dominant['value'] = 5
        self.root.update_idletasks()

        output_dir = self.create_output_directory("Genetic_dominant")

        sample_names = []
        for name in df.iloc[:, 0]:
            if pd.notna(name):
                sample_names.append(str(name).strip())
            else:
                sample_names.append(None)

        data = df.iloc[:, 1:]
        valid_indices = [i for i, name in enumerate(sample_names) if name and name != '']
        sample_names = [sample_names[i] for i in valid_indices]
        data = data.iloc[valid_indices]

        if len(sample_names) == 0:
            messagebox.showerror("Error", "No valid samples found.")
            self.update_status("No valid samples")
            return

        self.progress_bar_dominant['value'] = 10
        self.root.update_idletasks()

        self.results_text_dominant.delete(1.0, tk.END)
        self.results_text_dominant.insert(tk.END, "🧬 COMPLETE QUERY ANALYSIS - DOMINANT MARKERS\n")
        self.results_text_dominant.insert(tk.END, "=" * 60 + "\n\n")

        file_basename = os.path.splitext(os.path.basename(filename))[0]
        self.results_text_dominant.insert(tk.END, f"📁 File: {file_basename}\n")
        self.results_text_dominant.insert(tk.END, f"🔢 Samples: {len(sample_names)}\n\n")

        total_missing = 0
        missing_by_sample = {}
        missing_by_locus = {}

        for i, sample in enumerate(sample_names):
            sample_data = data.iloc[i]
            n_missing = sum(self.is_missing_dominant(val) for val in sample_data)
            if n_missing > 0:
                missing_by_sample[sample] = n_missing
                total_missing += n_missing

        for j in range(data.shape[1]):
            locus_data = data.iloc[:, j]
            n_missing = sum(self.is_missing_dominant(val) for val in locus_data)
            if n_missing > 0:
                missing_by_locus[f"Locus_{j+1}"] = n_missing

        if total_missing > 0:
            self.results_text_dominant.insert(tk.END,
                f"\n⚠️ **ATTENTION: MISSING DATA DETECTED** ⚠️\n")
            self.results_text_dominant.insert(tk.END,
                f"   Total empty cells: {total_missing}\n\n")
            self.results_text_dominant.insert(tk.END, "   📊 Distribution by sample:\n")
            for sample, n_missing in missing_by_sample.items():
                self.results_text_dominant.insert(tk.END,
                    f"      • {sample}: {n_missing} empty cells\n")
            self.results_text_dominant.insert(tk.END, "\n   📊 Distribution by locus:\n")
            for locus, n_missing in missing_by_locus.items():
                self.results_text_dominant.insert(tk.END,
                    f"      • {locus}: {n_missing} empty cells\n")
            self.results_text_dominant.insert(tk.END,
                "\n   🔧 For Jaccard similarity: using only valid loci in each comparison\n\n")
            self.root.update_idletasks()

        self.progress_bar_dominant['value'] = 20
        self.root.update_idletasks()

        self.results_text_dominant.insert(tk.END, "📊 CALCULATING GENETIC DIVERSITY...\n")
        self.root.update_idletasks()

        results = []
        for locus_idx in range(data.shape[1]):
            locus_data = data.iloc[:, locus_idx]
            locus_data_clean = []
            for x in locus_data:
                if not self.is_missing_dominant(x):
                    try:
                        x_str = str(x).upper().strip()
                        val = 1 if x_str in ['1', 'TRUE', 'YES', 'PRESENT', 'P', '+'] else 0
                        locus_data_clean.append(val)
                    except:
                        locus_data_clean.append(0)

            if len(locus_data_clean) == 0:
                continue

            count_zeros = sum(1 for x in locus_data_clean if x == 0)
            total_samples = len(locus_data_clean)
            q = count_zeros / total_samples
            p = 1 - q
            H = 2 * p * q

            results.append({
                'Locus': f"Locus_{locus_idx + 1}",
                'N_valid': total_samples,
                'Recessive_allele_freq': round(q, 4),
                'Dominant_allele_freq': round(p, 4),
                'Genetic_diversity_H': round(H, 4)
            })

        if results:
            mean_result = {
                'Locus': 'MEAN',
                'N_valid': round(np.mean([r['N_valid'] for r in results]), 1),
                'Recessive_allele_freq': round(
                    np.mean([r['Recessive_allele_freq'] for r in results]), 4),
                'Dominant_allele_freq': round(
                    np.mean([r['Dominant_allele_freq'] for r in results]), 4),
                'Genetic_diversity_H': round(
                    np.mean([r['Genetic_diversity_H'] for r in results]), 4)
            }
            results_with_mean = [mean_result] + results
        else:
            results_with_mean = results

        self.results_text_dominant.insert(tk.END, "📈 GENETIC DIVERSITY RESULTS:\n")
        for result in results_with_mean:
            if result['Locus'] == 'MEAN':
                self.results_text_dominant.insert(tk.END, f"📊 MEAN VALUES (all loci):\n")
                self.results_text_dominant.insert(tk.END,
                    f"   • Mean valid N per locus: {result['N_valid']:.1f}\n")
            else:
                self.results_text_dominant.insert(tk.END,
                    f"📍 {result['Locus']} (N={result['N_valid']}):\n")
            self.results_text_dominant.insert(tk.END,
                f"   • Recessive allele frequency (q): {result['Recessive_allele_freq']:.4f}\n")
            self.results_text_dominant.insert(tk.END,
                f"   • Dominant allele frequency (p): {result['Dominant_allele_freq']:.4f}\n")
            self.results_text_dominant.insert(tk.END,
                f"   • Genetic diversity (H): {result['Genetic_diversity_H']:.4f}\n\n")

        try:
            df_results = pd.DataFrame(results_with_mean)
            output_path = os.path.join(output_dir, "dominant_genetic_diversity.xlsx")
            with pd.ExcelWriter(output_path, engine='openpyxl') as writer:
                title_df = pd.DataFrame(
                    {"": ["GMDA - Genetic Diversity Analysis (Dominant Markers)"]})
                title_df.to_excel(writer, sheet_name='Diversity', index=False, header=False)
                df_results.to_excel(writer, sheet_name='Diversity', startrow=2, index=False)
            self.results_text_dominant.insert(tk.END,
                f"💾 Diversity results saved as: {output_path}\n")
        except Exception as e:
            self.results_text_dominant.insert(tk.END,
                f"❌ Error saving diversity results: {e}\n")

        self.progress_bar_dominant['value'] = 35
        self.root.update_idletasks()

        if len(sample_names) >= 2:
            self.results_text_dominant.insert(tk.END,
                "\n🔗 GENERATING SIMILARITY MATRIX AND DENDROGRAM...\n")
            self.root.update_idletasks()

            try:
                matrix = pd.DataFrame(index=sample_names, columns=sample_names, dtype=float)
                np.fill_diagonal(matrix.values, 1.0)

                total_pairs = len(data) * (len(data) - 1) / 2
                pairs_done = 0

                for i in range(len(data)):
                    for j in range(i + 1, len(data)):
                        sim = self.calcular_similaridade_dominante(
                            data.iloc[i].values, data.iloc[j].values)
                        matrix.iat[i, j] = sim
                        matrix.iat[j, i] = sim
                        pairs_done += 1
                        progress = 35 + (pairs_done / total_pairs * 25)
                        self.progress_bar_dominant['value'] = progress
                        self.root.update_idletasks()

                matrix_filename = os.path.join(
                    output_dir, "dominant_query_similarity_matrix.xlsx")
                with pd.ExcelWriter(matrix_filename, engine='openpyxl') as writer:
                    title_df = pd.DataFrame(
                        {"": ["GMDA - Dominant Query Similarity Matrix (Jaccard)"]})
                    title_df.to_excel(writer, sheet_name='Similarity Matrix',
                                      index=False, header=False)
                    matrix.to_excel(writer, sheet_name='Similarity Matrix', startrow=2)
                self.results_text_dominant.insert(tk.END,
                    f"💾 Similarity matrix saved as: {matrix_filename}\n")

                distance_matrix = 1 - matrix.astype(float).values
                np.fill_diagonal(distance_matrix, 0)
                if np.any(np.isnan(distance_matrix)):
                    distance_matrix = np.nan_to_num(distance_matrix, nan=1.0)

                condensed_dist = squareform(distance_matrix, checks=False)
                linkage_matrix = linkage(condensed_dist, method='average')

                self.dominant_query_distance_matrix = condensed_dist
                self.dominant_query_labels = sample_names.copy()
                
                # Calculate cophenetic correlation for query dendrogram
                query_cophenetic = self.calculate_cophenetic_correlation(condensed_dist, linkage_matrix)
                self.dominant_query_cophenetic_correlation = query_cophenetic
                
                self.results_text_dominant.insert(tk.END,
                    f"\n📊 Query Dendrogram Cophenetic Correlation: {query_cophenetic:.4f}\n")
                if query_cophenetic > 0.8:
                    self.results_text_dominant.insert(tk.END,
                        f"   ✅ Excellent representation (c > 0.8)\n")
                elif query_cophenetic > 0.7:
                    self.results_text_dominant.insert(tk.END,
                        f"   ✓ Good representation (c > 0.7)\n")
                elif query_cophenetic > 0.6:
                    self.results_text_dominant.insert(tk.END,
                        f"   ⚠️ Moderate representation (c > 0.6)\n")
                else:
                    self.results_text_dominant.insert(tk.END,
                        f"   ❌ Poor representation (c ≤ 0.6)\n")

                self.progress_bar_dominant['value'] = 70
                self.root.update_idletasks()

                if self.bootstrap_var_dominant.get():
                    self.results_text_dominant.insert(tk.END,
                        "🔄 Calculating bootstrap supports...\n")
                    self.root.update_idletasks()
                    supports = self.calculate_bootstrap_support_dominant(
                        data, linkage_matrix)
                    self.plot_dominant_dendrogram_query_only(
                        linkage_matrix, sample_names, supports, output_dir, missing_by_sample, query_cophenetic)
                else:
                    self.plot_dominant_dendrogram_query_only(
                        linkage_matrix, sample_names, output_dir=output_dir,
                        missing_info=missing_by_sample, cophenetic_corr=query_cophenetic)

                self.progress_bar_dominant['value'] = 80
                self.root.update_idletasks()

            except Exception as e:
                self.results_text_dominant.insert(tk.END,
                    f"❌ Error generating similarity matrix/dendrogram: {e}\n")

        if len(sample_names) >= 3:
            self.results_text_dominant.insert(tk.END, "\n📈 GENERATING PCoA...\n")
            self.root.update_idletasks()

            try:
                distance_matrix = np.zeros((len(sample_names), len(sample_names)))
                total_pairs_pcoa = len(sample_names) ** 2
                pairs_done_pcoa = 0

                for i in range(len(sample_names)):
                    for j in range(len(sample_names)):
                        if i != j:
                            sim = self.calcular_similaridade_dominante(
                                data.iloc[i].values, data.iloc[j].values)
                            distance_matrix[i, j] = 1 - sim
                        pairs_done_pcoa += 1
                        progress = 80 + (pairs_done_pcoa / total_pairs_pcoa * 15)
                        self.progress_bar_dominant['value'] = progress
                        self.root.update_idletasks()

                if np.any(np.isnan(distance_matrix)) or np.any(np.isinf(distance_matrix)):
                    distance_matrix = np.nan_to_num(
                        distance_matrix, nan=1.0, posinf=1.0, neginf=0.0)

                distance_matrix = (distance_matrix + distance_matrix.T) / 2
                np.fill_diagonal(distance_matrix, 0)
                distance_matrix += np.eye(len(distance_matrix)) * 1e-10

                mds = MDS(n_components=2, dissimilarity='precomputed',
                          random_state=42, max_iter=300)
                coords = mds.fit_transform(distance_matrix)

                self.dominant_query_pcoa_coords = coords

                total_inertia = np.sum(distance_matrix ** 2) / (2 * len(sample_names))
                inertia_per_component = [
                    np.sum((coords[:, i] - coords[:, i].mean()) ** 2)
                    for i in range(2)
                ]
                var_explained = [
                    inertia / total_inertia * 100
                    for inertia in inertia_per_component
                ]

                plt.figure(figsize=(10, 8))
                for i, (x, y) in enumerate(coords):
                    plt.scatter(x, y, color='steelblue', s=100, alpha=0.7)
                    plt.text(x, y, str(sample_names[i]), fontsize=10,
                             ha='center', va='bottom')
                plt.xlabel(f"PCo1 ({var_explained[0]:.1f}%)", fontsize=14)
                plt.ylabel(f"PCo2 ({var_explained[1]:.1f}%)", fontsize=14)
                plt.title("GMDA - Dominant Query PCoA (Jaccard)", fontsize=14)
                plt.grid(True, alpha=0.3, linestyle='--')
                plt.axhline(0, color='gray', linewidth=0.5)
                plt.axvline(0, color='gray', linewidth=0.5)
                plt.tight_layout()

                pcoa_path = os.path.join(output_dir, "dominant_query_PCoA.png")
                plt.savefig(pcoa_path, dpi=300, bbox_inches='tight')
                plt.close()

                self.results_text_dominant.insert(tk.END,
                    f"📈 PCoA plot saved as: {pcoa_path}\n")
                self.results_text_dominant.insert(tk.END,
                    f"📊 Variance explained: PCo1 = {var_explained[0]:.1f}%, "
                    f"PCo2 = {var_explained[1]:.1f}%\n")

                pcoa_data = pd.DataFrame({
                    'Sample': sample_names,
                    'PCo1': coords[:, 0],
                    'PCo2': coords[:, 1],
                    'Has_Missing': [
                        sample in missing_by_sample for sample in sample_names
                    ] if missing_by_sample else [False] * len(sample_names)
                })
                pcoa_excel_path = os.path.join(
                    output_dir, "dominant_query_PCoA_coordinates.xlsx")
                pcoa_data.to_excel(pcoa_excel_path, index=False)
                self.results_text_dominant.insert(tk.END,
                    f"💾 PCoA coordinates saved as: {pcoa_excel_path}\n")

                self.progress_bar_dominant['value'] = 98
                self.root.update_idletasks()

            except Exception as e:
                self.results_text_dominant.insert(tk.END,
                    f"❌ Error generating PCoA: {e}\n")

        self.progress_bar_dominant['value'] = 100
        self.root.update_idletasks()
        self.results_text_dominant.insert(tk.END,
            "\n✅ COMPLETE QUERY ANALYSIS FINISHED!\n")
        self.update_status("Complete query analysis for dominant markers completed")

    def plot_dominant_dendrogram_query_only(self, linkage_matrix, labels, supports=None,
                                            output_dir=None, missing_info=None, cophenetic_corr=None):
        if output_dir is None:
            output_dir = self.create_output_directory("Genetic_dominant")

        plt.figure(figsize=(12, 8))
        dendro = dendrogram(linkage_matrix, labels=labels,
                            leaf_rotation=90, leaf_font_size=10)

        if supports:
            icoord = np.array(dendro['icoord'])
            dcoord = np.array(dendro['dcoord'])
            for i, d in zip(icoord, dcoord):
                x = 0.5 * (i[1] + i[2])
                y = d[1]
                for cluster_idx, row in enumerate(linkage_matrix):
                    if abs(row[2] - y) < 1e-5:
                        support = supports.get(cluster_idx + len(labels), 0)
                        if support >= 50:
                            plt.text(x, y, f"{int(support)}%",
                                     ha='center', va='bottom',
                                     fontsize=8, color='red')
                        break

        title = "GMDA - Dominant Query UPGMA Dendrogram (Jaccard)"
        if supports:
            title += " (with bootstrap)"
        if cophenetic_corr is not None:
            title += f"\nCophenetic Correlation: {cophenetic_corr:.4f}"
        plt.title(title, fontsize=14)
        plt.ylabel("Distance (1 - Jaccard Similarity)", fontsize=14)
        plt.tight_layout()

        filename = "dominant_query_dendrogram"
        if supports:
            filename += "_with_bootstrap"
        plot_path = os.path.join(output_dir, f"{filename}.png")
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        plt.close()
        self.results_text_dominant.insert(tk.END,
            f"📊 Query dendrogram saved as: {plot_path}\n")

    # ==================== HAPLOID MARKER FUNCTIONS ====================

    def setup_genetic_haploid_tab(self):
        main_container = ModernFrame(self.genetic_haploid_frame)
        main_container.pack(fill="both", expand=True, padx=10, pady=10)

        input_card = CardFrame(main_container)
        input_card.pack(fill="x", pady=(0, 15))

        input_title = ModernLabel(input_card,
                                  text="Input Files (Haploid Markers - cpDNA, mtDNA, Y-chromosome, etc.)",
                                  font=("Segoe UI", 12, "bold"))
        input_title.pack(anchor="w", pady=(0, 15))

        info_label = ModernLabel(input_card,
                                 text="💡 Format: Row1=locus names, Col1=sample names. Each individual has one allele per locus. Empty cells = missing data.",
                                 font=("Segoe UI", 9),
                                 fg=MODERN_THEME["text_light"])
        info_label.pack(anchor="w", pady=(0, 10))

        db_frame = ModernFrame(input_card)
        db_frame.pack(fill="x", pady=8)

        label_database = ModernLabel(db_frame, text="Reference Database:", font=("Segoe UI", 10, "bold"))
        label_database.grid(row=0, column=0, sticky="w", padx=(0, 10))

        self.entry_database_haploid = ModernEntry(db_frame, width=50)
        self.entry_database_haploid.grid(row=0, column=1, padx=5, sticky="ew")

        btn_browse_db = ModernButton(db_frame, text="Browse",
                                     command=lambda: self.load_file(self.entry_database_haploid),
                                     bg=MODERN_THEME["accent"], fg="white")
        btn_browse_db.grid(row=0, column=2, padx=5)

        query_frame = ModernFrame(input_card)
        query_frame.pack(fill="x", pady=8)

        label_query = ModernLabel(query_frame, text="Query Samples:", font=("Segoe UI", 10, "bold"))
        label_query.grid(row=0, column=0, sticky="w", padx=(0, 10))

        self.entry_query_haploid = ModernEntry(query_frame, width=50)
        self.entry_query_haploid.grid(row=0, column=1, padx=5, sticky="ew")

        btn_browse_query = ModernButton(query_frame, text="Browse",
                                        command=lambda: self.load_file(self.entry_query_haploid),
                                        bg=MODERN_THEME["accent"], fg="white")
        btn_browse_query.grid(row=0, column=2, padx=5)

        options_frame = ModernFrame(input_card)
        options_frame.pack(fill="x", pady=10)

        self.bootstrap_var_haploid = tk.BooleanVar()
        check_bootstrap = tk.Checkbutton(options_frame, text="Use Bootstrap",
                                         variable=self.bootstrap_var_haploid,
                                         bg=MODERN_THEME["card_bg"],
                                         font=("Segoe UI", 10),
                                         selectcolor=MODERN_THEME["light"])
        check_bootstrap.pack(anchor="w")

        db_frame.columnconfigure(1, weight=1)
        query_frame.columnconfigure(1, weight=1)

        action_card = CardFrame(main_container)
        action_card.pack(fill="x", pady=(0, 15))

        action_title = ModernLabel(action_card, text="Analysis Tools", font=("Segoe UI", 12, "bold"))
        action_title.pack(anchor="w", pady=(0, 15))

        btn_frame = ModernFrame(action_card)
        btn_frame.pack(fill="x")

        btn_top3 = ModernButton(btn_frame,
                                text="🧬 Generate Similarity Matrix & Dendrogram",
                                command=self.calculate_haploid_similarity,
                                bg=MODERN_THEME["success"], fg="white")
        btn_top3.grid(row=0, column=0, padx=5, pady=5, sticky="ew")

        btn_save = ModernButton(btn_frame,
                                text="💾 Save Top 3 Results",
                                command=self.save_haploid_results,
                                bg=MODERN_THEME["accent"], fg="white")
        btn_save.grid(row=0, column=1, padx=5, pady=5, sticky="ew")

        btn_query_analysis = ModernButton(btn_frame,
                                          text="📈 Query Data Analysis",
                                          command=self.calculate_haploid_query_analysis,
                                          bg=MODERN_THEME["primary"], fg="white")
        btn_query_analysis.grid(row=1, column=0, columnspan=2, padx=5, pady=5, sticky="ew")

        btn_frame.columnconfigure(0, weight=1)
        btn_frame.columnconfigure(1, weight=1)

        progress_frame = ModernFrame(action_card)
        progress_frame.pack(fill="x", pady=(10, 0))

        self.progress_label_haploid = ModernLabel(progress_frame, text="Progress:", font=("Segoe UI", 9, "bold"))
        self.progress_label_haploid.pack(anchor="w")

        self.progress_bar_haploid = ttk.Progressbar(progress_frame, mode='determinate', length=100)
        self.progress_bar_haploid.pack(fill="x", pady=(5, 0))

        results_card = CardFrame(main_container)
        results_card.pack(fill="both", expand=True)

        results_title = ModernLabel(results_card, text="Results", font=("Segoe UI", 12, "bold"))
        results_title.pack(anchor="w", pady=(0, 10))

        text_frame = ModernFrame(results_card)
        text_frame.pack(fill="both", expand=True)

        self.results_text_haploid = scrolledtext.ScrolledText(text_frame,
                                                               bg="white",
                                                               fg=MODERN_THEME["text"],
                                                               font=("Consolas", 9),
                                                               relief="flat",
                                                               padx=10,
                                                               pady=10)
        self.results_text_haploid.pack(fill="both", expand=True)

    def calculate_haploid_similarity(self):
        self.update_status("Calculating genetic similarity for haploid markers...")
        self.progress_bar_haploid['value'] = 0
        self.progress_bar_haploid['maximum'] = 100
        self.root.update_idletasks()

        db_filename = self.entry_database_haploid.get().strip()
        query_filename = self.entry_query_haploid.get().strip()

        if not db_filename or not query_filename:
            messagebox.showwarning("Warning", "Please select both files.")
            self.update_status("Ready")
            return

        try:
            db_df, db_name_col = self.read_genetic_file(db_filename)
            query_df, query_name_col = self.read_genetic_file(query_filename)
        except Exception as e:
            messagebox.showerror("Error", f"Error reading files:\n{e}")
            self.update_status("Error reading files")
            return

        if db_df.empty or query_df.empty:
            messagebox.showerror("Error", "One or both files are empty.")
            self.update_status("Empty files")
            return

        self.progress_bar_haploid['value'] = 10
        self.root.update_idletasks()

        output_dir = self.create_output_directory("Genetic_haploid")

        self.results_text_haploid.delete(1.0, tk.END)
        self.results_text_haploid.insert(tk.END, "🔗 GENETIC SIMILARITY - HAPLOID MARKERS\n")
        self.results_text_haploid.insert(tk.END, "=" * 60 + "\n\n")

        self.haploid_top3_results = []
        total_queries = len(query_df)
        queries_processed = 0

        for idx_query, query_row in query_df.iterrows():
            query_name = query_row[query_name_col]
            if pd.isna(query_name) or str(query_name).strip() == '':
                continue

            query_data = query_row.drop(query_name_col).values
            similarities = {}

            for idx_ref, ref_row in db_df.iterrows():
                ref_name = ref_row[db_name_col]
                if pd.isna(ref_name) or str(ref_name).strip() == '':
                    continue
                ref_data = ref_row.drop(db_name_col).values
                similarity = self.calcular_similaridade_haploid(query_data, ref_data)
                similarities[ref_name] = similarity

            if not similarities:
                continue

            top3 = sorted(similarities.items(), key=lambda x: x[1], reverse=True)[:3]
            self.results_text_haploid.insert(tk.END, f"🔍 Query sample: {query_name}\n")
            for ref_name, score in top3:
                self.results_text_haploid.insert(tk.END, f"   📌 {ref_name}: {score * 100:.2f}%\n")
                self.haploid_top3_results.append((query_name, ref_name, score))

            self.haploid_top3_results.append(("", "", ""))
            queries_processed += 1
            progress = 10 + (queries_processed / total_queries * 40)
            self.progress_bar_haploid['value'] = progress
            self.root.update_idletasks()

        try:
            self.progress_bar_haploid['value'] = 60
            self.root.update_idletasks()
            self.generate_haploid_matrix_and_dendrogram(
                query_df, db_df, query_name_col, db_name_col,
                self.bootstrap_var_haploid.get(), output_dir
            )
        except Exception as e:
            self.results_text_haploid.insert(tk.END, f"\n❌ Error generating matrix/dendrogram: {e}\n")

        self.progress_bar_haploid['value'] = 100
        self.root.update_idletasks()
        self.update_status("Haploid genetic similarity analysis completed")

    def generate_haploid_matrix_and_dendrogram(self, query_df, db_df, query_name_col, db_name_col,
                                                bootstrap_enabled, output_dir):
        combined_df = pd.concat([query_df, db_df], ignore_index=True)
        name_col = query_name_col if query_name_col in combined_df.columns else db_name_col
        combined_df = combined_df.dropna(subset=[name_col])

        if combined_df.empty:
            self.results_text_haploid.insert(tk.END, "❌ Error: No valid samples found.\n")
            return

        names = []
        for name in combined_df.iloc[:, 0]:
            if pd.notna(name):
                names.append(str(name).strip())
            else:
                names.append(None)

        data = combined_df.iloc[:, 1:]
        valid_indices = [i for i, name in enumerate(names) if name and name != '']
        names = [names[i] for i in valid_indices]
        data = data.iloc[valid_indices]

        if len(names) < 2:
            self.results_text_haploid.insert(tk.END, "❌ Error: At least 2 samples required.\n")
            return

        missing_by_sample = {}
        for i, sample in enumerate(names):
            sample_data = data.iloc[i]
            n_missing = sum(self.is_missing_haploid(val) for val in sample_data)
            if n_missing > 0:
                missing_by_sample[sample] = n_missing

        if missing_by_sample:
            self.results_text_haploid.insert(tk.END,
                f"\n⚠️ Missing data detected in {len(missing_by_sample)} samples\n")
            for sample, count in missing_by_sample.items():
                self.results_text_haploid.insert(tk.END, f"   • {sample}: {count} empty cells\n")

        matrix = pd.DataFrame(index=names, columns=names, dtype=float)
        np.fill_diagonal(matrix.values, 1.0)

        total_pairs = len(data) * (len(data) - 1) / 2
        pairs_done = 0

        for i in range(len(data)):
            for j in range(i + 1, len(data)):
                sim = self.calcular_similaridade_haploid(data.iloc[i].values, data.iloc[j].values)
                matrix.iat[i, j] = sim
                matrix.iat[j, i] = sim
                pairs_done += 1
                progress = 60 + (pairs_done / total_pairs * 20)
                self.progress_bar_haploid['value'] = progress
                self.root.update_idletasks()

        matrix_filename = os.path.join(output_dir, "haploid_genetic_similarity_matrix.xlsx")
        try:
            with pd.ExcelWriter(matrix_filename, engine='openpyxl') as writer:
                title_df = pd.DataFrame({"": ["GMDA - Haploid Genetic Similarity Matrix"]})
                title_df.to_excel(writer, sheet_name='Similarity Matrix', index=False, header=False)
                matrix.to_excel(writer, sheet_name='Similarity Matrix', startrow=2)
            self.results_text_haploid.insert(tk.END,
                f"\n💾 Similarity matrix saved as: {matrix_filename}\n")
        except Exception as e:
            self.results_text_haploid.insert(tk.END, f"\n❌ Error saving matrix: {e}\n")

        distance_matrix = 1 - matrix.astype(float).values
        np.fill_diagonal(distance_matrix, 0)

        if np.any(np.isnan(distance_matrix)):
            distance_matrix = np.nan_to_num(distance_matrix, nan=1.0)

        try:
            condensed_dist = squareform(distance_matrix, checks=False)
            linkage_matrix = linkage(condensed_dist, method='average')

            self.haploid_linkage_matrix = linkage_matrix
            self.haploid_distance_matrix = condensed_dist
            self.haploid_labels = names.copy()
            
            # Calculate cophenetic correlation
            cophenetic_corr = self.calculate_cophenetic_correlation(condensed_dist, linkage_matrix)
            self.haploid_cophenetic_correlation = cophenetic_corr
            
            self.results_text_haploid.insert(tk.END,
                f"\n📊 Cophenetic Correlation Coefficient: {cophenetic_corr:.4f}\n")
            if cophenetic_corr > 0.8:
                self.results_text_haploid.insert(tk.END,
                    f"   ✅ Excellent representation (c > 0.8)\n")
            elif cophenetic_corr > 0.7:
                self.results_text_haploid.insert(tk.END,
                    f"   ✓ Good representation (c > 0.7)\n")
            elif cophenetic_corr > 0.6:
                self.results_text_haploid.insert(tk.END,
                    f"   ⚠️ Moderate representation (c > 0.6)\n")
            else:
                self.results_text_haploid.insert(tk.END,
                    f"   ❌ Poor representation (c ≤ 0.6)\n")

            self.progress_bar_haploid['value'] = 85
            self.root.update_idletasks()

            if bootstrap_enabled:
                self.results_text_haploid.insert(tk.END, "🔄 Calculating bootstrap supports...\n")
                self.root.update_idletasks()
                supports = self.calculate_bootstrap_support_haploid(data, linkage_matrix)
                self.plot_haploid_dendrogram(linkage_matrix, names, supports, output_dir, cophenetic_corr)
            else:
                self.plot_haploid_dendrogram(linkage_matrix, names, output_dir=output_dir, cophenetic_corr=cophenetic_corr)

            self.progress_bar_haploid['value'] = 95
            self.root.update_idletasks()

        except Exception as e:
            self.results_text_haploid.insert(tk.END, f"❌ Error computing dendrogram: {e}\n")

    def calculate_bootstrap_support_haploid(self, data, linkage_matrix, replicates=100):
        n = data.shape[0]
        clusters_suporte = defaultdict(int)
        original_clusters = self.get_cluster_members(linkage_matrix, n)

        for rep in range(replicates):
            try:
                sampled_cols = np.random.choice(
                    data.shape[1], size=data.shape[1], replace=True)
                data_boot = data.iloc[:, sampled_cols]

                matrix_boot = pd.DataFrame(
                    index=range(len(data_boot)),
                    columns=range(len(data_boot)), dtype=float)
                np.fill_diagonal(matrix_boot.values, 1.0)

                for i1 in range(len(data_boot)):
                    for i2 in range(i1 + 1, len(data_boot)):
                        sim = self.calcular_similaridade_haploid(
                            data_boot.iloc[i1].values, data_boot.iloc[i2].values)
                        matrix_boot.iat[i1, i2] = sim
                        matrix_boot.iat[i2, i1] = sim

                dist = 1 - matrix_boot.values
                if np.any(np.isnan(dist)):
                    dist = np.nan_to_num(dist, nan=1.0)

                condensed_dist = squareform(dist, checks=False)
                linkage_boot = linkage(condensed_dist, method='average')
                bootstrap_clusters = self.get_cluster_members(linkage_boot, n)

                for cluster_id, members in original_clusters.items():
                    if any(members == b_members
                           for b_members in bootstrap_clusters.values()):
                        clusters_suporte[cluster_id] += 1

                progress = 70 + (rep / replicates * 10)
                self.progress_bar_haploid['value'] = progress
                self.root.update_idletasks()

            except Exception:
                continue

        for c in clusters_suporte:
            clusters_suporte[c] = clusters_suporte[c] / replicates * 100

        return clusters_suporte

    def plot_haploid_dendrogram(self, linkage_matrix, labels, supports=None,
                                output_dir=None, cophenetic_corr=None):
        if output_dir is None:
            output_dir = self.create_output_directory("Genetic_haploid")

        plt.figure(figsize=(12, 8))
        dendro = dendrogram(linkage_matrix, labels=labels,
                            leaf_rotation=90, leaf_font_size=10)

        if supports:
            icoord = np.array(dendro['icoord'])
            dcoord = np.array(dendro['dcoord'])
            for i, d in zip(icoord, dcoord):
                x = 0.5 * (i[1] + i[2])
                y = d[1]
                for cluster_idx, row in enumerate(linkage_matrix):
                    if abs(row[2] - y) < 1e-5:
                        support = supports.get(cluster_idx + len(labels), 0)
                        if support >= 50:
                            plt.text(x, y, f"{int(support)}%",
                                     ha='center', va='bottom',
                                     fontsize=8, color='red')
                        break

        title = "GMDA - Haploid Markers UPGMA Dendrogram"
        if supports:
            title += " (with bootstrap)"
        if cophenetic_corr is not None:
            title += f"\nCophenetic Correlation: {cophenetic_corr:.4f}"
        plt.title(title, fontsize=14)
        plt.ylabel("Dissimilarity (1-S)", fontsize=14)
        plt.tight_layout()

        filename = "haploid_markers_dendrogram"
        if supports:
            filename += "_with_bootstrap"
        plot_path = os.path.join(output_dir, f"{filename}.png")
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        plt.close()
        self.results_text_haploid.insert(tk.END, f"📊 Dendrogram saved as: {plot_path}\n")

    def calculate_haploid_query_analysis(self):
        self.update_status("Performing complete query analysis for haploid markers...")
        self.progress_bar_haploid['value'] = 0
        self.progress_bar_haploid['maximum'] = 100
        self.root.update_idletasks()

        filename = self.entry_query_haploid.get().strip()
        if not filename:
            messagebox.showwarning("Warning", "Please select a query data file.")
            self.update_status("Ready")
            return

        try:
            df, name_col = self.read_genetic_file(filename)
        except Exception as e:
            messagebox.showerror("Error", f"Error reading file:\n{e}")
            self.update_status("Error reading file")
            return

        if df.empty:
            messagebox.showerror("Error", "The file is empty.")
            self.update_status("Empty file")
            return

        self.progress_bar_haploid['value'] = 5
        self.root.update_idletasks()

        output_dir = self.create_output_directory("Genetic_haploid")

        sample_names = []
        for name in df.iloc[:, 0]:
            if pd.notna(name):
                sample_names.append(str(name).strip())
            else:
                sample_names.append(None)

        data = df.iloc[:, 1:]
        valid_indices = [i for i, name in enumerate(sample_names) if name and name != '']
        sample_names = [sample_names[i] for i in valid_indices]
        data = data.iloc[valid_indices]

        if len(sample_names) == 0:
            messagebox.showerror("Error", "No valid samples found.")
            self.update_status("No valid samples")
            return

        self.progress_bar_haploid['value'] = 10
        self.root.update_idletasks()

        self.results_text_haploid.delete(1.0, tk.END)
        self.results_text_haploid.insert(tk.END, "🧬 COMPLETE QUERY ANALYSIS - HAPLOID MARKERS\n")
        self.results_text_haploid.insert(tk.END, "=" * 60 + "\n\n")

        file_basename = os.path.splitext(os.path.basename(filename))[0]
        self.results_text_haploid.insert(tk.END, f"📁 File: {file_basename}\n")
        self.results_text_haploid.insert(tk.END, f"🔢 Samples: {len(sample_names)}\n\n")

        total_missing = 0
        missing_by_sample = {}
        missing_by_locus = {}

        for i, sample in enumerate(sample_names):
            sample_data = data.iloc[i]
            n_missing = sum(self.is_missing_haploid(val) for val in sample_data)
            if n_missing > 0:
                missing_by_sample[sample] = n_missing
                total_missing += n_missing

        for j in range(data.shape[1]):
            locus_data = data.iloc[:, j]
            n_missing = sum(self.is_missing_haploid(val) for val in locus_data)
            if n_missing > 0:
                missing_by_locus[f"Locus_{j+1}"] = n_missing

        if total_missing > 0:
            self.results_text_haploid.insert(tk.END,
                f"\n⚠️ **WARNING: MISSING DATA DETECTED** ⚠️\n")
            self.results_text_haploid.insert(tk.END,
                f"   Total empty cells: {total_missing}\n\n")
            self.results_text_haploid.insert(tk.END, "   📊 Distribution by sample:\n")
            for sample, n_missing in missing_by_sample.items():
                self.results_text_haploid.insert(tk.END,
                    f"      • {sample}: {n_missing} empty cells\n")
            self.results_text_haploid.insert(tk.END, "\n   📊 Distribution by locus:\n")
            for locus, n_missing in missing_by_locus.items():
                self.results_text_haploid.insert(tk.END,
                    f"      • {locus}: {n_missing} empty cells\n")
            self.results_text_haploid.insert(tk.END,
                "\n   🔧 Only valid loci used in each comparison.\n\n")
            self.root.update_idletasks()

        self.progress_bar_haploid['value'] = 20
        self.root.update_idletasks()

        self.results_text_haploid.insert(tk.END, "📊 CALCULATING GENETIC DIVERSITY (HAPLOID)...\n")
        self.root.update_idletasks()

        diversity_results = []
        pop_A, pop_Ae, pop_h = [], [], []

        for locus_idx in range(data.shape[1]):
            locus_data = data.iloc[:, locus_idx]
            locus_data_clean = [x for x in locus_data if not self.is_missing_haploid(x)]

            if len(locus_data_clean) == 0:
                continue

            n = len(locus_data_clean)
            from collections import Counter
            allele_counts = Counter([str(x).strip() for x in locus_data_clean])
            A = len(allele_counts)

            freqs = np.array([count / n for count in allele_counts.values()])
            sum_p2 = np.sum(freqs ** 2)
            Ae = 1 / sum_p2 if sum_p2 > 0 else 0
            h = 1 - sum_p2

            diversity_results.append({
                'Locus': f"Locus_{locus_idx + 1}",
                'N_valid': n,
                'A': A,
                'Ae': round(Ae, 4),
                'h': round(h, 4),
                'Allele_frequencies': {str(k): round(v, 4) for k, v in
                                       sorted(allele_counts.items())}
            })
            pop_A.append(A)
            pop_Ae.append(Ae)
            pop_h.append(h)

        if diversity_results:
            mean_result = {
                'Locus': 'MEAN',
                'N_valid': round(np.mean([r['N_valid'] for r in diversity_results]), 1),
                'A': round(np.mean(pop_A), 4),
                'Ae': round(np.mean(pop_Ae), 4),
                'h': round(np.mean(pop_h), 4),
                'Allele_frequencies': {}
            }
            results_with_mean = [mean_result] + diversity_results
        else:
            results_with_mean = diversity_results

        self.results_text_haploid.insert(tk.END, "📈 GENETIC DIVERSITY RESULTS (HAPLOID):\n\n")
        for result in results_with_mean:
            if result['Locus'] == 'MEAN':
                self.results_text_haploid.insert(tk.END, "📊 MEAN VALUES (all loci):\n")
                self.results_text_haploid.insert(tk.END,
                    f"   • Mean valid N per locus: {result['N_valid']:.1f}\n")
                self.results_text_haploid.insert(tk.END,
                    f"   • Alleles per locus: {result['A']:.2f}\n")
                self.results_text_haploid.insert(tk.END,
                    f"   • Allelic richness (Ae): {result['Ae']:.4f}\n")
                self.results_text_haploid.insert(tk.END,
                    f"   • Gene diversity (h): {result['h']:.4f}\n\n")
            else:
                self.results_text_haploid.insert(tk.END,
                    f"📍 {result['Locus']} (N={result['N_valid']}):\n")
                if result['Allele_frequencies']:
                    freq_str = ", ".join(
                        [f"{k}:{v}" for k, v in result['Allele_frequencies'].items()])
                    self.results_text_haploid.insert(tk.END,
                        f"   • Allele frequencies: {freq_str}\n")
                self.results_text_haploid.insert(tk.END,
                    f"   • Alleles per locus: {result['A']:.2f}\n")
                self.results_text_haploid.insert(tk.END,
                    f"   • Allelic richness (Ae): {result['Ae']:.4f}\n")
                self.results_text_haploid.insert(tk.END,
                    f"   • Gene diversity (h): {result['h']:.4f}\n\n")

        try:
            df_export = pd.DataFrame([
                {'Locus': r['Locus'], 'N_valid': r['N_valid'],
                 'A': r['A'], 'Ae': r['Ae'], 'h': r['h']}
                for r in results_with_mean
            ])
            output_path = os.path.join(output_dir, "haploid_genetic_diversity.xlsx")
            with pd.ExcelWriter(output_path, engine='openpyxl') as writer:
                title_df = pd.DataFrame(
                    {"": ["GMDA - Genetic Diversity Analysis (Haploid Markers)"]})
                title_df.to_excel(writer, sheet_name='Diversity', index=False, header=False)
                df_export.to_excel(writer, sheet_name='Diversity', startrow=2, index=False)
            self.results_text_haploid.insert(tk.END,
                f"💾 Diversity results saved as: {output_path}\n")
        except Exception as e:
            self.results_text_haploid.insert(tk.END,
                f"❌ Error saving diversity results: {e}\n")

        self.progress_bar_haploid['value'] = 35
        self.root.update_idletasks()

        if len(sample_names) >= 2:
            self.results_text_haploid.insert(tk.END,
                "\n🔗 GENERATING SIMILARITY MATRIX AND DENDROGRAM...\n")
            self.root.update_idletasks()

            try:
                matrix = pd.DataFrame(index=sample_names, columns=sample_names, dtype=float)
                np.fill_diagonal(matrix.values, 1.0)

                total_pairs = len(data) * (len(data) - 1) / 2
                pairs_done = 0

                for i in range(len(data)):
                    for j in range(i + 1, len(data)):
                        sim = self.calcular_similaridade_haploid(
                            data.iloc[i].values, data.iloc[j].values)
                        matrix.iat[i, j] = sim
                        matrix.iat[j, i] = sim
                        pairs_done += 1
                        progress = 35 + (pairs_done / total_pairs * 25)
                        self.progress_bar_haploid['value'] = progress
                        self.root.update_idletasks()

                matrix_filename = os.path.join(
                    output_dir, "haploid_query_similarity_matrix.xlsx")
                with pd.ExcelWriter(matrix_filename, engine='openpyxl') as writer:
                    title_df = pd.DataFrame(
                        {"": ["GMDA - Haploid Query Similarity Matrix"]})
                    title_df.to_excel(writer, sheet_name='Similarity Matrix',
                                      index=False, header=False)
                    matrix.to_excel(writer, sheet_name='Similarity Matrix', startrow=2)
                self.results_text_haploid.insert(tk.END,
                    f"💾 Similarity matrix saved as: {matrix_filename}\n")

                distance_matrix = 1 - matrix.astype(float).values
                np.fill_diagonal(distance_matrix, 0)
                if np.any(np.isnan(distance_matrix)):
                    distance_matrix = np.nan_to_num(distance_matrix, nan=1.0)

                condensed_dist = squareform(distance_matrix, checks=False)
                linkage_matrix = linkage(condensed_dist, method='average')

                self.haploid_query_distance_matrix = condensed_dist
                self.haploid_query_labels = sample_names.copy()
                
                # Calculate cophenetic correlation for query dendrogram
                query_cophenetic = self.calculate_cophenetic_correlation(condensed_dist, linkage_matrix)
                self.haploid_query_cophenetic_correlation = query_cophenetic
                
                self.results_text_haploid.insert(tk.END,
                    f"\n📊 Query Dendrogram Cophenetic Correlation: {query_cophenetic:.4f}\n")
                if query_cophenetic > 0.8:
                    self.results_text_haploid.insert(tk.END,
                        f"   ✅ Excellent representation (c > 0.8)\n")
                elif query_cophenetic > 0.7:
                    self.results_text_haploid.insert(tk.END,
                        f"   ✓ Good representation (c > 0.7)\n")
                elif query_cophenetic > 0.6:
                    self.results_text_haploid.insert(tk.END,
                        f"   ⚠️ Moderate representation (c > 0.6)\n")
                else:
                    self.results_text_haploid.insert(tk.END,
                        f"   ❌ Poor representation (c ≤ 0.6)\n")

                self.progress_bar_haploid['value'] = 70
                self.root.update_idletasks()

                if self.bootstrap_var_haploid.get():
                    self.results_text_haploid.insert(tk.END,
                        "🔄 Calculating bootstrap supports...\n")
                    self.root.update_idletasks()
                    supports = self.calculate_bootstrap_support_haploid(data, linkage_matrix)
                    self.plot_haploid_query_dendrogram(linkage_matrix, sample_names,
                                                       supports, output_dir, query_cophenetic)
                else:
                    self.plot_haploid_query_dendrogram(linkage_matrix, sample_names,
                                                       output_dir=output_dir, cophenetic_corr=query_cophenetic)

                self.progress_bar_haploid['value'] = 80
                self.root.update_idletasks()

            except Exception as e:
                self.results_text_haploid.insert(tk.END,
                    f"❌ Error generating similarity matrix/dendrogram: {e}\n")

        if len(sample_names) >= 3:
            self.results_text_haploid.insert(tk.END, "\n📈 GENERATING PCoA...\n")
            self.root.update_idletasks()

            try:
                distance_matrix_full = np.zeros((len(sample_names), len(sample_names)))
                total_pairs_pcoa = len(sample_names) ** 2
                pairs_done_pcoa = 0

                for i in range(len(sample_names)):
                    for j in range(len(sample_names)):
                        if i != j:
                            sim = self.calcular_similaridade_haploid(
                                data.iloc[i].values, data.iloc[j].values)
                            distance_matrix_full[i, j] = 1 - sim
                        pairs_done_pcoa += 1
                        progress = 80 + (pairs_done_pcoa / total_pairs_pcoa * 15)
                        self.progress_bar_haploid['value'] = progress
                        self.root.update_idletasks()

                if np.any(np.isnan(distance_matrix_full)) or \
                        np.any(np.isinf(distance_matrix_full)):
                    distance_matrix_full = np.nan_to_num(
                        distance_matrix_full, nan=1.0, posinf=1.0, neginf=0.0)

                distance_matrix_full = (
                    distance_matrix_full + distance_matrix_full.T) / 2
                np.fill_diagonal(distance_matrix_full, 0)
                distance_matrix_full += np.eye(len(distance_matrix_full)) * 1e-10

                mds = MDS(n_components=2, dissimilarity='precomputed',
                          random_state=42, max_iter=300)
                coords = mds.fit_transform(distance_matrix_full)

                self.haploid_query_pcoa_coords = coords

                total_inertia = np.sum(distance_matrix_full ** 2) / (
                    2 * len(sample_names))
                inertia_per_component = [
                    np.sum((coords[:, i] - coords[:, i].mean()) ** 2)
                    for i in range(2)
                ]
                var_explained = [
                    inertia / total_inertia * 100
                    for inertia in inertia_per_component
                ]

                plt.figure(figsize=(10, 8))
                for i, (x, y) in enumerate(coords):
                    plt.scatter(x, y, color='steelblue', s=100, alpha=0.7)
                    plt.text(x, y, str(sample_names[i]), fontsize=10,
                             ha='center', va='bottom')
                plt.xlabel(f"PCo1 ({var_explained[0]:.1f}%)", fontsize=14)
                plt.ylabel(f"PCo2 ({var_explained[1]:.1f}%)", fontsize=14)
                plt.title("GMDA - Haploid Query PCoA", fontsize=14)
                plt.grid(True, alpha=0.3, linestyle='--')
                plt.axhline(0, color='gray', linewidth=0.5)
                plt.axvline(0, color='gray', linewidth=0.5)
                plt.tight_layout()

                pcoa_path = os.path.join(output_dir, "haploid_query_PCoA.png")
                plt.savefig(pcoa_path, dpi=300, bbox_inches='tight')
                plt.close()

                self.results_text_haploid.insert(tk.END,
                    f"📈 PCoA plot saved as: {pcoa_path}\n")
                self.results_text_haploid.insert(tk.END,
                    f"📊 Variance explained: PCo1 = {var_explained[0]:.1f}%, "
                    f"PCo2 = {var_explained[1]:.1f}%\n")

                pcoa_data = pd.DataFrame({
                    'Sample': sample_names,
                    'PCo1': coords[:, 0],
                    'PCo2': coords[:, 1],
                    'Has_Missing': [
                        sample in missing_by_sample for sample in sample_names
                    ] if missing_by_sample else [False] * len(sample_names)
                })
                pcoa_excel_path = os.path.join(
                    output_dir, "haploid_query_PCoA_coordinates.xlsx")
                pcoa_data.to_excel(pcoa_excel_path, index=False)
                self.results_text_haploid.insert(tk.END,
                    f"💾 PCoA coordinates saved as: {pcoa_excel_path}\n")

                self.progress_bar_haploid['value'] = 98
                self.root.update_idletasks()

            except Exception as e:
                self.results_text_haploid.insert(tk.END,
                    f"❌ Error generating PCoA: {e}\n")

        self.progress_bar_haploid['value'] = 100
        self.root.update_idletasks()
        self.results_text_haploid.insert(tk.END,
            "\n✅ COMPLETE QUERY ANALYSIS (HAPLOID) FINISHED!\n")
        self.update_status("Complete query analysis for haploid markers completed")

    def plot_haploid_query_dendrogram(self, linkage_matrix, labels, supports=None,
                                      output_dir=None, cophenetic_corr=None):
        if output_dir is None:
            output_dir = self.create_output_directory("Genetic_haploid")

        plt.figure(figsize=(12, 8))
        dendro = dendrogram(linkage_matrix, labels=labels,
                            leaf_rotation=90, leaf_font_size=10)

        if supports:
            icoord = np.array(dendro['icoord'])
            dcoord = np.array(dendro['dcoord'])
            for i, d in zip(icoord, dcoord):
                x = 0.5 * (i[1] + i[2])
                y = d[1]
                for cluster_idx, row in enumerate(linkage_matrix):
                    if abs(row[2] - y) < 1e-5:
                        support = supports.get(cluster_idx + len(labels), 0)
                        if support >= 50:
                            plt.text(x, y, f"{int(support)}%",
                                     ha='center', va='bottom',
                                     fontsize=8, color='red')
                        break

        title = "GMDA - Haploid Query UPGMA Dendrogram"
        if supports:
            title += " (with bootstrap)"
        if cophenetic_corr is not None:
            title += f"\nCophenetic Correlation: {cophenetic_corr:.4f}"
        plt.title(title, fontsize=14)
        plt.ylabel("Dissimilarity (1-S)", fontsize=14)
        plt.tight_layout()

        filename = "haploid_query_dendrogram"
        if supports:
            filename += "_with_bootstrap"
        plot_path = os.path.join(output_dir, f"{filename}.png")
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        plt.close()
        self.results_text_haploid.insert(tk.END,
            f"📊 Query dendrogram saved as: {plot_path}\n")

    def save_haploid_results(self):
        if not self.haploid_top3_results:
            messagebox.showinfo("No data", "No results available to save.")
            return

        output_dir = self.create_output_directory("Genetic_haploid")

        df_results = pd.DataFrame(
            self.haploid_top3_results, columns=["Query", "Reference", "Similarity"])
        df_results["Similarity (%)"] = pd.to_numeric(
            df_results["Similarity"], errors='coerce') * 100
        df_results.drop(columns=["Similarity"], inplace=True)

        save_path = os.path.join(output_dir, "Top3_haploid.tsv")

        try:
            with open(save_path, 'w', encoding='utf-8') as f:
                f.write("GMDA - Haploid Genetic Similarity Results\n\n")
                df_results.to_csv(f, sep="\t", index=False)
            self.results_text_haploid.insert(tk.END,
                f"\n💾 Results saved as TSV: {save_path}\n")
        except Exception as e:
            self.results_text_haploid.insert(tk.END, f"❌ Error saving TSV: {e}\n")

    # ==================== MORPHO TAB ====================

    def setup_morpho_tab(self):
        main_container = ModernFrame(self.morpho_frame)
        main_container.pack(fill="both", expand=True, padx=10, pady=10)

        input_card = CardFrame(main_container)
        input_card.pack(fill="x", pady=(0, 15))

        input_title = ModernLabel(input_card, text="Morphometric Data Input", font=("Segoe UI", 12, "bold"))
        input_title.pack(anchor="w", pady=(0, 15))

        file_frame = ModernFrame(input_card)
        file_frame.pack(fill="x", pady=8)

        label_file = ModernLabel(file_frame, text="Data file:", font=("Segoe UI", 10, "bold"))
        label_file.grid(row=0, column=0, sticky="w", padx=(0, 10))

        self.entry_file_morpho = ModernEntry(file_frame, width=50)
        self.entry_file_morpho.grid(row=0, column=1, padx=5, sticky="ew")

        btn_browse = ModernButton(file_frame, text="Browse",
                                command=lambda: self.load_file(self.entry_file_morpho),
                                bg=MODERN_THEME["accent"], fg="white")
        btn_browse.grid(row=0, column=2, padx=5)

        params_frame = ModernFrame(input_card)
        params_frame.pack(fill="x", pady=10)

        label_replicates = ModernLabel(params_frame, text="Replicates per sample:", font=("Segoe UI", 10))
        label_replicates.grid(row=0, column=0, sticky="w", padx=(0, 10))

        self.entry_replicates = ModernEntry(params_frame, textvariable=self.replicates_var, width=10)
        self.entry_replicates.grid(row=0, column=1, sticky="w", padx=5)

        self.bootstrap_var_morpho = tk.BooleanVar()
        check_bootstrap = tk.Checkbutton(params_frame, text="Use Bootstrap for Statistics",
                                       variable=self.bootstrap_var_morpho,
                                       bg=MODERN_THEME["card_bg"],
                                       font=("Segoe UI", 10),
                                       selectcolor=MODERN_THEME["light"])
        check_bootstrap.grid(row=0, column=2, sticky="w", padx=20)

        file_frame.columnconfigure(1, weight=1)

        action_card = CardFrame(main_container)
        action_card.pack(fill="x", pady=(0, 15))

        action_title = ModernLabel(action_card, text="Analysis Tools", font=("Segoe UI", 12, "bold"))
        action_title.pack(anchor="w", pady=(0, 15))

        btn_frame = ModernFrame(action_card)
        btn_frame.pack(fill="x")

        btn_analyze = ModernButton(btn_frame,
                                 text="📊 Analyze Data",
                                 command=self.calculate_morpho,
                                 bg=MODERN_THEME["success"], fg="white")
        btn_analyze.grid(row=0, column=0, padx=5, pady=5, sticky="ew")

        btn_save = ModernButton(btn_frame,
                              text="💾 Save Results",
                              command=self.save_morpho_results,
                              bg=MODERN_THEME["accent"], fg="white")
        btn_save.grid(row=0, column=1, padx=5, pady=5, sticky="ew")

        btn_frame.columnconfigure(0, weight=1)
        btn_frame.columnconfigure(1, weight=1)

        progress_frame = ModernFrame(action_card)
        progress_frame.pack(fill="x", pady=(10, 0))

        self.progress_label_morpho = ModernLabel(progress_frame, text="Progress:", font=("Segoe UI", 9, "bold"))
        self.progress_label_morpho.pack(anchor="w")

        self.progress_bar_morpho = ttk.Progressbar(progress_frame, mode='determinate', length=100)
        self.progress_bar_morpho.pack(fill="x", pady=(5, 0))

        results_card = CardFrame(main_container)
        results_card.pack(fill="both", expand=True)

        results_title = ModernLabel(results_card, text="Analysis Results", font=("Segoe UI", 12, "bold"))
        results_title.pack(anchor="w", pady=(0, 10))

        text_frame = ModernFrame(results_card)
        text_frame.pack(fill="both", expand=True)

        self.results_text_morpho = scrolledtext.ScrolledText(text_frame,
                                         bg="white",
                                         fg=MODERN_THEME["text"],
                                         font=("Consolas", 9),
                                         relief="flat",
                                         padx=10,
                                         pady=10)
        self.results_text_morpho.pack(fill="both", expand=True)

    def calculate_morpho(self):
        self.update_status("Performing morphometric analysis...")
        self.progress_bar_morpho['value'] = 0
        self.progress_bar_morpho['maximum'] = 100
        self.root.update_idletasks()

        filename = self.entry_file_morpho.get().strip()
        if not filename:
            messagebox.showwarning("Warning", "Please select a data file.")
            self.update_status("Ready")
            return

        try:
            df = pd.read_excel(filename, engine='openpyxl', header=0)
        except Exception as e:
            messagebox.showerror("Error", f"Error reading file:\n{e}")
            self.update_status("Error reading file")
            return

        if df.empty:
            messagebox.showerror("Error", "The file is empty.")
            self.update_status("Empty file")
            return

        sample_col_name = df.columns[0]

        if pd.isna(sample_col_name) or str(sample_col_name).strip() == '':
            messagebox.showerror("Error", "The first column header cannot be empty")
            self.update_status("Invalid file format")
            return

        try:
            n_replicates = int(self.replicates_var.get())
            if n_replicates < 1:
                raise ValueError("Number of replicates must be at least 1")
        except ValueError as e:
            messagebox.showerror("Error", f"Invalid number of replicates: {e}")
            self.update_status("Invalid replicates")
            return

        self.progress_bar_morpho['value'] = 5
        self.root.update_idletasks()

        output_dir = self.create_output_directory("Morphometric")
        variable_names = df.columns[1:]

        sample_names = []
        for i in range(0, len(df), n_replicates):
            sample_name = df.iloc[i, 0]
            if pd.isna(sample_name):
                continue
            sample_names.append(str(sample_name).strip())

        measurements = defaultdict(list)
        for i in range(0, len(df), n_replicates):
            sample_name = str(df.iloc[i, 0]).strip()
            if pd.isna(sample_name):
                continue
            for var in variable_names:
                replicates = []
                for j in range(n_replicates):
                    if i + j < len(df):
                        value = df.iloc[i + j][var]
                        if self.isanumber(value):
                            replicates.append(float(value))
                measurements[var].append(replicates)

        self.results_text_morpho.delete(1.0, tk.END)
        self.results_text_morpho.insert(tk.END, "📊 MORPHOMETRIC ANALYSIS RESULTS\n")
        self.results_text_morpho.insert(tk.END, "=" * 50 + "\n\n")

        file_basename = os.path.splitext(os.path.basename(filename))[0]
        self.results_text_morpho.insert(tk.END, f"📁 File: {file_basename}\n")
        self.results_text_morpho.insert(tk.END, f"🔢 Samples: {len(sample_names)}\n")
        self.results_text_morpho.insert(tk.END, f"🔄 Replicates: {n_replicates}\n")
        self.results_text_morpho.insert(tk.END, f"📋 Sample column: '{sample_col_name}'\n\n")

        self.morpho_results = []
        self.anova_results = []
        self.tukey_results = []
        self.effect_size_results = []

        mean_values = {}
        for var in variable_names:
            mean_values[var] = {}
            for i, sample in enumerate(sample_names):
                if measurements[var][i]:
                    mean_values[var][sample] = round(np.mean(measurements[var][i]), 2)
                else:
                    mean_values[var][sample] = np.nan

        # Salvar médias para exportação automática
        means_list = []
        for sample in sample_names:
            row = {'Sample': sample}
            for var in variable_names:
                row[var] = mean_values[var].get(sample, np.nan)
            means_list.append(row)
        self.morpho_means = means_list

        total_traits = len(variable_names)
        traits_processed = 0

        for var in variable_names:
            data_for_anova = []
            group_labels = []

            for i, sample_measures in enumerate(measurements[var]):
                if sample_measures:
                    data_for_anova.extend(sample_measures)
                    group_labels.extend([sample_names[i]] * len(sample_measures))

            if not data_for_anova or len(set(group_labels)) < 2:
                continue

            stats_result = self.calculate_statistics(data_for_anova)

            effect_sizes = {}
            ci_lower = {}
            ci_upper = {}

            if len(measurements[var]) >= 2:
                for i in range(len(measurements[var])):
                    for j in range(i + 1, len(measurements[var])):
                        if measurements[var][i] and measurements[var][j]:
                            group1 = measurements[var][i]
                            group2 = measurements[var][j]
                            n1, n2 = len(group1), len(group2)
                            mean1, mean2 = np.mean(group1), np.mean(group2)
                            pooled_std = np.sqrt(
                                ((n1-1)*np.std(group1)**2 + (n2-1)*np.std(group2)**2)
                                / (n1 + n2 - 2))
                            d = (mean1 - mean2) / pooled_std if pooled_std > 0 else 0.0
                            se = np.sqrt((n1 + n2)/(n1*n2) + d**2/(2*(n1+n2)))
                            ci_low = d - 1.96 * se
                            ci_high = d + 1.96 * se
                            pair_name = f"{sample_names[i]}-{sample_names[j]}"
                            effect_sizes[pair_name] = round(d, 2)
                            ci_lower[pair_name] = round(ci_low, 2)
                            ci_upper[pair_name] = round(ci_high, 2)

            if effect_sizes:
                for pair, d in effect_sizes.items():
                    self.effect_size_results.append({
                        'Trait': var,
                        'Comparison': pair,
                        "Cohen's d": d,
                        'CI 95% Lower': ci_lower[pair],
                        'CI 95% Upper': ci_upper[pair]
                    })

            self.results_text_morpho.insert(tk.END, f"📏 Trait: {var}\n")
            self.results_text_morpho.insert(tk.END,
                f"  • Mean: {stats_result['mean']:.2f} ± {stats_result['mean_std']:.2f}\n")
            self.results_text_morpho.insert(tk.END,
                f"  • Median: {stats_result['median']:.2f} ± {stats_result['median_std']:.2f}\n")

            if isinstance(stats_result['mode'], str):
                self.results_text_morpho.insert(tk.END,
                    f"  • Mode: {stats_result['mode']}\n")
            else:
                self.results_text_morpho.insert(tk.END,
                    f"  • Mode: {stats_result['mode']:.2f} ± {stats_result['mode_std']:.2f}\n")

            if effect_sizes:
                self.results_text_morpho.insert(tk.END,
                    "  • Effect sizes (Cohen's d) between groups:\n")
                for pair, d in effect_sizes.items():
                    self.results_text_morpho.insert(tk.END,
                        f"    ↳ {pair}: d = {d:.2f} "
                        f"(95% CI: {ci_lower[pair]:.2f} to {ci_upper[pair]:.2f})\n")

            p_value = None
            try:
                if len(measurements[var]) > 1:
                    groups = [m for m in measurements[var] if m]
                    if len(groups) >= 2:
                        f_val, p_val = stats.f_oneway(*groups)
                        anova_result = {'F': f_val, 'p': p_val}
                    else:
                        anova_result = None
                else:
                    anova_result = None

                if anova_result:
                    p_value = anova_result['p']
                    self.results_text_morpho.insert(tk.END,
                        f"  • ANOVA p-value: {p_value:.4f}\n")

                    if p_value < 0.05:
                        self.results_text_morpho.insert(tk.END,
                            "    ↳ ✅ Significant differences found (p < 0.05)\n")

                        tukey = pairwise_tukeyhsd(
                            endog=data_for_anova,
                            groups=group_labels,
                            alpha=0.05)

                        tukey_df = pd.DataFrame(
                            data=tukey._results_table.data[1:],
                            columns=tukey._results_table.data[0])

                        groups_tukey = self._assign_tukey_groups(
                            tukey_df, sample_names)

                        tukey_formatted = {
                            'Trait': var,
                            'Results': [
                                {'Sample': sample,
                                 'Mean': mean_values[var][sample],
                                 'Group': groups_tukey.get(sample, '?')}
                                for sample in sample_names
                            ]
                        }
                        self.tukey_results.append(tukey_formatted)

                        self.results_text_morpho.insert(tk.END,
                            "    ↳ Tukey HSD grouping:\n")
                        for result in tukey_formatted['Results']:
                            self.results_text_morpho.insert(tk.END,
                                f"      {result['Sample']}: "
                                f"{result['Mean']:.2f}{result['Group']}\n")

                        self.anova_results.append({
                            'Trait': var,
                            'ANOVA_p_value': p_value,
                            'Tukey_results': tukey_df.to_string()
                        })
                    else:
                        self.results_text_morpho.insert(tk.END,
                            "    ↳ ❌ No significant differences found (p ≥ 0.05)\n")
                        self.anova_results.append({
                            'Trait': var,
                            'ANOVA_p_value': p_value,
                            'Tukey_results': "Not applicable (p ≥ 0.05)"
                        })
                else:
                    self.results_text_morpho.insert(tk.END,
                        "  • ANOVA: Not enough groups for analysis\n")
                    self.anova_results.append({
                        'Trait': var,
                        'ANOVA_p_value': "N/A",
                        'Tukey_results': "Not enough groups for analysis"
                    })
            except Exception as e:
                self.results_text_morpho.insert(tk.END,
                    f"  • ❌ Error in ANOVA: {str(e)}\n")
                self.anova_results.append({
                    'Trait': var,
                    'ANOVA_p_value': "Error",
                    'Tukey_results': f"Error: {str(e)}"
                })

            self.results_text_morpho.insert(tk.END,
                f"  • Count: {stats_result['count']}\n\n")

            self.morpho_results.append({
                'Trait': var,
                'Count': stats_result['count'],
                'Mean': f"{stats_result['mean']:.2f} ± {stats_result['mean_std']:.2f}",
                'Median': f"{stats_result['median']:.2f} ± {stats_result['median_std']:.2f}",
                'Mode': stats_result['mode'],
                'ANOVA_p_value': f"{p_value:.4f}" if p_value is not None else "N/A"
            })

            traits_processed += 1
            progress = 5 + (traits_processed / total_traits * 45)
            self.progress_bar_morpho['value'] = progress
            self.root.update_idletasks()

        try:
            data_for_clustering = []
            for var in variable_names:
                var_means = [
                    np.mean(m) if m else np.nan
                    for m in measurements[var]
                ]
                data_for_clustering.append(var_means)

            df_clustering = pd.DataFrame(
                data_for_clustering, index=variable_names, columns=sample_names).T

            total_missing = df_clustering.isnull().sum().sum()
            missing_by_sample = df_clustering.isnull().sum(axis=1).to_dict()
            missing_by_trait = df_clustering.isnull().sum(axis=0).to_dict()

            if total_missing > 0:
                self.results_text_morpho.insert(tk.END,
                    f"\n⚠️ **WARNING: MISSING DATA DETECTED** ⚠️\n")
                self.results_text_morpho.insert(tk.END,
                    f"   Total missing values: {total_missing}\n\n")
                self.results_text_morpho.insert(tk.END,
                    "   📊 Distribution by sample:\n")
                for sample, n_missing in missing_by_sample.items():
                    if n_missing > 0:
                        self.results_text_morpho.insert(tk.END,
                            f"      • {sample}: {n_missing} missing values\n")
                self.results_text_morpho.insert(tk.END,
                    "\n   📊 Distribution by trait:\n")
                for trait, n_missing in missing_by_trait.items():
                    if n_missing > 0:
                        self.results_text_morpho.insert(tk.END,
                            f"      • {trait}: {n_missing} missing values\n")
                self.results_text_morpho.insert(tk.END,
                    "\n   🔧 Applying mean imputation...\n")
                self.root.update_idletasks()

                df_clustering = df_clustering.copy()
                for col in df_clustering.columns:
                    if df_clustering[col].isnull().any():
                        mean_val = df_clustering[col].mean(skipna=True)
                        df_clustering[col] = df_clustering[col].fillna(mean_val)
                        self.results_text_morpho.insert(tk.END,
                            f"      • {col}: imputed with mean = {mean_val:.2f}\n")

            self.progress_bar_morpho['value'] = 55
            self.root.update_idletasks()

            if len(df_clustering) > 1:
                morpho_result = self.generate_morpho_dendrogram(
                    df_clustering.values,
                    df_clustering.index.tolist(),
                    bootstrap=self.bootstrap_var_morpho.get(),
                    output_dir=output_dir,
                    missing_info=missing_by_sample if total_missing > 0 else None)

                if morpho_result:
                    self.morpho_linkage_matrix = morpho_result['linkage_matrix']
                    self.morpho_distance_matrix = morpho_result['distance_matrix']
                    self.morpho_labels = df_clustering.index.tolist()
                    self.morpho_query_distance_matrix = morpho_result['distance_matrix']
                    self.morpho_query_labels = df_clustering.index.tolist()
                    self.morpho_cophenetic_correlation = morpho_result.get('cophenetic_correlation')
                    self.morpho_query_cophenetic_correlation = morpho_result.get('cophenetic_correlation')

                self.progress_bar_morpho['value'] = 75
                self.root.update_idletasks()

                pca_result = self.generate_pca(
                    df_clustering.copy(),
                    df_clustering.index.tolist(),
                    output_dir,
                    missing_info=missing_by_sample if total_missing > 0 else None)
                if pca_result:
                    self.morpho_pca_coords = pca_result['pca_result']
                    self.morpho_query_pca_coords = pca_result['pca_result']

                self.progress_bar_morpho['value'] = 90
                self.root.update_idletasks()
            else:
                self.results_text_morpho.insert(tk.END,
                    "\n⚠️ Not enough valid samples for dendrogram and PCA\n")

        except Exception as e:
            self.results_text_morpho.insert(tk.END, f"\n❌ Error in analysis: {e}\n")

        self.progress_bar_morpho['value'] = 100
        self.root.update_idletasks()
        self.update_status("Morphometric analysis completed")

    def _assign_tukey_groups(self, tukey_df, sample_names):
        groups = defaultdict(set)
        for _, row in tukey_df.iterrows():
            group1 = row['group1']
            group2 = row['group2']
            reject = row['reject']
            if not reject:
                found = False
                for g in groups.values():
                    if group1 in g or group2 in g:
                        g.add(group1)
                        g.add(group2)
                        found = True
                        break
                if not found:
                    groups[len(groups)] = {group1, group2}

        all_samples_in_groups = set()
        for g in groups.values():
            all_samples_in_groups.update(g)

        for sample in sample_names:
            if sample not in all_samples_in_groups:
                groups[len(groups)] = {sample}

        letters = string.ascii_lowercase
        group_letters = {}
        for i, (_, samples) in enumerate(groups.items()):
            letter = letters[i] if i < len(letters) else f"z{i - len(letters) + 1}"
            for sample in samples:
                group_letters[sample] = letter

        return group_letters

    def save_morpho_results(self):
        if not self.morpho_results:
            messagebox.showinfo("No data", "No results available to save.")
            return

        output_dir = self.create_output_directory("Morphometric")
        df_stats = pd.DataFrame(self.morpho_results)
        df_stats = df_stats[['Trait', 'Count', 'Mean', 'Median', 'Mode', 'ANOVA_p_value']]

        if self.tukey_results:
            sample_data = {}
            for result in self.tukey_results:
                trait = result['Trait']
                for res in result['Results']:
                    sample = res['Sample']
                    if sample not in sample_data:
                        sample_data[sample] = {'Sample': sample}
                    sample_data[sample][f"{trait}_Mean"] = \
                        f"{res['Mean']:.2f}{res['Group']}"

            df_tukey = pd.DataFrame.from_dict(sample_data, orient='index')
            df_tukey = df_tukey.sort_index().reset_index(drop=True)
            cols = ['Sample'] + [c for c in df_tukey.columns if c != 'Sample']
            df_tukey = df_tukey[cols]
        else:
            df_tukey = pd.DataFrame()

        if self.effect_size_results:
            df_effect_sizes = pd.DataFrame(self.effect_size_results)
            df_effect_sizes = df_effect_sizes[
                ['Trait', 'Comparison', "Cohen's d", 'CI 95% Lower', 'CI 95% Upper']]
        else:
            df_effect_sizes = pd.DataFrame(
                columns=['Trait', 'Comparison', "Cohen's d",
                         'CI 95% Lower', 'CI 95% Upper'])

        save_path = os.path.join(output_dir, "Morphometric_results.xlsx")
        try:
            with pd.ExcelWriter(save_path, engine='openpyxl') as writer:
                title_df = pd.DataFrame({"": ["GMDA - Morphometric Statistics"]})
                title_df.to_excel(writer, sheet_name='Statistics',
                                  index=False, header=False)
                df_stats.to_excel(writer, sheet_name='Statistics',
                                  startrow=2, index=False)

                if not df_tukey.empty:
                    title_df = pd.DataFrame(
                        {"": ["GMDA - Tukey HSD Grouping (Mean±Group)"]})
                    title_df.to_excel(writer, sheet_name='Tukey_HSD',
                                      index=False, header=False)
                    df_tukey.to_excel(writer, sheet_name='Tukey_HSD',
                                      startrow=2, index=False)

                title_df = pd.DataFrame(
                    {"": ["GMDA - Effect Sizes and Confidence Intervals"]})
                title_df.to_excel(writer, sheet_name='Effect_Sizes',
                                  index=False, header=False)
                df_effect_sizes.to_excel(writer, sheet_name='Effect_Sizes',
                                         startrow=2, index=False)

            self.results_text_morpho.insert(tk.END,
                f"\n💾 Results saved as Excel: {save_path}\n")
        except Exception as e:
            self.results_text_morpho.insert(tk.END, f"❌ Error saving Excel: {e}\n")

    def calculate_statistics(self, data):
        if not data:
            return {
                'count': 0, 'mean': "N/A", 'mean_std': "N/A",
                'median': "N/A", 'median_std': "N/A",
                'mode': "N/A", 'mode_std': "N/A"
            }

        data = np.array(data)
        stats_result = {
            'count': len(data),
            'mean': round(np.mean(data), 2),
            'mean_std': round(np.std(data), 2),
            'median': round(np.median(data), 2),
            'median_std': round(np.std(data - np.median(data)), 2)
        }

        try:
            modes = multimode(data)
            if len(modes) == 1:
                stats_result['mode'] = round(modes[0], 2)
                stats_result['mode_std'] = round(np.std(data - modes[0]), 2)
            else:
                stats_result['mode'] = "multiple"
                stats_result['mode_std'] = "N/A"
        except Exception:
            stats_result['mode'] = "N/A"
            stats_result['mode_std'] = "N/A"

        return stats_result

    def generate_morpho_dendrogram(self, data, labels, bootstrap=False,
                                   replicates=100, output_dir=None,
                                   missing_info=None):
        if output_dir is None:
            output_dir = self.create_output_directory("Morphometric")

        supports = None
        cophenetic_corr = None

        try:
            data = data.copy()

            if np.any(np.isnan(data)):
                self.results_text_morpho.insert(tk.END,
                    f"\n⚠️ WARNING: Applying final imputation...\n")
                col_means = np.nanmean(data, axis=0)
                inds = np.where(np.isnan(data))
                data[inds] = np.take(col_means, inds[1])

            scaler = StandardScaler()
            data_scaled = scaler.fit_transform(data)

            dist_matrix = pdist(data_scaled, metric='euclidean')
            linkage_matrix = linkage(dist_matrix, method='average')
            
            # Calculate cophenetic correlation
            cophenetic_corr = self.calculate_cophenetic_correlation(dist_matrix, linkage_matrix)
            self.morpho_cophenetic_correlation = cophenetic_corr

            plt.figure(figsize=(12, 8))
            dendro = dendrogram(linkage_matrix, labels=labels, orientation='top',
                                leaf_rotation=90, leaf_font_size=16,
                                color_threshold=0.7*np.max(linkage_matrix[:, 2]))

            if bootstrap:
                supports = self.calculate_bootstrap_support_morpho(
                    data, linkage_matrix, replicates=replicates)

                icoord = np.array(dendro['icoord'])
                dcoord = np.array(dendro['dcoord'])
                for i, d in zip(icoord, dcoord):
                    x = 0.5 * (i[1] + i[2])
                    y = d[1]
                    for cluster_idx, row in enumerate(linkage_matrix):
                        if abs(row[2] - y) < 1e-5:
                            support = supports.get(cluster_idx + len(labels), 0)
                            if support >= 50:
                                plt.text(x, y, f"{int(support)}%",
                                         ha='center', va='bottom',
                                         fontsize=8, color='red')
                            break

            title = "GMDA - Morphometric UPGMA Dendrogram"
            if bootstrap:
                title += f" (Bootstrap n={replicates})"
            title += f"\nCophenetic Correlation: {cophenetic_corr:.4f}"

            plt.title(title, fontsize=14, pad=20)
            plt.xlabel("Samples", fontsize=12)
            plt.ylabel("Euclidean Distance", fontsize=12)
            plt.grid(True, alpha=0.3, linestyle='--')
            plt.tight_layout()

            filename = "Morphometric_Dendrogram"
            if bootstrap:
                filename += "_with_bootstrap"
            plot_path = os.path.join(output_dir, f"{filename}.png")
            plt.savefig(plot_path, dpi=300, bbox_inches='tight')
            plt.close()

            self.results_text_morpho.insert(tk.END,
                f"\n📊 Dendrogram saved as: {plot_path}\n")
            self.results_text_morpho.insert(tk.END,
                f"📊 Cophenetic Correlation Coefficient: {cophenetic_corr:.4f}\n")
            if cophenetic_corr > 0.8:
                self.results_text_morpho.insert(tk.END,
                    f"   ✅ Excellent representation (c > 0.8)\n")
            elif cophenetic_corr > 0.7:
                self.results_text_morpho.insert(tk.END,
                    f"   ✓ Good representation (c > 0.7)\n")
            elif cophenetic_corr > 0.6:
                self.results_text_morpho.insert(tk.END,
                    f"   ⚠️ Moderate representation (c > 0.6)\n")
            else:
                self.results_text_morpho.insert(tk.END,
                    f"   ❌ Poor representation (c ≤ 0.6)\n")

            return {
                'linkage_matrix': linkage_matrix,
                'distance_matrix': dist_matrix,
                'bootstrap_supports': supports,
                'plot_path': plot_path,
                'cophenetic_correlation': cophenetic_corr
            }

        except Exception as e:
            plt.close()
            self.results_text_morpho.insert(tk.END,
                f"\n❌ Error generating dendrogram: {e}\n")
            return None

    def calculate_bootstrap_support_morpho(self, data, linkage_matrix, replicates=100):
        n = data.shape[0]
        clusters_suporte = defaultdict(int)
        original_clusters = self.get_cluster_members(linkage_matrix, n)

        for rep in range(replicates):
            try:
                sampled_cols = np.random.choice(
                    data.shape[1], size=data.shape[1], replace=True)
                data_boot = data[:, sampled_cols]

                scaler = StandardScaler()
                data_boot_scaled = scaler.fit_transform(data_boot)

                dist_boot = pdist(data_boot_scaled, metric='euclidean')
                linkage_boot = linkage(dist_boot, method='average')
                bootstrap_clusters = self.get_cluster_members(linkage_boot, n)

                for cluster_id, members in original_clusters.items():
                    if any(members == b_members
                           for b_members in bootstrap_clusters.values()):
                        clusters_suporte[cluster_id] += 1

                progress = 55 + (rep / replicates * 15)
                self.progress_bar_morpho['value'] = progress
                self.root.update_idletasks()

            except Exception:
                continue

        for c in clusters_suporte:
            clusters_suporte[c] = clusters_suporte[c] / replicates * 100

        return clusters_suporte

    def generate_pca(self, data, sample_names, output_dir=None, missing_info=None):
        if output_dir is None:
            output_dir = self.create_output_directory("Morphometric")

        try:
            data = data.copy()

            if data.isnull().any().any():
                self.results_text_morpho.insert(tk.END,
                    f"\n⚠️ WARNING: Applying final imputation for PCA...\n")
                for col in data.columns:
                    if data[col].isnull().any():
                        mean_val = data[col].mean()
                        data[col] = data[col].fillna(mean_val)

            scaler = StandardScaler()
            data_scaled = scaler.fit_transform(data)

            pca = PCA()
            pca_result = pca.fit_transform(data_scaled)
            explained_var = pca.explained_variance_ratio_

            plt.figure(figsize=(10, 8))
            plt.scatter(pca_result[:, 0], pca_result[:, 1],
                        alpha=0.7, s=100, c='steelblue')

            for i, name in enumerate(sample_names):
                plt.annotate(str(name), (pca_result[i, 0], pca_result[i, 1]),
                             xytext=(5, 5), textcoords='offset points', fontsize=16)

            loadings = pca.components_.T * np.sqrt(pca.explained_variance_)
            for i, var in enumerate(data.columns):
                plt.arrow(0, 0, loadings[i, 0]*3, loadings[i, 1]*3,
                          color='red', alpha=0.7, head_width=0.1)
                plt.text(loadings[i, 0]*3.2, loadings[i, 1]*3.2, var,
                         color='red', ha='center', va='center', fontsize=10)

            pc1_var = explained_var[0] * 100
            pc2_var = explained_var[1] * 100

            plt.xlabel(f'PC1 ({pc1_var:.1f}%)', fontsize=16)
            plt.ylabel(f'PC2 ({pc2_var:.1f}%)', fontsize=16)
            plt.title('GMDA - Principal Component Analysis', fontsize=14)
            plt.grid(True, alpha=0.3)
            plt.axhline(0, color='gray', linewidth=0.5)
            plt.axvline(0, color='gray', linewidth=0.5)
            plt.tight_layout()

            pca_path = os.path.join(output_dir, "Morphometric_PCA.png")
            plt.savefig(pca_path, dpi=300, bbox_inches='tight')
            plt.close()

            self.results_text_morpho.insert(tk.END,
                f"\n📈 PCA plot saved as: {pca_path}\n")
            self.results_text_morpho.insert(tk.END, "\n📊 PCA Variance Explained:\n")
            for i, var in enumerate(explained_var):
                self.results_text_morpho.insert(tk.END, f"  PC{i+1}: {var*100:.1f}%\n")

            return {
                'pca_result': pca_result,
                'explained_variance': explained_var,
                'loadings': loadings,
                'plot_path': pca_path
            }

        except Exception as e:
            plt.close()
            self.results_text_morpho.insert(tk.END,
                f"\n❌ Error generating PCA: {e}\n")
            return None

    # ==================== CORRELATION TAB ====================

    def setup_correlation_tab(self):
        outer_frame = ModernFrame(self.correlation_frame)
        outer_frame.pack(fill="both", expand=True)

        left_panel = ModernFrame(outer_frame)
        left_panel.pack(side="left", fill="y", padx=(10, 5), pady=10)

        canvas = tk.Canvas(left_panel, bg=MODERN_THEME["bg"],
                           width=460, highlightthickness=0)
        scrollbar = ttk.Scrollbar(left_panel, orient="vertical",
                                  command=canvas.yview)
        canvas.configure(yscrollcommand=scrollbar.set)

        scrollbar.pack(side="right", fill="y")
        canvas.pack(side="left", fill="both", expand=True)

        scroll_frame = ModernFrame(canvas)
        scroll_window = canvas.create_window((0, 0), window=scroll_frame, anchor="nw")

        def _on_frame_configure(event):
            canvas.configure(scrollregion=canvas.bbox("all"))

        def _on_canvas_configure(event):
            canvas.itemconfig(scroll_window, width=event.width)

        scroll_frame.bind("<Configure>", _on_frame_configure)
        canvas.bind("<Configure>", _on_canvas_configure)

        def _on_mousewheel(event):
            canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")

        canvas.bind_all("<MouseWheel>", _on_mousewheel)

        info_label = tk.Label(scroll_frame,
                              text="💡 Run genetic and morphometric analyses first",
                              font=("Segoe UI", 9),
                              bg=MODERN_THEME["bg"],
                              fg=MODERN_THEME["warning"],
                              wraplength=420, justify="left")
        info_label.pack(anchor="w", pady=(5, 8), padx=5)

        btn_cfg = dict(relief="flat", bd=0, padx=6, pady=4,
                       font=("Segoe UI", 8), cursor="hand2", fg="white",
                       anchor="w", wraplength=420)

        def make_btn(parent, text, command, color):
            b = tk.Button(parent, text=text, command=command,
                          bg=color, **btn_cfg)
            b.pack(fill="x", padx=4, pady=2)
            return b

        tk.Label(scroll_frame, text="🌳  Dendrogram Comparisons",
                 font=("Segoe UI", 9, "bold"),
                 bg=MODERN_THEME["bg"], fg=MODERN_THEME["text"]).pack(
            anchor="w", padx=5, pady=(8, 2))

        make_btn(scroll_frame,
                 "Dendrogram codominant genetic  vs  Dendrogram morphometric",
                 self.calculate_mantel_test_codominant_vs_morpho_dendrograms,
                 MODERN_THEME["primary"])
        make_btn(scroll_frame,
                 "Dendrogram dominant genetic  vs  Dendrogram morphometric",
                 self.calculate_mantel_test_dominant_vs_morpho_dendrograms,
                 MODERN_THEME["primary"])
        make_btn(scroll_frame,
                 "Dendrogram haploid genetic  vs  Dendrogram morphometric",
                 self.calculate_mantel_test_haploid_vs_morpho_dendrograms,
                 MODERN_THEME["primary"])
        make_btn(scroll_frame,
                 "Dendrogram codominant genetic  vs  Dendrogram dominant genetic",
                 self.calculate_mantel_test_codominant_vs_dominant_dendrograms,
                 MODERN_THEME["warning"])
        make_btn(scroll_frame,
                 "Dendrogram codominant genetic  vs  Dendrogram haploid genetic",
                 self.calculate_mantel_test_codominant_vs_haploid_dendrograms,
                 MODERN_THEME["warning"])
        make_btn(scroll_frame,
                 "Dendrogram dominant genetic  vs  Dendrogram haploid genetic",
                 self.calculate_mantel_test_dominant_vs_haploid_dendrograms,
                 MODERN_THEME["warning"])

        tk.Label(scroll_frame, text="📈  PCoA / PCA Comparisons",
                 font=("Segoe UI", 9, "bold"),
                 bg=MODERN_THEME["bg"], fg=MODERN_THEME["text"]).pack(
            anchor="w", padx=5, pady=(12, 2))

        make_btn(scroll_frame,
                 "PCoA codominant genetic  vs  PCA morphometric",
                 self.calculate_mantel_test_codominant_vs_morpho_pcoa,
                 MODERN_THEME["success"])
        make_btn(scroll_frame,
                 "PCoA dominant genetic  vs  PCA morphometric",
                 self.calculate_mantel_test_dominant_vs_morpho_pcoa,
                 MODERN_THEME["success"])
        make_btn(scroll_frame,
                 "PCoA haploid genetic  vs  PCA morphometric",
                 self.calculate_mantel_test_haploid_vs_morpho_pcoa,
                 MODERN_THEME["success"])
        make_btn(scroll_frame,
                 "PCoA codominant genetic  vs  PCoA dominant genetic",
                 self.calculate_mantel_test_codominant_vs_dominant_pcoa,
                 MODERN_THEME["accent"])
        make_btn(scroll_frame,
                 "PCoA codominant genetic  vs  PCoA haploid genetic",
                 self.calculate_mantel_test_codominant_vs_haploid_pcoa,
                 MODERN_THEME["accent"])
        make_btn(scroll_frame,
                 "PCoA dominant genetic  vs  PCoA haploid genetic",
                 self.calculate_mantel_test_dominant_vs_haploid_pcoa,
                 MODERN_THEME["accent"])

        self.progress_bar = ttk.Progressbar(scroll_frame, mode='indeterminate')
        self.progress_bar.pack(fill="x", padx=4, pady=(12, 6))

        right_panel = ModernFrame(outer_frame)
        right_panel.pack(side="left", fill="both", expand=True,
                         padx=(5, 10), pady=10)

        results_title = ModernLabel(right_panel, text="Correlation Results",
                                    font=("Segoe UI", 12, "bold"))
        results_title.pack(anchor="w", pady=(0, 6))

        self.results_text_corr = scrolledtext.ScrolledText(right_panel,
                                       bg="white",
                                       fg=MODERN_THEME["text"],
                                       font=("Consolas", 9),
                                       relief="flat",
                                       padx=10,
                                       pady=10)
        self.results_text_corr.pack(fill="both", expand=True)

    def perform_mantel_test(self, dist1, dist2, permutations=1000):
        matrix1 = squareform(dist1)
        matrix2 = squareform(dist2)
        result = mantel(matrix1, matrix2, method='pearson', permutations=permutations)
        return result[0], result[1], result[2]

    def _prepare_common_condensed(self, dist_mat1, labels1, dist_mat2, labels2):
        set1 = set(labels1)
        set2 = set(labels2)
        common = list(set1 & set2)

        if len(common) < 3:
            return None, None, None

        idx1 = [i for i, l in enumerate(labels1) if l in common]
        idx2 = [i for i, l in enumerate(labels2) if l in common]

        common_ordered = [labels1[i] for i in idx1]

        sq1 = squareform(dist_mat1)
        sq2 = squareform(dist_mat2)

        sub1 = sq1[np.ix_(idx1, idx1)]
        sub2 = sq2[np.ix_(idx2, idx2)]

        return squareform(sub1, checks=False), squareform(sub2, checks=False), common_ordered

    def _run_mantel_generic(self, dist1_condensed, dist2_condensed,
                            samples_ordered, label1, label2, output_dir):
        corr, p_value, z_score = self.perform_mantel_test(dist1_condensed, dist2_condensed)
        analysis_type = f"{label1}/{label2}"

        self.results_text_corr.delete(1.0, tk.END)
        header = f"🔗 MANTEL TEST RESULTS ({label1.upper()} vs {label2.upper()} - QUERY ONLY)\n"
        self.results_text_corr.insert(tk.END, header)
        self.results_text_corr.insert(tk.END, "=" * 70 + "\n\n")
        self.results_text_corr.insert(tk.END,
            f"📊 Mantel correlation coefficient: {corr:.4f}\n")
        self.results_text_corr.insert(tk.END, f"📈 P-value: {p_value:.4f}\n")
        self.results_text_corr.insert(tk.END, f"📋 Z-score: {z_score:.4f}\n")
        self.results_text_corr.insert(tk.END,
            f"🔢 Number of common Query samples: {len(samples_ordered)}\n")
        self.results_text_corr.insert(tk.END,
            f"📝 Common Query samples: {', '.join(samples_ordered)}\n\n")

        if p_value < 0.05:
            self.results_text_corr.insert(tk.END,
                "✅ The correlation is statistically significant (p < 0.05)\n")
        else:
            self.results_text_corr.insert(tk.END,
                "❌ The correlation is not statistically significant (p ≥ 0.05)\n")

    def save_mantel_results(self, dist1, dist2, corr, p_value, z_score,
                            sample_names, analysis_type, output_dir):
        default_filename = f"Mantel_Test_{analysis_type.replace('/', '_')}.pdf"
        save_path = os.path.join(output_dir, default_filename)

        try:
            with PdfPages(save_path) as pdf:
                fig, ax = plt.subplots(figsize=(8.5, 11))
                ax.axis('off')
                text_content = [
                    f"GMDA - Mantel Test Analysis ({analysis_type})",
                    "=" * 50, "",
                    f"Mantel correlation coefficient: {corr:.4f}",
                    f"P-value: {p_value:.4f}",
                    f"Z-score: {z_score:.4f}",
                    f"Number of Query samples: {len(sample_names)}", "",
                    "Query sample names:", ", ".join(sample_names)
                ]
                ax.text(0.1, 0.9, "\n".join(text_content), transform=ax.transAxes,
                        verticalalignment='top', fontsize=12)
                pdf.savefig(fig, bbox_inches='tight')
                plt.close(fig)

                fig, ax = plt.subplots(figsize=(8, 8))
                ax.scatter(dist1, dist2, alpha=0.7)
                z = np.polyfit(dist1, dist2, 1)
                p = np.poly1d(z)
                ax.plot(dist1, p(dist1), "r--", alpha=0.8)
                ax.set_xlabel(f"{analysis_type.split('/')[0]} distance", fontsize=12)
                ax.set_ylabel(f"{analysis_type.split('/')[-1]} distance", fontsize=12)
                ax.set_title(
                    f"Mantel Test ({analysis_type})\n(r = {corr:.3f}, p = {p_value:.3f})",
                    fontsize=14)
                ax.grid(True, alpha=0.3)
                pdf.savefig(fig, bbox_inches='tight')
                plt.close(fig)

            self.results_text_corr.insert(tk.END,
                f"\n💾 Results saved as PDF: {save_path}\n")
        except Exception as e:
            self.results_text_corr.insert(tk.END,
                f"❌ Error saving results: {str(e)}\n")

    def save_mantel_pca_results(self, dist1, dist2, corr, p_value, z_score,
                                sample_names, analysis_type, output_dir):
        save_path = os.path.join(output_dir,
                                 f"Mantel_Test_{analysis_type.replace('/', '_')}.pdf")
        try:
            with PdfPages(save_path) as pdf:
                fig, ax = plt.subplots(figsize=(8.5, 11))
                ax.axis('off')
                text_content = [
                    f"GMDA - Mantel Test Analysis ({analysis_type})",
                    "=" * 50, "",
                    f"Mantel correlation coefficient: {corr:.4f}",
                    f"P-value: {p_value:.4f}",
                    f"Z-score: {z_score:.4f}",
                    f"Number of Query samples: {len(sample_names)}", "",
                    "Query sample names:", ", ".join(sample_names)
                ]
                ax.text(0.1, 0.9, "\n".join(text_content), transform=ax.transAxes,
                        verticalalignment='top', fontsize=12)
                pdf.savefig(fig, bbox_inches='tight')
                plt.close(fig)

                fig, ax = plt.subplots(figsize=(8, 8))
                ax.scatter(dist1, dist2, alpha=0.7)
                z = np.polyfit(dist1, dist2, 1)
                p = np.poly1d(z)
                ax.plot(dist1, p(dist1), "r--", alpha=0.8)
                ax.set_xlabel(f"{analysis_type.split('/')[0]} distance", fontsize=12)
                ax.set_ylabel(f"{analysis_type.split('/')[-1]} distance", fontsize=12)
                ax.set_title(
                    f"Mantel Test ({analysis_type})\n(r = {corr:.3f}, p = {p_value:.3f})",
                    fontsize=14)
                ax.grid(True, alpha=0.3)
                pdf.savefig(fig, bbox_inches='tight')
                plt.close(fig)

            self.results_text_corr.insert(tk.END,
                f"\n💾 Results saved as PDF: {save_path}\n")
        except Exception as e:
            self.results_text_corr.insert(tk.END,
                f"❌ Error saving results: {str(e)}\n")

    # ==================== MANTEL TEST METHODS ====================

    def calculate_mantel_test_codominant_vs_morpho_dendrograms(self):
        self.update_status("Running Mantel test: Codominant vs Morphometric dendrograms...")
        output_dir = self.create_output_directory("Mantel_correlations")

        if self.genetic_query_distance_matrix is None or \
                self.morpho_query_distance_matrix is None:
            messagebox.showwarning("Warning",
                "Please run both codominant genetic and morphometric Query analyses first.")
            self.update_status("Ready")
            return

        self.progress_bar.start(10)
        labels_g = [str(l) for l in self.genetic_query_labels]
        labels_m = [str(l) for l in self.morpho_query_labels]

        c1, c2, common = self._prepare_common_condensed(
            self.genetic_query_distance_matrix, labels_g,
            self.morpho_query_distance_matrix, labels_m)

        if common is None:
            messagebox.showerror("Error", "Not enough common Query samples.")
            self.progress_bar.stop()
            self.update_status("Not enough common Query samples")
            return

        corr, p_value, z_score = self.perform_mantel_test(c1, c2)
        self.save_mantel_results(c1, c2, corr, p_value, z_score, common,
                                 "Codominant_Markers/Morphometric_Dendrograms_Query",
                                 output_dir)

        self.results_text_corr.delete(1.0, tk.END)
        self.results_text_corr.insert(tk.END,
            "🔗 MANTEL TEST RESULTS (CODOMINANT vs MORPHOMETRIC DENDROGRAMS - QUERY ONLY)\n")
        self.results_text_corr.insert(tk.END, "=" * 70 + "\n\n")
        self.results_text_corr.insert(tk.END,
            f"📊 Mantel correlation coefficient: {corr:.4f}\n")
        self.results_text_corr.insert(tk.END, f"📈 P-value: {p_value:.4f}\n")
        self.results_text_corr.insert(tk.END, f"📋 Z-score: {z_score:.4f}\n")
        self.results_text_corr.insert(tk.END,
            f"🔢 Number of common Query samples: {len(common)}\n")
        self.results_text_corr.insert(tk.END,
            f"📝 Common Query samples: {', '.join(common)}\n\n")
        if p_value < 0.05:
            self.results_text_corr.insert(tk.END,
                "✅ The correlation is statistically significant (p < 0.05)\n")
        else:
            self.results_text_corr.insert(tk.END,
                "❌ The correlation is not statistically significant (p ≥ 0.05)\n")

        self.progress_bar.stop()
        self.update_status("Mantel test completed")

    def calculate_mantel_test_dominant_vs_morpho_dendrograms(self):
        self.update_status("Running Mantel test: Dominant vs Morphometric dendrograms...")
        output_dir = self.create_output_directory("Mantel_correlations")

        if self.dominant_query_distance_matrix is None or \
                self.morpho_query_distance_matrix is None:
            messagebox.showwarning("Warning",
                "Please run both dominant genetic and morphometric Query analyses first.")
            self.update_status("Ready")
            return

        self.progress_bar.start(10)
        labels_d = [str(l) for l in self.dominant_query_labels]
        labels_m = [str(l) for l in self.morpho_query_labels]

        c1, c2, common = self._prepare_common_condensed(
            self.dominant_query_distance_matrix, labels_d,
            self.morpho_query_distance_matrix, labels_m)

        if common is None:
            messagebox.showerror("Error", "Not enough common Query samples.")
            self.progress_bar.stop()
            self.update_status("Not enough common Query samples")
            return

        corr, p_value, z_score = self.perform_mantel_test(c1, c2)
        self.save_mantel_results(c1, c2, corr, p_value, z_score, common,
                                 "Dominant_Markers/Morphometric_Dendrograms_Query",
                                 output_dir)

        self.results_text_corr.delete(1.0, tk.END)
        self.results_text_corr.insert(tk.END,
            "🔗 MANTEL TEST RESULTS (DOMINANT vs MORPHOMETRIC DENDROGRAMS - QUERY ONLY)\n")
        self.results_text_corr.insert(tk.END, "=" * 70 + "\n\n")
        self.results_text_corr.insert(tk.END,
            f"📊 Mantel correlation coefficient: {corr:.4f}\n")
        self.results_text_corr.insert(tk.END, f"📈 P-value: {p_value:.4f}\n")
        self.results_text_corr.insert(tk.END, f"📋 Z-score: {z_score:.4f}\n")
        self.results_text_corr.insert(tk.END,
            f"🔢 Number of common Query samples: {len(common)}\n")
        self.results_text_corr.insert(tk.END,
            f"📝 Common Query samples: {', '.join(common)}\n\n")
        if p_value < 0.05:
            self.results_text_corr.insert(tk.END,
                "✅ The correlation is statistically significant (p < 0.05)\n")
        else:
            self.results_text_corr.insert(tk.END,
                "❌ The correlation is not statistically significant (p ≥ 0.05)\n")

        self.progress_bar.stop()
        self.update_status("Mantel test completed")

    def calculate_mantel_test_haploid_vs_morpho_dendrograms(self):
        self.update_status("Running Mantel test: Haploid vs Morphometric dendrograms...")
        output_dir = self.create_output_directory("Mantel_correlations")

        if self.haploid_query_distance_matrix is None or \
                self.morpho_query_distance_matrix is None:
            messagebox.showwarning("Warning",
                "Please run both haploid genetic and morphometric Query analyses first.")
            self.update_status("Ready")
            return

        self.progress_bar.start(10)
        labels_h = [str(l) for l in self.haploid_query_labels]
        labels_m = [str(l) for l in self.morpho_query_labels]

        c1, c2, common = self._prepare_common_condensed(
            self.haploid_query_distance_matrix, labels_h,
            self.morpho_query_distance_matrix, labels_m)

        if common is None:
            messagebox.showerror("Error", "Not enough common Query samples.")
            self.progress_bar.stop()
            self.update_status("Not enough common Query samples")
            return

        corr, p_value, z_score = self.perform_mantel_test(c1, c2)
        self.save_mantel_results(c1, c2, corr, p_value, z_score, common,
                                 "Haploid_Markers/Morphometric_Dendrograms_Query",
                                 output_dir)

        self.results_text_corr.delete(1.0, tk.END)
        self.results_text_corr.insert(tk.END,
            "🔗 MANTEL TEST RESULTS (HAPLOID vs MORPHOMETRIC DENDROGRAMS - QUERY ONLY)\n")
        self.results_text_corr.insert(tk.END, "=" * 70 + "\n\n")
        self.results_text_corr.insert(tk.END,
            f"📊 Mantel correlation coefficient: {corr:.4f}\n")
        self.results_text_corr.insert(tk.END, f"📈 P-value: {p_value:.4f}\n")
        self.results_text_corr.insert(tk.END, f"📋 Z-score: {z_score:.4f}\n")
        self.results_text_corr.insert(tk.END,
            f"🔢 Number of common Query samples: {len(common)}\n")
        self.results_text_corr.insert(tk.END,
            f"📝 Common Query samples: {', '.join(common)}\n\n")
        if p_value < 0.05:
            self.results_text_corr.insert(tk.END,
                "✅ The correlation is statistically significant (p < 0.05)\n")
        else:
            self.results_text_corr.insert(tk.END,
                "❌ The correlation is not statistically significant (p ≥ 0.05)\n")

        self.progress_bar.stop()
        self.update_status("Mantel test completed")

    def calculate_mantel_test_codominant_vs_dominant_dendrograms(self):
        self.update_status("Running Mantel test: Codominant vs Dominant dendrograms...")
        output_dir = self.create_output_directory("Mantel_correlations")

        if self.genetic_query_distance_matrix is None or \
                self.dominant_query_distance_matrix is None:
            messagebox.showwarning("Warning",
                "Please run both codominant and dominant genetic Query analyses first.")
            self.update_status("Ready")
            return

        self.progress_bar.start(10)
        labels_g = [str(l) for l in self.genetic_query_labels]
        labels_d = [str(l) for l in self.dominant_query_labels]

        c1, c2, common = self._prepare_common_condensed(
            self.genetic_query_distance_matrix, labels_g,
            self.dominant_query_distance_matrix, labels_d)

        if common is None:
            messagebox.showerror("Error", "Not enough common Query samples.")
            self.progress_bar.stop()
            self.update_status("Not enough common Query samples")
            return

        corr, p_value, z_score = self.perform_mantel_test(c1, c2)
        self.save_mantel_results(c1, c2, corr, p_value, z_score, common,
                                 "Codominant_Markers/Dominant_Markers_Dendrograms_Query",
                                 output_dir)

        self.results_text_corr.delete(1.0, tk.END)
        self.results_text_corr.insert(tk.END,
            "🔗 MANTEL TEST RESULTS (CODOMINANT vs DOMINANT DENDROGRAMS - QUERY ONLY)\n")
        self.results_text_corr.insert(tk.END, "=" * 70 + "\n\n")
        self.results_text_corr.insert(tk.END,
            f"📊 Mantel correlation coefficient: {corr:.4f}\n")
        self.results_text_corr.insert(tk.END, f"📈 P-value: {p_value:.4f}\n")
        self.results_text_corr.insert(tk.END, f"📋 Z-score: {z_score:.4f}\n")
        self.results_text_corr.insert(tk.END,
            f"🔢 Number of common Query samples: {len(common)}\n")
        self.results_text_corr.insert(tk.END,
            f"📝 Common Query samples: {', '.join(common)}\n\n")
        if p_value < 0.05:
            self.results_text_corr.insert(tk.END,
                "✅ The correlation is statistically significant (p < 0.05)\n")
        else:
            self.results_text_corr.insert(tk.END,
                "❌ The correlation is not statistically significant (p ≥ 0.05)\n")

        self.progress_bar.stop()
        self.update_status("Mantel test completed")

    def calculate_mantel_test_codominant_vs_haploid_dendrograms(self):
        self.update_status("Running Mantel test: Codominant vs Haploid dendrograms...")
        output_dir = self.create_output_directory("Mantel_correlations")

        if self.genetic_query_distance_matrix is None or \
                self.haploid_query_distance_matrix is None:
            messagebox.showwarning("Warning",
                "Please run both codominant and haploid genetic Query analyses first.")
            self.update_status("Ready")
            return

        self.progress_bar.start(10)
        labels_c = [str(l) for l in self.genetic_query_labels]
        labels_h = [str(l) for l in self.haploid_query_labels]

        c1, c2, common = self._prepare_common_condensed(
            self.genetic_query_distance_matrix, labels_c,
            self.haploid_query_distance_matrix, labels_h)

        if common is None:
            messagebox.showerror("Error", "Not enough common Query samples.")
            self.progress_bar.stop()
            self.update_status("Not enough common Query samples")
            return

        corr, p_value, z_score = self.perform_mantel_test(c1, c2)
        self.save_mantel_results(c1, c2, corr, p_value, z_score, common,
                                 "Codominant_Markers/Haploid_Markers_Dendrograms_Query",
                                 output_dir)

        self.results_text_corr.delete(1.0, tk.END)
        self.results_text_corr.insert(tk.END,
            "🔗 MANTEL TEST RESULTS (CODOMINANT vs HAPLOID DENDROGRAMS - QUERY ONLY)\n")
        self.results_text_corr.insert(tk.END, "=" * 70 + "\n\n")
        self.results_text_corr.insert(tk.END,
            f"📊 Mantel correlation coefficient: {corr:.4f}\n")
        self.results_text_corr.insert(tk.END, f"📈 P-value: {p_value:.4f}\n")
        self.results_text_corr.insert(tk.END, f"📋 Z-score: {z_score:.4f}\n")
        self.results_text_corr.insert(tk.END,
            f"🔢 Number of common Query samples: {len(common)}\n")
        self.results_text_corr.insert(tk.END,
            f"📝 Common Query samples: {', '.join(common)}\n\n")
        if p_value < 0.05:
            self.results_text_corr.insert(tk.END,
                "✅ The correlation is statistically significant (p < 0.05)\n")
        else:
            self.results_text_corr.insert(tk.END,
                "❌ The correlation is not statistically significant (p ≥ 0.05)\n")

        self.progress_bar.stop()
        self.update_status("Mantel test completed")

    def calculate_mantel_test_dominant_vs_haploid_dendrograms(self):
        self.update_status("Running Mantel test: Dominant vs Haploid dendrograms...")
        output_dir = self.create_output_directory("Mantel_correlations")

        if self.dominant_query_distance_matrix is None or \
                self.haploid_query_distance_matrix is None:
            messagebox.showwarning("Warning",
                "Please run both dominant and haploid genetic Query analyses first.")
            self.update_status("Ready")
            return

        self.progress_bar.start(10)
        labels_d = [str(l) for l in self.dominant_query_labels]
        labels_h = [str(l) for l in self.haploid_query_labels]

        c1, c2, common = self._prepare_common_condensed(
            self.dominant_query_distance_matrix, labels_d,
            self.haploid_query_distance_matrix, labels_h)

        if common is None:
            messagebox.showerror("Error", "Not enough common Query samples.")
            self.progress_bar.stop()
            self.update_status("Not enough common Query samples")
            return

        corr, p_value, z_score = self.perform_mantel_test(c1, c2)
        self.save_mantel_results(c1, c2, corr, p_value, z_score, common,
                                 "Dominant_Markers/Haploid_Markers_Dendrograms_Query",
                                 output_dir)

        self.results_text_corr.delete(1.0, tk.END)
        self.results_text_corr.insert(tk.END,
            "🔗 MANTEL TEST RESULTS (DOMINANT vs HAPLOID DENDROGRAMS - QUERY ONLY)\n")
        self.results_text_corr.insert(tk.END, "=" * 70 + "\n\n")
        self.results_text_corr.insert(tk.END,
            f"📊 Mantel correlation coefficient: {corr:.4f}\n")
        self.results_text_corr.insert(tk.END, f"📈 P-value: {p_value:.4f}\n")
        self.results_text_corr.insert(tk.END, f"📋 Z-score: {z_score:.4f}\n")
        self.results_text_corr.insert(tk.END,
            f"🔢 Number of common Query samples: {len(common)}\n")
        self.results_text_corr.insert(tk.END,
            f"📝 Common Query samples: {', '.join(common)}\n\n")
        if p_value < 0.05:
            self.results_text_corr.insert(tk.END,
                "✅ The correlation is statistically significant (p < 0.05)\n")
        else:
            self.results_text_corr.insert(tk.END,
                "❌ The correlation is not statistically significant (p ≥ 0.05)\n")

        self.progress_bar.stop()
        self.update_status("Mantel test completed")

    def calculate_mantel_test_codominant_vs_morpho_pcoa(self):
        self.update_status("Running Mantel test: Codominant PCoA vs Morphometric PCA...")
        output_dir = self.create_output_directory("Mantel_correlations")

        if self.genetic_query_pcoa_coords is None or self.morpho_query_pca_coords is None:
            messagebox.showwarning("Warning",
                "Please run both codominant PCoA and morphometric PCA Query analyses first.")
            self.update_status("Ready")
            return

        self.progress_bar.start(10)
        labels_g = [str(l) for l in self.genetic_query_labels]
        labels_m = [str(l) for l in self.morpho_query_labels]

        common = list(set(labels_g) & set(labels_m))
        if len(common) < 3:
            messagebox.showerror("Error", "Not enough common Query samples.")
            self.progress_bar.stop()
            self.update_status("Not enough common Query samples")
            return

        idx_g = [i for i, l in enumerate(labels_g) if l in common]
        idx_m = [i for i, l in enumerate(labels_m) if l in common]
        common_ordered = [labels_g[i] for i in idx_g]

        dist_g = pdist(self.genetic_query_pcoa_coords[idx_g], metric='euclidean')
        dist_m = pdist(self.morpho_query_pca_coords[idx_m], metric='euclidean')

        corr, p_value, z_score = self.perform_mantel_test(dist_g, dist_m)
        self.save_mantel_pca_results(dist_g, dist_m, corr, p_value, z_score,
                                     common_ordered,
                                     "Codominant_PCoA/Morphometric_PCA_Query", output_dir)

        self.results_text_corr.delete(1.0, tk.END)
        self.results_text_corr.insert(tk.END,
            "🔗 MANTEL TEST RESULTS (CODOMINANT PCoA vs MORPHOMETRIC PCA - QUERY ONLY)\n")
        self.results_text_corr.insert(tk.END, "=" * 70 + "\n\n")
        self.results_text_corr.insert(tk.END,
            f"📊 Mantel correlation coefficient: {corr:.4f}\n")
        self.results_text_corr.insert(tk.END, f"📈 P-value: {p_value:.4f}\n")
        self.results_text_corr.insert(tk.END, f"📋 Z-score: {z_score:.4f}\n")
        self.results_text_corr.insert(tk.END,
            f"🔢 Number of common Query samples: {len(common_ordered)}\n")
        self.results_text_corr.insert(tk.END,
            f"📝 Common Query samples: {', '.join(common_ordered)}\n\n")
        if p_value < 0.05:
            self.results_text_corr.insert(tk.END,
                "✅ The correlation is statistically significant (p < 0.05)\n")
        else:
            self.results_text_corr.insert(tk.END,
                "❌ The correlation is not statistically significant (p ≥ 0.05)\n")

        self.progress_bar.stop()
        self.update_status("Mantel test completed")

    def calculate_mantel_test_dominant_vs_morpho_pcoa(self):
        self.update_status("Running Mantel test: Dominant PCoA vs Morphometric PCA...")
        output_dir = self.create_output_directory("Mantel_correlations")

        if self.dominant_query_pcoa_coords is None or self.morpho_query_pca_coords is None:
            messagebox.showwarning("Warning",
                "Please run both dominant PCoA and morphometric PCA Query analyses first.")
            self.update_status("Ready")
            return

        self.progress_bar.start(10)
        labels_d = [str(l) for l in self.dominant_query_labels]
        labels_m = [str(l) for l in self.morpho_query_labels]

        common = list(set(labels_d) & set(labels_m))
        if len(common) < 3:
            messagebox.showerror("Error", "Not enough common Query samples.")
            self.progress_bar.stop()
            self.update_status("Not enough common Query samples")
            return

        idx_d = [i for i, l in enumerate(labels_d) if l in common]
        idx_m = [i for i, l in enumerate(labels_m) if l in common]
        common_ordered = [labels_d[i] for i in idx_d]

        dist_d = pdist(self.dominant_query_pcoa_coords[idx_d], metric='euclidean')
        dist_m = pdist(self.morpho_query_pca_coords[idx_m], metric='euclidean')

        corr, p_value, z_score = self.perform_mantel_test(dist_d, dist_m)
        self.save_mantel_pca_results(dist_d, dist_m, corr, p_value, z_score,
                                     common_ordered,
                                     "Dominant_PCoA/Morphometric_PCA_Query", output_dir)

        self.results_text_corr.delete(1.0, tk.END)
        self.results_text_corr.insert(tk.END,
            "🔗 MANTEL TEST RESULTS (DOMINANT PCoA vs MORPHOMETRIC PCA - QUERY ONLY)\n")
        self.results_text_corr.insert(tk.END, "=" * 70 + "\n\n")
        self.results_text_corr.insert(tk.END,
            f"📊 Mantel correlation coefficient: {corr:.4f}\n")
        self.results_text_corr.insert(tk.END, f"📈 P-value: {p_value:.4f}\n")
        self.results_text_corr.insert(tk.END, f"📋 Z-score: {z_score:.4f}\n")
        self.results_text_corr.insert(tk.END,
            f"🔢 Number of common Query samples: {len(common_ordered)}\n")
        self.results_text_corr.insert(tk.END,
            f"📝 Common Query samples: {', '.join(common_ordered)}\n\n")
        if p_value < 0.05:
            self.results_text_corr.insert(tk.END,
                "✅ The correlation is statistically significant (p < 0.05)\n")
        else:
            self.results_text_corr.insert(tk.END,
                "❌ The correlation is not statistically significant (p ≥ 0.05)\n")

        self.progress_bar.stop()
        self.update_status("Mantel test completed")

    def calculate_mantel_test_haploid_vs_morpho_pcoa(self):
        self.update_status("Running Mantel test: Haploid PCoA vs Morphometric PCA...")
        output_dir = self.create_output_directory("Mantel_correlations")

        if self.haploid_query_pcoa_coords is None or \
                self.morpho_query_pca_coords is None:
            messagebox.showwarning("Warning",
                "Please run both haploid PCoA and morphometric PCA Query analyses first.")
            self.update_status("Ready")
            return

        self.progress_bar.start(10)
        labels_h = [str(l) for l in self.haploid_query_labels]
        labels_m = [str(l) for l in self.morpho_query_labels]

        common = list(set(labels_h) & set(labels_m))
        if len(common) < 3:
            messagebox.showerror("Error", "Not enough common Query samples.")
            self.progress_bar.stop()
            self.update_status("Not enough common Query samples")
            return

        idx_h = [i for i, l in enumerate(labels_h) if l in common]
        idx_m = [i for i, l in enumerate(labels_m) if l in common]
        common_ordered = [labels_h[i] for i in idx_h]

        coords_h = self.haploid_query_pcoa_coords[idx_h]
        coords_m = self.morpho_query_pca_coords[idx_m]

        dist_h = pdist(coords_h, metric='euclidean')
        dist_m = pdist(coords_m, metric='euclidean')

        corr, p_value, z_score = self.perform_mantel_test(dist_h, dist_m)
        self.save_mantel_pca_results(dist_h, dist_m, corr, p_value, z_score,
                                     common_ordered,
                                     "Haploid_PCoA/Morphometric_PCA_Query", output_dir)

        self.results_text_corr.delete(1.0, tk.END)
        self.results_text_corr.insert(tk.END,
            "🔗 MANTEL TEST RESULTS (HAPLOID PCoA vs MORPHOMETRIC PCA - QUERY ONLY)\n")
        self.results_text_corr.insert(tk.END, "=" * 70 + "\n\n")
        self.results_text_corr.insert(tk.END,
            f"📊 Mantel correlation coefficient: {corr:.4f}\n")
        self.results_text_corr.insert(tk.END, f"📈 P-value: {p_value:.4f}\n")
        self.results_text_corr.insert(tk.END, f"📋 Z-score: {z_score:.4f}\n")
        self.results_text_corr.insert(tk.END,
            f"🔢 Number of common Query samples: {len(common_ordered)}\n")
        self.results_text_corr.insert(tk.END,
            f"📝 Common Query samples: {', '.join(common_ordered)}\n\n")
        if p_value < 0.05:
            self.results_text_corr.insert(tk.END,
                "✅ The correlation is statistically significant (p < 0.05)\n")
        else:
            self.results_text_corr.insert(tk.END,
                "❌ The correlation is not statistically significant (p ≥ 0.05)\n")

        self.progress_bar.stop()
        self.update_status("Mantel test completed")

    def calculate_mantel_test_codominant_vs_dominant_pcoa(self):
        self.update_status("Running Mantel test: Codominant PCoA vs Dominant PCoA...")
        output_dir = self.create_output_directory("Mantel_correlations")

        if self.genetic_query_pcoa_coords is None or \
                self.dominant_query_pcoa_coords is None:
            messagebox.showwarning("Warning",
                "Please run both codominant PCoA and dominant PCoA Query analyses first.")
            self.update_status("Ready")
            return

        self.progress_bar.start(10)
        labels_g = [str(l) for l in self.genetic_query_labels]
        labels_d = [str(l) for l in self.dominant_query_labels]

        common = list(set(labels_g) & set(labels_d))
        if len(common) < 3:
            messagebox.showerror("Error", "Not enough common Query samples.")
            self.progress_bar.stop()
            self.update_status("Not enough common Query samples")
            return

        idx_g = [i for i, l in enumerate(labels_g) if l in common]
        idx_d = [i for i, l in enumerate(labels_d) if l in common]
        common_ordered = [labels_g[i] for i in idx_g]

        dist_g = pdist(self.genetic_query_pcoa_coords[idx_g], metric='euclidean')
        dist_d = pdist(self.dominant_query_pcoa_coords[idx_d], metric='euclidean')

        corr, p_value, z_score = self.perform_mantel_test(dist_g, dist_d)
        self.save_mantel_pca_results(dist_g, dist_d, corr, p_value, z_score,
                                     common_ordered,
                                     "Codominant_PCoA/Dominant_PCoA_Query", output_dir)

        self.results_text_corr.delete(1.0, tk.END)
        self.results_text_corr.insert(tk.END,
            "🔗 MANTEL TEST RESULTS (CODOMINANT PCoA vs DOMINANT PCoA - QUERY ONLY)\n")
        self.results_text_corr.insert(tk.END, "=" * 70 + "\n\n")
        self.results_text_corr.insert(tk.END,
            f"📊 Mantel correlation coefficient: {corr:.4f}\n")
        self.results_text_corr.insert(tk.END, f"📈 P-value: {p_value:.4f}\n")
        self.results_text_corr.insert(tk.END, f"📋 Z-score: {z_score:.4f}\n")
        self.results_text_corr.insert(tk.END,
            f"🔢 Number of common Query samples: {len(common_ordered)}\n")
        self.results_text_corr.insert(tk.END,
            f"📝 Common Query samples: {', '.join(common_ordered)}\n\n")
        if p_value < 0.05:
            self.results_text_corr.insert(tk.END,
                "✅ The correlation is statistically significant (p < 0.05)\n")
        else:
            self.results_text_corr.insert(tk.END,
                "❌ The correlation is not statistically significant (p ≥ 0.05)\n")

        self.progress_bar.stop()
        self.update_status("Mantel test completed")

    def calculate_mantel_test_codominant_vs_haploid_pcoa(self):
        self.update_status("Running Mantel test: Codominant PCoA vs Haploid PCoA...")
        output_dir = self.create_output_directory("Mantel_correlations")

        if self.genetic_query_pcoa_coords is None or \
                self.haploid_query_pcoa_coords is None:
            messagebox.showwarning("Warning",
                "Please run both codominant PCoA and haploid PCoA Query analyses first.")
            self.update_status("Ready")
            return

        self.progress_bar.start(10)
        labels_c = [str(l) for l in self.genetic_query_labels]
        labels_h = [str(l) for l in self.haploid_query_labels]

        common = list(set(labels_c) & set(labels_h))
        if len(common) < 3:
            messagebox.showerror("Error", "Not enough common Query samples.")
            self.progress_bar.stop()
            self.update_status("Not enough common Query samples")
            return

        idx_c = [i for i, l in enumerate(labels_c) if l in common]
        idx_h = [i for i, l in enumerate(labels_h) if l in common]
        common_ordered = [labels_c[i] for i in idx_c]

        dist_c = pdist(self.genetic_query_pcoa_coords[idx_c], metric='euclidean')
        dist_h = pdist(self.haploid_query_pcoa_coords[idx_h], metric='euclidean')

        corr, p_value, z_score = self.perform_mantel_test(dist_c, dist_h)
        self.save_mantel_pca_results(dist_c, dist_h, corr, p_value, z_score,
                                     common_ordered,
                                     "Codominant_PCoA/Haploid_PCoA_Query", output_dir)

        self.results_text_corr.delete(1.0, tk.END)
        self.results_text_corr.insert(tk.END,
            "🔗 MANTEL TEST RESULTS (CODOMINANT PCoA vs HAPLOID PCoA - QUERY ONLY)\n")
        self.results_text_corr.insert(tk.END, "=" * 70 + "\n\n")
        self.results_text_corr.insert(tk.END,
            f"📊 Mantel correlation coefficient: {corr:.4f}\n")
        self.results_text_corr.insert(tk.END, f"📈 P-value: {p_value:.4f}\n")
        self.results_text_corr.insert(tk.END, f"📋 Z-score: {z_score:.4f}\n")
        self.results_text_corr.insert(tk.END,
            f"🔢 Number of common Query samples: {len(common_ordered)}\n")
        self.results_text_corr.insert(tk.END,
            f"📝 Common Query samples: {', '.join(common_ordered)}\n\n")
        if p_value < 0.05:
            self.results_text_corr.insert(tk.END,
                "✅ The correlation is statistically significant (p < 0.05)\n")
        else:
            self.results_text_corr.insert(tk.END,
                "❌ The correlation is not statistically significant (p ≥ 0.05)\n")

        self.progress_bar.stop()
        self.update_status("Mantel test completed")

    def calculate_mantel_test_dominant_vs_haploid_pcoa(self):
        self.update_status("Running Mantel test: Dominant PCoA vs Haploid PCoA...")
        output_dir = self.create_output_directory("Mantel_correlations")

        if self.dominant_query_pcoa_coords is None or \
                self.haploid_query_pcoa_coords is None:
            messagebox.showwarning("Warning",
                "Please run both dominant PCoA and haploid PCoA Query analyses first.")
            self.update_status("Ready")
            return

        self.progress_bar.start(10)
        labels_d = [str(l) for l in self.dominant_query_labels]
        labels_h = [str(l) for l in self.haploid_query_labels]

        common = list(set(labels_d) & set(labels_h))
        if len(common) < 3:
            messagebox.showerror("Error", "Not enough common Query samples.")
            self.progress_bar.stop()
            self.update_status("Not enough common Query samples")
            return

        idx_d = [i for i, l in enumerate(labels_d) if l in common]
        idx_h = [i for i, l in enumerate(labels_h) if l in common]
        common_ordered = [labels_d[i] for i in idx_d]

        dist_d = pdist(self.dominant_query_pcoa_coords[idx_d], metric='euclidean')
        dist_h = pdist(self.haploid_query_pcoa_coords[idx_h], metric='euclidean')

        corr, p_value, z_score = self.perform_mantel_test(dist_d, dist_h)
        self.save_mantel_pca_results(dist_d, dist_h, corr, p_value, z_score,
                                     common_ordered,
                                     "Dominant_PCoA/Haploid_PCoA_Query", output_dir)

        self.results_text_corr.delete(1.0, tk.END)
        self.results_text_corr.insert(tk.END,
            "🔗 MANTEL TEST RESULTS (DOMINANT PCoA vs HAPLOID PCoA - QUERY ONLY)\n")
        self.results_text_corr.insert(tk.END, "=" * 70 + "\n\n")
        self.results_text_corr.insert(tk.END,
            f"📊 Mantel correlation coefficient: {corr:.4f}\n")
        self.results_text_corr.insert(tk.END, f"📈 P-value: {p_value:.4f}\n")
        self.results_text_corr.insert(tk.END, f"📋 Z-score: {z_score:.4f}\n")
        self.results_text_corr.insert(tk.END,
            f"🔢 Number of common Query samples: {len(common_ordered)}\n")
        self.results_text_corr.insert(tk.END,
            f"📝 Common Query samples: {', '.join(common_ordered)}\n\n")
        if p_value < 0.05:
            self.results_text_corr.insert(tk.END,
                "✅ The correlation is statistically significant (p < 0.05)\n")
        else:
            self.results_text_corr.insert(tk.END,
                "❌ The correlation is not statistically significant (p ≥ 0.05)\n")

        self.progress_bar.stop()
        self.update_status("Mantel test completed")

    # ==================== OUTLIER DETECTION TAB ====================

    def setup_outlier_tab(self):
        main_container = ModernFrame(self.outlier_frame)
        main_container.pack(fill="both", expand=True, padx=10, pady=10)

        info_card = CardFrame(main_container)
        info_card.pack(fill="x", pady=(0, 15))

        info_title = ModernLabel(info_card, text="PCA-Based Outlier Detection", font=("Segoe UI", 12, "bold"))
        info_title.pack(anchor="w", pady=(0, 10))

        info_text = ModernLabel(info_card, 
            text="Detects outliers using PCA (Principal Component Analysis) combined with Elliptic Envelope.\n"
                 "Outliers are samples that deviate significantly from the main structure of the data.\n\n"
                 "💡 This analysis does NOT require predefined populations!",
            font=("Segoe UI", 9),
            fg=MODERN_THEME["text_light"],
            wraplength=800)
        info_text.pack(anchor="w", pady=(0, 5))

        data_card = CardFrame(main_container)
        data_card.pack(fill="x", pady=(0, 15))

        data_title = ModernLabel(data_card, text="Select Data for Outlier Detection", font=("Segoe UI", 11, "bold"))
        data_title.pack(anchor="w", pady=(0, 10))

        type_frame = ModernFrame(data_card)
        type_frame.pack(fill="x", pady=5)

        self.outlier_data_type = tk.StringVar(value="codominant")
        
        radio_cod = tk.Radiobutton(type_frame, text="Codominant Markers", variable=self.outlier_data_type,
                                   value="codominant", bg=MODERN_THEME["card_bg"], font=("Segoe UI", 9))
        radio_cod.pack(side="left", padx=10)
        
        radio_dom = tk.Radiobutton(type_frame, text="Dominant Markers", variable=self.outlier_data_type,
                                   value="dominant", bg=MODERN_THEME["card_bg"], font=("Segoe UI", 9))
        radio_dom.pack(side="left", padx=10)
        
        radio_hap = tk.Radiobutton(type_frame, text="Haploid Markers", variable=self.outlier_data_type,
                                   value="haploid", bg=MODERN_THEME["card_bg"], font=("Segoe UI", 9))
        radio_hap.pack(side="left", padx=10)
        
        radio_morph = tk.Radiobutton(type_frame, text="Morphometric Data", variable=self.outlier_data_type,
                                     value="morphometric", bg=MODERN_THEME["card_bg"], font=("Segoe UI", 9))
        radio_morph.pack(side="left", padx=10)

        file_frame = ModernFrame(data_card)
        file_frame.pack(fill="x", pady=10)

        label_file = ModernLabel(file_frame, text="Data File:", font=("Segoe UI", 10, "bold"))
        label_file.grid(row=0, column=0, sticky="w", padx=(0, 10))

        self.entry_outlier_file = ModernEntry(file_frame, width=50)
        self.entry_outlier_file.grid(row=0, column=1, padx=5, sticky="ew")

        btn_browse = ModernButton(file_frame, text="Browse",
                                  command=lambda: self.load_file(self.entry_outlier_file),
                                  bg=MODERN_THEME["accent"], fg="white")
        btn_browse.grid(row=0, column=2, padx=5)

        file_frame.columnconfigure(1, weight=1)

        params_frame = ModernFrame(data_card)
        params_frame.pack(fill="x", pady=10)

        label_contamination = ModernLabel(params_frame, text="Contamination Rate (0.01-0.2):", font=("Segoe UI", 9))
        label_contamination.grid(row=0, column=0, sticky="w", padx=(0, 10))
        
        self.contamination_var = tk.DoubleVar(value=0.05)
        self.entry_contamination = ModernEntry(params_frame, textvariable=self.contamination_var, width=10)
        self.entry_contamination.grid(row=0, column=1, sticky="w", padx=5)
        
        label_threshold = ModernLabel(params_frame, text="Mahalanobis Threshold:", font=("Segoe UI", 9))
        label_threshold.grid(row=1, column=0, sticky="w", padx=(0, 10), pady=(5, 0))
        
        self.threshold_var = tk.DoubleVar(value=3.0)
        self.entry_threshold = ModernEntry(params_frame, textvariable=self.threshold_var, width=10)
        self.entry_threshold.grid(row=1, column=1, sticky="w", padx=5, pady=(5, 0))

        action_card = CardFrame(main_container)
        action_card.pack(fill="x", pady=(0, 15))

        btn_analyze = ModernButton(action_card,
                                   text="🔍 Detect Outliers",
                                   command=self.detect_outliers_pca,
                                   bg=MODERN_THEME["success"], fg="white",
                                   font=("Segoe UI", 11, "bold"))
        btn_analyze.pack(pady=10)

        self.outlier_progress_bar = ttk.Progressbar(action_card, mode='indeterminate', length=100)
        self.outlier_progress_bar.pack(fill="x", pady=(5, 0))

        results_card = CardFrame(main_container)
        results_card.pack(fill="both", expand=True)

        results_title = ModernLabel(results_card, text="Outlier Detection Results", font=("Segoe UI", 12, "bold"))
        results_title.pack(anchor="w", pady=(0, 10))

        self.results_text_outlier = scrolledtext.ScrolledText(results_card,
                                                               bg="white",
                                                               fg=MODERN_THEME["text"],
                                                               font=("Consolas", 9),
                                                               relief="flat",
                                                               padx=10,
                                                               pady=10)
        self.results_text_outlier.pack(fill="both", expand=True)

    def detect_outliers_pca(self):
        self.update_status("Detecting outliers using PCA...")
        self.outlier_progress_bar.start(10)
        self.root.update_idletasks()

        filename = self.entry_outlier_file.get().strip()
        if not filename:
            messagebox.showwarning("Warning", "Please select a data file.")
            self.outlier_progress_bar.stop()
            self.update_status("Ready")
            return

        data_type = self.outlier_data_type.get()
        contamination = self.contamination_var.get()
        mahalanobis_threshold = self.threshold_var.get()

        if contamination < 0.01 or contamination > 0.2:
            messagebox.showwarning("Warning", "Contamination rate should be between 0.01 and 0.2")
            contamination = max(0.01, min(0.2, contamination))

        try:
            if data_type == 'morphometric':
                df, sample_names, data_matrix = self._prepare_morpho_data_for_outlier(filename)
            else:
                df, sample_names, data_matrix = self._prepare_genetic_data_for_outlier(filename, data_type)

            if data_matrix is None or len(sample_names) < 3:
                messagebox.showerror("Error", "Not enough valid samples for outlier detection (minimum 3 required)")
                self.outlier_progress_bar.stop()
                self.update_status("Error")
                return

            self.results_text_outlier.delete(1.0, tk.END)
            self.results_text_outlier.insert(tk.END, "=" * 70 + "\n")
            self.results_text_outlier.insert(tk.END, f"📈 PCA-BASED OUTLIER DETECTION - {data_type.upper()}\n")
            self.results_text_outlier.insert(tk.END, "=" * 70 + "\n\n")
            self.results_text_outlier.insert(tk.END, f"📁 File: {os.path.basename(filename)}\n")
            self.results_text_outlier.insert(tk.END, f"🔢 Samples: {len(sample_names)}\n")
            self.results_text_outlier.insert(tk.END, f"📊 Contamination rate: {contamination}\n")
            self.results_text_outlier.insert(tk.END, f"📏 Mahalanobis threshold: {mahalanobis_threshold}\n\n")

            self.outlier_progress_bar['value'] = 30
            self.root.update_idletasks()

            scaler = StandardScaler()
            data_scaled = scaler.fit_transform(data_matrix)

            n_components = min(5, data_scaled.shape[1], data_scaled.shape[0] - 1)
            pca = PCA(n_components=n_components)
            pca_result = pca.fit_transform(data_scaled)
            
            explained_variance = pca.explained_variance_ratio_
            cumulative_variance = np.cumsum(explained_variance)

            self.results_text_outlier.insert(tk.END, "📊 PCA Variance Explained:\n")
            for i, var in enumerate(explained_variance[:5]):
                self.results_text_outlier.insert(tk.END, f"   PC{i+1}: {var*100:.2f}%\n")
            self.results_text_outlier.insert(tk.END, f"   Cumulative (PC1-PC{min(5, n_components)}): {cumulative_variance[min(4, n_components-1)]*100:.2f}%\n\n")

            self.outlier_progress_bar['value'] = 50
            self.root.update_idletasks()

            from sklearn.covariance import EllipticEnvelope
            detector = EllipticEnvelope(contamination=contamination, random_state=42)
            outliers_ee = detector.fit_predict(pca_result)
            
            outlier_indices_ee = [i for i, val in enumerate(outliers_ee) if val == -1]
            inlier_indices_ee = [i for i, val in enumerate(outliers_ee) if val == 1]

            self.outlier_progress_bar['value'] = 70
            self.root.update_idletasks()

            from scipy.spatial.distance import mahalanobis
            from scipy.stats import chi2
            
            mean_pca = np.mean(pca_result, axis=0)
            cov_pca = np.cov(pca_result.T)
            inv_cov_pca = np.linalg.pinv(cov_pca)
            
            mahalanobis_distances = []
            for i in range(len(pca_result)):
                dist = mahalanobis(pca_result[i], mean_pca, inv_cov_pca)
                mahalanobis_distances.append(dist)
            
            chi2_threshold = np.sqrt(chi2.ppf(0.95, n_components))
            outlier_indices_maha = [i for i, d in enumerate(mahalanobis_distances) if d > mahalanobis_threshold]

            consensus_outliers = list(set(outlier_indices_ee) & set(outlier_indices_maha))
            
            self.outlier_progress_bar['value'] = 85
            self.root.update_idletasks()

            self.results_text_outlier.insert(tk.END, "🔍 OUTLIER DETECTION RESULTS\n")
            self.results_text_outlier.insert(tk.END, "-" * 50 + "\n\n")

            self.results_text_outlier.insert(tk.END, f"📌 Method 1 - Elliptic Envelope (contamination={contamination}):\n")
            if outlier_indices_ee:
                self.results_text_outlier.insert(tk.END, f"   ⚠️ {len(outlier_indices_ee)} potential outliers detected:\n")
                for idx in outlier_indices_ee:
                    self.results_text_outlier.insert(tk.END, f"      • {sample_names[idx]}\n")
            else:
                self.results_text_outlier.insert(tk.END, "   ✅ No outliers detected.\n")
            self.results_text_outlier.insert(tk.END, "\n")

            self.results_text_outlier.insert(tk.END, f"📌 Method 2 - Mahalanobis Distance (threshold={mahalanobis_threshold}):\n")
            if outlier_indices_maha:
                self.results_text_outlier.insert(tk.END, f"   ⚠️ {len(outlier_indices_maha)} potential outliers detected:\n")
                for idx in outlier_indices_maha:
                    self.results_text_outlier.insert(tk.END, f"      • {sample_names[idx]} (distance={mahalanobis_distances[idx]:.2f})\n")
            else:
                self.results_text_outlier.insert(tk.END, "   ✅ No outliers detected.\n")
            self.results_text_outlier.insert(tk.END, "\n")

            self.results_text_outlier.insert(tk.END, f"📌 Consensus (detected by both methods):\n")
            if consensus_outliers:
                self.results_text_outlier.insert(tk.END, f"   ⚠️ {len(consensus_outliers)} high-confidence outliers:\n")
                for idx in consensus_outliers:
                    self.results_text_outlier.insert(tk.END, f"      • {sample_names[idx]}\n")
            else:
                self.results_text_outlier.insert(tk.END, "   ✅ No high-confidence outliers detected.\n")
            self.results_text_outlier.insert(tk.END, "\n")

            self.results_text_outlier.insert(tk.END, "📋 SUMMARY TABLE\n")
            self.results_text_outlier.insert(tk.END, "-" * 50 + "\n")
            self.results_text_outlier.insert(tk.END, f"{'Sample':<25} {'Elliptic':<12} {'Mahalanobis':<12} {'Consensus':<10}\n")
            self.results_text_outlier.insert(tk.END, "-" * 60 + "\n")
            
            for i, name in enumerate(sample_names):
                is_ee = "OUTLIER" if i in outlier_indices_ee else "Normal"
                is_maha = "OUTLIER" if i in outlier_indices_maha else "Normal"
                is_consensus = "⚠️ YES" if i in consensus_outliers else "No"
                self.results_text_outlier.insert(tk.END, f"{name:<25} {is_ee:<12} {is_maha:<12} {is_consensus:<10}\n")
            
            self.results_text_outlier.insert(tk.END, "\n")

            self.outlier_progress_bar['value'] = 90
            self.root.update_idletasks()

            output_dir = self.create_output_directory("Outlier_Detection")
            
            fig, axes = plt.subplots(1, 2, figsize=(14, 6))
            
            colors_ee = ['red' if i in outlier_indices_ee else 'steelblue' for i in range(len(sample_names))]
            axes[0].scatter(pca_result[:, 0], pca_result[:, 1], 
                               c=colors_ee, s=100, alpha=0.7, edgecolors='black', linewidth=0.5)
            
            for i, name in enumerate(sample_names):
                axes[0].annotate(name, (pca_result[i, 0], pca_result[i, 1]),
                                xytext=(5, 5), textcoords='offset points', fontsize=8, alpha=0.8)
            
            axes[0].set_xlabel(f'PC1 ({explained_variance[0]*100:.1f}%)', fontsize=12)
            axes[0].set_ylabel(f'PC2 ({explained_variance[1]*100:.1f}%)', fontsize=12)
            axes[0].set_title(f'PCA - Elliptic Envelope Outliers\n(contamination={contamination})', fontsize=12)
            axes[0].grid(True, alpha=0.3)
            
            from matplotlib.patches import Patch
            legend_elements = [Patch(facecolor='steelblue', label='Normal', alpha=0.7),
                              Patch(facecolor='red', label='Outlier', alpha=0.7)]
            axes[0].legend(handles=legend_elements, loc='best')
            
            colors_maha = ['red' if i in outlier_indices_maha else 'steelblue' for i in range(len(sample_names))]
            axes[1].scatter(pca_result[:, 0], pca_result[:, 1], 
                               c=colors_maha, s=100, alpha=0.7, edgecolors='black', linewidth=0.5)
            
            for i, name in enumerate(sample_names):
                axes[1].annotate(name, (pca_result[i, 0], pca_result[i, 1]),
                                xytext=(5, 5), textcoords='offset points', fontsize=8, alpha=0.8)
            
            axes[1].set_xlabel(f'PC1 ({explained_variance[0]*100:.1f}%)', fontsize=12)
            axes[1].set_ylabel(f'PC2 ({explained_variance[1]*100:.1f}%)', fontsize=12)
            axes[1].set_title(f'PCA - Mahalanobis Distance Outliers\n(threshold={mahalanobis_threshold})', fontsize=12)
            axes[1].grid(True, alpha=0.3)
            axes[1].legend(handles=legend_elements, loc='best')
            
            plt.tight_layout()
            plot_path1 = os.path.join(output_dir, f"pca_outliers_{data_type}.png")
            plt.savefig(plot_path1, dpi=300, bbox_inches='tight')
            plt.close()
            
            self.results_text_outlier.insert(tk.END, f"📊 PCA plots saved: {plot_path1}\n\n")

            fig, ax = plt.subplots(figsize=(10, 6))
            
            colors_dist = ['red' if i in outlier_indices_maha else 'steelblue' for i in range(len(sample_names))]
            ax.bar(range(len(sample_names)), mahalanobis_distances, color=colors_dist, alpha=0.7)
            ax.axhline(y=mahalanobis_threshold, color='red', linestyle='--', linewidth=2, label=f'Threshold ({mahalanobis_threshold})')
            
            ax.set_xlabel('Sample Index', fontsize=12)
            ax.set_ylabel('Mahalanobis Distance', fontsize=12)
            ax.set_title('Mahalanobis Distance Distribution', fontsize=14)
            ax.set_xticks(range(len(sample_names)))
            ax.set_xticklabels(sample_names, rotation=45, ha='right', fontsize=8)
            ax.legend()
            ax.grid(True, alpha=0.3)
            
            plt.tight_layout()
            plot_path2 = os.path.join(output_dir, f"mahalanobis_distances_{data_type}.png")
            plt.savefig(plot_path2, dpi=300, bbox_inches='tight')
            plt.close()
            
            self.results_text_outlier.insert(tk.END, f"📊 Mahalanobis distances plot saved: {plot_path2}\n")

            results_df = pd.DataFrame({
                'Sample': sample_names,
                'Elliptic_Envelope_Outlier': ['Yes' if i in outlier_indices_ee else 'No' for i in range(len(sample_names))],
                'Mahalanobis_Distance': [round(d, 3) for d in mahalanobis_distances],
                'Mahalanobis_Outlier': ['Yes' if i in outlier_indices_maha else 'No' for i in range(len(sample_names))],
                'Consensus_Outlier': ['Yes' if i in consensus_outliers else 'No' for i in range(len(sample_names))]
            })
            
            excel_path = os.path.join(output_dir, f"outlier_results_{data_type}.xlsx")
            with pd.ExcelWriter(excel_path, engine='openpyxl') as writer:
                title_df = pd.DataFrame({"": [f"GMDA - PCA Outlier Detection ({data_type.upper()})"]})
                title_df.to_excel(writer, sheet_name='Outlier_Results', index=False, header=False)
                results_df.to_excel(writer, sheet_name='Outlier_Results', startrow=2, index=False)
                
                pca_stats = pd.DataFrame({
                    'PC': [f'PC{i+1}' for i in range(len(explained_variance))],
                    'Explained_Variance_%': [round(v*100, 2) for v in explained_variance],
                    'Cumulative_%': [round(c*100, 2) for c in cumulative_variance]
                })
                pca_stats.to_excel(writer, sheet_name='PCA_Statistics', startrow=0, index=False)
            
            self.results_text_outlier.insert(tk.END, f"\n💾 Results saved to: {excel_path}\n")
            
            self.results_text_outlier.insert(tk.END, "\n✅ Outlier detection completed!\n")
            self.results_text_outlier.insert(tk.END, "=" * 70 + "\n")
            
            self.outlier_progress_bar.stop()
            self.update_status("Outlier detection completed")
            
        except Exception as e:
            messagebox.showerror("Error", f"Outlier detection failed:\n{str(e)}")
            self.results_text_outlier.insert(tk.END, f"\n❌ Error: {str(e)}\n")
            self.update_status("Error in outlier detection")
            self.outlier_progress_bar.stop()

    def _prepare_genetic_data_for_outlier(self, filename, data_type):
        try:
            df, name_col = self.read_genetic_file(filename)
            
            sample_names = []
            for name in df[name_col]:
                if pd.notna(name):
                    sample_names.append(str(name).strip())
            
            data = df.drop(columns=[name_col])
            
            valid_indices = [i for i, name in enumerate(sample_names) if name and name != '']
            sample_names = [sample_names[i] for i in valid_indices]
            data = data.iloc[valid_indices]
            
            if len(sample_names) < 3:
                return df, sample_names, None
            
            if data_type == 'codominant':
                data_matrix = []
                for i in range(len(data)):
                    row = []
                    for col in data.columns:
                        val = data.iloc[i][col]
                        if pd.isna(val):
                            row.append(np.nan)
                        else:
                            try:
                                row.append(float(val))
                            except (ValueError, TypeError):
                                row.append(hash(str(val)) % 1000)
                    data_matrix.append(row)
                data_matrix = np.array(data_matrix, dtype=float)
                
            elif data_type == 'dominant':
                data_matrix = []
                for i in range(len(data)):
                    row = []
                    for col in data.columns:
                        val = data.iloc[i][col]
                        if pd.isna(val):
                            row.append(np.nan)
                        else:
                            val_str = str(val).upper().strip()
                            if val_str in ['1', 'TRUE', 'YES', 'PRESENT', 'P', '+']:
                                row.append(1.0)
                            else:
                                row.append(0.0)
                    data_matrix.append(row)
                data_matrix = np.array(data_matrix, dtype=float)
                
            else:
                data_matrix = []
                for i in range(len(data)):
                    row = []
                    for col in data.columns:
                        val = data.iloc[i][col]
                        if pd.isna(val):
                            row.append(np.nan)
                        else:
                            row.append(hash(str(val)) % 1000)
                    data_matrix.append(row)
                data_matrix = np.array(data_matrix, dtype=float)
            
            from sklearn.impute import SimpleImputer
            imputer = SimpleImputer(strategy='mean')
            data_matrix = imputer.fit_transform(data_matrix)
            
            return df, sample_names, data_matrix
            
        except Exception as e:
            raise Exception(f"Error preparing genetic data: {str(e)}")

    def _prepare_morpho_data_for_outlier(self, filename):
        try:
            df = pd.read_excel(filename, engine='openpyxl', header=0)
            
            if df.empty:
                raise ValueError("File is empty")
            
            sample_col = df.columns[0]
            sample_names = []
            for name in df[sample_col]:
                if pd.notna(name):
                    sample_names.append(str(name).strip())
            
            data = df.drop(columns=[sample_col])
            
            valid_indices = [i for i, name in enumerate(sample_names) if name and name != '']
            sample_names = [sample_names[i] for i in valid_indices]
            data = data.iloc[valid_indices]
            
            if len(sample_names) < 3:
                return df, sample_names, None
            
            data_matrix = []
            for i in range(len(data)):
                row = []
                for col in data.columns:
                    val = data.iloc[i][col]
                    if pd.isna(val):
                        row.append(np.nan)
                    else:
                        try:
                            row.append(float(val))
                        except (ValueError, TypeError):
                            row.append(np.nan)
                data_matrix.append(row)
            
            data_matrix = np.array(data_matrix, dtype=float)
            
            from sklearn.impute import SimpleImputer
            imputer = SimpleImputer(strategy='mean')
            data_matrix = imputer.fit_transform(data_matrix)
            
            return df, sample_names, data_matrix
            
        except Exception as e:
            raise Exception(f"Error preparing morphometric data: {str(e)}")


if __name__ == "__main__":
    root = tk.Tk()
    app = GMDA(root)
    root.mainloop()
