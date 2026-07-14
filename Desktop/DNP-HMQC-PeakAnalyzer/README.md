# DNP-HMQC-PeakAnalyzer 🔬

An automated Python tool designed to parse, filter, analyze, and visualize 2D NMR ($^1H$ - $^{13}C$ HMQC) peak lists enhanced by **Dynamic Nuclear Polarization (DNP)**. 

This repository was built as a clean, modular, and fully tested final project, addressing all previous feedback regarding hardcoded paths, lack of documentation, and minimal analysis.

---

## 🌟 Key Features
* **Zero Hardcoded Paths:** Completely dynamic file path resolution using Python's `os` library—works out-of-the-box on any machine.
* **Expanded Comparative Analysis:** Automatically processes and compares multiple dataset files (e.g., `protein_1` and `protein_2`), computing comparative statistical metrics (mean, median, standard deviation) on DNP enhancement factors.
* **Exceptional Signal Filtering:** Identifies and extracts specific residue peaks exhibiting high DNP enhancements (Enhancement > 2.5).
* **Automated Scientific Visualizations:** Generates and exports high-resolution figures to the `/outputs` folder:
  1. A 2D NMR Chemical Shift Correlation Map scaled by enhancement values.
  2. A DNP Enhancement distribution comparison (Boxplot + Strip-plot).
* **Automated Unit Testing:** Fully-equipped suite of tests to verify directory structure and data integrity.

---

## 📁 Repository Structure
```text
DNP-HMQC-PeakAnalyzer/
│
├── data/                    # Raw input NMR/DNP peak lists (CSVs)
│   ├── protein_1.csv
│   └── protein_2.csv
│
├── outputs/                 # Auto-generated analytical plots & statistics
│
├── tests/                   # Automated unit tests
│   └── test_analyzer.py
│
├── requirements.txt         # Package dependencies
├── README.md                # Project documentation
└── main.py                  # Main execution script