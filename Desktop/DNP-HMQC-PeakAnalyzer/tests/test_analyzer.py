import os
import pytest
import pandas as pd
from main import get_project_paths

def test_project_paths():
    """בדיקה שהקוד מזהה את התיקיות בצורה דינמית ומייצר את תיקיית outputs"""
    data_dir, output_dir = get_project_paths()
    
    assert os.path.exists(data_dir) is True
    assert os.path.exists(output_dir) is True
    assert "data" in data_dir
    assert "outputs" in output_dir

def test_csv_data_structure():
    """בדיקה שקובצי הנתונים קיימים ומכילים את העמודות הנכונות לפרויקט"""
    data_dir, _ = get_project_paths()
    path_p1 = os.path.join(data_dir, "protein_1.csv")
    path_p2 = os.path.join(data_dir, "protein_2.csv")
    
    expected_columns = ['residue', 'h', 'c', 'enhancement']
    
    for path in [path_p1, path_p2]:
        assert os.path.exists(path) is True, f"Missing data file: {path}"
        df = pd.read_csv(path)
        # בדיקה שכל העמודות הנדרשות קיימות בקובץ
        for col in expected_columns:
            assert col in df.columns, f"Column '{col}' is missing in {path}"