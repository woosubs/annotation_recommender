# constants.py
"""
Constants for modlues
"""

import os

BASE_DIR = '/Users/woosubs/Desktop/AutomateAnnotation/'
DATA_DIR = os.path.join(BASE_DIR, "DATA")
ALGO_DIR = os.path.join(DATA_DIR, "algo")
CHEBI_DIR = os.path.join(DATA_DIR, "chebi")
RHEA_DIR = os.path.join(DATA_DIR, "rhea")

CHEBI = "chebi"
RHEA = "rhea"
KEGG_REACTION = "kegg.reaction"