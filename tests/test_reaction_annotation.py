# test_reaction_annotation.py
# Test for ReactionAnnotation class

import libsbml
import os
import sys
import unittest

from annotation_recommender import species_annotation as sa
from annotation_recommender import reaction_annotation as ra
from annotation_recommender import constants as cn
from annotation_recommender import tools


E_COLI_PATH = os.path.join(os.getcwd(), 'e_coli_core.xml')
BIOMD_248_PATH = os.path.join(os.getcwd(), 'BIOMD0000000248.xml')
# Below ??
# M_FDP_C = 'M_fdp_c'
# M_ATP_C = 'M_atp_c'
# ONESET_SPECIES_IDS = [M_FDP_C, M_ATP_C]

# ID of a reaction
R_PFK = 'R_PFK'

#############################
# Tests
#############################
class TestReactionAnnotation(unittest.TestCase):

  def setUp(self):
    self.spec_cl = sa.SpeciesAnnotation(libsbml_fpath = E_COLI_PATH)
    self.reac_cl = ra.ReactionAnnotation(libsbml_fpath = E_COLI_PATH)

  def testGetMatchScore(self):
    one_dict = {'R1': {'M1': 0.7}}
    one_match_score = self.reac_cl.getMatchScore(score_dict=one_dict)
    self.assertEqual(one_match_score, 0.7)