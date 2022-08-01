# test_species_annotation.py
# Testing SpeciesAnnotation class


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


#############################
# Tests
#############################
class TestSpeciesAnnotation(unittest.TestCase):

  def setUp(self):
    self.spec_cl = sa.SpeciesAnnotation(libsbml_path = E_COLI_PATH)
