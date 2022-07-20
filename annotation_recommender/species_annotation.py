"""Annotation for Species."""

from annotation_recommender import constants as cn
from annotation_recommender import tools

import editdistance
import libsbml
import numpy as np
import os
import pickle

# below might be in constants or main script
with open(os.path.join(cn.CHEBI_DIR, 'chebi_shortened_formula_30apr2022.pickle'), 'rb') as f:
  ref_shortened_chebi_to_formula = pickle.load(f)
with open(os.path.join(cn.CHEBI_DIR, 'chebi_synonyms.pickle'), 'rb') as f:
  chebi_synonyms = pickle.load(f)
#
chebi_low_synonyms = dict()
for one_k in chebi_synonyms.keys():
  chebi_low_synonyms[one_k] = list(set([val.lower() for val in chebi_synonyms[one_k]]))


class SpeciesAnnotation(object):

  def __init__(self, libsbml_fpath=None):
    # self.exist_annotation stores existing CHEBI annotations in the model
    # If none exists, set None
    if libsbml_fpath is not None:
      reader = libsbml.SBMLReader()
      document = reader.readSBML(libsbml_fpath)
      self.model = document.getModel()
      exist_annotation_raw = {val.getId():tools.getQualifierFromString(val.getAnnotationString(), cn.CHEBI) \
                        for val in self.model.getListOfSpecies()}
      exist_annotation_filt = {val:exist_annotation_raw[val] for val in exist_annotation_raw.keys() \
                               if exist_annotation_raw[val] is not None}
      self.exist_annotation = {k:tools.transformCHEBIToFormula(exist_annotation_filt[k], ref_shortened_chebi_to_formula) \
                               for k in exist_annotation_filt.keys()}
    else:
      self.exist_annotation = None
    # self.pred_annotation stores predicted species annotations
    self.pred_annotation = None
      

  def predictSpeciesCandidatesByName(self, inp_spec_list=None):
    """
    Predict list of species annotations
    using species names/IDs.
    Rule is 1) use species name, 
    2) if not provided, use species ID.
  
    Parameters
    ----------
    inp_spec_list: str-list (or iterable list of strings)
        List of species IDs to extract names

    Returns
    -------
    result:  
        Should be {species_id: 
                      {'chebi':[CHEBI terms]},
                      {'score': match_score},
                      {'formula': [chemical formulas in string]}
                      {'formula2chebi': [CHEHBI terms]}
                  }
        match_score is expected to be between 0.0-1.0
    """
    result = dict()
    if inp_spec_list is None:
      spec_list = [val.getId() for val in self.model.getListOfSpecies()]
    else:
      spec_list = inp_spec_list
    for one_spec_id in spec_list:
      one_result = dict()
      one_spec_name = self.model.getSpecies(one_spec_id).name.lower()
      if len(one_spec_name) == 0:
        one_spec_name = one_spec_id.lower()
      # For now, choose the terms that are included in the CHEBI-formula mapping reference
      dist_dict_min = {one_k:np.min([editdistance.eval(one_spec_name, val) for val in chebi_low_synonyms[one_k]]) \
                       for one_k in chebi_low_synonyms.keys() if one_k in ref_shortened_chebi_to_formula.keys()}
      min_min_dist = np.min([dist_dict_min[val] for val in dist_dict_min.keys()])
      one_match_score = 1 - min_min_dist/len(one_spec_name)
      one_result['match_score'] = one_match_score
      min_min_chebis = [one_k for one_k in dist_dict_min.keys() \
                        if dist_dict_min[one_k]==min_min_dist and one_k in ref_shortened_chebi_to_formula.keys()]
      # predicted formula of the species
      one_result['chebi'] = min_min_chebis


      # formul2chebi should be provided as an independent reference dictionary
      formula2chebi = dict()
      for one_chebi in min_min_chebis:
        one_itm = ref_shortened_chebi_to_formula[one_chebi]
        if one_itm in formula2chebi.keys():
          formula2chebi[one_itm].append(one_chebi)
        else:
          formula2chebi[one_itm] = [one_chebi]
      one_result['formula2chebi'] = formula2chebi
    
      min_min_formula = list(formula2chebi.keys())
      one_result['formula'] = min_min_formula
      result[one_spec_id] = one_result
    return result

