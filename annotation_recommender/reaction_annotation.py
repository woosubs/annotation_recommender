"""Annotation for Reactions."""

from annotation_recommender import constants as cn
from annotation_recommender import tools

import libsbml
import numpy as np
import os
import pickle
import pandas as pd


with open(os.path.join(cn.RHEA_DIR, 'rhea_all2bi.pkl'), 'rb') as f:
  ref_rhea2bi = pickle.load(f)
with open(os.path.join(cn.RHEA_DIR, 'kegg2rhea_bi.pickle'), 'rb') as handle:
  ref_kegg2rhea_bi = pickle.load(handle)
with open(os.path.join(cn.RHEA_DIR, 'rhea2chebi_reference.pkl'), 'rb') as f:
  ref_rhea_to_chebi = pickle.load(f)
with open(os.path.join(cn.CHEBI_DIR, 'chebi_shortened_formula_30apr2022.pickle'), 'rb') as f:
  ref_shortened_chebi_to_formula = pickle.load(f)
with open(os.path.join(cn.ALGO_DIR, 'binary_ref_df.pickle'), 'rb') as handle:
    ref_mat = pickle.load(handle)


class ReactionAnnotation(object):

  def __init__(self, libsbml_fpath=None, exist_qualifier=cn.RHEA):
    # self.exist_annotation stores 
    # existing KEGG Reaction or Rhea annotations in the model.
    # If none exists, set None.
    if libsbml_fpath is not None:
      reader = libsbml.SBMLReader()
      document = reader.readSBML(libsbml_fpath)
      self.model = document.getModel()
      # Annotation of Rhea
      reac_dict_raw_rhea = {r.getId():tools.getQualifierFromString(r.getAnnotationString(), cn.RHEA) \
                           for r in self.model.getListOfReactions()}
      reac_dict_raw_filt_rhea = {k:reac_dict_raw_rhea[k] \
                                 for k in reac_dict_raw_rhea.keys() \
                                 if reac_dict_raw_rhea[k] is not None}
      reac_dict_format_rhea = {k:['RHEA:'+val for val in reac_dict_raw_filt_rhea[k]] \
                                 for k in reac_dict_raw_filt_rhea.keys()}
      reac_dict_rhea = dict()
      for one_id in reac_dict_format_rhea.keys():
        one_itm = list(set([ref_rhea2bi[val] for val in reac_dict_format_rhea[one_id] \
                   if val in ref_rhea2bi.keys()]))
        if len(one_itm) > 0:
          reac_dict_rhea[one_id] = one_itm
      # Annotation of KEGG (mapped to corresponding Rhea BI term) 
      reac_dict_raw_kegg = {r.getId():tools.getQualifierFromString(r.getAnnotationString(), cn.KEGG_REACTION) \
                           for r in self.model.getListOfReactions()}
      reac_dict_raw_filt_kegg = {k:reac_dict_raw_kegg[k] \
                                 for k in reac_dict_raw_kegg.keys() \
                                 if reac_dict_raw_kegg[k] is not None}
      reac_dict_kegg = {k:[ref_kegg2rhea_bi[val] \
                           for val in reac_dict_raw_filt_kegg[k] if val in ref_kegg2rhea_bi.keys()] \
                        for k in reac_dict_raw_filt_kegg.keys()}
      exist_annotation = reac_dict_rhea
      for one_id in reac_dict_kegg.keys():
        if one_id in exist_annotation.keys():
          exist_annotation[one_id] = list(set(exist_annotation[one_id] + reac_dict_kegg[one_id]))
      self.exist_annotation = exist_annotation
    else:
      self.model = None
      self.exist_annotation = None
    # Attributes after prediction
    self.candidates = None
    self.match_score = None
    self.sum_match_score = None
    self.query_df = None

  def getMatchScore(self, inp_match_score):
    """
    Calculate match scores using 
    the given dictionary
  
    Parameters
    ----------
    inp_match_score: dict {reaction_id: {candidate_id: float}}
  
    Returns
    -------
    match_score: float
    """
    return np.sum([np.max([inp_match_score[val][k] \
                           for k in inp_match_score[val].keys()]) \
                   for val in inp_match_score.keys()])

  def getReactionComponents(self,
                            inp_reaction):
    """
    Get component of reactions in species IDs
    (both reactants and products)
    of a reaction. 

    Parameters
    ----------
    reaction_id: str/libsbml.Reaction


    Returns
    -------
    r_components: None/str-list (list of species IDs)
    """
    if isinstance(inp_reaction, libsbml.Reaction):
      one_reaction = inp_reaction
    elif isinstance(inp_reaction, str):
      one_reaction = self.model.getReaction(inp_reaction)
    else:
      return None
    reactants = [val.species for val in one_reaction.getListOfReactants()]
    products = [val.species for val in one_reaction.getListOfProducts()]
    r_components = list(set(reactants + products))
    return r_components

  def predictAnnotation(self,
                        inp_spec_dict,
                        inp_reac_list=None,
                        inp_ref_mat=ref_mat):
    """
    Predict 1) reaction annotation candidates 
    and 2) match score of them
    using species dict (argument) etc.
    Also updates several class attributes as a result. 
  
    Parameters
    ----------
    dict: inp_spec_dict
        Dictionoary, {species id: formula(str-list)}
    inp_reac_list: str-list
        IDs of reactions to predict. If default, will do all reactions
      
    Returns
    -------
    pred_match_score: dict
        Confidence score of each prediction
        {reaction ID: {Rhea ID: float between 0.0-1.0}}
    """
    # get libsbml.reaction and their IDs
    if inp_reac_list is not None:
      reactions = [self.model.getReaction(val) for val in inp_reac_list]
      reaction_ids = inp_reac_list
    else:
      reactions = self.model.getListOfReactions()
      reaction_ids = [val.getId() for val in reactions]
    # get dictionary of reaction ID: species component
    r2pred_spec_formulas = dict()
    for one_reaction in reactions:
      r2pred_spec_formulas[one_reaction.getId()] = {one_spec:inp_spec_dict[one_spec] \
                                                    for one_spec in self.getReactionComponents(one_reaction)}
    # prepare query df for prediction
    query_df = pd.DataFrame(0, 
                            index=inp_ref_mat.columns,
                            columns=reaction_ids)
    for one_rid in reaction_ids:
      one_set_species = r2pred_spec_formulas[one_rid]
      # for each species element of the select reaction
      for one_spec_key in one_set_species.keys():
        one_spec = one_set_species[one_spec_key]
        # For each one_rid, set the values 1.0
        query_df.loc[[val for val in one_spec if val in query_df.index], one_rid] = 1
    multi_mat = inp_ref_mat.dot(query_df)
    maxes = multi_mat.max()
    #
    # Collect candidates and calculate confidence score
    pred_cands = dict()
    pred_match_score = dict()
    for one_rid in maxes.index:
      one_multi = multi_mat.loc[:,one_rid]
      candidates = one_multi[one_multi==maxes[one_rid]].index
      # cand_data; (number of element matches, candidates)
      pred_cands[one_rid] = candidates
      # Now, confidence (calculated per each candidate)
      match_score_per_cand = dict()
      for one_cand in candidates:
        if one_cand in ref_rhea2bi.keys():
          num_matches = maxes[one_rid]
          num_maxpos_matches = len(inp_ref_mat.loc[one_cand, :].to_numpy().nonzero()[0])
          match_score_per_cand[one_cand] = num_matches / num_maxpos_matches
      pred_match_score[one_rid] = match_score_per_cand
    self.candidates = pred_cands
    self.match_score = pred_match_score
    self.sum_match_score = self.getMatchScore(pred_match_score)
    self.query_df = query_df
    return pred_match_score


  def updateSpeciesByAReaction(self, 
                               inp_rid, inp_spec_dict,
                               inp_rhea, inp_ref_mat):
    """
    Update predicted species annotation
    using predicted rhea annotation candidate
    for a single reaction.
    Current version works for when there is just
    one candidate (i.e., one RHEA term).  
  
    Parameters
    ----------
    inp_rid: str
        Reactino ID to match the candidate
    inp_spec_dict: dict
        {species ID: chemical formula}
    inp_rhea: str
        A RHEA term predicted 
    inp_ref_mat: pandas.DataFrame
  
    Returns
    -------
    dict: {species ID: [list of CHEBI terms]}
          Suggested mapping from species ID to a chebi term
    """
    # Chebi terms (reference) associated with the given Rhea term
    # If there is no such Rhea term in the reference, return None
    if inp_rhea in ref_rhea2bi.keys():
      if ref_rhea2bi[inp_rhea] in ref_rhea_to_chebi.keys():
        rhea_term_to_chebi_elements = [val for val in ref_rhea_to_chebi[ref_rhea2bi[inp_rhea]] \
                                       if val in ref_shortened_chebi_to_formula.keys()]
      else:
        return None
    else:
      return None
    # Similarly, species formula (reference) associated with the Rhea term. 
    one_r_elements_row = inp_ref_mat.loc[inp_rhea, :]
    one_r_elements = one_r_elements_row[one_r_elements_row!=0].index
    # Next, species formulas (query) associated with the predicted Rhea term.
    one_r_query_elements_row = self.query_df.loc[:, inp_rid]
    one_r_query_elements = one_r_query_elements_row[one_r_query_elements_row!=0].index
    #
    # Identify species not included vs. species included (in query), compared to ref_mat
    # species_not_included = [val for val in one_r_elements if val not in one_r_query_elements]
    formula_included = [val for val in one_r_elements if val in one_r_query_elements]
    # {predicted species ID: formula} for all elements in the reaction 
    all_species_in_a_reaction = self.getReactionComponents(inp_rid)
    spec2predicted_formula = {one_spec:inp_spec_dict[one_spec] \
                                       for one_spec in all_species_in_a_reaction}  
    #
    possibly_correct_species_ids = []
    chebi_term_already_used_in_prediction = []
    for one_incl_formula in formula_included:
      for one_specid in spec2predicted_formula.keys():
        # If one included formula can be found in one of the predicted spec_dict,
        # it will be possibly correct (we are looking for incorrect ones)
        if one_incl_formula in spec2predicted_formula[one_specid]:
          possibly_correct_species_ids.append(one_specid)
          for one_chebi_term in rhea_term_to_chebi_elements:
            if one_incl_formula == ref_shortened_chebi_to_formula[one_chebi_term]:
              chebi_term_already_used_in_prediction.append(one_chebi_term)
    remaining_chebi = [val for val in rhea_term_to_chebi_elements \
                       if val not in chebi_term_already_used_in_prediction]
    remaining_specid = [val for val in all_species_in_a_reaction \
                        if val not in possibly_correct_species_ids]
    if len(remaining_specid) == 1 and len(remaining_chebi) == 1:    
      return {remaining_specid[0]: [remaining_chebi[0]]}
    else:
      None
