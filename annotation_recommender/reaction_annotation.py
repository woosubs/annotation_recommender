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
    self.query_df = None

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
    and 2) confide`nce of them
    using species dict (argument) etc.
  
    Parameters
    ----------
    dict: inp_spec_dict
        Dictionoary, {species id: formula(str-list)}
    inp_reac_list: str-list
        IDs of reactions to predict. If default, will do all reactions
      
    Returns
    -------
    pred_cands: dict
        Predicted candidates
        {reaction ID: list of RHEA IDs}
    pred_conf: dict
        Confidence score of each prediction
        {reaction ID: {Rhea ID: float between 0.0-1.0}}
    query_df: pandas.DataFrame
        Query Dataframe for reference 
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
    self.query_df = query_df
    return pred_match_score






