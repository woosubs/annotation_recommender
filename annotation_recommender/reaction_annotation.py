"""Annotation for Reactions."""

from annotation_recommender import constants as cn
from annotation_recommender import tools

import libsbml
import numpy as np
import os
import pickle


with open(os.path.join(cn.RHEA_DIR, 'rhea_all2bi.pkl'), 'rb') as f:
  ref_rhea2bi = pickle.load(f)
with open(os.path.join(cn.RHEA_DIR, 'kegg2rhea_bi.pickle'), 'rb') as handle:
  ref_kegg2rhea_bi = pickle.load(handle)


class ReactionAnnotation(object):

  def __init__(self, libsbml_fpath=None, exist_qualifier='rhea'):
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
      self.exist_annotation = None