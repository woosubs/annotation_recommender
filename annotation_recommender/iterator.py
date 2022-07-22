# iterator.py
"""
Script for iteration;
can be refactored into class in future..
"""

from annotation_recommender import constants as cn
from annotation_recommender import tools
from annotation_recommender import species_annotation as sa
from annotation_recommender import reaction_annotation as ra


def iterateAndGetUpdatedResults(cur_candidates_dict,
                                cur_spec_formula_dict,
                                cur_reac_match_score,
                                inp_one_cands = None,
                                num_iter=10,
                                iter_func=updateSpeciesByAReaction,
                                show_message=False):
  """
  Using current species and reaction annotations
  (for species, get formulafied one),
  update both species and reaction annotations
  
  Parameters
  ----------
  cur_canadidates_dict: dict
  cur_spec_formula_dict: dict 
  cur_reac_match_score: float
  inp_one_cands: dict
      If not provided, create its own
  num_iter: int
  iter_func: function
  show_message: bool
      Prints out iteration number and match score
      
  Returns
  -------
  ?:
  rep: int 
      Last iteration position (actual run - 1) when algorithm quit
  """
  flag = False
  if show_message:
    print("Initial match score: %.02f" % cur_reac_match_score)
    print("*************************")
  for rep in range(0, num_iter):
    if flag:
      break
    if show_message:
      print("Iteration %d" % (rep+1))
    if inp_one_cands is None:
      one_cands = {val:cur_candidates_dict[val] for val in cur_candidates_dict.keys() if len(cur_candidates_dict[val])==1}
    else: 
      one_cands = inp_one_cands
    all_upd_spec = dict()
    for one_k in one_cands.keys():
      one_upd_spec = iter_func(inp_rid=one_k,
                                                 inp_model=model,
                                                 inp_spec_dict=cur_spec_formula_dict,
                                                 inp_rhea=one_cands[one_k][0],
                                                 inp_ref_mat=ref_mat,
                                                 inp_query_df=query_df)
      all_upd_spec = ip.updateDictKeyToList(all_upd_spec, one_upd_spec)
    #
    # update species dictionary to use
    upd_spec_formula_dict = dict()
    for one_k in cur_spec_formula_dict.keys():
      if one_k in all_upd_spec.keys():
        upd_spec_formula_dict[one_k] = list(set([ref_shortened_chebi_to_formula[val] for val in all_upd_spec[one_k]]))
      else:
        upd_spec_formula_dict[one_k] = cur_spec_formula_dict[one_k]
    
    #
    # Using upd_spec_formula_dict, predict reaction again
    new_candidates_dict, new_candidates_conf, new_query_df = ip.getReactionPredResults(inp_spec_dict=upd_spec_formula_dict,
                                                                           inp_model=model,
                                                                           inp_list_reactions=None,
                                                                           inp_ref_mat=ref_mat)
    #
    upd_reac_match_score = ip.getMatchScore(new_candidates_conf)
    # Check wheter to continue;
    if upd_reac_match_score > cur_reac_match_score:
      cur_candidates_dict = new_candidates_dict
      cur_spec_formula_dict = upd_spec_formula_dict
      cur_reac_match_score = upd_reac_match_score
      if show_message:
        print("Updated match score: %.02f" % cur_reac_match_score)
        print("*************************")
    else:
      flag = True
      if show_message:
        print("Updated match score: %.02f" % cur_reac_match_score)
        print("Score not increasing. Quitting iteration...")
  if show_message:
    print("\nCalculation finished.")
  return cur_candidates_dict, cur_spec_formula_dict, cur_reac_match_score, new_query_df, rep



def updateSpeciesByAReaction(inp_rid, inp_model,
                             inp_spec_dict, inp_rhea,
                             inp_ref_mat, inp_query_df):
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
  inp_model: libsbml.Model
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
  # Chebi terms associated with the given Rhea term
  # if there is no such Rhea term in the reference, return None
  if inp_rhea in ref_rhea2bi.keys():
    if ref_rhea2bi[inp_rhea] in ref_rhea_to_chebi.keys():
      rhea_term_to_chebi_elements = [val for val in ref_rhea_to_chebi[ref_rhea2bi[inp_rhea]] \
                                     if val in ref_shortened_chebi_to_formula.keys()]

    else:
      return None
  else:
  	return None
  # From ref_mat, checks reference-level components in species formula; 
  one_r_elements_row = inp_ref_mat.loc[inp_rhea, :]
  one_r_elements = one_r_elements_row[one_r_elements_row!=0].index
  # From query_df (represents what species formula was used to predict the Rhea term)
  one_r_query_elements_row = inp_query_df.loc[:, inp_rid]
  one_r_query_elements = one_r_query_elements_row[one_r_query_elements_row!=0].index
  # Identify species not included vs. species included (in query), compared to ref_mat
  species_not_included = [val for val in one_r_elements if val not in one_r_query_elements]
  species_included = [val for val in one_r_elements if val in one_r_query_elements]
  # Using species_included
  one_reaction = inp_model.getReaction(inp_rid)
  reactants = [val.species for val in one_reaction.getListOfReactants()]
  products = [val.species for val in one_reaction.getListOfProducts()]
  all_species_in_a_reaction = list(set(reactants + products)) 
  # {predicted species ID: formula} for all elements in the reaction 
  spec2predicted_formula = {one_spec:inp_spec_dict[one_spec] \
                                     for one_spec in all_species_in_a_reaction}  
  #
  possibly_correct_species_ids = []
  chebi_term_already_used_by_original = []
  for one_incl_spec in species_included:
    for one_pos_incor_spec in spec2predicted_formula.keys():
      # if one included species formula is included in one of the predicted spec_dict,
      # it will be possibly correct (we are looking for incorrect ones)
      if one_incl_spec in spec2predicted_formula[one_pos_incor_spec]:
        possibly_correct_species_ids.append(one_pos_incor_spec)
        for one_chebi_term in rhea_term_to_chebi_elements:
          if one_incl_spec == ref_shortened_chebi_to_formula[one_chebi_term]:
            chebi_term_already_used_by_original.append(one_chebi_term)
  remaining_chebi = [val for val in rhea_term_to_chebi_elements \
                     if val not in chebi_term_already_used_by_original]
  remaining_specid = [val for val in all_species_in_a_reaction \
                      if val not in possibly_correct_species_ids]
  if len(remaining_specid) == 1 and len(remaining_chebi) == 1:    
    return {remaining_specid[0]: [remaining_chebi[0]]}
  else:
    None