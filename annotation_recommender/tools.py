# tools.py

import re


def getOntologyFromString(string_annotation):
  """
  Parse string and return string annotation,
  marked as <bqbiol:is> or <bqbiol:isVersionOf>.
  If neither exists, return None.

  Parameters
  ----------
  str: string_annotation

  Returns
  -------
  list-tuple (ontology type, ontology id)
       Return [] if none is provided
  """
  # first, extracts strings tagged as bqbiol:is or bqbiol:isVersionOf.
  is_str = ''
  isVersionOf_str = ''
  is_str_match = re.findall('<bqbiol:is[^a-zA-Z].*?<\/bqbiol:is>',
                            string_annotation,
                            flags=re.DOTALL)
  if len(is_str_match)>0:
    is_str_match_filt = [s.replace("      ", "") for s in is_str_match]
    is_str = '\n'.join(is_str_match_filt)  
  #
  is_VersionOf_str_match = re.findall('<bqbiol:isVersionOf[^a-zA-Z].*?<\/bqbiol:isVersionOf>',
                                      string_annotation,
                                      flags=re.DOTALL)  
  #
  if len(is_VersionOf_str_match) > 0:
    is_VersionOf_str_match_filt = [s.replace("      ", "") for s in is_VersionOf_str_match]
    isVersionOf_str = '\n'.join(is_VersionOf_str_match_filt) 
  #
  combined_str = is_str + isVersionOf_str
  if combined_str == '':
    return []
  identifiers_list = re.findall('identifiers\.org/.*/', combined_str)
  return [(r.split('/')[1],r.split('/')[2].replace('\"', '')) \
          for r in identifiers_list]

def getQualifierFromString(input_str, qualifier):
  """
  Parses string and returns an identifier. 
  If not, return None

  Parameters
  ----------
  str: string_annotation

  Returns
  -------
  str (ontology Id)
      Return None if none is provided
  """
  ontologies = getOntologyFromString(input_str)
  # To make sure it works, make it lower
  qualifier_list = [val for val in ontologies if val[0]==qualifier.lower()]
  if qualifier_list:
    return [val[1] for val in qualifier_list]
  else:
    return None


def transformCHEBIToFormula(inp_list, ref_to_formula_dict):
  """
  transform input list of CHEBI terms
  to list of annotations. 
  
  Parameters
  ----------
  inp_list: str-list
  
  Returns
  -------
  res: str-list
  """
  inp_formulas = [ref_to_formula_dict[val] for val in inp_list \
                  if val in ref_to_formula_dict.keys()]
  res = list(set([val for val in inp_formulas if val is not None]))
  return res