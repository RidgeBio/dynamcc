import functools
import re
from tornado.web import RequestHandler

from dynamcc_handlers import *

class SequenceHandler(RequestHandler):
    def get(self):
        self.render("sequence.html")

    def post(self):
        if "table" in self.request.files:
            sorted_dict = util.BuildCustomUsageDict(self.request.files["table"][0])
            organism_name = "user uploaded usage table"
        else:
            selected_organism = self.get_argument("usage_table")
            if selected_organism in organism_mapping:
                sorted_dict = util.BuildUsageDict(organism_mapping[selected_organism])
                organism_name = organism_names[selected_organism]
            else:
                raise ValueError("unknown organism: ", selected_organism)

        compression_method = self.get_argument("compression_method")
        if compression_method == "rank":
            threshold = self.get_argument("input_rank")
        else:
            threshold = self.get_argument("input_usage")
        backbone = self.get_argument("backbone")
        edits = self.get_argument("edits")

        sequences = optimise_codons_seq(backbone, edits, sorted_dict, compression_method, threshold)
        self.render("sequence_results.html", organism=organism_name, backbone=backbone, edits=edits, sequences=sequences, rank=compression_method == "rank", threshold_value=threshold)

_EDIT_PATTERN = re.compile(r"([A-Z])(\d+)(-?)([A-Z]+)")

def _parse_edit(edit):
    '''Parses an edit.'''
    match = _EDIT_PATTERN.match(edit)
    if not match:
        raise ValueError("Invalid edit: {}".format(edit))
    orig = match.group(1)
    pos = int(match.group(2)) - 1
    invert = match.group(3) == '-'
    new = set(match.group(4))
    if invert:
        new = aa.difference(new)
    return (orig, pos, new)

def optimise_codons_seq(aa_seq, edits, usage_dict, compression_method, threshold):
    edited_positions = dict()
    for edit in edits.split(','):
        if not edit:
            continue
        prev_aa, pos, new_aa = _parse_edit(edit)
        if pos >= len(aa_seq) or aa_seq[pos] != prev_aa:
            raise ValueError("Edit '{}' is invalid".format(edit))
        if pos in edited_positions:
            edited_positions[pos] = edited_positions[pos] + new_aa
        else:
            edited_positions[pos] = new_aa
    res = None
    for pos in range(0, len(aa_seq)):
        if pos in edited_positions:
            codons = optimise_codons(edited_positions.get(pos), usage_dict, compression_method, threshold)
        else:
            codons = optimise_codons(aa_seq[pos], usage_dict, compression_method, threshold)
        if res:
            res = itertools.product(res, codons)
            res = ["".join(r) for r in res]
        else:
            res = codons
    return res


def memoize(func):
    cache = {}
    
    def clear_cache():
        cache.clear()

    @functools.wraps(func)
    def wrapper(*args):
        key = tuple(str(a) for a in args)
        if key not in cache:
            cache[key] = func(*args)
        return cache[key]
    wrapper.clear_cache = clear_cache
    
    return wrapper

@memoize
def optimise_codons(selected_aa, usage_dict, compression_method, threshold):
    rules_dict, _ = util.BuildRulesDict('rules.txt')
    remove_aa = list(aa.difference(selected_aa).union({'X'}))
    filtered_dict = util.EditUsageDict(remove_aa, usage_dict.copy())
    if compression_method == 'rank':
        selection = 'R'
        new_dict = RemoveCodonByRank(int(threshold), filtered_dict)
    else:
        selection = 'U'
        new_dict = RemoveLowCodons(float(threshold), filtered_dict)
    
    codon_order = new_dict.keys()
    codon_count = BuildCodonCount(new_dict, codon_order)
    redundancy = 0
    best_result = start_multiprocessing(new_dict, rules_dict, selection, codon_count, redundancy, processes = 3)
    return best_result['BestReducedList']
