import ctypes
import re
import argparse
import copy
import time
from bidict import bidict
from scipy import sparse
import numpy as np


NEWATOM_PREFIX = "x_"
KW_FACT = "FACT"
KW_AND = "AND"
KW_OR = "OR"


C_LPMAT = ctypes.CDLL('../build/liblpmatrixgrounder.so')


def make_clist(lst):
    return (ctypes.c_char_p * len(lst))(*[x.encode() for x in lst])

# filename = "../samples/example.lp".encode()
# C_LPMAT.ground_single.argtypes = [ctypes.c_char_p]
# C_LPMAT.ground_single.restype = ctypes.c_char_p
# lp = C_LPMAT.ground_single(filename)

# print(lp)

# listof_files = ["../samples/model_graphcoloring.lp", "../samples/graph_6nodes.lp"]
# C_LPMAT.ground_list.argtypes = [ctypes.POINTER(ctypes.c_char_p), ctypes.c_size_t]
# C_LPMAT.ground_list.restype = ctypes.c_char_p
# lp = C_LPMAT.ground_list(make_clist(listof_files), len(listof_files))

# print(lp)


def remapAtoms(atoms_mapping, atoms_mapping_negations, atoms_mapping_newatom):
    new_atoms_mapping = bidict()

    all_symbols = list(set(atoms_mapping.keys()).difference(
        atoms_mapping_negations).difference(atoms_mapping_newatom))
    all_symbols = sorted([x for x in all_symbols])

    all_symbols_negated = sorted([x for x in atoms_mapping_negations])
    all_symbols_introduced = sorted([x for x in atoms_mapping_newatom])

    atoms_mapping_curidx = 0
    for sym in all_symbols + all_symbols_introduced + all_symbols_negated:
        new_atoms_mapping[sym] = atoms_mapping_curidx
        atoms_mapping_curidx += 1

    return new_atoms_mapping


class grounder_ret(ctypes.Structure):
    _fields_ = [
        ("rawdata", ctypes.c_char_p),
        ("duration_clingo", ctypes.c_int),
        ("duration_internal", ctypes.c_int),
    ]


def groundUsingClingo(filenames):
    C_LPMAT.ground_list.argtypes = [
        ctypes.POINTER(ctypes.c_char_p), ctypes.c_size_t]
    C_LPMAT.ground_list.restype = grounder_ret
    ret_data = C_LPMAT.ground_list(make_clist(filenames), len(filenames))
    aspif = ret_data.rawdata
    # print(aspif)
    rawdata = aspif.decode("utf-8").split("\n")
    return hande_aspif(rawdata), \
        float(ret_data.duration_clingo * 10**-9), \
        float(ret_data.duration_internal * 10**-9)


def hande_aspif(rawdata):
    # rules = dict()
    rules_mapped = dict()
    # {0: ('OR', [...]), 1: ('FACT', [...]), 4: ('AND', [...]), 5: ('AND', [...])}

    introduced_atoms_counter = 1
    atoms_negation = set()
    atoms_introduced = set()
    atoms_max = 0
    for line in rawdata:
        # print(line)
        if "<=" not in line or "." not in line:
            continue
        parts = re.split('<=', line.strip().replace(".", ""))
        head = int(parts[0])
        if head > atoms_max:
            atoms_max = head
        if len(parts) > 1:
            body = [int(x) for x in parts[1].strip().split(',')]
            for b in body:
                if b < 0:
                    atoms_negation.add(b)
                    b = -b
                if b > atoms_max:
                    atoms_max = head
        else:   # handle fact
            body = list()

        # update rule
        if head not in rules_mapped:
            if len(body) > 0:
                rules_mapped[head] = (KW_AND, list())
            else:
                rules_mapped[head] = (KW_FACT, list())
        rules_mapped[head][1].append(body)

    # identify mapping for negations
    atoms_negation = sorted(list(atoms_negation), reverse=True)
    # print(atoms_max, atoms_negation)
    introduced_atoms_counter = atoms_max + 1
    # rules_orginal = copy.deepcopy(rules_mapped)

    updated_rules = dict()
    for head, (ruletype, rulebodies) in rules_mapped.items():
        # print(head, rulebodies)
        if len(rulebodies) > 1:
            for idx, body in enumerate(rulebodies):
                if len(body) > 1:
                    # print("introduce new atom for", body)
                    natom = introduced_atoms_counter
                    atoms_introduced.add(natom)
                    introduced_atoms_counter += 1
                    rules_mapped[head][1][idx] = [natom]
                    updated_rules[natom] = (KW_AND, list())
                    updated_rules[natom][1].append(body)
            nbody = list()
            for rb in rules_mapped[head][1]:
                nbody += rb
            rules_mapped[head] = (KW_OR, [nbody])
        elif len(rulebodies) > 0:
            if len(rulebodies[0]) == 0:
                rules_mapped[head] = (KW_FACT, rulebodies)
    rules_mapped.update(updated_rules)

    # mapping for negations
    atoms_negation_mapped = list(atoms_negation)
    for idx, atom in enumerate(atoms_negation):
        atoms_negation_mapped[idx] = introduced_atoms_counter + idx

    # update mapping for negations
    for head, (ruletype, rulebodies) in rules_mapped.items():
        for idx, body in enumerate(rulebodies):
            newbody = [-b if b < 0 else b for b in body]
            if body != newbody:
                rules_mapped[head][1][idx] = newbody

    return rules_mapped, introduced_atoms_counter + len(atoms_negation), atoms_negation_mapped


def buildMatrixMapped(rules, n_atoms, atoms_negation):
    mp = np.zeros((n_atoms, n_atoms), dtype=float)
    for head, (ruletype, rulebodies) in rules.items():
        for body in rulebodies:
            if ruletype == KW_AND:
                val = 1.0 / len(body)
                for b in body:
                    mp[head, b] = val
            elif ruletype == KW_OR:
                for b in body:
                    mp[head, b] = 1.0
            elif ruletype == KW_FACT:
                mp[head, head] = 1.0
    for neg in atoms_negation:
        mp[neg, neg] = 1.0
    return mp


pre_row = np.empty(100_000_000, dtype=np.int32)
pre_col = np.empty(100_000_000, dtype=np.int32)
pre_val = np.empty(100_000_000, dtype=float)


def buildMatrixMapped_sparsefast(rules, n_atoms, atoms_negation):
    nnz = 0
    for h, (ruletype, rulebodies) in rules.items():
        for body in rulebodies:
            if ruletype == KW_AND:
                v = 1.0 / len(body)
                for b in body:
                    pre_row[nnz] = h
                    pre_col[nnz] = b
                    pre_val[nnz] = v
                    nnz = nnz + 1
            elif ruletype == KW_OR:
                for b in body:
                    pre_row[nnz] = h
                    pre_col[nnz] = b
                    pre_val[nnz] = 1.0
                    nnz = nnz + 1
            elif ruletype == KW_FACT:
                pre_row[nnz] = h
                pre_col[nnz] = h
                pre_val[nnz] = 1.0
                nnz = nnz + 1
    for neg in atoms_negation:
        pre_row[nnz] = neg
        pre_col[nnz] = neg
        pre_val[nnz] = 1.0
        nnz = nnz + 1
    ms = sparse.coo_matrix((pre_val[: nnz], (pre_row[: nnz], pre_col[: nnz])), shape=(
        n_atoms, n_atoms)).tocsr()
    return ms


def buildMatrixMapped_sparse(rules, n_atoms, atoms_negation):
    row, col, val = [], [], []
    for h, (ruletype, rulebodies) in rules.items():
        for body in rulebodies:
            if ruletype == KW_AND:
                v = 1.0 / len(body)
                for b in body:
                    row.append(h)
                    col.append(b)
                    val.append(v)
            elif ruletype == KW_OR:
                for b in body:
                    row.append(h)
                    col.append(b)
                    val.append(1.0)
            elif ruletype == KW_FACT:
                row.append(h)
                col.append(h)
                val.append(1.0)
    for neg in atoms_negation:
        row.append(neg)
        col.append(neg)
        val.append(1.0)
    ms = sparse.coo_matrix((val, (row, col)), shape=(
        n_atoms, n_atoms)).tocsr()
    return ms


def main():
    parser = argparse.ArgumentParser(description='LP Matrix Grounder')
    parser.add_argument('--input', help='Input file',
                        type=argparse.FileType('r'), nargs='+')
    parser.add_argument('--output', help='Output file')

    args = parser.parse_args()

    inputfiles = args.input
    outputfile = args.output

    start_time = time.time()
    (rules, n_atoms, atoms_negation), duration_clingo, duration_internal = groundUsingClingo(
        [f.name for f in inputfiles])
    executation_time = time.time() - start_time
    print("Grounding time              :", executation_time)
    print("---- duration_clingo        :", duration_clingo)
    print("---- duration_internal      :", duration_internal)
    print("---- overhead               :", executation_time -
          duration_clingo - duration_internal)
    # print(rules)

    start_time = time.time()
    mp = buildMatrixMapped(rules, n_atoms, atoms_negation)
    executation_time = time.time() - start_time
    print("Building dense matrix time  :", executation_time)
    # print(mp)

    start_time = time.time()
    ms = buildMatrixMapped_sparse(rules, n_atoms, atoms_negation)
    executation_time = time.time() - start_time
    print("Building sparse matrix time :", executation_time)
    # print(ms)

    start_time = time.time()
    ms = buildMatrixMapped_sparsefast(rules, n_atoms, atoms_negation)
    executation_time = time.time() - start_time
    print("Building fast sparse time   :", executation_time)
    # print(ms)


if __name__ == "__main__":
    main()
