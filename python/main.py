import ctypes
import re
import argparse
import copy
import time
import os
from bidict import bidict
from scipy import sparse
import numpy as np
import pandas as pd


NEWATOM_PREFIX = "x_"
KW_FACT = "FACT"
KW_AND = "AND"
KW_OR = "OR"


C_LPMAT = ctypes.CDLL('../build/liblpmatrixgrounder.so')


def make_clist(lst):
    return (ctypes.c_char_p * len(lst))(*[x.encode() for x in lst])


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
        ("duration_clingo", ctypes.c_float),
        ("duration_internal", ctypes.c_float),
    ]


def groundUsingClingo(filenames):
    C_LPMAT.ground_list.argtypes = [
        ctypes.POINTER(ctypes.c_char_p), ctypes.c_size_t]
    C_LPMAT.ground_list.restype = grounder_ret
    ret_data = C_LPMAT.ground_list(make_clist(filenames), len(filenames))
    aspif = ret_data.rawdata
    rawdata = aspif.decode("utf-8").split("\n")
    start_time = time.time()
    ret = handle_aspif(rawdata)
    duration = time.time() - start_time
    return ret, duration, \
        ret_data.duration_clingo, \
        ret_data.duration_internal


def handle_aspif(rawdata):
    rules_mapped = dict()
    introduced_atoms_counter = 1
    atoms_negation = set()
    atoms_introduced = set()
    atoms_max = 0
    original_rule_count = 0
    for line in rawdata:
        if "." not in line:
            continue
        original_rule_count += 1
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
    introduced_atoms_counter = atoms_max + 1

    updated_rules = dict()
    for head, (ruletype, rulebodies) in rules_mapped.items():
        if len(rulebodies) > 1:
            for idx, body in enumerate(rulebodies):
                if len(body) > 1:
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

    return rules_mapped, atoms_max, introduced_atoms_counter + len(atoms_negation), \
        atoms_negation_mapped, original_rule_count


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


ARRAY_SIZE = 1_000_000_000
pre_row = np.empty(ARRAY_SIZE, dtype=np.int32)
pre_col = np.empty(ARRAY_SIZE, dtype=np.int32)
pre_val = np.empty(ARRAY_SIZE, dtype=float)


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
    return ms, nnz


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


def main2():
    parser = argparse.ArgumentParser(description='LP Matrix Grounder')
    parser.add_argument('--input', help='Input file',
                        type=argparse.FileType('r'), nargs='+')
    parser.add_argument('--output', help='Output file')
    args = parser.parse_args()
    inputfiles = args.input
    outputfile = args.output
    conductExperiment([f.name for f in inputfiles])


def conductExperiment(inputfiles):
    start_time = time.time()
    (rules, n_atoms_original, n_atoms, atoms_negation, original_rule_count), \
        duration_handle, duration_clingo, duration_internal = groundUsingClingo(inputfiles)
    executation_groundtime = time.time() - start_time
    
    print("Grounding time              :", executation_groundtime)
    print("---- duration_clingo        :", duration_clingo)
    print("---- duration_internal      :", duration_internal)
    print("---- duration_handle        :", duration_handle)
    print("---- overhead               :", executation_groundtime -
          duration_clingo - duration_internal - duration_handle)

    start_time = time.time()
    try:
        mp = buildMatrixMapped(rules, n_atoms, atoms_negation)
        executation_densematrix_time = time.time() - start_time
        print("Building dense matrix time  :", executation_densematrix_time)
        print(mp)
    except np.core._exceptions._ArrayMemoryError:
        print("Unable to allocate an array with shape (%d, %d) and data type float64" % (
            n_atoms, n_atoms))
        executation_densematrix_time = "NA"

    start_time = time.time()
    ms = buildMatrixMapped_sparse(rules, n_atoms, atoms_negation)
    executation_sparsematrix_time = time.time() - start_time
    print("Building sparse matrix time :", executation_sparsematrix_time)

    start_time = time.time()
    ms, nnz = buildMatrixMapped_sparsefast(rules, n_atoms, atoms_negation)
    executation_sparsematrixfast_time = time.time() - start_time
    print("Building fast sparse time   :", executation_sparsematrixfast_time)

    avg_array = []
    for head, (ruletype, rulebodies) in rules.items():
        for idx, body in enumerate(rulebodies):
            avg_array.append(len(body))

    dataset = "name"    # will be updated
    no_rules_original = original_rule_count
    no_atoms = n_atoms_original
    no_negations = len(atoms_negation)
    no_rules = len(rules)
    matrix_size = n_atoms
    avg_rulelength = np.average(avg_array)
    nnz = nnz
    sparsity = 1.0 - nnz / (matrix_size * matrix_size)
    dense_size = matrix_size * matrix_size * 4
    sparse_size = (len(ms.indices) + len(ms.indptr) + len(ms.data)) * 4
    clingo_time = duration_clingo
    internal_time = duration_internal
    handling_time = duration_handle
    grounding_time = executation_groundtime
    dense_matrix_time = executation_densematrix_time
    sparse_matrix_time = executation_sparsematrix_time
    sparse_matrix_fast_time = executation_sparsematrixfast_time

    return dataset, no_rules_original, no_atoms, no_negations, \
        no_rules, matrix_size, avg_rulelength,  \
        nnz, \
        clingo_time, internal_time, handling_time, \
        dense_matrix_time, sparse_matrix_fast_time


def format_tex(number):
    return "$\\num{%d}$" % number


def main():
    df = pd.DataFrame(columns=['Data instance', '$m$', '$n$', '$k$', '$m\'$', '$n\'$', '$\overline{l}$',
                               '$\eta_z$',
                               '$t_{clingo}$', '$t_{int}$', '$t_{alg}$',
                               '$t_{asg}^{(d)}$', '$t_{asg}^{(s)}$'])
    datasets = dict({"3-Coloring": ("../experiments/datasets/", "3Coloring-Clasp.txt", "instances"),
                    #  "HamiltonianCycle": ("../experiments/datasets/", "HamiltonianCycle-Clasp.txt", "instances"),
                     #  "LatinSquares": ("../experiments/datasets/", "LatinSquares-Clasp.txt", "instances"),
                     #  "NQueens": ("../experiments/datasets/", "queens-Clasp.txt", "instances"),
                     #  "TransitiveClosure": ("../experiments/datasets/", "Reachability-Clasp.txt", "instances")
                     })
    for dataname, (basedir, modelfile, instancedir) in datasets.items():
        modelfile = os.path.join(basedir, dataname, modelfile)
        instancedir = os.path.join(basedir, dataname, instancedir)
        for (dirpath, dirnames, filenames) in os.walk(instancedir):
            for file in sorted(filenames):
                instancefile = os.path.join(dirpath, file)
                print("--- Processing:", modelfile, instancefile, "...")
                idx = len(df.index)
                df.loc[idx] = conductExperiment([modelfile, instancefile])
                # range(8, 13):
                for col in ["$t_{clingo}$", "$t_{int}$", "$t_{alg}$", "$t_{asg}^{(d)}$", "$t_{asg}^{(s)}$"]:
                    if df.iloc[idx, df.columns.get_loc(col)] != "NA":
                        df.iloc[idx, df.columns.get_loc(col)] *= 1000
                df.iloc[idx, df.columns.get_loc(
                    'Data instance')] = os.path.join(dataname, file)

                df.to_csv("results.csv", float_format='%.12f')

                df_copy = df.copy()
                s = df_copy.style.format(subset=['Data instance'], escape="latex").\
                    format(subset=["$\overline{l}$", "$t_{clingo}$", "$t_{int}$", "$t_{alg}$", "$t_{asg}^{(d)}$", "$t_{asg}^{(s)}$"], thousands=",", precision=3).\
                    format(subset=['Data instance'], escape="latex").format(subset=[
                        '$m$', '$n$', '$k$', '$m\'$', '$n\'$', '$\eta_z$'], thousands=",").hide_index()

                s.to_latex(
                    "table.tex",
                    column_format="l|rrrrrrr|rr|r|rr",
                    hrules=True,
                    multirow_align="t", multicol_align="r",
                )


if __name__ == "__main__":
    main()
