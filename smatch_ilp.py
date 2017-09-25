#!/usr/bin/env python

"""
This python module provides an implementation of SMATCH using ILP powered by Gurobi

author: Thamme Gowda (tg@isi.edu)
date: August 30, 2017
"""

from amr import AMR
from gurobipy import Model as GRBModel
from gurobipy import GRB
import argparse

import logging as log
log.basicConfig(level=log.WARNING)


class SmatchILP(object):

    """
    This class provides an implementation of Smatch using 0-1 Integer Linear programming powered by Gurobi
    """
    def __init__(self, amr1, amr2):
        """
        Two AMR graphs
        :param amr1:  the first AMR
        :param amr2: the second AMR
        """
        self.arg1 = amr1
        self.arg2 = amr2
        self.arg1vars = {}
        self.arg2vars = {}
        # initialized later in build_model
        self.vars, self.trpls = None, None
        self.model = None
        self.arg1size, self.arg2size = None, None
        self.build_model()

    def build_model(self):
        """
        Constructs GUROBI ILP model
        :return: None
        """
        insts1, attrs1, rels1 = self.arg1.get_triples()
        insts2, attrs2, rels2 = self.arg2.get_triples()
        for items, shld_norm in [(insts1, True), (insts2, True), (attrs1, True),
                                 (attrs2, True), (rels1, False), (rels2, False)]:
            for i in range(len(items)):
                # GUROBI cant handle Unicode so step down to ASCII
                items[i] = [items[i][0].encode('ascii', 'ignore').lower(),
                            items[i][1].encode('ascii', 'ignore'),
                            items[i][2].encode('ascii', 'ignore')]
                # normalize concept names -- instances and attributes
                if shld_norm:
                    items[i][2] = SmatchILP.normalize(items[i][2])

        # Attributes are same as relations
        rels1.extend(attrs1)
        rels2.extend(attrs2)

        log.debug("AMR 1 Instances:\n  %s" % insts1)
        log.debug("AMR 1 Relations:\n  %s" % rels1)
        log.debug("AMR 2 Instances:\n  %s" % insts2)
        log.debug("AMR 2 Relations:\n  %s" % rels2)

        for index, items in [(self.arg1vars, insts1), (self.arg2vars, insts2)]:
            for name, var, concept in items:
                assert name == 'instance'  # relation name is instance ==> variable definition
                assert var not in index  # variable name is unique
                index[var] = concept

        var_choices = set()  # possible variable matches
        for v1 in self.arg1vars.keys():
            for v2 in self.arg2vars.keys():
                var_choices.add((v1, v2))

        # instances are relations too
        rels1.extend(insts1)
        rels2.extend(insts2)

        self.arg1size = len(rels1)
        self.arg2size = len(rels2)

        trpl_choices = set()
        trpl_var_consts = {}
        for name1, var11, var12 in rels1:
            id1 = "%s:%s:%s" % (name1, var11, var12)
            for name2, var21, var22 in rels2:
                possible = 0
                id2 = "%s:%s:%s" % (name2, var21, var22)
                # triple name matches && first argument to triples can be matched
                if name1 == name2 and (var11, var21) in var_choices:
                    # second argument to triple can also be matched OR
                    possible += 1
                    if (var12, var22) in var_choices or (
                            # they are the same concepts
                            # var12 not in self.arg1vars and var22 not in self.arg2vars and
                            var12 == var22):
                        possible += 1
                        trpl_choices.add((id1, id2))
                        # constrains between variables and triples
                        trpl_var_consts[id1, id2] = [(var11, var21)]
                        # if second argument is also variable

                        if (var12, var22) in var_choices:
                            trpl_var_consts[id1, id2].append((var12, var22))
                log.debug('\t %s <--> %s ? %s ' % (id1, id2, possible))

        # Add variables to ILP model
        model = GRBModel('Smatch ILP')
        if log.getLogger().getEffectiveLevel() >= log.INFO:
            model.Params.OutputFlag = 0         # disable output
        log.info("Number of possible variable matches %s" % len(var_choices))
        log.info("Number of possible triple matches %s" % len(trpl_choices))

        self.vars = model.addVars(var_choices, vtype=GRB.BINARY, name="v")
        self.trpls = model.addVars(trpl_choices, vtype=GRB.BINARY, name="t")

        # constraints
        for v1 in self.arg1vars:
            model.addConstr(self.vars.sum(v1, '*') <= 1, name='to max 1 var')
        for v2 in self.arg2vars:
            model.addConstr(self.vars.sum('*', v2) <= 1, name='from max 1 var')

        for trpl_idx, var_idxs in trpl_var_consts.items():
            for var_idx in var_idxs:
                model.addConstr(self.trpls[trpl_idx] <= self.vars[var_idx], name="%s::%s" % (trpl_idx, var_idx))

        # objective
        model.setObjective(self.trpls.sum(), GRB.MAXIMIZE)
        self.model = model

        # stats for how big the problem is
        var_trpl_consts_count = sum(len(x) for x in trpl_var_consts.values())
        num_constr = len(var_choices) + len(trpl_choices) + var_trpl_consts_count
        num_vars = len(var_choices) + len(trpl_choices)
        log.info("ILP SIZE: %d binary variables (%d vars + %d triple vars)" % (num_vars, len(var_choices), len(trpl_choices)))
        log.info("ILP SIZE: %d constraints (%d b/w arg vars and triples)" % (num_constr, var_trpl_consts_count))


    @staticmethod
    def normalize(item):
        """
        lowercase and remove quote signifiers from items that are about to be compared
        """
        item = item.lower().strip().rstrip('_')
        return item

    def solve(self):
        """
        :return: smatch score, triple_match count
        """
        self.model.optimize()
        assert self.model.status == GRB.OPTIMAL

        num_trpl_match = self.model.objVal
        precision = num_trpl_match / self.arg1size
        recall = num_trpl_match / self.arg2size
        smatch_score = SmatchILP.f_mneasure(precision, recall)

        if log.getLogger().getEffectiveLevel() <= log.INFO:
            variables = self.model.getAttr('x', self.vars)
            triples = self.model.getAttr('x', self.trpls)
            var_matchings = [pair for (pair, assignment) in variables.items() if assignment]
            triple_matchings = [pair for (pair, assignment) in triples.items() if assignment]
            log.info("Number of Triple Matched: %d" % num_trpl_match)
            log.info("Triples in first AMR : %d" % self.arg1size)
            log.info("Triples in second AMR : %d" % self.arg2size)
            log.info("Matched Variables:")
            log.info('\n'.join('  %s -> %s' % (v1, v2) for v1, v2 in var_matchings))
            log.info("Matched Triples:")
            log.info('\n'.join('  %s -> %s' % (v1, v2) for v1, v2 in triple_matchings))
        return smatch_score, num_trpl_match

    @staticmethod
    def f_mneasure(prec, recall):
        denominator = (prec + recall)
        if abs(0.0 - denominator) < 1e-6:
            f_score = 0.0
        else:
            f_score = 2 * prec * recall / denominator
        return  f_score

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Smatch ILP')
    parser.add_argument('amrfile', nargs=2, help='Path to file having AMR. The second AMR should be a gold.')
    parser.add_argument('-v', help='Verbose (log level = INFO)', action='store_true')
    parser.add_argument('-vv', help='Verbose (log level = DEBUG)', action='store_true')
    parser.add_argument('-s', '--significant', type=int, default=2, help='significant digits to output (default: 2)')
    parser.add_argument('--ms', action='store_true', default=False,
                        help='Output multiple scores (one AMR pair a score)'
                             'instead of a single document-level smatch score (Default: false)')
    args = vars(parser.parse_args())
    if args['v']:
        log.getLogger().setLevel(level=log.INFO)

    if args['vv']:
        log.getLogger().setLevel(level=log.DEBUG)

    file1, file2 = args['amrfile']
    float_fmt = '%%.%df' % args['significant']
    # Note: instead of computing overage, we are summing all AMRs
    total_match, file1_count, file2_count = 0, 0, 0
    for amr1, amr2 in AMR.read_amrs(file1, file2):
        smatch = SmatchILP(amr1, amr2)
        score, match_count = smatch.solve()
        total_match += match_count
        file1_count += smatch.arg1size
        file2_count += smatch.arg2size
        if args['ms']:
            out = float_fmt % score
            print('F-score: %s' % out)
    if total_match > 0:
        prec = total_match / file1_count
        recall = total_match / file2_count
        smatch_score = SmatchILP.f_mneasure(prec, recall)
        out = float_fmt % smatch_score
        print('\nAggregated F-score: %s' % out)
    else:
        print("No AMRs are found or matched")
