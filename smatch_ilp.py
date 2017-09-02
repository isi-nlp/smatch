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
        insts1, rels1 = self.arg1.get_triples2()
        insts2, rels2 = self.arg2.get_triples2()
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

        # FIXME: double check :: ignoring TOP relations
        '''
        rels1 = [t for t in rels1 if t[0] != 'TOP']
        rels2 = [t for t in rels2 if t[0] != 'TOP']
        '''

        self.arg1size = len(rels1)
        self.arg2size = len(rels2)

        trpl_choices = set()
        trpl_var_consts = {}
        for name1, var11, var12 in rels1:
            id1 = "%s:%s:%s" % (name1, var11, var12)
            for name2, var21, var22 in rels2:
                id2 = "%s:%s:%s" % (name2, var21, var22)

                # triple name matches && first argument to triples can be matched
                if name1 == name2 and (var11, var21) in var_choices:
                    # second argument to triple can also be matched OR
                    if (var12, var22) in var_choices or (
                        # they are the same concepts
                            var12 not in self.arg1vars and
                            var22 not in self.arg2vars and
                            var12 == var22):
                        trpl_choices.add((id1, id2))
                        # constrains between variables and triples
                        trpl_var_consts[id1, id2] = [(var11, var21)]
                        if (var12, var22) in var_choices:
                            trpl_var_consts[id1, id2].append((var11, var21))

        # Add variables to ILP model
        model = GRBModel('Smatch ILP')
        if log.getLogger().getEffectiveLevel() >= log.INFO:
            model.Params.OutputFlag = 0         # disable output
        log.info("Number of possible variable Matches %s" % len(var_choices))
        log.info("Number of possible triple matches %s" % len(trpl_var_consts))

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

    def solve(self):
        """
        :return: smatch score
        """
        self.model.optimize()
        assert self.model.status == GRB.OPTIMAL
        num_trpl_match = self.model.objVal
        precision = num_trpl_match / min(self.arg1size, self.arg2size)
        recall = num_trpl_match / max(self.arg1size, self.arg2size)
        denominator = (precision + recall)
        if abs(0.0 - denominator) < 1e-6:
            smatch_score = 0.0
        else:
            smatch_score = 2 * precision * recall / denominator

        if log.getLogger().getEffectiveLevel() <= log.INFO:
            variables = self.model.getAttr('x', self.vars)
            triples = self.model.getAttr('x', self.trpls)
            var_matchings = [pair for (pair, assignment) in variables.items() if assignment]
            triple_matchings = [pair for (pair, assignment) in triples.items() if assignment]
            log.info("Variables:")
            log.info('\n'.join('  %s -> %s' % (v1, v2) for v1, v2 in var_matchings))
            log.info("Triples:")
            log.info('\n'.join('  %s -> %s' % (v1, v2) for v1, v2 in triple_matchings))
        return smatch_score

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Smatch ILP')
    parser.add_argument('amrfile', nargs=2, help='Path to file having AMR')
    parser.add_argument('-v', help='Verbose (log level = INFO)', action='store_true')
    parser.add_argument('-vv', help='Verbose (log level = DEBUG)', action='store_true')
    args = vars(parser.parse_args())
    if args['v']:
        log.getLogger().setLevel(level=log.INFO)

    if args['vv']:
        log.getLogger().setLevel(level=log.DEBUG)

    file1, file2 = args['amrfile']
    for amr1, amr2 in AMR.read_amrs(file1, file2):
        smatch = SmatchILP(amr1, amr2)
        score = smatch.solve()
        print('%.3f' % score)
