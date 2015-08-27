#!/usr/bin/env python

from __future__ import print_function

from psamm.command import Command, main
from psamm.formula import Formula
from psamm.rpair import predict_rpair

from six import iteritems, integer_types


class ReactionPairPredictCommand(Command):
    """Predict reaction pairs from compound formulas."""

    name = "rpairpredict"
    title = "Predict reaction pairs from compound formulas"

    def run(self):
        compound_formula = {}
        for compound in self._model.parse_compounds():
            f = Formula.parse(compound.formula).flattened()
            compound_formula[compound.id] = f

        for entry in self._model.parse_reactions():
            reaction = entry.equation

            print('===== {} ====='.format(entry.id))
            print(entry.name)
            print(reaction)
            print()

            if any(not isinstance(value, integer_types) for _, value
                    in reaction.compounds):
                print('Non-integer stoichiometry! Skipping reaction...')
                continue

            transfer, balance = predict_rpair(reaction, compound_formula)
            for (inst1, inst2), f in iteritems(transfer):
                print('{}, {}: {}'.format(inst1, inst2, f))
            if len(balance) != 0:
                print('Not balanced! {}'.format(balance))
            print()


if __name__ == '__main__':
    main(ReactionPairPredictCommand)
