#!/usr/bin/env python

import argparse
import os
import logging
from collections import Counter
from itertools import product

from psamm.formula import Formula, Atom
from psamm.rpair import predict_rpair
from psamm.datasource import kegg
from psamm.lpsolver import generic
from psamm import rpair_milp

from six import iteritems, integer_types


logger = logging.getLogger(__name__)


def formula_is_only_h(form):
    form_dict = dict(form.items())
    return len(form_dict) == 1 and Atom('H') in form_dict


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)

    parser = argparse.ArgumentParser(
        description='Evaluate compound pair prediction on KEGG')
    parser.add_argument('kegg', help='KEGG database directory',
                        type=str)
    parser.add_argument('--method', help='Prediction method',
                        default='greedy', choices=['greedy', 'mapmaker'])
    args = parser.parse_args()

    kegg_dir = args.kegg

    compound_formula = {}
    compound_name = {}
    compound_not_parsed = 0
    compound_is_variable = 0
    with open(os.path.join(kegg_dir, 'compound', 'compound'), 'r') as f:
        for entry in kegg.parse_compound_file(f):
            compound_name[entry.id] = entry.name
            if entry.formula is not None:
                try:
                    form = Formula.parse(entry.formula)
                except ValueError:
                    logger.debug('Unable to parse formula of {}: {}'.format(
                        entry.id, entry.formula))
                    compound_not_parsed += 1
                else:
                    if form.is_variable():
                        compound_is_variable += 1
                        continue
                    compound_formula[entry.id] = form

    if args.method == 'mapmaker':
        solver = generic.Solver(integer=True)

    skipped_non_integer = 0
    skipped_unknown_formula = 0
    skipped_empty_rpairs = 0
    skipped_invalid_rpairs = 0

    unbalanced_reactions = set()
    trivial_reaction_count = 0
    trivial_incorrect = set()

    tp, fp, fn = 0, 0, 0
    reaction_count = 0
    reaction_matched = 0
    kegg_pairs = 0
    with open(os.path.join(kegg_dir, 'reaction', 'reaction'), 'r') as f:
        for entry in kegg.parse_reaction_file(f):
            if entry.equation is None:
                continue

            reaction = kegg.parse_reaction(entry.equation)

            if any(not isinstance(value, integer_types) for _, value
                    in reaction.compounds):
                logger.info(
                    '{}: Non-integer stoichiometry!'
                    ' Skipping reaction...'.format(entry.id))
                skipped_non_integer += 1
                continue

            if any(compound.name not in compound_formula for compound, _
                    in reaction.compounds):
                logger.info(
                    '{}: Unknown formula for reaction compounds!'
                    ' Skipping reaction...'.format(entry.id))
                skipped_unknown_formula += 1
                continue

            entry_rpairs = list(entry.rpairs)
            if len(entry_rpairs) == 0:
                logger.info(
                    '{}: Skipping reaction with no reaction pairs!'.format(
                        entry.id))
                skipped_empty_rpairs += 1
                continue

            rpair_compounds = set()
            invalid_compounds = set()
            for _, transfer_names, _ in entry_rpairs:
                rpair_compounds.update(transfer_names)
            for c, _ in reaction.compounds:
                if (c.name not in rpair_compounds and
                        Atom('C') in compound_formula[c.name]):
                    invalid_compounds.add(c)

            if len(invalid_compounds) > 0:
                logger.info(
                    '{}: Skipping reaction where reaction pair carbon-transfers'
                    ' are invalid! {}'.format(entry.id, invalid_compounds))
                skipped_invalid_rpairs += 1
                continue

            # Count carbons on either side
            c_compounds_left, c_compounds_right = 0, 0
            for c, v in reaction.compounds:
                if Atom('C') in compound_formula[c.name]:
                    if v < 0:
                        c_compounds_left += 1
                    else:
                        c_compounds_right += 1

            logger.debug(entry.id)

            if args.method == 'greedy':
                instance_transfer, balance = predict_rpair(
                    reaction, compound_formula)
                if len(balance) > 0:
                    if any(not formula_is_only_h(form) for _, form
                           in iteritems(balance)):
                        logger.info('{}: Skipping unbalanced equation!'.format(
                            entry.id))
                        logger.info('Balance: {}'.format(balance))
                        unbalanced_reactions.add(entry.id)
                        continue
                    else:
                        logger.info('Ignoring unbalanced H in equation...')

                transfer = {}
                for ((c1, n1), (c2, n2)), form in iteritems(instance_transfer):
                    pair = c1, c2
                    if pair not in transfer:
                        transfer[pair] = form
                    else:
                        transfer[pair] |= form
            else:
                try:
                    transfer = rpair_milp.predict_rpair(
                        reaction, compound_formula, solver)
                except rpair_milp.UnbalancedReactionError:
                    logger.info('{}: Skipping unbalanced equation!'.format(
                        entry.id))
                    unbalanced_reactions.add(entry.id)
                    continue

            predicted = {}
            for (c1, c2), form in iteritems(transfer):
                transfer_names = tuple(sorted([c1.name, c2.name]))
                form_dict = dict(form.items())
                if Atom('C') in form_dict:
                    rp_type = 'main'
                elif not formula_is_only_h(form):
                    rp_type = 'leave'
                else:
                    continue

                if (transfer_names not in predicted or
                        predicted[transfer_names] == 'leave'):
                    predicted[transfer_names] = rp_type

            actual = {}
            for rpair in entry_rpairs:
                _, transfer_names, rp_type = rpair
                if transfer_names not in actual:
                    actual[transfer_names] = rp_type
                    kegg_pairs += 1

            entry_tp, entry_fp, entry_fn = 0, 0, 0
            c_pair_incorrect = 0
            for pair, rp_type in iteritems(predicted):
                if pair in actual:
                    entry_tp += 1
                else:
                    entry_fp += 1
                    print('{}: {} ({}): {} not in KEGG'.format(
                        entry.id, entry.name, entry.definition,
                        tuple(compound_name[x] for x in pair)))

                    if all(Atom('C') in compound_formula[c] for c in pair):
                        c_pair_incorrect += 1
            for pair, rp_type in iteritems(actual):
                if pair not in predicted:
                    entry_fn += 1
                    print('{}: {} ({}): {} not predicted'.format(
                        entry.id, entry.name, entry.definition,
                        tuple(compound_name[x] for x in pair)))

                    if all(Atom('C') in compound_formula[c] for c in pair):
                        c_pair_incorrect += 1

            reaction_count += 1
            false = entry_fp + entry_fn
            if false > 0:
                logger.info('===== {} ====='.format(entry.id))
                logger.info(entry.name)
                logger.info(entry.definition)
                logger.info(reaction)

                logger.info('Predicted: {}'.format(predicted))
                logger.info('Actual: {}'.format(actual))
                logger.info('TP: {}, FP: {}, FN: {}'.format(
                    entry_tp, entry_fp, entry_fn))
            else:
                reaction_matched += 1

            if c_compounds_left <= 1 or c_compounds_right <= 1:
                trivial_reaction_count += 1
                if c_pair_incorrect > 0:
                    trivial_incorrect.add(entry.id)

            tp += entry_tp
            fp += entry_fp
            fn += entry_fn

    logger.info(
        'Compounds: {}, (skipped: unable to parse: {}; variable: {})'.format(
            len(compound_formula), compound_not_parsed, compound_is_variable))
    logger.info(
        'Reactions: {}, (skipped: non-integer stoichiometry: {};'
        ' unknown formula: {}; unbalanced reaction: {};'
        ' missing rpair annotation: {}; invalid rpairs: {})'.format(
            reaction_count, skipped_non_integer, skipped_unknown_formula,
            len(unbalanced_reactions), skipped_empty_rpairs,
            skipped_invalid_rpairs))
    logger.info('Reactions perfectly matched: {}'.format(reaction_matched))
    logger.info('Trivial reactions: {} (incorrect: {})'.format(
        trivial_reaction_count, len(trivial_incorrect)))
    logger.info('Total KEGG reaction pairs: {}'.format(kegg_pairs))
    logger.info('TP: {}, FP: {}, FN: {}'.format(tp, fp, fn))
