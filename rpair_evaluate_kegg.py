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
                        continue
                    compound_formula[entry.id] = form

    if args.method == 'mapmaker':
        solver = generic.Solver(integer=True)

    skipped_non_integer = 0
    skipped_unknown_formula = 0
    skipped_unbalanced = 0
    skipped_empty_rpairs = 0
    skipped_invalid_rpairs = 0

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
                logger.info('Non-integer stoichiometry! Skipping reaction...')
                skipped_non_integer += 1
                continue

            if any(compound.name not in compound_formula for compound, _
                    in reaction.compounds):
                logger.info('Unknown formula for reaction compounds!'
                            ' Skipping reaction...')
                skipped_unknown_formula += 1
                continue

            entry_rpairs = list(entry.rpairs)
            if len(entry_rpairs) == 0:
                logger.info('Skipping reaction with no reaction pairs!')
                skipped_empty_rpairs += 1
                continue

            rpair_compounds = set()
            invalid_compounds = set()
            for _, transfer_names, _ in entry_rpairs:
                rpair_compounds.update(transfer_names)
            for c, _ in reaction.compounds:
                if c.name not in rpair_compounds:
                    invalid_compounds.add(c)

            if len(invalid_compounds) > 0:
                logger.info(
                    'Skipping reaction where a carbon-containing has no'
                    ' reaction pairs! {}'.format(invalid_compounds))
                skipped_invalid_rpairs += 1
                continue

            logger.debug(entry.id)

            if args.method == 'greedy':
                instance_transfer, balance = predict_rpair(
                    reaction, compound_formula)
                if len(balance) > 0:
                    if any(not formula_is_only_h(form) for _, form
                           in iteritems(balance)):
                        logger.info('Skipping unbalanced equation!')
                        logger.info('Balance: {}'.format(balance))
                        skipped_unbalanced += 1
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
                    logger.info('Skipping unbalanced equation!')
                    skipped_unbalanced += 1
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
            for rpair in entry.rpairs:
                _, transfer_names, _ = rpair
                if transfer_names not in actual:
                    actual[transfer_names] = 'main'
                    kegg_pairs += 1

            entry_tp, entry_fp, entry_fn = 0, 0, 0
            for pair, rp_type in iteritems(predicted):
                if pair in actual:
                        entry_tp += 1
                else:
                        entry_fp += 1
                        print('{}: {} ({}): {} not in KEGG'.format(
                            entry.id, entry.name, entry.definition,
                            tuple(compound_name[x] for x in pair)))
            for pair, rp_type in iteritems(actual):
                if pair not in actual:
                    entry_fn += 1
                    print('{}: {} ({}): {} not predicted'.format(
                        entry.id, entry.name, entry.definition,
                        tuple(compound_name[x] for x in pair)))

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

            tp += entry_tp
            fp += entry_fp
            fn += entry_fn

    logger.info('Compounds: {}, (not parsed: {})'.format(
        len(compound_formula), compound_not_parsed))
    logger.info('Reactions: {} (matched: {})'.format(
        reaction_count, reaction_matched))
    logger.info('Total KEGG reaction pairs: {}'.format(kegg_pairs))
    logger.info(
        'Skipped: {} (non-integer), {} (unknown formula), {} (unbalanced),'
        ' {} (empty rpairs), {} (invalid pairs)'.format(
            skipped_non_integer, skipped_unknown_formula, skipped_unbalanced,
            skipped_empty_rpairs, skipped_invalid_rpairs))
    logger.info('TP: {}, FP: {}, FN: {}'.format(tp, fp, fn))
