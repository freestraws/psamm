#!/usr/bin/env python

import sys
import logging

from psamm.formula import Formula, Atom
from psamm.rpair import predict_rpair
from psamm.datasource import kegg

from six import iteritems, integer_types


logger = logging.getLogger(__name__)


def formula_is_only_h(form):
    form_dict = dict(form.items())
    return len(form_dict) == 1 and Atom('H') in form_dict


if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)

    compound_formula = {}
    compound_not_parsed = 0
    with open(sys.argv[1], 'r') as f:
        for entry in kegg.parse_compound_file(f):
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

    logger.info('Compounds: {}, (not parsed: {})'.format(
        len(compound_formula), compound_not_parsed))

    skipped_non_integer, skipped_unknown_formula, skipped_unbalanced = 0, 0, 0
    matched, mismatched, missing = 0, 0, 0
    reaction_count = 0
    with open(sys.argv[2], 'r') as f:
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

            logger.debug(entry.id)

            transfer, balance = predict_rpair(reaction, compound_formula)
            if len(balance) > 0:
                if any(not formula_is_only_h(form) for _, form
                       in iteritems(balance)):
                    logger.info('Skipping unbalanced equation!')
                    logger.info('Balance: {}'.format(balance))
                    skipped_unbalanced += 1
                    continue
                else:
                    logger.info('Ignoring unbalanced H in equation...')

            predicted = {}
            for (inst1, inst2), form in iteritems(transfer):
                c1, n1 = inst1
                c2, n2 = inst2
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
                _, transfer_names, rp_type = rpair
                actual[transfer_names] = rp_type

            entry_matched, entry_mismatched, entry_missing = 0, 0, 0
            for pair, rp_type in iteritems(predicted):
                if rp_type == 'main':
                    if pair in actual:
                        entry_matched += 1
                    else:
                        entry_mismatched += 1
            for pair, rp_type in iteritems(actual):
                if rp_type == 'main':
                    if pair not in predicted:
                        entry_missing += 1

            if entry_mismatched > 0 or entry_missing > 0:
                logger.info('===== {} ====='.format(entry.id))
                logger.info(entry.name)
                logger.info(entry.definition)
                logger.info(reaction)

                logger.info('Predicted: {}'.format(predicted))
                logger.info('Actual: {}'.format(actual))

            matched += entry_matched
            mismatched += entry_mismatched
            missing += entry_missing

            reaction_count += 1

    logger.info('Reactions: {}'.format(reaction_count))
    logger.info('Skipped: {} (non-integer), {} (unknown formula),'
                ' {} (unbalanced)'.format(
        skipped_non_integer, skipped_unknown_formula, skipped_unbalanced))
    logger.info('Matched: {}, mismatched: {}, missing: {}'.format(
        matched, mismatched, missing))
