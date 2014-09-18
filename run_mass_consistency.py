#!/usr/bin/env python

import argparse
import csv

from metnet.metabolicmodel import MetabolicDatabase
from metnet.massconsistency import MassConsistencyCheck

if __name__ == '__main__':
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Run FastGapFill on a metabolic model')
    parser.add_argument('--database', required=True, metavar='reactionfile', action='append',
                        type=argparse.FileType('r'), default=[],
                        help='Reaction definition list to use as database')
    parser.add_argument('--compounds', metavar='compoundfile', action='append',
                        type=argparse.FileType('r'), default=[],
                        help='Optional compound information table')
    parser.add_argument('reactionlist', nargs='?', type=argparse.FileType('r'),
                        help='Model definition')
    args = parser.parse_args()

    database = MetabolicDatabase.load_from_files(*args.database)

    # Load model from file if given, otherwise run on full database
    if args.reactionlist:
        model = database.load_model_from_file(args.reactionlist)
    else:
        model = database.get_model(database.reactions)

    # Load compound information
    compounds = {}
    for compound_table in args.compounds:
        compound_table.readline() # Skip header
        for row in csv.reader(compound_table, delimiter='\t'):
            cpdid, names = row[:2]
            synonyms = names.split(',<br>')
            name = synonyms.pop()
            compounds[cpdid] = name

    # Create a set of known mass-inconsistent reactions
    exchange = set()

    for rxnid in model.reaction_set:
        rx = database.get_reaction(rxnid)
        if len(rx.left) == 0 or len(rx.right) == 0:
            exchange.add(rxnid)

    # Create set of compounds allowed to have mass zero
    zeromass = set()
    zeromass.add('cpd11632') # Photon
    zeromass.add('cpd12713') # Electron

    mass_consistency = MassConsistencyCheck()

    print 'Mass consistency on model (database)...'
    epsilon = 1e-5
    compound_iter = mass_consistency.check_compound_consistency(model, exchange, zeromass)

    print 'Compound consistency...'
    good = 0
    total = 0
    for compound, mass in sorted(compound_iter, key=lambda x: x[1], reverse=True):
        compound_id, comp = compound
        if mass >= 1-epsilon or compound_id in zeromass:
            good += 1
        total += 1
        print '{}: {}'.format(compounds.get(compound_id, compound_id), mass)
    print 'Consistent compounds: {}/{}'.format(good, total)

    print 'Is consistent? {}'.format(mass_consistency.is_consistent(model, exchange, zeromass))

    print 'Reaction consistency...'
    reaction_iter, compound_iter = mass_consistency.check_reaction_consistency(model, exchange, zeromass)
    for reaction_id, residual in sorted(reaction_iter, key=lambda x: abs(x[1]), reverse=True):
        if abs(residual) >= epsilon:
            reaction = database.get_reaction(reaction_id).translated_compounds(lambda x: compounds.get(x, x))
            print '{}\t{}\t{}'.format(reaction_id, residual, reaction)