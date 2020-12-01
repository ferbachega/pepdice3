def save_XYZ_to_file(molecule, filename):
    with open(filename, "a") as output_file:

        text = ''
        n = 0

        for residue_i in molecule.residues:
            for atom_i in residue_i.atoms:
                text += ("{}\t{}\n".format(atom_i.name,
                                           "\t".join([str(round(c, 2)) for c in atom_i.pos])))
                n = n + 1

        output_file.write(str(n) + "\n\n")
        output_file.write(text)


def load_XYZ_to_system (molecule, filename):
    """ Function doc """
    pass
