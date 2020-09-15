package ru.spb.arcadia.jnj.aam.check;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.smsd.AtomAtomMapping;
import org.openscience.smsd.Isomorphism;
import org.openscience.smsd.interfaces.Algorithm;
import ru.spb.arcadia.jnj.aam.io.RDFReader;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.*;
import java.util.stream.Collectors;

public class CheckIsomorphism {
    public static void main(String[] args) throws Exception {
        String fn = "inputs/externalExperts_135.rdf";
        IReaction reaction = new RDFReader(fn, "Reaction_ID").read_reactions()[0];
        IAtomContainer reactant = get_biggest(reaction.getReactants());
        IAtomContainer product = get_biggest(reaction.getProducts());
        System.out.println("# of atoms in reactant: " + reactant.getAtomCount());
        System.out.println("# of atoms in product: " + product.getAtomCount());
        int count = 0;
        Algorithm algorithm = Algorithm.CDKMCS;
        boolean bond_type_flag = true;
        boolean match_atom_type = true;
        int nt = 5;
        String out_fn = "algorithm_" + algorithm + "_bond_" + bond_type_flag + "_atom_" + match_atom_type + ".txt";
        BufferedWriter bw = new BufferedWriter(new FileWriter(out_fn));
        for (int t = 0; t < nt; t++) {
            Isomorphism iso = new Isomorphism(reactant, product, algorithm, bond_type_flag, false, match_atom_type);
            iso.setChemFilters(true, true, true);
            bw.write("# of mappings: " + iso.getMappingCount());
            bw.newLine();
            int num_correct = 0;
            for (int i = 0; i < iso.getMappingCount(); i++) {
                AtomAtomMapping aam = iso.getAllAtomMapping().get(i);
                bw.write("mapping #" + i);
                bw.newLine();
                Map<String, String> id_pairs = new HashMap<>();
                for (Map.Entry<IAtom, IAtom> pair : aam.getMappingsByAtoms().entrySet()) {
                    String id1 = pair.getKey().getID();
                    String id2 = pair.getValue().getID();
                    id_pairs.put(id1, id2);
                    bw.write("\t" + id1 + " : " + id2);
                    bw.newLine();
                }
                if (id_pairs.get("2").equals("2") && id_pairs.get("6").equals("6"))
                    num_correct++;
                else if (id_pairs.get("6").equals("2") && id_pairs.get("2").equals("6"))
                    num_correct++;
                else if (id_pairs.get("6").equals("7") && id_pairs.get("2").equals("3"))
                    num_correct++;
                else if (id_pairs.get("2").equals("7") && id_pairs.get("6").equals("3"))
                    num_correct++;
            }
            bw.write("----------");
            bw.newLine();
            if (num_correct == iso.getMappingCount())
                count++;
        }
        bw.close();
        System.out.println("Algorithm: " + algorithm);
        System.out.println("\tbond type flag: " + bond_type_flag);
        System.out.println("\tatom match type: " + match_atom_type);
        System.out.println("\tcount: " + count + "/" + nt);
    }

    private static IAtomContainer get_biggest(IAtomContainerSet molset) {
        IAtomContainer res = null;
        int count = 0;
        for (IAtomContainer mol : molset.atomContainers())
            if (mol.getAtomCount() > count) {
                count = mol.getAtomCount();
                res = mol;
            }
        return res;
    }
}
