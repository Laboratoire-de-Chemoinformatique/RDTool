package ru.spb.arcadia.jnj.aam;

import org.openscience.cdk.interfaces.*;
import org.openscience.smsd.tools.ExtAtomContainerManipulator;
import static org.openscience.cdk.aromaticity.Kekulization.kekulize;

import java.util.HashMap;
import java.util.Map;

import static org.openscience.cdk.CDKConstants.MAPPED;
import static org.openscience.smsd.tools.ExtAtomContainerManipulator.cloneWithIDs;

public class MappingFixer {
    private final IReaction mapped_reaction;

    public MappingFixer(IReaction mapped_reaction) throws Exception {
        this.mapped_reaction = mapped_reaction;
        for (IAtomContainer ac : mapped_reaction.getReactants().atomContainers()) {
            ExtAtomContainerManipulator.setNullHCountToZero(ac);
        }
        for (IAtomContainer ac : mapped_reaction.getProducts().atomContainers()) {
            ExtAtomContainerManipulator.setNullHCountToZero(ac);
        }
        for (IAtomContainer ac : mapped_reaction.getReactants().atomContainers()) kekulize(ac);
        for (IAtomContainer ac : mapped_reaction.getProducts().atomContainers()) kekulize(ac);
//        System.out.println("MappingFixer: # of mappings = " + mapped_reaction.getMappingCount());
        assign_ids();
    }

    public IReaction get_result() throws Exception {
        return clone_mapped();
    }

    private void assign_ids() throws Exception {
        // drop
        for (IAtomContainer ac : mapped_reaction.getReactants().atomContainers())
            for (IAtom a : ac.atoms())
                a.setID(Integer.toString(-1));
        for (IAtomContainer ac : mapped_reaction.getProducts().atomContainers())
            for (IAtom a : ac.atoms())
                a.setID(Integer.toString(-1));
        // mapped atoms
        int counter = 1;
        for (IMapping mapping : mapped_reaction.mappings()) {
            IAtom fst_atom = (IAtom) mapping.getChemObject(0);
            IAtom snd_atom = (IAtom) mapping.getChemObject(1);
            fst_atom.setID(Integer.toString(counter));
            fst_atom.setFlag(MAPPED, true);
            snd_atom.setID(Integer.toString(counter++));
            snd_atom.setFlag(MAPPED, true);
        }
    }

    private IReaction clone_mapped() throws Exception {
        // insert reactants and products
        IReaction reaction = mapped_reaction.getBuilder().newInstance(IReaction.class);
        for (int i = 0; i < mapped_reaction.getReactantCount(); i++) {
            IAtomContainer mol = mapped_reaction.getReactants().getAtomContainer(i);
            mol = cloneWithIDs(mol);
            reaction.addReactant(mol);
        }
        for (int i = 0; i < mapped_reaction.getProductCount(); i++) {
            IAtomContainer mol = mapped_reaction.getProducts().getAtomContainer(i);
            mol = cloneWithIDs(mol);
            reaction.addProduct(mol);
        }
        reaction.setID(mapped_reaction.getID());
        reaction.setDirection(mapped_reaction.getDirection());
        // add mapping
        Map<IAtom, IAtom> mappings = new HashMap<>();
        for (IAtomContainer ac : reaction.getReactants().atomContainers())
            for (IAtom atom : ac.atoms()) {
                IAtom atom2 = get_container_atom_by_id(reaction.getProducts(), atom.getID());
                if (atom2 != null)
                    mappings.put(atom, atom2);
            }
        mappings.keySet().stream().filter((key) -> (key != null && mappings.get(key) != null)).map(
                (key) -> reaction.getBuilder().newInstance(IMapping.class, key, mappings.get(key))).forEachOrdered(
                (mappingObject) -> {reaction.addMapping(mappingObject);});
        return reaction;
    }

    private static IAtom get_container_atom_by_id(IAtomContainerSet products, String mappingID) {
        for (IAtomContainer ac : products.atomContainers()) {
            for (IAtom atom : ac.atoms()) {
                if (atom.getID().equals(mappingID)) {
                    return atom;
                }
            }
        }
        return null;
    }
}
