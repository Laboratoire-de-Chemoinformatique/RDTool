/*
 * Copyright (C) 2003-2018 Syed Asad Rahman <asad @ ebi.ac.uk>.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301  USA
 */
package uk.ac.ebi.reactionblast.mapping.helper;

import static org.openscience.cdk.CDKConstants.MAPPED;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IMapping;
import org.openscience.cdk.interfaces.IReaction;
import ru.spb.arcadia.jnj.aam.io.Utils;
import uk.ac.ebi.reactionblast.tools.BasicDebugger;

/**
 * @Author: Syed Asad Rahman <asad @ ebi.ac.uk>
 * @Date: 2009/06/3
 * @Revision: 1.10
 */
public class MappingHandler extends BasicDebugger {

    /**
     *
     * @param MappedReaction
     */
    public static void cleanMapping(IReaction MappedReaction) {
        int count = MappedReaction.getMappingCount();
        for (int i = count; i > 0; i--) {
            MappedReaction.removeMapping(i);
        }

        for (int eMol = 0; eMol < MappedReaction.getReactantCount(); eMol++) {
            IAtomContainer eMolecule = MappedReaction.getReactants().getAtomContainer(eMol);
            for (int eAtom = 0; eAtom < eMolecule.getAtomCount(); eAtom++) {

//                IAtom atomE = eMolecule.getAtom(eAtom);
//                System.out.println("Atom: " + atomE.getSymbol());
                IAtom atomEMap = MappedReaction.getReactants().getAtomContainer(eMol).getAtom(eAtom);
//                System.out.println("AtomCopy: " + atomEMap.getSymbol());
                String atomLabel = Integer.toString(-1);
                atomEMap.setFlag(MAPPED, false);
                atomEMap.setID(atomLabel);
            }
        }

        for (int pMol = 0; pMol < MappedReaction.getProductCount(); pMol++) {
            IAtomContainer pMolecule = MappedReaction.getProducts().getAtomContainer(pMol);
            for (int pAtom = 0; pAtom < pMolecule.getAtomCount(); pAtom++) {

                IAtom atomPMap = MappedReaction.getProducts().getAtomContainer(pMol).getAtom(pAtom);
                String atomLabel = Integer.toString(-1);
                atomPMap.setFlag(MAPPED, false);
                atomPMap.setID(atomLabel);
            }
        }
    }

    /**
     *
     * @param expLabReaction
     * @param MappedReaction
     * @param counter
     * @return
     */
    protected static synchronized int setMappingFlags(IReaction expLabReaction, IReaction MappedReaction, int counter) {
        IAtomContainerSet expEductSet = expLabReaction.getReactants();
        IAtomContainerSet expProductSet = expLabReaction.getProducts();

        for (IMapping map : expLabReaction.mappings()) {
            IAtom I_Atom = (IAtom) map.getChemObject(0);
            IAtom J_Atom = (IAtom) map.getChemObject(1);
            if (I_Atom != null && J_Atom != null) {

                /*
                *******************************
                * Mapping the Reactants ******************************
                 */
                boolean eFlag = false;
                IAtom firstAtom = null;
                IAtom secondAtom = null;
                for (int eMol = 0; eMol < expEductSet.getAtomContainerCount(); eMol++) {
                    IAtomContainer eMolecule = expEductSet.getAtomContainer(eMol);
                    for (int eAtom = 0; eAtom < eMolecule.getAtomCount(); eAtom++) {
                        if (I_Atom.getID().trim().equalsIgnoreCase(eMolecule.getAtom(eAtom).getID().trim())) {
                            String atomLabel = Integer.toString(counter);
                            firstAtom = MappedReaction.getReactants().getAtomContainer(eMol).getAtom(eAtom);
                            firstAtom.setID(atomLabel);
                            firstAtom.setFlag(MAPPED, true);
                            eFlag = true;
                            break;
                        }
                    }

                    if (eFlag) {
                        break;
                    }

                }
                /*
                *******************************
                * Mapping the Products ******************************
                 */
                boolean pFlag = false;
                for (int pMol = 0; pMol < expProductSet.getAtomContainerCount(); pMol++) {
                    IAtomContainer pMolecule = expProductSet.getAtomContainer(pMol);
                    for (int pAtom = 0; pAtom < pMolecule.getAtomCount(); pAtom++) {

                        if (J_Atom.getID().trim().equalsIgnoreCase(pMolecule.getAtom(pAtom).getID().trim())) {
//                            System.out.println("Hi Matched product");
//                            System.out.println("ID:" + J_Atom.getID().trim());

                            String atomLabel = Integer.toString(counter);

                            secondAtom = MappedReaction.getProducts().getAtomContainer(pMol).getAtom(pAtom);
                            secondAtom.setID(atomLabel);
                            secondAtom.setFlag(MAPPED, true);
                            IMapping mappingObject = MappedReaction.getBuilder().newInstance(IMapping.class, firstAtom, secondAtom);
                            MappedReaction.addMapping(mappingObject);
                            counter++;
                            pFlag = true;
                            break;
                        }
                    }

                    if (pFlag) {
                        break;
                    }
                }
            }
        }
        return counter;
    }

    /**
     *
     * @param MappedReaction
     * @param ReactionWithUniqueSTOICHIOMETRY
     * @param coreMappedReaction
     * @param counter
     * @return
     */
    protected static int setMappingFlags(IReaction MappedReaction, IReaction ReactionWithUniqueSTOICHIOMETRY, IReaction coreMappedReaction, int counter, boolean log) {
        IAtomContainerSet expEductSet = ReactionWithUniqueSTOICHIOMETRY.getReactants();
        IAtomContainerSet expProductSet = ReactionWithUniqueSTOICHIOMETRY.getProducts();
        for (IMapping map : coreMappedReaction.mappings()) {
            IAtom I_Atom = (IAtom) map.getChemObject(0);
            IAtom J_Atom = (IAtom) map.getChemObject(1);
            if (I_Atom != null && J_Atom != null) {
                // reactants
                boolean eFlag = false;
                for (int eMol = 0; eMol < expEductSet.getAtomContainerCount(); eMol++) {
                    IAtomContainer eMolecule = expEductSet.getAtomContainer(eMol);
                    for (int eAtom = 0; eAtom < eMolecule.getAtomCount(); eAtom++) {
                        String mapping_atom_id = I_Atom.getID().trim();
                        String other_atom_id = eMolecule.getAtom(eAtom).getID().trim();
                        if (mapping_atom_id.equalsIgnoreCase(other_atom_id)) {
                            String atomLabel = Integer.toString(counter);
                            IAtomContainer target_ac = MappedReaction.getReactants().getAtomContainer(eMol);
                            IAtom sourceAtom = eMolecule.getAtom(eAtom);
                            IAtom firstAtom = target_ac.getAtom(eAtom);
                            firstAtom.setID(atomLabel);
                            firstAtom.setFlag(MAPPED, true);
                            eFlag = true;
                            if (false) {
                                System.out.println("Update in reactants:");
                                System.out.println("\tMapping atom: " + I_Atom.getSymbol() + ", " + I_Atom.getID());
                                System.out.println("\tSource atom: " + sourceAtom.getSymbol() + ", " + sourceAtom.getID());
                                System.out.println("\tTarget atom: " + firstAtom.getSymbol() + ", " + firstAtom.getID());
                            }
                            break;
                        }
                    }
                    if (eFlag) break;
                }
                // products
                boolean pFlag = false;
                for (int pMol = 0; pMol < expProductSet.getAtomContainerCount(); pMol++) {
                    IAtomContainer pMolecule = expProductSet.getAtomContainer(pMol);
                    for (int pAtom = 0; pAtom < pMolecule.getAtomCount(); pAtom++)
                        if (J_Atom.getID().trim().equalsIgnoreCase(pMolecule.getAtom(pAtom).getID().trim())) {
                            String atomLabel = Integer.toString(counter);
                            IAtom sourceAtom = pMolecule.getAtom(pAtom);
                            IAtom secondAtom = MappedReaction.getProducts().getAtomContainer(pMol).getAtom(pAtom);
                            secondAtom.setID(atomLabel);
                            secondAtom.setFlag(MAPPED, true);
                            counter++;
                            pFlag = true;
                            if (false) {
                                System.out.println("Update in products:");
                                System.out.println("\tMapping atom: " + J_Atom.getSymbol() + ", " + J_Atom.getID());
                                System.out.println("\tSource atom: " + sourceAtom.getSymbol() + ", " + sourceAtom.getID());
                                System.out.println("\tTarget atom: " + secondAtom.getSymbol() + ", " + secondAtom.getID());
                            }
                            break;
                        }
                    if (pFlag) break;
                }
            }
        }
        return counter;
    }

    protected static int emulateSetmappingFlags(IReaction mapped_reaction, IReaction ReactionWithUniqueSTOICHIOMETRY, IReaction core_mapped_reaction, int counter) {
        for (IAtomContainer mol : core_mapped_reaction.getReactants().atomContainers())
            Utils.show_bonds("INITIAL", mol);
        for (IAtomContainer mol : core_mapped_reaction.getProducts().atomContainers())
            Utils.show_bonds("INITIAL", mol);
        IAtomContainerSet educt_set = ReactionWithUniqueSTOICHIOMETRY.getReactants();
        IAtomContainerSet product_set = ReactionWithUniqueSTOICHIOMETRY.getProducts();
        System.out.println("Start writing detailed remapping");
        for (IMapping mapping : core_mapped_reaction.mappings()) {
            String re_id = mapping.getChemObject(0).getID().trim();
            String pr_id = mapping.getChemObject(1).getID().trim();
            System.out.println("\t" + re_id + " -> " + pr_id);
            for (int i = 0; i < educt_set.getAtomContainerCount(); i++) {
                IAtomContainer mol = educt_set.getAtomContainer(i);
                int ind = find_atom_index(mol, re_id);
                if (ind != -1) {
                    String label = Integer.toString(counter);
                    IAtom source_atom = mol.getAtom(ind);
                    IAtom target_atom = mapped_reaction.getReactants().getAtomContainer(i).getAtom(ind);
                    target_atom.setID(label);
                    target_atom.setFlag(MAPPED, true);
//                    System.out.println("\t[reagent] source: " + format_atom(source_atom) + "; target: " + format_atom(target_atom));
                    break;
                }
            }
            for (int i = 0; i < product_set.getAtomContainerCount(); i++) {
                IAtomContainer mol = product_set.getAtomContainer(i);
                int ind = find_atom_index(mol, pr_id);
                if (ind != -1) {
                    String label = Integer.toString(counter++);
                    IAtom source_atom = mol.getAtom(ind);
                    IAtom target_atom = mapped_reaction.getProducts().getAtomContainer(i).getAtom(ind);
                    target_atom.setID(label);
                    target_atom.setFlag(MAPPED, true);
//                    System.out.println("\t[product] source: " + format_atom(source_atom) + "; target: " + format_atom(target_atom));
                    break;
                }
            }
        }
//        for (IAtomContainer mol : mapped_reaction.getReactants().atomContainers())
//            Utils.show_bonds("FINAL", mol);
//        for (IAtomContainer mol : mapped_reaction.getProducts().atomContainers())
//            Utils.show_bonds("FINAL", mol);
        return counter;
    }

    private static int find_atom_index(IAtomContainer mol, String other_id) {
        for (int i = 0; i < mol.getAtomCount(); i++) {
            String id = mol.getAtom(i).getID().trim();
            if (id.equalsIgnoreCase(other_id)) return i;
        }
        return -1;
    }

    private static String format_atom(IAtom atom) {
        return atom.getSymbol() + " (" + atom.getID() + ")";
    }
}
