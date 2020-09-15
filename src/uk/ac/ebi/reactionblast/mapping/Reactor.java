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
package uk.ac.ebi.reactionblast.mapping;

import java.io.IOException;
import java.io.Serializable;
import static java.lang.System.out;
import java.util.ArrayList;
import static java.util.Collections.synchronizedMap;
import static java.util.Collections.synchronizedSortedMap;
import static java.util.Collections.unmodifiableList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.TreeMap;
import static java.util.logging.Level.SEVERE;
import static org.openscience.cdk.CDKConstants.ATOM_ATOM_MAPPING;
import static org.openscience.cdk.CDKConstants.MAPPED;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IMapping;
import org.openscience.cdk.interfaces.IReaction;
import static org.openscience.cdk.interfaces.IReaction.Direction.BIDIRECTIONAL;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.smiles.SmilesGenerator;
import static org.openscience.cdk.smiles.SmilesGenerator.generic;
import static org.openscience.cdk.smiles.SmilesGenerator.unique;
import static org.openscience.cdk.tools.manipulator.AtomContainerSetManipulator.getTotalFormalCharge;

import ru.spb.arcadia.jnj.aam.Mapping;
import ru.spb.arcadia.jnj.aam.io.Utils;
import uk.ac.ebi.reactionblast.mapping.algorithm.CalculationProcess;
import uk.ac.ebi.reactionblast.mapping.container.MoleculeMoleculeMapping;
import uk.ac.ebi.reactionblast.mapping.helper.AbstractReactor;
import uk.ac.ebi.reactionblast.mapping.interfaces.IMappingAlgorithm;
import static uk.ac.ebi.reactionblast.tools.ExtAtomContainerManipulator.aromatizeMolecule;
import static uk.ac.ebi.reactionblast.tools.ExtAtomContainerManipulator.cloneWithIDs;
import static uk.ac.ebi.reactionblast.tools.ExtAtomContainerManipulator.percieveAtomTypesAndConfigureAtoms;
import static uk.ac.ebi.reactionblast.tools.ExtAtomContainerManipulator.removeHydrogensExceptSingleAndPreserveAtomID;
import static uk.ac.ebi.reactionblast.tools.ExtReactionManipulatorTool.deepClone;
import static java.lang.Integer.parseInt;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.valueOf;
import static java.lang.System.getProperty;
import static java.util.Arrays.sort;
import static java.util.Collections.synchronizedList;

import static org.openscience.cdk.geometry.GeometryTools.has2DCoordinates;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.tools.ILoggingTool;
import static org.openscience.cdk.tools.LoggingToolFactory.createLoggingTool;
import static org.openscience.cdk.tools.manipulator.AtomContainerManipulator.getBondArray;

/**
 *
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class Reactor extends AbstractReactor implements Serializable {

    private static final boolean DEBUG = false;
    private static final long serialVersionUID = 197816786981017L;
    private final static ILoggingTool LOGGER
            = createLoggingTool(Reactor.class);
    private final Map<Integer, Integer> rLabelledAtoms;
    private final Map<Integer, Integer> pLabelledAtoms;
    private final Map<Integer, Integer> inputRankLabelledAtomsReactant;
    private final Map<Integer, Integer> inputRankLabelledAtomsProduct;
    private final Map<Integer, IAtomContainer> educts;
    private final Map<Integer, IAtomContainer> products;
    private final List<IBond> rBonds;
    private final List<IBond> pBonds;
    private final IReaction reactionWithSTOICHIOMETRY;
    private final boolean partialMapping;
    private final IMappingAlgorithm algorithm;
    private MoleculeMoleculeMapping reactionBlastMolMapping;
    private Integer substrateAtomCounter;
    private Integer productAtomCounter;
    private int delta;
    private boolean balanceFlag;
    private IReaction reactionWithUniqueSTOICHIOMETRY;
    static final String NEW_LINE = getProperty("line.separator");

    //~--- constructors -------------------------------------------------------
    /**
     *
     * @param reaction Reaction
     * @param balanceReaction eg. balance hydrogens in the reaction if its not
     * balanced
     * @param partialMapping (without Hydrogens is set True, its faster)
     * @throws Exception
     */
    public Reactor(IReaction reaction,
            boolean partialMapping,
            IMappingAlgorithm algorithm)
            throws Exception {
//        System.LOGGER.debug("In Reaction");
//        SmilesGenerator withAtomClasses = SmilesGenerator.unique().aromatic().withAtomClasses();
//        System.LOGGER.debug("Input reaction to be mapped " + withAtomClasses.createReactionSMILES(reaction));
        this.partialMapping = partialMapping;
        this.algorithm = algorithm;
        this.reactionWithSTOICHIOMETRY = reaction.getBuilder().newInstance(IReaction.class);
        this.reactionWithUniqueSTOICHIOMETRY = reaction.getBuilder().newInstance(IReaction.class);
        this.balanceFlag = true;

        this.inputRankLabelledAtomsReactant = synchronizedMap(new HashMap<>());
        this.inputRankLabelledAtomsProduct = synchronizedMap(new HashMap<>());
        this.rLabelledAtoms = synchronizedMap(new HashMap<>());
        this.pLabelledAtoms = synchronizedMap(new HashMap<>());
        this.rBonds = synchronizedList(new ArrayList<>());
        this.pBonds = synchronizedList(new ArrayList<>());

        this.educts = synchronizedSortedMap(new TreeMap<>());
        this.products = synchronizedSortedMap(new TreeMap<>());

        this.substrateAtomCounter = 1;
        this.productAtomCounter = 1;
        if (DEBUG) {
            out.println("|++++++++++++++++++++++++++++|");
            out.println("|i. Reactor Initialized");
        }
        cleanMapping(reaction);
        if (DEBUG) {
            out.println("|++++++++++++++++++++++++++++|");
            super.printReaction(reaction);
            out.println("|ii. Create Mapping Objects");
        }
        copyReferenceReaction(reaction);
        expandReaction();
        checkReactionBalance();
        if (DEBUG) {
            out.println("|iii. Compute atom-atom Mappings");
        }
        calculateAtomAtomMapping();
        if (DEBUG) {
            super.printReaction(reactionWithUniqueSTOICHIOMETRY);
            out.println("|++++++++++++++++++++++++++++|");
            out.println("|iv. Done|");
            out.println("|++++++++++++++++++++++++++++|" + NEW_LINE + NEW_LINE);
        }
    }

    @Override
    public String toString() {
        //SmilesGenerator smiles = new SmilesGenerator(SmiFlavor.Unique | SmiFlavor.UseAromaticSymbols | SmiFlavor.AtomAtomMap);
        SmilesGenerator smiles;
        if (partialMapping) {
            //else CDKToBeam throws an error "Aromatic bond connects non-aromatic atomic atoms"
            smiles = new SmilesGenerator(SmiFlavor.Unique
                    | SmiFlavor.AtomAtomMap
                    | SmiFlavor.Stereo);
        } else {
            smiles = new SmilesGenerator(SmiFlavor.Unique
                    | SmiFlavor.UseAromaticSymbols
                    | SmiFlavor.AtomAtomMap
                    | SmiFlavor.Stereo);
        }

        String createReactionSMILES = "";
        try {
            createReactionSMILES = smiles.create(reactionWithUniqueSTOICHIOMETRY);
        } catch (CDKException ex) {
            LOGGER.error(SEVERE, null, ex);
        }
        return "Reactor{" + "partialMapping=" + partialMapping + ", algorithm=" + algorithm
                + ", mapping=" + createReactionSMILES + '}';
    }

    private synchronized void copyReferenceReaction(IReaction referenceReaction) throws CDKException, IOException, Exception {
        try {
            for (int i = 0; i < referenceReaction.getReactantCount(); i++) {
                IAtomContainer refMol = referenceReaction.getReactants().getAtomContainer(i);
                IAtomContainer mol = cloneWithIDs(refMol);//refMol.clone();
                mol = canonicalLabelling(mol);

                mol.setID(referenceReaction.getReactants().getAtomContainer(i).getID());
                Double st = referenceReaction.getReactantCoefficient(refMol);
                aromatizeMolecule(mol);
                reactionWithSTOICHIOMETRY.addReactant(mol, st);
            }
        } catch (CloneNotSupportedException | CDKException e) {
            LOGGER.error(SEVERE, null, e);
        }
        try {
            for (int i = 0; i < referenceReaction.getProductCount(); i++) {
                IAtomContainer refMol = referenceReaction.getProducts().getAtomContainer(i);
                IAtomContainer mol = cloneWithIDs(refMol);//refMol.clone();
                mol = canonicalLabelling(mol);

                mol.setID(referenceReaction.getProducts().getAtomContainer(i).getID());
                Double st = referenceReaction.getProductCoefficient(refMol);
                aromatizeMolecule(mol);
                reactionWithSTOICHIOMETRY.addProduct(mol, st);
            }
            reactionWithSTOICHIOMETRY.setID(referenceReaction.getID());
            reactionWithSTOICHIOMETRY.setDirection(referenceReaction.getDirection());
        } catch (CloneNotSupportedException | CDKException e) {
            LOGGER.error(SEVERE, null, e);
        }
    }

    private synchronized void expandReaction() throws CloneNotSupportedException {
        for (int i = 0; i < reactionWithSTOICHIOMETRY.getReactantCount(); i++) {
            IAtomContainer _react = reactionWithSTOICHIOMETRY.getReactants().getAtomContainer(i);
            Double stoichiometry = reactionWithSTOICHIOMETRY.getReactantCoefficient(_react);
            while (stoichiometry > 0.0) {
                stoichiometry -= 1;
                IAtomContainer _reactDup = cloneWithIDs(_react);
                _reactDup.setID(_react.getID());
                _reactDup.setProperty("STOICHIOMETRY", 1.0);
                reactionWithUniqueSTOICHIOMETRY.addReactant(_reactDup, 1.0);
            }
        }

        for (int j = 0; j < reactionWithSTOICHIOMETRY.getProductCount(); j++) {
            IAtomContainer _prod = reactionWithSTOICHIOMETRY.getProducts().getAtomContainer(j);
            Double stoichiometry = reactionWithSTOICHIOMETRY.getProductCoefficient(_prod);
            while (stoichiometry > 0.0) {
                stoichiometry -= 1;
                IAtomContainer prodDup = cloneWithIDs(_prod);
                prodDup.setID(_prod.getID());
                prodDup.setProperty("STOICHIOMETRY", 1.0);
                reactionWithUniqueSTOICHIOMETRY.addProduct(prodDup, 1.0);
            }
        }
        reactionWithUniqueSTOICHIOMETRY.setID(
                reactionWithSTOICHIOMETRY.getID() == null
                ? "MappedReaction (ecBLAST)"
                : reactionWithSTOICHIOMETRY.getID());
        reactionWithUniqueSTOICHIOMETRY.setDirection(reactionWithSTOICHIOMETRY.getDirection() == null
                ? BIDIRECTIONAL
                : reactionWithSTOICHIOMETRY.getDirection());
        LabelAtoms();
        BondCollection();
    }

    private synchronized void LabelAtoms() {
        int new_atom_rank_index_reactant = 1;
        int new_atom_rank_index_product = 1;
        for (int i = 0; i < reactionWithUniqueSTOICHIOMETRY.getReactantCount(); i++) {
            IAtomContainer container = reactionWithUniqueSTOICHIOMETRY.getReactants().getAtomContainer(i);
            for (int k = 0; k < container.getAtomCount(); k++) {
                String counter = (substrateAtomCounter).toString();
                substrateAtomCounter += 1;
                IAtom atom = container.getAtom(k);
                atom.setID(counter);
                rLabelledAtoms.put(atom.hashCode(), i);
                if (atom.getProperty("OLD_RANK") != null) {
                    getInputRankLabelledAtomsReactant().put((int) atom.getProperty("OLD_RANK"), (new_atom_rank_index_reactant++));
                }
            }
            educts.put(i, container);
        }
        for (int j = 0; j < reactionWithUniqueSTOICHIOMETRY.getProductCount(); j++) {
            IAtomContainer container = reactionWithUniqueSTOICHIOMETRY.getProducts().getAtomContainer(j);
            for (int k = 0; k < container.getAtomCount(); k++) {
                String counter = (productAtomCounter).toString();
                productAtomCounter += 1;
                IAtom atom = container.getAtom(k);
                atom.setID(counter);
                pLabelledAtoms.put(atom.hashCode(), j);
                if (atom.getProperty("OLD_RANK") != null) {
                    getInputRankLabelledAtomsProduct().put((int) atom.getProperty("OLD_RANK"), (new_atom_rank_index_product++));
                }
            }
            products.put(j, container);
        }
    }

    private synchronized void BondCollection() {

        for (int i = 0; i < reactionWithUniqueSTOICHIOMETRY.getReactantCount(); i++) {
            IAtomContainer mol = reactionWithUniqueSTOICHIOMETRY.getReactants().getAtomContainer(i);
            for (int j = 0; j < mol.getBondCount(); j++) {
                IBond bond = mol.getBond(j);
                rBonds.add(bond);
            }
        }

        for (int i = 0; i < reactionWithUniqueSTOICHIOMETRY.getProductCount(); i++) {
            IAtomContainer mol = reactionWithUniqueSTOICHIOMETRY.getProducts().getAtomContainer(i);
            for (int j = 0; j < mol.getBondCount(); j++) {
                IBond bond = mol.getBond(j);
                pBonds.add(bond);
            }
        }

    }

    private synchronized void checkReactionBalance() throws IOException {
        IAtomContainerSet reactantSet = getExpandedReactants();
        IAtomContainerSet productSet = getExpandedProducts();
        HashMap<String, Integer> AtomMap = new HashMap<>();
        for (int i = 0; i < reactantSet.getAtomContainerCount(); i++) {
            IAtomContainer rMol = reactantSet.getAtomContainer(i);
            Iterator<IAtom> rAtomIterator = rMol.atoms().iterator();
            while (rAtomIterator.hasNext()) {
                IAtom rAtom = rAtomIterator.next();
                if (!rAtom.getSymbol().equals("H")) {
//                    System.out.println("E: " + rAtom.getSymbol());
                    if (AtomMap.containsKey(rAtom.getSymbol())) {
                        int count = AtomMap.get(rAtom.getSymbol()) + 1;
                        AtomMap.put(rAtom.getSymbol(), count);
                    } else {
                        AtomMap.put(rAtom.getSymbol(), 1);
                    }
                }
            }
        }

        for (int j = 0; j < productSet.getAtomContainerCount(); j++) {
            IAtomContainer pMol = productSet.getAtomContainer(j);
            Iterator<IAtom> pAtomIterator = pMol.atoms().iterator();
            while (pAtomIterator.hasNext()) {
                IAtom pAtom = pAtomIterator.next();
                if (!pAtom.getSymbol().equals("H")) {
//                    System.out.println("P: " + atomIndex.getSymbol());
                    if (AtomMap.containsKey(pAtom.getSymbol())) {
                        int count = AtomMap.get(pAtom.getSymbol()) - 1;
                        AtomMap.put(pAtom.getSymbol(), count);
                    } else if (!AtomMap.containsKey(pAtom.getSymbol())) {
                        AtomMap.put(pAtom.getSymbol(), 1);
                        this.balanceFlag = false;
                        break;
                    }
                }
            }
        }

        for (Map.Entry<String, Integer> I : AtomMap.entrySet()) {
//            System.out.println("A: " + I.getKey() + " V: " + I.getValue());
            if (I.getValue() != 0) {
                this.balanceFlag = false;
                break;
            }
        }
//        System.out.println("Bal: " + balanceFlag);
    }

    private synchronized void calculateAtomAtomMapping() throws IOException, Exception {

        try {
            IReaction reactionCopy = copyReaction(reactionWithUniqueSTOICHIOMETRY, partialMapping);
            CalculationProcess calP
                    = new CalculationProcess(partialMapping, reactionCopy, getAlgorithm());
            delta = calP.getDelta();
            IReaction mappedReaction = calP.getMappedReaction();
            reactionWithUniqueSTOICHIOMETRY = getMapping(mappedReaction);
            setReactionBlastMolMapping(calP.getReactionBlastMolMapping());
        } catch (Exception ex) {
            LOGGER.error(SEVERE, null, ex);
        }
    }

    private static void analyze_mapping(IReaction mapped_reaction) {
        boolean good_mapping = false;
        IMapping res = null;
        String re_id = find_source_id(mapped_reaction.getReactants().atomContainers());
        for (IMapping mapping : mapped_reaction.mappings()) {
            String fst_id = mapping.getChemObject(0).getID();
            String snd_id = mapping.getChemObject(1).getID();
            if (re_id.equals(fst_id)) {
                res = mapping;
                for (IAtomContainer mol : mapped_reaction.getProducts().atomContainers())
                    if (has_3_nbr(mol, snd_id)) {
                        good_mapping = true;
                        break;
                    }
                break;
            }
        }
        out.println("\tIs mapping good: " + good_mapping);
        IAtom source = (IAtom) res.getChemObject(0);
        IAtom target = (IAtom) res.getChemObject(1);
//        out.println("\t\t" + source.getID() + " <---> " + target.getID());
    }

    private static String find_source_id(Iterable<IAtomContainer> reactants) {
        for (IAtomContainer mol : reactants)
            if (mol.getAtomCount() == 3)
                for (IAtom a : mol.atoms())
                    if (a.getSymbol().equals("C")) return a.getID();
        return null;
    }

    private static boolean has_3_nbr(IAtomContainer mol, String id) {
        for (IAtom a : mol.atoms())
            if (a.getID().equals(id))
                for (IBond b : mol.bonds())
                    if (b.getAtom(0).getID().equals(id) || b.getAtom(1).getID().equals(id)) {
                        IAtom a2 = b.getOther(a);
                        int count = 0;
                        for (IBond b2 : mol.bonds())
                            if (b2.getAtom(0).getID().equals(a2.getID()) || b2.getAtom(1).getID().equals(a2.getID()))
                                count++;
                        if (count == 3) return true;
                    }
        return false;
    }

    private synchronized IReaction getMapping(IReaction coreMappedReaction) throws IOException, CDKException, CloneNotSupportedException {
        IReaction mappedReaction = deepClone(reactionWithUniqueSTOICHIOMETRY);
        cleanMapping(mappedReaction);

//        printReaction(mappedReaction);

        /*
        * This section set the mappingMap ID for the mapped atoms
         */
        int counter = 1;

        counter = setMappingFlags(mappedReaction, reactionWithUniqueSTOICHIOMETRY, coreMappedReaction, counter, algorithm == IMappingAlgorithm.MIN);

        /*
        * This section set the mappingMap ID for the unmapped atoms
        *
         */
        for (int eMol = 0; eMol < mappedReaction.getReactantCount(); eMol++) {
            IAtomContainer eMolecule = mappedReaction.getReactants().getAtomContainer(eMol);
            for (int eAtom = 0; eAtom < eMolecule.getAtomCount(); eAtom++) {
                IAtom atom = mappedReaction.getReactants().getAtomContainer(eMol).getAtom(eAtom);
                if (!atom.getSymbol().equalsIgnoreCase("H") && atom.getID().equalsIgnoreCase("-1")) {
                    String atomLabel = Integer.toString(counter);
                    atom.setID(atomLabel);
                    atom.setFlag(MAPPED, false);
                }
                counter += 1;
            }
        }

        for (int pMol = 0; pMol < mappedReaction.getProductCount(); pMol++) {
            IAtomContainer pMolecule = mappedReaction.getProducts().getAtomContainer(pMol);
            for (int pAtom = 0; pAtom < pMolecule.getAtomCount(); pAtom++) {
                IAtom atom = mappedReaction.getProducts().getAtomContainer(pMol).getAtom(pAtom);
                if (!atom.getSymbol().equalsIgnoreCase("H") && atom.getID().equalsIgnoreCase("-1")) {
                    String atomLabel = Integer.toString(counter);
                    atom.setID(atomLabel);
                    atom.setFlag(MAPPED, false);
                    counter += 1;
                }
            }
        }

        /*
        * This section will mark map common H atoms. example H-R + H <=> R-H + H Here R-H will be mapped to the R-H.
         */
        for (int eMol = 0; eMol < mappedReaction.getReactantCount(); eMol++) {
            IAtomContainer eMolecule = mappedReaction.getReactants().getAtomContainer(eMol);
            for (int eAtom = 0; eAtom < eMolecule.getAtomCount(); eAtom++) {
                IAtom atom = mappedReaction.getReactants().getAtomContainer(eMol).getAtom(eAtom);
                if (!atom.getSymbol().equalsIgnoreCase("H") && !atom.getID().equalsIgnoreCase("-1")) {
                    List<IAtom> eductConnAtoms = eMolecule.getConnectedAtomsList(atom);
                    List<IAtom> productHAtoms = markHAroundCoreAtoms(atom.getID(), mappedReaction.getProducts());
                    for (IAtom eAtomH : eductConnAtoms) {
                        //Collect ummmarked H and map common ones
                        if (eAtomH.getID().equalsIgnoreCase("-1") && eAtomH.getSymbol().equalsIgnoreCase("H")) {
                            if (!productHAtoms.isEmpty()) {
                                String atomLabel = Integer.toString(counter);
                                eAtomH.setID(atomLabel);
                                eAtomH.setFlag(MAPPED, true);
                                IAtom pAtomH = productHAtoms.iterator().next();
                                pAtomH.setID(atomLabel);
                                pAtomH.setFlag(MAPPED, true);
                                productHAtoms.remove(pAtomH);
                                counter += 1;
                            } else {
                                break;
                            }
                        }
                    }
                }
            }
        }

        /*
        *
        * This section will mark single unmapped H atoms on both sides. example R-H + H <=> H + R-H Here H and H will
        * be marked and matched
        *
         */
        List<IAtom> unMappedSingleHAtEduct = collectUnMappedSingleHAtoms(mappedReaction.getReactants());
        List<IAtom> unMappedSingleHAtProduct = collectUnMappedSingleHAtoms(mappedReaction.getProducts());

        /*
        * Mark single unmapped Hs on both sides
         */
        for (IAtom eAtomH : unMappedSingleHAtEduct) {
            if (!unMappedSingleHAtProduct.isEmpty()) {
                String atomLabel = Integer.toString(counter);
                eAtomH.setID(atomLabel);
                eAtomH.setFlag(MAPPED, true);
                IAtom pAtomH = unMappedSingleHAtProduct.iterator().next();
                pAtomH.setID(atomLabel);
                pAtomH.setFlag(MAPPED, true);
                unMappedSingleHAtProduct.remove(pAtomH);
                counter += 1;
            } else {
                break;
            }
        }

        /*
        *
        * This section will mark unmapped H atoms on both sides. example R-H + H <=> H-H + R Here H and H will be
        * marked and matched
        *
         */
        List<IAtom> unMappedHAtEduct = collectUnMappedHAtoms(mappedReaction.getReactants());
        List<IAtom> unMappedHAtProduct = collectUnMappedHAtoms(mappedReaction.getProducts());

        /*
        * Mark single unmapped Hs on both sides
         */
        for (IAtom eAtomH : unMappedHAtEduct) {
            if (!unMappedHAtProduct.isEmpty()) {
                String atomLabel = Integer.toString(counter);
                eAtomH.setID(atomLabel);
                eAtomH.setFlag(MAPPED, true);
                IAtom pAtomH = unMappedHAtProduct.iterator().next();
                pAtomH.setID(atomLabel);
                pAtomH.setFlag(MAPPED, true);
                unMappedHAtProduct.remove(pAtomH);
                counter += 1;
            } else {
                break;
            }
        }

        /*
        * Mark unmapped H atoms i.e. protonations
         */
        counter = markUnMappedHAtoms(mappedReaction, counter);

        /*
        * Create atom-atom mappingMap objects to be stored in a map
        *
         */
//        Map<IAtom, IAtom> mappings = new HashMap<>();
        for (IAtomContainer ac1 : mappedReaction.getReactants().atomContainers()) {
            for (IAtom atom1 : ac1.atoms()) {
                IAtom atom2 = getContainerAtomByID(mappedReaction.getProducts(), atom1.getID());
                if (atom2 != null) {
//                    mappings.put(atom1, atom2);
//                    String sym1 = atom1.getSymbol(), sym2 = atom2.getSymbol();
//                    int ind1 = atom1.getIndex(), ind2 = atom2.getIndex();
//                    out.println(sym1 + " : " + sym2 + ", " + ind1 + " : " + ind2);
                    IMapping m = mappedReaction.getBuilder().newInstance(IMapping.class, atom1, atom2);
                    mappedReaction.addMapping(m);
                }
            }
        }

        /*
         * Store atom-atom mappingMap objects in the reaction
         *
         */
//        mappings.keySet().stream().filter((key)
//                -> (key != null && mappings.get(key) != null)).map((key)
//                -> mappedReaction.getBuilder().newInstance(IMapping.class, key, mappings.get(key))).forEachOrdered((mappingObject) -> {
////                    IAtom fst = (IAtom) mappingObject.getChemObject(0);
////                    IAtom snd = (IAtom) mappingObject.getChemObject(1);
////                    out.println(fst.getSymbol() + " : " + snd.getSymbol());
//            mappedReaction.addMapping(mappingObject);
//        });

        /*
        * Canonical labelling of each molecule is done and mappingMap number corresponds to the lables
        *
         */
//        System.out.println("Counter Before " + counter);
        counter = setCanonicalMappingLabels(mappedReaction);
//        System.out.println("Counter After " + counter);
//        System.out.println("mappedReaction After setMapping: " + mappedReaction.getMappingCount());
//        printReaction(mappedReaction);
        return mappedReaction;
    }

//~--- get methods --------------------------------------------------------
    /**
     *
     * @return reactantSet expanded STOICHIOMETRY
     * @throws java.io.IOException
     */
    @Override
    public synchronized IAtomContainerSet getExpandedReactants() throws IOException {
        return reactionWithUniqueSTOICHIOMETRY.getReactants();
    }

    /**
     *
     * @return productSet expanded STOICHIOMETRY
     * @throws java.io.IOException
     */
    @Override
    public synchronized IAtomContainerSet getExpandedProducts() throws IOException {
        return reactionWithUniqueSTOICHIOMETRY.getProducts();
    }

    /**
     *
     * @return IReaction object with unique atom labeling
     * @throws Exception
     */
    @Override
    public synchronized IReaction getReactionWithAtomAtomMapping() throws Exception {
        return reactionWithUniqueSTOICHIOMETRY;
    }

    /**
     *
     * @param i Index Ith position
     * @return Stoichiometry weight of the reactant molecule at ith Position
     */
    @Override
    public synchronized Double getExpandedReactantStoichiometry(
            int i) {
        IAtomContainer Mol = reactionWithUniqueSTOICHIOMETRY.getReactants().getAtomContainer(i);
        return reactionWithUniqueSTOICHIOMETRY.getReactantCoefficient(Mol);

    }

    /**
     *
     * @param i Index at I th position
     * @return Stoichiometry weight of the product molecule at i th Position
     *
     */
    @Override
    public synchronized Double getExpandedProductStoichiometry(int i) {
        IAtomContainer Mol = reactionWithUniqueSTOICHIOMETRY.getProducts().getAtomContainer(i);
        return reactionWithUniqueSTOICHIOMETRY.getProductCoefficient(Mol);
    }

    /**
     *
     * @return true if its a balanced reaction else false
     * @throws IOException
     *
     */
    @Override
    public synchronized boolean getReactionBalanceFlag() throws IOException {

        boolean flag = true;
        if (!Objects.equals(this.getLabledReactantAtomsCount(), this.getLabledProductAtomsCount())) {
            flag = false;
        }
        if (!getReactionBalanceFlagWithoutHydrogen()) {
            flag = false;
        }
        return flag;
    }

    /**
     *
     * @return @throws IOException
     */
    @Override
    public synchronized boolean getReactionBalanceFlagWithChargeBalance() throws IOException {
        boolean flag = true;
        if (!Objects.equals(this.getLabledReactantAtomsCount(), this.getLabledProductAtomsCount())) {
            flag = false;
        }
        if (getTotalFormalCharge(this.getExpandedReactants()) != getTotalFormalCharge(this.getExpandedProducts())) {
            flag = false;
        }
        if (!getReactionBalanceFlagWithoutHydrogen()) {
            flag = false;
        }
        return flag;
    }

    /**
     *
     * @return true if its a balanced reaction else false Note: This does not
     * consider whether Hydrogens are balanced or not
     *
     *
     */
    @Override
    public synchronized boolean getReactionBalanceFlagWithoutHydrogen() {
        return this.balanceFlag;
    }

    /**
     *
     * @return this will return IAtom Vector of total Reactant atom count with
     * unique labeling
     *
     *
     */
    private synchronized List<IAtom> getLabledReactantAtoms() {
        List<IAtom> reactantAtoms = new ArrayList<>();
        IAtomContainerSet MSet = reactionWithUniqueSTOICHIOMETRY.getReactants();
        for (int j = 0; j
                < MSet.getAtomContainerCount(); j++) {
            IAtomContainer M = MSet.getAtomContainer(j);
            for (int k = 0; k
                    < M.getAtomCount(); k++) {
                reactantAtoms.add(M.getAtom(k));
            }
        }
        return unmodifiableList(reactantAtoms);
    }

    /**
     *
     * @return this will return IAtom Vector of total Product atom count with
     * unique labelling
     *
     */
    private synchronized List<IAtom> getLabledProductAtoms() {
        List<IAtom> productAtoms = new ArrayList<>();
        IAtomContainerSet MSet = reactionWithUniqueSTOICHIOMETRY.getProducts();
        for (int j = 0; j
                < MSet.getAtomContainerCount(); j++) {
            IAtomContainer M = MSet.getAtomContainer(j);
            for (int k = 0; k
                    < M.getAtomCount(); k++) {
                productAtoms.add(M.getAtom(k));
            }
        }
        return unmodifiableList(productAtoms);
    }

    /**
     *
     * @return this will return the total reactant
     *
     * atom count with unique labeling
     *
     *
     */
    private synchronized Integer getLabledReactantAtomsCount() {
        return getLabledReactantAtoms().size();
    }

    /**
     *
     * @return this will return the total product
     *
     * atom count with unique labeling
     *
     */
    private synchronized Integer getLabledProductAtomsCount() {
        return getLabledProductAtoms().size();

    }

    /**
     *
     * @return bonds of reactantSet
     */
    @Override
    public synchronized List<IBond> getEductBonds() {
        return unmodifiableList(rBonds);
    }

    /**
     * @return bonds of productSet
     */
    @Override
    public synchronized List<IBond> getProductBonds() {
        return unmodifiableList(pBonds);
    }

    /**
     *
     * @return
     */
    @Override
    public synchronized int getMappingCount() {
        return reactionWithUniqueSTOICHIOMETRY.getMappingCount();
    }

    private synchronized IReaction copyReaction(IReaction orignalReaction, boolean removeHydrogen) throws Exception {
//        System.out.println("R size: " + orignalReaction.getReactantCount() + " , " + orignalReaction.getProductCount());
        IReaction copiedReaction = reactionWithUniqueSTOICHIOMETRY.getBuilder().newInstance(IReaction.class);

        for (int i = 0; i < orignalReaction.getReactantCount(); i++) {
            IAtomContainer mol = orignalReaction.getReactants().getAtomContainer(i);
            Double st = orignalReaction.getReactantCoefficient(mol);
            IAtomContainer newMol = cloneWithIDs(mol);
            for (int index = 0; index < mol.getAtomCount(); index++) {
                mol.getAtom(index).setProperty("index", index);
                IAtom atom = newMol.getAtom(index);
                atom.setProperty("index", index);
            }
//            System.out.println("Hydrogen Before" + newMol.getAtomCount());

            percieveAtomTypesAndConfigureAtoms(newMol);
            if (removeHydrogen) {
                newMol = removeHydrogensExceptSingleAndPreserveAtomID(newMol);
            }
//            System.out.println("Hydrogen After" + newMol.getAtomCount());
            copiedReaction.addReactant(newMol, st);
        }
        for (int i = 0; i < orignalReaction.getProductCount(); i++) {
            IAtomContainer mol = orignalReaction.getProducts().getAtomContainer(i);
            Double st = orignalReaction.getProductCoefficient(mol);
            IAtomContainer newMol = cloneWithIDs(mol);
            for (int index = 0; index < mol.getAtomCount(); index++) {
                mol.getAtom(index).setProperty("index", index);
                IAtom atom = newMol.getAtom(index);
                atom.setProperty("index", index);
            }
//            System.out.println("Hydrogen Before " + newMol.getAtomCount());

            percieveAtomTypesAndConfigureAtoms(newMol);
            if (removeHydrogen) {
                newMol = removeHydrogensExceptSingleAndPreserveAtomID(newMol);
            }
//            System.out.println("Hydrogen After " + newMol.getAtomCount());
            copiedReaction.addProduct(newMol, st);
        }
        copiedReaction.setFlags(orignalReaction.getFlags());
        copiedReaction.setID(orignalReaction.getID());
        copiedReaction.setDirection(orignalReaction.getDirection());
        copiedReaction.notifyChanged();
        return copiedReaction;
    }

    /**
     *
     * @param id
     * @param molSet
     * @return
     */
    private synchronized List<IAtom> markHAroundCoreAtoms(String id, IAtomContainerSet molSet) {

        List<IAtom> list = new ArrayList<>();
        for (int pMol = 0; pMol < molSet.getAtomContainerCount(); pMol++) {
            IAtomContainer pMolecule = molSet.getAtomContainer(pMol);
            for (int pAtom = 0; pAtom < pMolecule.getAtomCount(); pAtom++) {
                IAtom atom = molSet.getAtomContainer(pMol).getAtom(pAtom);
                if (!atom.getSymbol().equalsIgnoreCase("H") && !atom.getID().equalsIgnoreCase("-1")) {
                    if (atom.getID().equalsIgnoreCase(id)) {
                        List<IAtom> conAtoms = pMolecule.getConnectedAtomsList(atom);
                        conAtoms.stream().filter((atomH) -> (atomH.getID().equalsIgnoreCase("-1") && atomH.getSymbol().equalsIgnoreCase("H"))).forEach((atomH) -> {
                            list.add(atomH);
                        });
                    }
                }
            }
        }
        return list;
    }

    /**
     * @param molSet
     * @return
     */
    private synchronized List<IAtom> collectUnMappedSingleHAtoms(IAtomContainerSet molSet) {

        List<IAtom> list = new ArrayList<>();
        for (int index = 0; index < molSet.getAtomContainerCount(); index++) {
            IAtomContainer mol = molSet.getAtomContainer(index);
            if (mol.getAtomCount() == 1) {
                for (int atomIndex = 0; atomIndex < mol.getAtomCount(); atomIndex++) {
                    IAtom atom = molSet.getAtomContainer(index).getAtom(atomIndex);
                    if (atom.getSymbol().equalsIgnoreCase("H")
                            && !atom.getFlag(MAPPED)
                            && atom.getID().equalsIgnoreCase("-1")) {
                        list.add(atom);
                    }
                }
            }
        }
        return list;
    }

    /**
     * @param molSet
     * @return
     */
    private synchronized List<IAtom> collectUnMappedHAtoms(IAtomContainerSet molSet) {

        List<IAtom> list = new ArrayList<>();
        for (int index = 0; index < molSet.getAtomContainerCount(); index++) {
            IAtomContainer mol = molSet.getAtomContainer(index);
            for (int atomIndex = 0; atomIndex < mol.getAtomCount(); atomIndex++) {
                IAtom atom = molSet.getAtomContainer(index).getAtom(atomIndex);
                if (atom.getSymbol().equalsIgnoreCase("H")
                        && !atom.getFlag(MAPPED)
                        && atom.getID().equalsIgnoreCase("-1")) {
                    list.add(atom);
                }
            }
        }
        return list;
    }

    /**
     *
     * @param mappedReaction
     * @param counter
     * @return updated Counter
     */
    private synchronized int markUnMappedHAtoms(IReaction mappedReaction, int counter) {

        int localCounter = counter;

        /*
        * Mark unmapped H atoms
         */
        for (int eMol = 0; eMol < mappedReaction.getReactantCount(); eMol++) {
            IAtomContainer eMolecule = mappedReaction.getReactants().getAtomContainer(eMol);
            for (int eAtom = 0; eAtom < eMolecule.getAtomCount(); eAtom++) {
                IAtom atom = mappedReaction.getReactants().getAtomContainer(eMol).getAtom(eAtom);
                if (atom.getSymbol().equalsIgnoreCase("H") && !atom.getFlag(MAPPED)
                        && atom.getID().equalsIgnoreCase("-1")) {
                    String atomLabel = Integer.toString(localCounter);
                    atom.setFlag(MAPPED, false);
                    atom.setID(atomLabel);
                    localCounter += 1;
                }
            }
        }

        for (int pMol = 0; pMol < mappedReaction.getProductCount(); pMol++) {
            IAtomContainer pMolecule = mappedReaction.getProducts().getAtomContainer(pMol);
            for (int pAtom = 0; pAtom < pMolecule.getAtomCount(); pAtom++) {
                IAtom atom = mappedReaction.getProducts().getAtomContainer(pMol).getAtom(pAtom);
                if (atom.getSymbol().equalsIgnoreCase("H") && !atom.getFlag(MAPPED)
                        && atom.getID().equalsIgnoreCase("-1")) {
                    String atomLabel = Integer.toString(localCounter);
                    atom.setID(atomLabel);
                    atom.setFlag(MAPPED, false);
                    localCounter += 1;
                }
            }
        }
        return localCounter;
    }

    /**
     * @return the delta
     */
    public synchronized int getDelta() {
        return delta;
    }

    /**
     * @return the reactionBlastMolMapping
     */
    public synchronized MoleculeMoleculeMapping getReactionBlastMolMapping() {
        return reactionBlastMolMapping;
    }

    /**
     * @param reactionBlastMolMapping the reactionBlastMolMapping to set
     */
    private synchronized void setReactionBlastMolMapping(MoleculeMoleculeMapping reactionBlastMolMapping) {
        this.reactionBlastMolMapping = reactionBlastMolMapping;
    }

    private synchronized IAtom getContainerAtomByID(IAtomContainerSet products, String mappingID) {
        for (IAtomContainer ac : products.atomContainers()) {
            for (IAtom atom : ac.atoms()) {
                if (atom.getID().equals(mappingID)) {
                    return atom;
                }
            }
        }
        return null;
    }

    private synchronized int setCanonicalMappingLabels(IReaction mappedReaction) throws CDKException {
//        ICanonicalMoleculeLabeller cng = new SignatureMoleculeLabeller();

        IAtomContainerSet rMolSet = mappedReaction.getReactants();
        IAtomContainerSet pMolSet = mappedReaction.getProducts();

        Map<IAtom, IAtom> mappingMap = new HashMap<>();

        for (IMapping aaMapping : mappedReaction.mappings()) {
            aaMapping.getChemObject(0).removeProperty(ATOM_ATOM_MAPPING);
            aaMapping.getChemObject(1).removeProperty(ATOM_ATOM_MAPPING);
            mappingMap.put((IAtom) aaMapping.getChemObject(0), (IAtom) aaMapping.getChemObject(1));
        }

        /*
        Re-arrange the molecule index for mappings
         */
        for (IAtomContainer mol : rMolSet.atomContainers()) {
            List<Integer> atom_index = new ArrayList<>();
            try {
                int[] p = new int[mol.getAtomCount()];
                String smiles = unique().create(mol, p);
            } catch (CDKException e) {
                LOGGER.error(SEVERE, null, e);
            }
//            int[] canonicalPermutation = cng.getCanonicalPermutation(mol);
//            permuteWithoutClone(canonicalPermutation, mol);

            for (IAtom a : mol.atoms()) {
                if (!a.getSymbol().equalsIgnoreCase("H")) {
                    atom_index.add(mol.indexOf(a));
                }
            }
            for (IAtom a : mol.atoms()) {
                if (a.getSymbol().equalsIgnoreCase("H")) {
                    atom_index.add(mol.indexOf(a));
                }
            }
            int[] array = new int[atom_index.size()];

            int index = 0;
            for (int c : atom_index) {
                array[index] = c;
                index++;
            }
            permuteWithoutClone(array, mol);
        }

        for (IAtomContainer mol : pMolSet.atomContainers()) {

            try {
                int[] p = new int[mol.getAtomCount()];
                String smiles = unique().create(mol, p);
            } catch (CDKException e) {
                LOGGER.error(SEVERE, null, e);
            }

//            int[] canonicalPermutation = cng.getCanonicalPermutation(mol);
//            permuteWithoutClone(canonicalPermutation, mol);
            List<Integer> atom_index = new ArrayList<>();
            for (IAtom a : mol.atoms()) {
                if (!a.getSymbol().equalsIgnoreCase("H")) {
                    atom_index.add(mol.indexOf(a));
                }
            }
            for (IAtom a : mol.atoms()) {
                if (a.getSymbol().equalsIgnoreCase("H")) {
                    atom_index.add(mol.indexOf(a));
                }
            }
            int[] array = new int[atom_index.size()];

            int index = 0;
            for (int c : atom_index) {
                array[index] = c;
                index++;
            }
            permuteWithoutClone(array, mol);
        }

        int counter = 1;
        for (IAtomContainer mol : rMolSet.atomContainers()) {
            /*
            * Assign mappingMap to non H atoms in reactant and product
             */
            for (IAtom qAtom : mol.atoms()) {
                if (mappingMap.containsKey(qAtom) && !qAtom.getSymbol().equalsIgnoreCase("H")) {
//                    System.out.println("Atom " + qAtom.getSymbol() + " new Rank: " + counter);
                    String id = valueOf(counter);
                    qAtom.setID(id);
                    mappingMap.get(qAtom).setID(id);
                    qAtom.setProperty(ATOM_ATOM_MAPPING, parseInt(qAtom.getID()));
                    mappingMap.get(qAtom).setProperty(ATOM_ATOM_MAPPING, parseInt(mappingMap.get(qAtom).getID()));
                    counter++;
                }
            }
        }

        for (IAtomContainer mol : rMolSet.atomContainers()) {
            /*
            * Assign mappingMap to non H atoms in reactant and product
             */
            for (IAtom qAtom : mol.atoms()) {
                if (mappingMap.containsKey(qAtom) && qAtom.getSymbol().equalsIgnoreCase("H")) {
//                    System.out.println("Atom " + qAtom.getSymbol() + " new Rank: " + counter);
                    String id = valueOf(counter);
                    qAtom.setID(id);
                    mappingMap.get(qAtom).setID(id);
                    qAtom.setProperty(ATOM_ATOM_MAPPING, parseInt(qAtom.getID()));
                    mappingMap.get(qAtom).setProperty(ATOM_ATOM_MAPPING, parseInt(mappingMap.get(qAtom).getID()));
                    counter++;
                }
            }
        }

        /*
        * Assign mappingMap to atoms which are not mapped in the reactant
         */
        for (IAtomContainer mol : rMolSet.atomContainers()) {
            for (IAtom qAtom : mol.atoms()) {
                if (!mappingMap.containsKey(qAtom)) {
                    String id = valueOf(counter);
                    qAtom.setID(id);
                    qAtom.setProperty(ATOM_ATOM_MAPPING, parseInt(qAtom.getID()));
                    counter++;
                }
            }
        }

        /*
        * Assign mappingMap to atoms which are not mapped in the product
         */
        for (IAtomContainer mol : pMolSet.atomContainers()) {
            for (IAtom tAtom : mol.atoms()) {
                if (!mappingMap.containsValue(tAtom)) {
                    String id = valueOf(counter);
                    tAtom.setID(id);
                    tAtom.setProperty(ATOM_ATOM_MAPPING, parseInt(tAtom.getID()));
                    counter++;
                }
            }
        }
        /*
        Finally permute molecules based on the atom mapping rank
         */
        for (IAtomContainer mol : pMolSet.atomContainers()) {
            TreeMap<Integer, Integer> mapping_rank = new TreeMap<>();
            for (IAtom a : mol.atoms()) {
                mapping_rank.put((Integer) a.getProperty(ATOM_ATOM_MAPPING), mol.indexOf(a));
            }
            int[] mappingIndexPermutation = new int[mapping_rank.size()];
            int index = 0;
            for (int i : mapping_rank.values()) {
                mappingIndexPermutation[index] = i;
                index++;
            }
            permuteWithoutClone(mappingIndexPermutation, mol);
        }

        mappingMap.clear();
        return counter;
    }

    /**
     * @return the algorithm
     */
    public synchronized IMappingAlgorithm getAlgorithm() {
        return algorithm;
    }

    private IAtomContainer canonicalLabelling(IAtomContainer org_mol) throws CloneNotSupportedException, CDKException {

        IAtomContainer cloneMolecule = cloneWithIDs(org_mol);

        if (DEBUG) {
            LOGGER.debug("Orignal");
            printAtoms(cloneMolecule);
        }

        if (DEBUG) {
            LOGGER.debug(NEW_LINE + "mol before: ");
            printAtoms(cloneMolecule);
        }
        /*
        Use the Canonical labelling from the SMILES
        IMP: Suggested by John May
         */
        int[] p = new int[cloneMolecule.getAtomCount()];

        try {
            String smiles = unique().create(cloneMolecule, p);
            if (DEBUG) {
                LOGGER.debug("smiles " + smiles);
            }
        } catch (CDKException e) {
            LOGGER.error(SEVERE, null, e);
        }

        permuteWithoutClone(p, cloneMolecule);

        if (DEBUG) {
            LOGGER.debug("mol after: ");
            printAtoms(cloneMolecule);
        }

        /*
        Generate 2D Diagram without cloning
         */
        if (!has2DCoordinates(cloneMolecule)) {
            try {
                /*
                Clone it else it will loose mol ID
                 */
                StructureDiagramGenerator sdg = new StructureDiagramGenerator();
                sdg.setMolecule(cloneMolecule, false);
                sdg.generateCoordinates();
            } catch (CDKException e) {
            }
        }

        /*
        Set the IDs to -1 very IMP
         */
        for (IAtom atom : cloneMolecule.atoms()) {
            atom.setID("-1");
        }

        /*
        Set the IDs to container
         */
        if (org_mol.getID() != null) {
            cloneMolecule.setID(org_mol.getID());
        }

        if (DEBUG) {
            LOGGER.debug("Processed");
            printAtoms(cloneMolecule);
            LOGGER.debug("canonicalMolecule: "
                    + generic().create(cloneMolecule)
                    + NEW_LINE + NEW_LINE);
        }

        return cloneMolecule;
    }

    /*
    This is a very imp code modified by John May
    The idea is to canonicalise the atoms and bonds
     */
    private void permuteWithoutClone(int[] p, IAtomContainer atomContainer) {
        int n = atomContainer.getAtomCount();
        if (DEBUG) {
            LOGGER.debug("permuting " + java.util.Arrays.toString(p));
        }
        IAtom[] permutedAtoms = new IAtom[n];

        for (int i = 0; i < n; i++) {
            IAtom atom = atomContainer.getAtom(i);
            permutedAtoms[p[i]] = atom;
            atom.setProperty("label", p[i]);
        }
        atomContainer.setAtoms(permutedAtoms);

        IBond[] bonds = getBondArray(atomContainer);
        sort(bonds, (IBond o1, IBond o2) -> {
            int u = o1.getAtom(0).getProperty("label");
            int v = o1.getAtom(1).getProperty("label");
            int x = o2.getAtom(0).getProperty("label");
            int y = o2.getAtom(1).getProperty("label");
            int min1 = min(u, v);
            int min2 = min(x, y);
            int max1 = max(u, v);
            int max2 = max(x, y);

            int minCmp = Integer.compare(min1, min2);
            if (minCmp != 0) {
                return minCmp;
            }
            int maxCmp = Integer.compare(max1, max2);
            if (maxCmp != 0) {
                return maxCmp;
            }
            LOGGER.debug("pokemon!");
            throw new InternalError();
        });
        atomContainer.setBonds(bonds);
    }

    /**
     * Old Atom Rank in the reactant mapped to new Rank
     *
     * @return the inputRankLabelledAtomsReactant
     */
    public Map<Integer, Integer> getInputRankLabelledAtomsReactant() {
        return inputRankLabelledAtomsReactant;
    }

    /**
     * Old Atom Rank in the product mapped to new Rank
     *
     * @return the inputRankLabelledAtomsProduct
     */
    public Map<Integer, Integer> getInputRankLabelledAtomsProduct() {
        return inputRankLabelledAtomsProduct;
    }
}
