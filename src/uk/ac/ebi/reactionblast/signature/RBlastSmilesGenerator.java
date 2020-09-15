package uk.ac.ebi.reactionblast.signature;

/*
 * $Revision$ $Author$ $Date$
 *
 * Copyright (C) 2002-2007 Oliver Horlacher
 *
 * Contact: cdk-devel@lists.sourceforge.net
 *
 * This program is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either version 2.1 of the License, or (at your option)
 * any later version. All we ask is that proper credit is given for our work, which includes - but is not limited to -
 * adding the above copyright notice to the beginning of your source code files, and to any copyright notice that you
 * may distribute with programs based on this work.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License along with this program; if not, write to
 * the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
import java.io.IOException;
import static java.lang.Math.PI;
import static java.lang.Math.abs;
import java.util.ArrayList;
import java.util.Collection;
import static java.util.Collections.sort;
import java.util.Iterator;
import java.util.List;
import java.util.Objects;
import java.util.TreeMap;
import java.util.Vector;

import static org.openscience.cdk.CDKConstants.BONDORDER_DOUBLE;
import static org.openscience.cdk.CDKConstants.ISAROMATIC;
import static org.openscience.cdk.CDKConstants.UNSET;
import static org.openscience.cdk.CDKConstants.VISITED;
import org.openscience.cdk.PseudoAtom;

import org.openscience.cdk.config.IsotopeFactory;
import static org.openscience.cdk.config.Isotopes.getInstance;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.geometry.BondTools;
import static org.openscience.cdk.geometry.BondTools.giveAngle;
import static org.openscience.cdk.geometry.BondTools.giveAngleBothMethods;
import static org.openscience.cdk.geometry.BondTools.giveAngleFromMiddle;
import static org.openscience.cdk.geometry.BondTools.isCisTrans;
import static org.openscience.cdk.geometry.BondTools.isLeft;
import static org.openscience.cdk.geometry.BondTools.isSquarePlanar;
import static org.openscience.cdk.geometry.BondTools.isStereo;
import static org.openscience.cdk.geometry.BondTools.isTetrahedral;
import static org.openscience.cdk.geometry.BondTools.isTrigonalBipyramidalOrOctahedral;
import static org.openscience.cdk.graph.ConnectivityChecker.partitionIntoMolecules;
import static org.openscience.cdk.graph.invariant.MorganNumbersTools.getMorganNumbersWithElementSymbol;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import static org.openscience.cdk.interfaces.IAtomType.Hybridization.PLANAR3;
import org.openscience.cdk.interfaces.IBond;
import static org.openscience.cdk.interfaces.IBond.Order.DOUBLE;
import static org.openscience.cdk.interfaces.IBond.Order.SINGLE;
import static org.openscience.cdk.interfaces.IBond.Order.TRIPLE;
import static org.openscience.cdk.interfaces.IBond.Stereo.NONE;
import org.openscience.cdk.interfaces.IChemObject;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.interfaces.IIsotope;
import org.openscience.cdk.interfaces.IPseudoAtom;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.ringsearch.AllRingsFinder;
import static org.openscience.cdk.ringsearch.RingPartitioner.partitionRings;
import static org.openscience.cdk.tools.manipulator.RingSetManipulator.getAllAtomContainers;
import static uk.ac.ebi.reactionblast.tools.ExtAtomContainerManipulator.aromatizeDayLight;
import static uk.ac.ebi.reactionblast.tools.ExtAtomContainerManipulator.percieveAtomTypesAndConfigureAtoms;
import uk.ac.ebi.reactionblast.tools.labelling.ICanonicalMoleculeLabeller;

/**
 * Generates SMILES strings {
 *
 * @cdk.cite WEI88, WEI89}. It takes into account the isotope and formal charge
 * information of the atoms. In addition to this it takes stereochemistry in
 * account for both Bond's and Atom's. Via the flag useAromaticity it can be set
 * if only SP2-hybridized atoms shall be set to lower case (default) or atoms,
 * which are SP2 or aromatic.
 *
 * <p>
 * Some example code:  <pre>
 * IAtomContainer benzene; // single/aromatic bonds between 6 carbons
 * SmilesGenerator sg = new SmilesGenerator();
 * String smiles = sg.createSMILES(benzene); // C1CCCCC1
 * sg.setUseAromaticityFlag(true);
 * smiles = sg.createSMILES(benzene); // c1ccccc1
 * IAtomContainer benzene2; // one of the two kekule structures with explicit double bond orders
 * String smiles2 = sg.createSMILES(benzene2); // C1=CC=CC=C1
 * </pre> <b>Note</b>Due to the way the initial atom labeling is constructed,
 * ensure that the input molecule is appropriately configured. In absence of
 * such configuration it is possible that different forms of the same molecule
 * will not result in the same canonical SMILES.
 *
 * @author Oliver Horlacher
 * @author Stefan Kuhn (chiral smiles) @cdk.created 2002-02-26 @cdk.keyword
 * SMILES, generator @cdk.module smiles
 * @cdk.bug 1793446
 */
public class RBlastSmilesGenerator {

    //private final static boolean debug = false;
    /**
     * The number of rings that have been opened
     */
    private int ringMarker = 0;
    /**
     * Collection of all the bonds that were broken
     */
    private List<BrokenBond> brokenBonds = new ArrayList<>();
    /**
     * The isotope factory which is used to write the mass is needed
     */
    private IsotopeFactory isotopeFactory;
    AllRingsFinder ringFinder;
    /**
     * RingSet that holds all rings of the molecule
     */
    private IRingSet rings = null;
    /**
     * The canonical labler
     */
//    private CanonicalLabeler canLabler = new CanonicalLabeler();
    private ICanonicalMoleculeLabeller labeller;
    private final String RING_CONFIG = "stereoconfig";
    private final String UP = "up";
    private final String DOWN = "down";
    private boolean useAromaticityFlag = false;

    /**
     * Create the SMILES generator.
     */
    public RBlastSmilesGenerator() {
    }

    /**
     * Create the SMILES generator.
     *
     * @param useAromaticityFlag if false only SP2-hybridized atoms will be
     * lower case (default), true=SP2 or aromaticity trigger lower case (same as
     * using setUseAromaticityFlag later)
     */
    public RBlastSmilesGenerator(boolean useAromaticityFlag) {
        this.useAromaticityFlag = useAromaticityFlag;
    }

    /**
     *
     * @param useAromaticityFlag
     * @param labeller
     */
    public RBlastSmilesGenerator(boolean useAromaticityFlag, ICanonicalMoleculeLabeller labeller) {
        this.useAromaticityFlag = useAromaticityFlag;
        this.labeller = labeller;
    }

    /**
     * Tells if a certain bond is center of a valid double bond configuration.
     *
     * @param container The atomcontainer.
     * @param bond The bond.
     * @return true=is a potential configuration, false=is not.
     */
    public boolean isValidDoubleBondConfiguration(IAtomContainer container, IBond bond) {
        IAtom atom0 = bond.getAtom(0);
        IAtom atom1 = bond.getAtom(1);
        List<IAtom> connectedAtoms = container.getConnectedAtomsList(atom0);
        IAtom from = null;
        for (IAtom connectedAtom : connectedAtoms) {
            if (connectedAtom != atom1) {
                from = connectedAtom;
            }
        }
        boolean[] array = new boolean[container.getBondCount()];
        for (int i = 0; i < array.length; i++) {
            array[i] = true;
        }
        if (isStartOfDoubleBond(container, atom0, from, array) && isEndOfDoubleBond(container, atom1, atom0, array) && !bond.getFlag(ISAROMATIC)) {
            return (true);
        } else {
            return (false);
        }
    }

    /**
     * Provide a reference to a RingSet that holds ALL rings of the
     * molecule.<BR> During creation of a SMILES the aromaticity of the molecule
     * has to be detected. This, in turn, requires the determination of all
     * rings of the molecule. If this computationally expensive calculation has
     * been done beforehand, a RingSet can be handed over to the SmilesGenerator
     * to save the effort of another all-rings- calculation.
     *
     * @param rings RingSet that holds ALL rings of the molecule
     * @return reference to the SmilesGenerator object this method was called
     * for
     */
    public RBlastSmilesGenerator setRings(IRingSet rings) {
        this.rings = rings;
        return this;
    }

    /**
     * Generate canonical SMILES from the <code>molecule</code>. This method
     * canonically labels the molecule but does not perform any checks on the
     * chemical validity of the molecule. IMPORTANT: A precomputed Set of All
     * Rings (SAR) can be passed to this SmilesGenerator in order to avoid
     * recomputing it. Use setRings() to assign the SAR.
     *
     * @param molecule The molecule to evaluate
     * @see
     * org.openscience.cdk.graph.invariant.CanonicalLabeler#canonLabel(IAtomContainer)
     * @return the SMILES representation of the molecule
     */
    public synchronized String createSMILES(IAtomContainer molecule) {
        try {
            return (createSMILES(molecule, false, new boolean[molecule.getBondCount()]));
        } catch (CDKException exception) {
            // This exception can only happen if a chiral smiles is requested
            return ("");
        }
    }

    /**
     * Generate a SMILES for the given <code>Reaction</code>.
     *
     * @param reaction the reaction in question
     * @return the SMILES representation of the reaction
     * @throws org.openscience.cdk.exception.CDKException if there is an error
     * during SMILES generation
     */
    public synchronized String createSMILES(IReaction reaction) throws CDKException {
        StringBuilder reactionSMILES = new StringBuilder();
        IAtomContainerSet reactants = reaction.getReactants();
        for (int i = 0; i < reactants.getAtomContainerCount(); i++) {
            reactionSMILES.append(createSMILES(reactants.getAtomContainer(i)));
            if (i + 1 < reactants.getAtomContainerCount()) {
                reactionSMILES.append('.');
            }
        }
        reactionSMILES.append('>');
        IAtomContainerSet agents = reaction.getAgents();
        for (int i = 0; i < agents.getAtomContainerCount(); i++) {
            reactionSMILES.append(createSMILES(agents.getAtomContainer(i)));
            if (i + 1 < agents.getAtomContainerCount()) {
                reactionSMILES.append('.');
            }
        }
        reactionSMILES.append('>');
        IAtomContainerSet products = reaction.getProducts();
        for (int i = 0; i < products.getAtomContainerCount(); i++) {
            reactionSMILES.append(createSMILES(products.getAtomContainer(i)));
            if (i + 1 < products.getAtomContainerCount()) {
                reactionSMILES.append('.');
            }
        }
        return reactionSMILES.toString();
    }

    /**
     * Generate canonical and chiral SMILES from the <code>molecule</code>. This
     * method canonicaly lables the molecule but dose not perform any checks on
     * the chemical validity of the molecule. The chiral smiles is done like in
     * the
     * <a href="http://www.daylight.com/dayhtml/doc/theory/theory.smiles.html">
     * daylight theory manual</a> . I did not find rules for canonical and
     * chiral smiles, therefore there is no guarantee that the smiles complies
     * to any externeal rules, but it is canonical compared to other smiles
     * produced by this method. The method checks if there are 2D coordinates
     * but does not check if coordinates make sense. Invalid stereo
     * configurations are ignored; if there are no valid stereo configuration
     * the smiles will be the same as the non-chiral one. Note that often stereo
     * configurations are only complete and can be converted to a smiles if
     * explicit Hs are given. IMPORTANT: A precomputed Set of All Rings (SAR)
     * can be passed to this SmilesGenerator in order to avoid recomputing it.
     * Use setRings() to assign the SAR.
     *
     * @param molecule The molecule to evaluate.
     * @param doubleBondConfiguration Should E/Z configurations be read at these
     * positions? If the flag at position X is set to true, an E/Z configuration
     * will be written from coordinates around bond X, if false, it will be
     * ignored. If flag is true for a bond which does not constitute a valid
     * double bond configuration, it will be ignored (meaning setting all to
     * true will create E/Z indication will be pu in the smiles wherever
     * possible, but note the coordinates might be arbitrary).
     * @exception CDKException At least one atom has no Point2D; coordinates are
     * needed for creating the chiral smiles.
     * @see
     * org.openscience.cdk.graph.invariant.CanonicalLabeler#canonLabel(IAtomContainer)
     * @return the SMILES representation of the molecule
     */
    public synchronized String createChiralSMILES(IAtomContainer molecule, boolean[] doubleBondConfiguration) throws CDKException {
        return (createSMILES(molecule, true, doubleBondConfiguration));
    }

    /**
     * Generate canonical SMILES from the <code>molecule</code>. This method
     * canonicaly lables the molecule but dose not perform any checks on the
     * chemical validity of the molecule. This method also takes care of
     * multiple molecules. IMPORTANT: A precomputed Set of All Rings (SAR) can
     * be passed to this SmilesGenerator in order to avoid recomputing it. Use
     * setRings() to assign the SAR.
     *
     * @param molecule The molecule to evaluate.
     * @param chiral true=SMILES will be chiral, false=SMILES. will not be
     * chiral.
     * @param doubleBondConfiguration Should E/Z configurations be read at these
     * positions? If the flag at position X is set to true, an E/Z configuration
     * will be written from coordinates around bond X, if false, it will be
     * ignored. If flag is true for a bond which does not constitute a valid
     * double bond configuration, it will be ignored (meaning setting all to
     * true will create E/Z indication will be pu in the smiles wherever
     * possible, but note the coordinates might be arbitrary).
     * @exception CDKException At least one atom has no Point2D; coordinates are
     * needed for crating the chiral smiles. This excpetion can only be thrown
     * if chiral smiles is created, ignore it if you want a non-chiral smiles
     * (createSMILES(GraphAtomContainer) does not throw an exception).
     * @see
     * org.openscience.cdk.graph.invariant.CanonicalLabeler#canonLabel(IAtomContainer)
     * @return the SMILES representation of the molecule
     */
    public synchronized String createSMILES(IAtomContainer molecule, boolean chiral, boolean doubleBondConfiguration[]) throws CDKException {
        IAtomContainerSet moleculeSet = partitionIntoMolecules(molecule);
        if (moleculeSet.getAtomContainerCount() > 1) {
            StringBuilder fullSMILES = new StringBuilder();
            for (int i = 0; i < moleculeSet.getAtomContainerCount(); i++) {
                IAtomContainer molPart = moleculeSet.getAtomContainer(i);
                fullSMILES.append(createSMILESWithoutCheckForMultipleMolecules(
                        molPart, chiral, doubleBondConfiguration));
                if (i < (moleculeSet.getAtomContainerCount() - 1)) {
                    // are there more molecules?
                    fullSMILES.append('.');
                }
            }
            return fullSMILES.toString();
        } else {
            return (createSMILESWithoutCheckForMultipleMolecules(molecule,
                    chiral, doubleBondConfiguration));
        }
    }

    /**
     * Generate canonical SMILES from the <code>molecule</code>. This method
     * canonicaly lables the molecule but dose not perform any checks on the
     * chemical validity of the molecule. Does not care about multiple
     * molecules. IMPORTANT: A precomputed Set of All Rings (SAR) can be passed
     * to this SmilesGenerator in order to avoid recomputing it. Use setRings()
     * to assign the SAR.
     *
     * @param molecule The molecule to evaluate.
     * @param chiral true=SMILES will be chiral, false=SMILES will not be
     * chiral.
     * @param doubleBondConfiguration Should E/Z configurations be read at these
     * positions? If the flag at position X is set to true, an E/Z configuration
     * will be written from coordinates around bond X, if false, it will be
     * ignored. If flag is true for a bond which does not constitute a valid
     * double bond configuration, it will be ignored (meaning setting all to
     * true will create E/Z indication will be pu in the smiles wherever
     * possible, but note the coordinates might be arbitrary).
     * @exception CDKException At least one atom has no Point2D; coordinates are
     * needed for creating the chiral smiles. This excpetion can only be thrown
     * if chiral smiles is created, ignore it if you want a non-chiral smiles
     * (createSMILES(GraphAtomContainer) does not throw an exception).
     * @see
     * org.openscience.cdk.graph.invariant.CanonicalLabeler#canonLabel(IAtomContainer)
     * @return the SMILES representation of the molecule
     */
    public synchronized String createSMILESWithoutCheckForMultipleMolecules(IAtomContainer molecule, boolean chiral, boolean doubleBondConfiguration[]) throws CDKException {

        if (molecule.getAtomCount() == 0) {
            return "";
        }
        int[] canonicalLabels = labeller.getCanonicalPermutation(molecule);
        brokenBonds.clear();
        ringMarker = 0;
        IAtom start = null;
        for (int i = 0; i < molecule.getAtomCount(); i++) {
            IAtom atom = molecule.getAtom(i);
            if (chiral && atom.getPoint2d() == null) {
                throw new CDKException("Atom number " + i + " has no 2D coordinates, but 2D coordinates are needed for creating chiral smiles");
            }
            //LOGGER.debug("Setting all VISITED flags to false");
            atom.setFlag(VISITED, false);
            // set the start to the atom labelled '0'
            if (canonicalLabels[i] == 0) {
                start = atom;
            }
        }

        //detect aromaticity
        if (useAromaticityFlag || chiral) {
            if (rings == null) {
                if (ringFinder == null) {
                    ringFinder = new AllRingsFinder();
                }
                rings = ringFinder.findAllRings(molecule);
            }
            percieveAtomTypesAndConfigureAtoms(molecule);
            aromatizeDayLight(molecule);
        }
        if (chiral && rings.getAtomContainerCount() > 0) {
            List v = partitionRings(rings);
            //LOGGER.debug("RingSystems: " + v.size());
            for (int i = 0; i < v.size(); i++) {
                int counter = 0;
                Iterator<IAtomContainer> containers = getAllAtomContainers((IRingSet) v.get(i)).iterator();
                while (containers.hasNext()) {
                    IAtomContainer allrings = containers.next();
                    for (int k = 0; k < allrings.getAtomCount(); k++) {
                        if (!isStereo(molecule, allrings.getAtom(k)) && hasWedges(molecule, allrings.getAtom(k)) != null) {
                            IBond bond = molecule.getBond(allrings.getAtom(k), hasWedges(molecule, allrings.getAtom(k)));
                            if (bond.getStereo() == IBond.Stereo.UP) {
                                allrings.getAtom(k).setProperty(RING_CONFIG, UP);
                            } else {
                                allrings.getAtom(k).setProperty(RING_CONFIG, DOWN);
                            }
                            counter++;
                        }
                    }
                    if (counter == 1) {
                        for (int k = 0; k < allrings.getAtomCount(); k++) {
                            IBond bond = molecule.getBond(allrings.getAtom(k), hasWedges(molecule, allrings.getAtom(k)));
                            if (bond != null) {
                                if (bond.getStereo() == IBond.Stereo.UP) {
                                    allrings.getAtom(k).setProperty(RING_CONFIG, UP);
                                } else {
                                    allrings.getAtom(k).setProperty(RING_CONFIG, DOWN);
                                }
                            }
                        }
                    }
                }
            }
        }

        StringBuffer l = new StringBuffer();
        createSMILES(start, l, molecule, chiral, doubleBondConfiguration, canonicalLabels, useAromaticityFlag);
        rings = null;

        // remove all CanonicalLable/InvariancePair props
        for (int k = 0; k < molecule.getAtomCount(); k++) {
            molecule.getAtom(k).removeProperty("CanonicalLable");
            molecule.getAtom(k).removeProperty("InvariancePair");
        }

        return l.toString();
    }

    private IAtom hasWedges(IAtomContainer ac, IAtom a) {
        List<IAtom> atoms = ac.getConnectedAtomsList(a);
        //      for (int i = 0; i < atoms.size(); i++)
//      {
//          atomi = (IAtom)atoms.get(i);
//          if (ac.getBond(a, atomi).getStereo() != IBond.Stereo.NONE && !atomi.getSymbol().equals("H"))
//          {
//              return (atomi);
//          }
//      }
        for (IAtom atom : atoms) {
            if (ac.getBond(a, atom).getStereo() != NONE) {
                return (atom);
            }
        }
        return (null);
    }

    /**
     * Says if an atom is the end of a double bond configuration
     *
     * @param atom The atom which is the end of configuration
     * @param container The atomContainer the atom is in
     * @param parent The atom we came from
     * @param doubleBondConfiguration The array indicating where double bond
     * configurations are specified (this method ensures that there is actually
     * the possibility of a double bond configuration)
     * @return false=is not end of configuration, true=is
     */
    private boolean isEndOfDoubleBond(IAtomContainer container, IAtom atom, IAtom parent, boolean[] doubleBondConfiguration) {
        IBond bond = container.getBond(atom, parent);
        if (bond != null
                || doubleBondConfiguration.length <= container.indexOf(bond)
                || !doubleBondConfiguration[container.indexOf(bond)]) {
            return false;
        }
        // TO-DO: We make the silent assumption of unset hydrogen count equals zero hydrogen count here.
        int lengthAtom = container.getConnectedBondsCount(atom) + ((Objects.equals(atom.getImplicitHydrogenCount(), UNSET)) ? 0 : atom.getImplicitHydrogenCount());
        // TO-DO: We make the silent assumption of unset hydrogen count equals zero hydrogen count here.
        int lengthParent = container.getConnectedBondsCount(parent) + ((Objects.equals(parent.getImplicitHydrogenCount(), UNSET)) ? 0 : parent.getImplicitHydrogenCount());
        if (container.getBond(atom, parent) != null) {
            if (container.getBond(atom, parent).getOrder() == IBond.Order.DOUBLE
                    && (lengthAtom == 3 || (lengthAtom == 2 && atom.getSymbol().equals("N")))
                    && (lengthParent == 3 || (lengthParent == 2 && parent.getSymbol().equals("N")))) {
                List<IAtom> atoms = container.getConnectedAtomsList(atom);
                IAtom one = null;
                IAtom two = null;
                IAtom atomi = null;
                for (int i = 0; i < atoms.size(); i++) {
                    atomi = container.getAtom(i);
                    if (atomi != parent && one == null) {
                        one = atomi;
                    } else if (atomi != parent && one != null) {
                        two = atomi;
                    }
                }
                String[] morgannumbers = getMorganNumbersWithElementSymbol(container);
                if ((one != null && two == null && atom.getSymbol().equals("N") && abs(giveAngleBothMethods(parent, atom, one, true)) > PI / 10) || (!atom.getSymbol().equals("N") && one != null && two != null && !morgannumbers[container.indexOf(one)].equals(morgannumbers[container.indexOf(two)]))) {
                    return (true);
                } else {
                    return (false);
                }
            }
        }
        return (false);
    }

    /**
     * Says if an atom is the start of a double bond configuration
     *
     * @param a The atom which is the start of configuration
     * @param container The atomContainer the atom is in
     * @param parent The atom we came from
     * @param doubleBondConfiguration The array indicating where double bond
     * configurations are specified (this method ensures that there is actually
     * the possibility of a double bond configuration)
     * @return false=is not start of configuration, true=is
     */
    private boolean isStartOfDoubleBond(IAtomContainer container, IAtom a, IAtom parent, boolean[] doubleBondConfiguration) {
        // TO-DO: We make the silent assumption of unset hydrogen count equals zero hydrogen count here.
        int lengthAtom = container.getConnectedBondsCount(a) + ((Objects.equals(a.getImplicitHydrogenCount(), UNSET)) ? 0 : a.getImplicitHydrogenCount());
        if (lengthAtom != 3 && (lengthAtom != 2 && !a.getSymbol().equals("N"))) {
            return (false);
        }
        List<IAtom> atoms = container.getConnectedAtomsList(a);
        IAtom one = null;
        IAtom two = null;
        boolean doubleBond = false;
        IAtom nextAtom = null;
        for (IAtom atomi : atoms) {
            if (atomi != parent && container.getBond(atomi, a).getOrder() == IBond.Order.DOUBLE
                    && isEndOfDoubleBond(container, atomi, a, doubleBondConfiguration)) {
                doubleBond = true;
                nextAtom = atomi;
            }
            if (atomi != nextAtom && one == null) {
                one = atomi;
            } else if (atomi != nextAtom && one != null) {
                two = atomi;
            }
        }
        String[] morgannumbers = getMorganNumbersWithElementSymbol(container);

        IBond bond = container.getBond(a, nextAtom);
        if (bond == null) {
            return false;
        }
        if (one != null && ((!a.getSymbol().equals("N") && two != null
                && !morgannumbers[container.indexOf(one)].equals(morgannumbers[container.indexOf(two)])
                && doubleBond && doubleBondConfiguration[container.indexOf(bond)])
                || (doubleBond && a.getSymbol().equals("N") && abs(giveAngleBothMethods(nextAtom, a, parent, true)) > PI / 10))) {
            return (true);
        } else {
            return (false);
        }
    }

    /**
     * Gets the bondBroken attribute of the SmilesGenerator object
     */
    private boolean isBondBroken(IAtom a1, IAtom a2) {
        for (BrokenBond bond : brokenBonds) {
            if ((bond.getA1().equals(a1) || bond.getA1().equals(a2)) && (bond.getA2().equals(a1) || bond.getA2().equals(a2))) {
                return (true);
            }
        }
        return false;
    }

    /**
     * Determines if the atom <code>a</code> is a atom with a ring marker.
     *
     * @param a the atom to test
     * @return true if the atom participates in a bond that was broken in the
     * first pass.
     */
//  private boolean isRingOpening(IAtom a)
//  {
//      Iterator it = brokenBonds.iterator();
//      while (it.hasNext())
//      {
//          BrokenBond bond = (BrokenBond) it.next();
//          if (bond.getA1().equals(a) || bond.getA2().equals(a))
//          {
//              return true;
//          }
//      }
//      return false;
//  }
    /**
     * Determines if the atom <code>a</code> is a atom with a ring marker.
     *
     * @return true if the atom participates in a bond that was broken in the
     * first pass.
     */
    private boolean isRingOpening(IAtom a1, List v) {
        return brokenBonds.stream().anyMatch((BrokenBond bond) -> v.stream().anyMatch((aV) -> ((bond.getA1().equals(a1) && bond.getA2().equals(aV)) || (bond.getA1().equals(aV) && bond.getA2().equals(a1)))));
    }

    /**
     * Return the neighbours of atom <code>a</code> in canonical order with the
     * atoms that have high bond order at the front.
     *
     * @param a the atom whose neighbours are to be found.
     * @param container the GraphAtomContainer that is being parsed.
     * @return Vector of atoms in canonical oreder.
     */
    private List getCanNeigh(final IAtom a, final IAtomContainer container, final int[] canonicalLabels) {
        List<IAtom> v = container.getConnectedAtomsList(a);
        if (v.size() > 1) {
            sort(v, (IAtom a1, IAtom a2) -> {
                int l1 = canonicalLabels[container.indexOf(a1)];
                int l2 = canonicalLabels[container.indexOf(a2)];
                if (l1 < l2) {
                    return -1;
                } else if (l1 > l2) {
                    return 1;
                } else {
                    return 0;
                }
            });
        }
        return v;
    }

    /**
     * Gets the ringOpenings attribute of the SmilesGenerator object
     */
    private List getRingOpenings(IAtom a, List vbonds) {
        Iterator it = brokenBonds.iterator();
        List v = new Vector(10);
        while (it.hasNext()) {
            BrokenBond bond = (BrokenBond) it.next();
            if (bond.getA1().equals(a) || bond.getA2().equals(a)) {
                v.add(bond.getMarker());
                if (vbonds != null) {
                    vbonds.add(bond.getA1().equals(a) ? bond.getA2() : bond.getA1());
                }
            }
        }
        sort(v);
        return v;
    }

    /**
     * Returns true if the <code>atom</code> in the <code>container</code> has
     * been marked as a chiral center by the user.
     */
//  private boolean isChiralCenter(IAtom atom, IAtomContainer container)
//  {
//      IBond[] bonds = container.getConnectedBonds(atom);
//      for (int i = 0; i < bonds.length; i++)
//      {
//          IBond bond = bonds[i];
//          int stereo = bond.getStereo();
//          if (stereo == IBond.Stereo.DOWN ||
//                  stereo == IBond.Stereo.UP)
//          {
//              return true;
//          }
//      }
//      return false;
//  }
    /**
     * Gets the last atom object (not Vector) in a Vector as created by
     * createDSFTree.
     *
     * @param v The Vector
     * @param result The feature to be added to the Atoms attribute
     */
    private void addAtoms(List v, List result) {
        v.stream().forEach((aV) -> {
            if (aV instanceof IAtom) {
                result.add(aV);
            } else {
                addAtoms((List) aV, result);
            }
        });
    }

    /**
     * Performes a DFS search on the <code>atomContainer</code>. Then parses the
     * resulting tree to create the SMILES string.
     *
     * @param a the atom to start the search at.
     * @param line the StringBuffer that the SMILES is to be appended to.
     * @param chiral true=SMILES will be chiral, false=SMILES will not be
     * chiral.
     * @param doubleBondConfiguration Should E/Z configurations be read at these
     * positions? If the flag at position X is set to true, an E/Z configuration
     * will be written from coordinates around bond X, if false, it will be
     * ignored. If flag is true for a bond which does not constitute a valid
     * double bond configuration, it will be ignored (meaning setting all to
     * true will create E/Z indication will be pu in the smiles wherever
     * possible, but note the coordinates might be arbitrary).
     * @param atomContainer the GraphAtomContainer that the SMILES string is
     * generated for.
     * @param useAromaticity true=aromaticity or sp2 will trigger lower case
     * letters, wrong=only sp2
     */
    private void createSMILES(IAtom a, StringBuffer line, IAtomContainer atomContainer, boolean chiral, boolean[] doubleBondConfiguration, int[] canonicalLabels, boolean useAromaticity) {
        List tree = new Vector();

        // set all ISVISITED labels to FALSE
        Iterator atoms = atomContainer.atoms().iterator();
        while (atoms.hasNext()) {
            ((IChemObject) atoms.next()).setFlag(VISITED, false);
        }

        createDFSTree(a, tree, null, atomContainer, canonicalLabels);
        //LOGGER.debug("Done with tree");

        parseChain(tree, line, atomContainer, null, chiral, doubleBondConfiguration, new Vector(), useAromaticity);
    }

    /**
     * Recursively perform a DFS search on the <code>container</code> placing
     * atoms and branches in the vector <code>tree</code>.
     *
     * @param a the atom being visited.
     * @param tree vector holding the tree.
     * @param parent the atom we came from.
     * @param container the GraphAtomContainer that we are parsing.
     */
    private void createDFSTree(IAtom a, List tree, IAtom parent, IAtomContainer container, int[] canonicalLabels) {
        tree.add(a);
        List neighbours = new ArrayList(getCanNeigh(a, container, canonicalLabels));
        neighbours.remove(parent);
        IAtom next;
        a.setFlag(VISITED, true);
        //LOGGER.debug("Starting with DFSTree and GraphAtomContainer of size " + container.getAtomCount());
        //LOGGER.debug("Current Atom has " + neighbours.size() + " neighbours");
        Iterator iter = neighbours.iterator();
        while (iter.hasNext()) {
            next = (IAtom) iter.next();
            if (!next.getFlag(VISITED)) {
                if (!iter.hasNext()) {
                    //Last neighbour therefore in this chain
                    createDFSTree(next, tree, a, container, canonicalLabels);
                } else {
                    List branch = new Vector();
                    tree.add(branch);
                    //LOGGER.debug("adding branch");
                    createDFSTree(next, branch, a, container, canonicalLabels);
                }
            } else {
                //Found ring closure between next and a
                //LOGGER.debug("found ringclosure in DFTTreeCreation");
                ringMarker++;
                BrokenBond bond = new BrokenBond(a, next, ringMarker);
                if (!brokenBonds.contains(bond)) {
                    brokenBonds.add(bond);
                } else {
                    ringMarker--;
                }
            }
        }
    }

    /**
     * Parse a branch
     */
    private void parseChain(List v, StringBuffer buffer, IAtomContainer container, IAtom parent, boolean chiral, boolean[] doubleBondConfiguration, List atomsInOrderOfSmiles, boolean useAromaticity) {
        int positionInVector = 0;
        IAtom atom;
        //LOGGER.debug("in parse chain. Size of tree: " + v.size());
        for (int h = 0; h < v.size(); h++) {
            Object o = v.get(h);
            if (o instanceof IAtom) {
                atom = (IAtom) o;
                if (parent != null) {
                    parseBond(buffer, atom, parent, container, useAromaticity);
                } else if (chiral && isStereo(container, atom)) {
                    parent = (IAtom) ((List) v.get(1)).get(0);
                }
                parseAtom(atom, buffer, container, chiral, doubleBondConfiguration, parent, atomsInOrderOfSmiles, v, useAromaticity);
                //LOGGER.debug("in parseChain after parseAtom()");
                /*
                 * The principle of making chiral smiles is quite simple, although the code is pretty uggly. The Atoms
                 * connected to the chiral center are put in sorted[] in the order they have to appear in the smiles.
                 * Then the Vector v is rearranged according to sorted[]
                 */
                if (chiral && isStereo(container, atom) && container.getBond(parent, atom) != null) {
                    //LOGGER.debug("in parseChain in isChiral");
                    IAtom[] sorted = null;
                    List chiralNeighbours = container.getConnectedAtomsList(atom);
                    if (isTetrahedral(container, atom, false) > 0) {
                        sorted = new IAtom[3];
                    }
                    if (BondTools.isTetrahedral(container, atom, false) == 1) {
                        if (container.getBond(parent, atom).getStereo() == IBond.Stereo.DOWN) {
                            for (int i = 0; i < chiralNeighbours.size(); i++) {
                                if (chiralNeighbours.get(i) != parent) {
                                    if ((container.getBond((IAtom) chiralNeighbours.get(i), atom).getStereo() == null
                                            || container.getBond((IAtom) chiralNeighbours.get(i), atom).getStereo() == NONE)
                                            && isLeft(((IAtom) chiralNeighbours.get(i)), parent, atom) && !isBondBroken((IAtom) chiralNeighbours.get(i), atom)) {
                                        sorted[2] = (IAtom) chiralNeighbours.get(i);
                                    }
                                    if ((container.getBond((IAtom) chiralNeighbours.get(i), atom).getStereo() == null
                                            || container.getBond((IAtom) chiralNeighbours.get(i), atom).getStereo() == NONE)
                                            && !isLeft(((IAtom) chiralNeighbours.get(i)), parent, atom) && !isBondBroken((IAtom) chiralNeighbours.get(i), atom)) {
                                        sorted[1] = (IAtom) chiralNeighbours.get(i);
                                    }
                                    if (container.getBond((IAtom) chiralNeighbours.get(i), atom).getStereo() == IBond.Stereo.UP && !isBondBroken((IAtom) chiralNeighbours.get(i), atom)) {
                                        sorted[0] = (IAtom) chiralNeighbours.get(i);
                                    }
                                }
                            }
                        }
                        if (container.getBond(parent, atom).getStereo() == IBond.Stereo.UP) {
                            for (int i = 0; i < chiralNeighbours.size(); i++) {
                                if (chiralNeighbours.get(i) != parent) {
                                    if ((container.getBond((IAtom) chiralNeighbours.get(i), atom).getStereo() == null
                                            || container.getBond((IAtom) chiralNeighbours.get(i), atom).getStereo() == NONE)
                                            && isLeft(((IAtom) chiralNeighbours.get(i)), parent, atom) && !isBondBroken((IAtom) chiralNeighbours.get(i), atom)) {
                                        sorted[1] = (IAtom) chiralNeighbours.get(i);
                                    }
                                    if ((container.getBond((IAtom) chiralNeighbours.get(i), atom).getStereo() == null
                                            || container.getBond((IAtom) chiralNeighbours.get(i), atom).getStereo() == NONE)
                                            && !isLeft(((IAtom) chiralNeighbours.get(i)), parent, atom) && !isBondBroken((IAtom) chiralNeighbours.get(i), atom)) {
                                        sorted[2] = (IAtom) chiralNeighbours.get(i);
                                    }
                                    if (container.getBond((IAtom) chiralNeighbours.get(i), atom).getStereo() == IBond.Stereo.DOWN && !isBondBroken((IAtom) chiralNeighbours.get(i), atom)) {
                                        sorted[0] = (IAtom) chiralNeighbours.get(i);
                                    }
                                }
                            }
                        }
                        if (container.getBond(parent, atom).getStereo() == UNSET || container.getBond(parent, atom).getStereo() == NONE) {
                            boolean normalBindingIsLeft = false;
                            for (int i = 0; i < chiralNeighbours.size(); i++) {
                                if (chiralNeighbours.get(i) != parent) {
                                    if (container.getBond((IAtom) chiralNeighbours.get(i), atom).getStereo() == null
                                            || container.getBond((IAtom) chiralNeighbours.get(i), atom).getStereo() == NONE) {
                                        if (isLeft(((IAtom) chiralNeighbours.get(i)), parent, atom)) {
                                            normalBindingIsLeft = true;
                                            break;
                                        }
                                    }
                                }
                            }
                            for (int i = 0; i < chiralNeighbours.size(); i++) {
                                if (chiralNeighbours.get(i) != parent) {
                                    if (normalBindingIsLeft) {
                                        if (container.getBond((IAtom) chiralNeighbours.get(i), atom).getStereo() == null
                                                || container.getBond((IAtom) chiralNeighbours.get(i), atom).getStereo() == NONE) {
                                            sorted[0] = (IAtom) chiralNeighbours.get(i);
                                        }
                                        if (container.getBond((IAtom) chiralNeighbours.get(i), atom).getStereo() == IBond.Stereo.UP) {
                                            sorted[2] = (IAtom) chiralNeighbours.get(i);
                                        }
                                        if (container.getBond((IAtom) chiralNeighbours.get(i), atom).getStereo() == IBond.Stereo.DOWN) {
                                            sorted[1] = (IAtom) chiralNeighbours.get(i);
                                        }
                                    } else {
                                        if (container.getBond((IAtom) chiralNeighbours.get(i), atom).getStereo() == IBond.Stereo.UP) {
                                            sorted[1] = (IAtom) chiralNeighbours.get(i);
                                        }
                                        if (container.getBond((IAtom) chiralNeighbours.get(i), atom).getStereo() == null
                                                || container.getBond((IAtom) chiralNeighbours.get(i), atom).getStereo() == NONE) {
                                            sorted[0] = (IAtom) chiralNeighbours.get(i);
                                        }
                                        if (container.getBond((IAtom) chiralNeighbours.get(i), atom).getStereo() == IBond.Stereo.DOWN) {
                                            sorted[2] = (IAtom) chiralNeighbours.get(i);
                                        }
                                    }
                                }
                            }
                        }
                    }
                    if (BondTools.isTetrahedral(container, atom, false) == 2) {
                        if (container.getBond(parent, atom).getStereo() == IBond.Stereo.UP) {
                            for (int i = 0; i < chiralNeighbours.size(); i++) {
                                if (chiralNeighbours.get(i) != parent) {
                                    if (container.getBond((IAtom) chiralNeighbours.get(i), atom).getStereo() == IBond.Stereo.DOWN && BondTools.isLeft(((IAtom) chiralNeighbours.get(i)), parent, atom) && !isBondBroken((IAtom) chiralNeighbours.get(i), atom)) {
                                        sorted[1] = (IAtom) chiralNeighbours.get(i);
                                    }
                                    if (container.getBond((IAtom) chiralNeighbours.get(i), atom).getStereo() == IBond.Stereo.DOWN && !BondTools.isLeft(((IAtom) chiralNeighbours.get(i)), parent, atom) && !isBondBroken((IAtom) chiralNeighbours.get(i), atom)) {
                                        sorted[2] = (IAtom) chiralNeighbours.get(i);
                                    }
                                    if (container.getBond((IAtom) chiralNeighbours.get(i), atom).getStereo() == IBond.Stereo.UP && !isBondBroken((IAtom) chiralNeighbours.get(i), atom)) {
                                        sorted[0] = (IAtom) chiralNeighbours.get(i);
                                    }
                                }
                            }
                        }
                        if (container.getBond(parent, atom).getStereo() == IBond.Stereo.DOWN) {
                            double angle1 = 0;
                            double angle2 = 0;
                            IAtom atom1 = null;
                            IAtom atom2 = null;
                            for (int i = 0; i < chiralNeighbours.size(); i++) {
                                if (chiralNeighbours.get(i) != parent) {
                                    if (container.getBond((IAtom) chiralNeighbours.get(i), atom).getStereo() == IBond.Stereo.UP && !isBondBroken((IAtom) chiralNeighbours.get(i), atom)) {
                                        if (angle1 == 0) {
                                            angle1 = giveAngle(atom, parent, (IAtom) chiralNeighbours.get(i));
                                            atom1 = (IAtom) chiralNeighbours.get(i);
                                        } else {
                                            angle2 = giveAngle(atom, parent, (IAtom) chiralNeighbours.get(i));
                                            atom2 = (IAtom) chiralNeighbours.get(i);
                                        }
                                    }
                                    if (container.getBond((IAtom) chiralNeighbours.get(i), atom).getStereo() == IBond.Stereo.DOWN && !isBondBroken((IAtom) chiralNeighbours.get(i), atom)) {
                                        sorted[1] = (IAtom) chiralNeighbours.get(i);
                                    }
                                }
                            }
                            if (angle1 < angle2) {
                                sorted[0] = atom2;
                                sorted[2] = atom1;
                            } else {
                                sorted[0] = atom1;
                                sorted[2] = atom2;
                            }
                        }
                    }
                    if (BondTools.isTetrahedral(container, atom, false) == 3) {
                        if (container.getBond(parent, atom).getStereo() == IBond.Stereo.UP) {
                            TreeMap hm = new TreeMap();
                            for (int i = 0; i < chiralNeighbours.size(); i++) {
                                if (chiralNeighbours.get(i) != parent && !isBondBroken((IAtom) chiralNeighbours.get(i), atom)) {
                                    hm.put(giveAngle(atom, parent, ((IAtom) chiralNeighbours.get(i))), i);
                                }
                            }
                            Object[] ohere = hm.values().toArray();
                            for (int i = ohere.length - 1; i > -1; i--) {
                                sorted[i] = ((IAtom) chiralNeighbours.get(((Number) ohere[i]).intValue()));
                            }
                        }
                        if (container.getBond(parent, atom).getStereo() == null
                                || container.getBond(parent, atom).getStereo() == NONE) {
                            double angle1 = 0;
                            double angle2 = 0;
                            IAtom atom1 = null;
                            IAtom atom2 = null;
                            for (int i = 0; i < chiralNeighbours.size(); i++) {
                                if (chiralNeighbours.get(i) != parent) {
                                    if ((container.getBond((IAtom) chiralNeighbours.get(i), atom).getStereo() == null
                                            || container.getBond((IAtom) chiralNeighbours.get(i), atom).getStereo() == NONE)
                                            && !isBondBroken((IAtom) chiralNeighbours.get(i), atom)) {
                                        if (angle1 == 0) {
                                            angle1 = giveAngle(atom, parent, (IAtom) chiralNeighbours.get(i));
                                            atom1 = (IAtom) chiralNeighbours.get(i);
                                        } else {
                                            angle2 = giveAngle(atom, parent, (IAtom) chiralNeighbours.get(i));
                                            atom2 = (IAtom) chiralNeighbours.get(i);
                                        }
                                    }
                                    if (container.getBond((IAtom) chiralNeighbours.get(i), atom).getStereo() == IBond.Stereo.UP && !isBondBroken((IAtom) chiralNeighbours.get(i), atom)) {
                                        sorted[0] = (IAtom) chiralNeighbours.get(i);
                                    }
                                }
                            }
                            if (angle1 < angle2) {
                                sorted[1] = atom2;
                                sorted[2] = atom1;
                            } else {
                                sorted[1] = atom1;
                                sorted[2] = atom2;
                            }
                        }
                    }
                    if (BondTools.isTetrahedral(container, atom, false) == 4) {
                        if (container.getBond(parent, atom).getStereo() == IBond.Stereo.DOWN) {
                            TreeMap hm = new TreeMap();
                            for (int i = 0; i < chiralNeighbours.size(); i++) {
                                if (chiralNeighbours.get(i) != parent && !isBondBroken((IAtom) chiralNeighbours.get(i), atom)) {
                                    hm.put(giveAngle(atom, parent, ((IAtom) chiralNeighbours.get(i))), i);
                                }
                            }
                            Object[] ohere = hm.values().toArray();
                            for (int i = ohere.length - 1; i > -1; i--) {
                                sorted[i] = ((IAtom) chiralNeighbours.get(((Number) ohere[i]).intValue()));
                            }
                        }
                        if (container.getBond(parent, atom).getStereo() == null
                                || container.getBond(parent, atom).getStereo() == NONE) {
                            double angle1 = 0;
                            double angle2 = 0;
                            IAtom atom1 = null;
                            IAtom atom2 = null;
                            for (int i = 0; i < chiralNeighbours.size(); i++) {
                                if (chiralNeighbours.get(i) != parent) {
                                    if ((container.getBond((IAtom) chiralNeighbours.get(i), atom).getStereo() == null
                                            || container.getBond((IAtom) chiralNeighbours.get(i), atom).getStereo() == NONE)
                                            && !isBondBroken((IAtom) chiralNeighbours.get(i), atom)) {
                                        if (angle1 == 0) {
                                            angle1 = giveAngle(atom, parent, (IAtom) chiralNeighbours.get(i));
                                            atom1 = (IAtom) chiralNeighbours.get(i);
                                        } else {
                                            angle2 = giveAngle(atom, parent, (IAtom) chiralNeighbours.get(i));
                                            atom2 = (IAtom) chiralNeighbours.get(i);
                                        }
                                    }
                                    if (container.getBond((IAtom) chiralNeighbours.get(i), atom).getStereo() == IBond.Stereo.DOWN && !isBondBroken((IAtom) chiralNeighbours.get(i), atom)) {
                                        sorted[2] = (IAtom) chiralNeighbours.get(i);
                                    }
                                }
                            }
                            if (angle1 < angle2) {
                                sorted[1] = atom2;
                                sorted[0] = atom1;
                            } else {
                                sorted[1] = atom1;
                                sorted[0] = atom2;
                            }
                        }
                    }
                    if (BondTools.isTetrahedral(container, atom, false) == 5) {
                        if (container.getBond(parent, atom).getStereo() == IBond.Stereo.DOWN) {
                            for (int i = 0; i < chiralNeighbours.size(); i++) {
                                if (chiralNeighbours.get(i) != parent) {
                                    if (container.getBond((IAtom) chiralNeighbours.get(i), atom).getStereo() == IBond.Stereo.UP) {
                                        sorted[0] = (IAtom) chiralNeighbours.get(i);
                                    }
                                    if (container.getBond((IAtom) chiralNeighbours.get(i), atom).getStereo() == null
                                            || container.getBond((IAtom) chiralNeighbours.get(i), atom).getStereo() == NONE) {
                                        sorted[2] = (IAtom) chiralNeighbours.get(i);
                                    }
                                    if (container.getBond((IAtom) chiralNeighbours.get(i), atom).getStereo() == IBond.Stereo.DOWN) {
                                        sorted[1] = (IAtom) chiralNeighbours.get(i);
                                    }
                                }
                            }
                        }
                        if (container.getBond(parent, atom).getStereo() == IBond.Stereo.UP) {
                            for (int i = 0; i < chiralNeighbours.size(); i++) {
                                if (chiralNeighbours.get(i) != parent) {
                                    if (container.getBond((IAtom) chiralNeighbours.get(i), atom).getStereo() == IBond.Stereo.DOWN && BondTools.isLeft(((IAtom) chiralNeighbours.get(i)), parent, atom) && !isBondBroken((IAtom) chiralNeighbours.get(i), atom)) {
                                        sorted[0] = (IAtom) chiralNeighbours.get(i);
                                    }
                                    if (container.getBond((IAtom) chiralNeighbours.get(i), atom).getStereo() == IBond.Stereo.DOWN && !BondTools.isLeft(((IAtom) chiralNeighbours.get(i)), parent, atom) && !isBondBroken((IAtom) chiralNeighbours.get(i), atom)) {
                                        sorted[2] = (IAtom) chiralNeighbours.get(i);
                                    }
                                    if (container.getBond((IAtom) chiralNeighbours.get(i), atom).getStereo() == null
                                            || container.getBond((IAtom) chiralNeighbours.get(i), atom).getStereo() == NONE) {
                                        sorted[1] = (IAtom) chiralNeighbours.get(i);
                                    }
                                }
                            }
                        }
                        if (container.getBond(parent, atom).getStereo() == UNSET || container.getBond(parent, atom).getStereo() == NONE) {
                            for (int i = 0; i < chiralNeighbours.size(); i++) {
                                if (chiralNeighbours.get(i) != parent) {
                                    if (container.getBond((IAtom) chiralNeighbours.get(i), atom).getStereo() == IBond.Stereo.DOWN && BondTools.isLeft(((IAtom) chiralNeighbours.get(i)), parent, atom) && !isBondBroken((IAtom) chiralNeighbours.get(i), atom)) {
                                        sorted[0] = (IAtom) chiralNeighbours.get(i);
                                    }
                                    if (container.getBond((IAtom) chiralNeighbours.get(i), atom).getStereo() == IBond.Stereo.DOWN && !BondTools.isLeft(((IAtom) chiralNeighbours.get(i)), parent, atom) && !isBondBroken((IAtom) chiralNeighbours.get(i), atom)) {
                                        sorted[2] = (IAtom) chiralNeighbours.get(i);
                                    }
                                    if (container.getBond((IAtom) chiralNeighbours.get(i), atom).getStereo() == IBond.Stereo.UP) {
                                        sorted[1] = (IAtom) chiralNeighbours.get(i);
                                    }
                                }
                            }
                        }
                    }
                    if (BondTools.isTetrahedral(container, atom, false) == 6) {
                        if (container.getBond(parent, atom).getStereo() == IBond.Stereo.UP) {
                            for (int i = 0; i < chiralNeighbours.size(); i++) {
                                if (chiralNeighbours.get(i) != parent) {
                                    if (container.getBond((IAtom) chiralNeighbours.get(i), atom).getStereo() == IBond.Stereo.UP) {
                                        sorted[0] = (IAtom) chiralNeighbours.get(i);
                                    }
                                    if (container.getBond((IAtom) chiralNeighbours.get(i), atom).getStereo() == null
                                            || container.getBond((IAtom) chiralNeighbours.get(i), atom).getStereo() == NONE) {
                                        sorted[2] = (IAtom) chiralNeighbours.get(i);
                                    }
                                    if (container.getBond((IAtom) chiralNeighbours.get(i), atom).getStereo() == IBond.Stereo.DOWN) {
                                        sorted[1] = (IAtom) chiralNeighbours.get(i);
                                    }
                                }
                            }
                        }
                        if (container.getBond(parent, atom).getStereo() == IBond.Stereo.DOWN) {
                            for (int i = 0; i < chiralNeighbours.size(); i++) {
                                if (chiralNeighbours.get(i) != parent) {
                                    if (container.getBond((IAtom) chiralNeighbours.get(i), atom).getStereo() == IBond.Stereo.UP && BondTools.isLeft(((IAtom) chiralNeighbours.get(i)), parent, atom) && !isBondBroken((IAtom) chiralNeighbours.get(i), atom)) {
                                        sorted[2] = (IAtom) chiralNeighbours.get(i);
                                    }
                                    if (container.getBond((IAtom) chiralNeighbours.get(i), atom).getStereo() == IBond.Stereo.UP && !BondTools.isLeft(((IAtom) chiralNeighbours.get(i)), parent, atom) && !isBondBroken((IAtom) chiralNeighbours.get(i), atom)) {
                                        sorted[0] = (IAtom) chiralNeighbours.get(i);
                                    }
                                    if (container.getBond((IAtom) chiralNeighbours.get(i), atom).getStereo() == null
                                            || container.getBond((IAtom) chiralNeighbours.get(i), atom).getStereo() == NONE) {
                                        sorted[1] = (IAtom) chiralNeighbours.get(i);
                                    }
                                }
                            }
                        }
                        if (container.getBond(parent, atom).getStereo() == UNSET || container.getBond(parent, atom).getStereo() == NONE) {
                            for (int i = 0; i < chiralNeighbours.size(); i++) {
                                if (chiralNeighbours.get(i) != parent) {
                                    if (container.getBond((IAtom) chiralNeighbours.get(i), atom).getStereo() == IBond.Stereo.UP && BondTools.isLeft(((IAtom) chiralNeighbours.get(i)), parent, atom) && !isBondBroken((IAtom) chiralNeighbours.get(i), atom)) {
                                        sorted[2] = (IAtom) chiralNeighbours.get(i);
                                    }
                                    if (container.getBond((IAtom) chiralNeighbours.get(i), atom).getStereo() == IBond.Stereo.UP && !BondTools.isLeft(((IAtom) chiralNeighbours.get(i)), parent, atom) && !isBondBroken((IAtom) chiralNeighbours.get(i), atom)) {
                                        sorted[0] = (IAtom) chiralNeighbours.get(i);
                                    }
                                    if (container.getBond((IAtom) chiralNeighbours.get(i), atom).getStereo() == IBond.Stereo.DOWN) {
                                        sorted[1] = (IAtom) chiralNeighbours.get(i);
                                    }
                                }
                            }
                        }
                    }
                    if (isSquarePlanar(container, atom)) {
                        sorted = new IAtom[3];
                        //This produces a U=SP1 order in every case
                        TreeMap hm = new TreeMap();
                        for (int i = 0; i < chiralNeighbours.size(); i++) {
                            if (chiralNeighbours.get(i) != parent && !isBondBroken((IAtom) chiralNeighbours.get(i), atom)) {
                                hm.put(giveAngle(atom, parent, ((IAtom) chiralNeighbours.get(i))), i);
                            }
                        }
                        Object[] ohere = hm.values().toArray();
                        for (int i = 0; i < ohere.length; i++) {
                            sorted[i] = ((IAtom) chiralNeighbours.get(((Integer) ohere[i])));
                        }
                    }
                    if (isTrigonalBipyramidalOrOctahedral(container, atom) != 0) {
                        sorted = new IAtom[container.getConnectedBondsCount(atom) - 1];
                        TreeMap hm = new TreeMap();
                        if (container.getBond(parent, atom).getStereo() == IBond.Stereo.UP) {
                            for (int i = 0; i < chiralNeighbours.size(); i++) {
                                if (container.getBond(atom, (IAtom) chiralNeighbours.get(i)).getStereo() == null
                                        || container.getBond(atom, (IAtom) chiralNeighbours.get(i)).getStereo() == NONE) {
                                    hm.put(giveAngle(atom, parent, ((IAtom) chiralNeighbours.get(i))), i);
                                }
                                if (container.getBond(atom, (IAtom) chiralNeighbours.get(i)).getStereo() == IBond.Stereo.DOWN) {
                                    sorted[sorted.length - 1] = (IAtom) chiralNeighbours.get(i);
                                }
                            }
                            Object[] ohere = hm.values().toArray();
                            for (int i = 0; i < ohere.length; i++) {
                                sorted[i] = ((IAtom) chiralNeighbours.get(((Number) ohere[i]).intValue()));
                            }
                        }
                        if (container.getBond(parent, atom).getStereo() == IBond.Stereo.DOWN) {
                            for (int i = 0; i < chiralNeighbours.size(); i++) {
                                if (container.getBond(atom, (IAtom) chiralNeighbours.get(i)).getStereo() == null
                                        || container.getBond(atom, (IAtom) chiralNeighbours.get(i)).getStereo() == NONE) {
                                    hm.put(giveAngle(atom, parent, ((IAtom) chiralNeighbours.get(i))), i);
                                }
                                if (container.getBond(atom, (IAtom) chiralNeighbours.get(i)).getStereo() == IBond.Stereo.UP) {
                                    sorted[sorted.length - 1] = (IAtom) chiralNeighbours.get(i);
                                }
                            }
                            Object[] ohere = hm.values().toArray();
                            for (int i = 0; i < ohere.length; i++) {
                                sorted[i] = ((IAtom) chiralNeighbours.get(((Number) ohere[i]).intValue()));
                            }
                        }
                        if (container.getBond(parent, atom).getStereo() == null
                                || container.getBond(parent, atom).getStereo() == NONE) {
                            for (int i = 0; i < chiralNeighbours.size(); i++) {
                                if (chiralNeighbours.get(i) != parent) {
                                    if (container.getBond(atom, (IAtom) chiralNeighbours.get(i)).getStereo() == null
                                            || container.getBond(atom, (IAtom) chiralNeighbours.get(i)).getStereo() == NONE) {
                                        hm.put((giveAngleFromMiddle(atom, parent, (IAtom) chiralNeighbours.get(i))), i);
                                    }
                                    if (container.getBond(atom, (IAtom) chiralNeighbours.get(i)).getStereo() == IBond.Stereo.UP) {
                                        sorted[0] = (IAtom) chiralNeighbours.get(i);
                                    }
                                    if (container.getBond(atom, (IAtom) chiralNeighbours.get(i)).getStereo() == IBond.Stereo.DOWN) {
                                        sorted[sorted.length - 2] = (IAtom) chiralNeighbours.get(i);
                                    }
                                }
                            }
                            Object[] ohere = hm.values().toArray();
                            sorted[sorted.length - 1] = ((IAtom) chiralNeighbours.get(((Number) ohere[ohere.length - 1]).intValue()));
                            if (ohere.length == 2) {
                                sorted[sorted.length - 3] = ((IAtom) chiralNeighbours.get(((Number) ohere[0]).intValue()));
                                if (giveAngleFromMiddle(atom, parent, ((IAtom) chiralNeighbours.get(((Number) ohere[1]).intValue()))) < 0) {
                                    IAtom dummy = sorted[sorted.length - 2];
                                    sorted[sorted.length - 2] = sorted[0];
                                    sorted[0] = dummy;
                                }
                            }
                            if (ohere.length == 3) {
                                sorted[sorted.length - 3] = sorted[sorted.length - 2];
                                sorted[sorted.length - 2] = ((IAtom) chiralNeighbours.get(((Number) ohere[ohere.length - 2]).intValue()));
                                sorted[sorted.length - 4] = ((IAtom) chiralNeighbours.get(((Number) ohere[ohere.length - 3]).intValue()));
                            }
                        }
                    }
                    //This builds an onew[] containing the objects after the center of the chirality in the order given by sorted[]
                    if (sorted != null) {
                        int numberOfAtoms = 3;
                        if (isTrigonalBipyramidalOrOctahedral(container, atom) != 0) {
                            numberOfAtoms = container.getConnectedBondsCount(atom) - 1;
                        }
                        Object[] omy = new Object[numberOfAtoms];
                        Object[] onew = new Object[numberOfAtoms];
                        for (int k = getRingOpenings(atom, null).size(); k < numberOfAtoms; k++) {
                            if (positionInVector + 1 + k - getRingOpenings(atom, null).size() < v.size()) {
                                omy[k] = v.get(positionInVector + 1 + k - getRingOpenings(atom, null).size());
                            }
                        }
                        for (int k = 0; k < sorted.length; k++) {
                            if (sorted[k] != null) {
                                for (Object omy1 : omy) {
                                    if (omy1 instanceof IAtom) {
                                        if (omy1 == sorted[k]) {
                                            onew[k] = omy1;
                                        }
                                    } else if (omy1 == null) {
                                        onew[k] = null;
                                    } else if (((List) omy1).get(0) == sorted[k]) {
                                        onew[k] = omy1;
                                    }
                                }
                            } else {
                                onew[k] = null;
                            }
                        }
                        //This is a workaround for 3624.MOL.2 I don't have a better solution currently
                        boolean doubleentry = false;
                        for (int m = 0; m < onew.length; m++) {
                            for (int k = 0; k < onew.length; k++) {
                                if (m != k && onew[k] == onew[m]) {
                                    doubleentry = true;
                                }
                            }
                        }
                        if (!doubleentry) {
                            //Make sure that the first atom in onew is the first one in the original smiles order. This is important to have a canonical smiles.
                            if (positionInVector + 1 < v.size()) {
                                Object atomAfterCenterInOriginalSmiles = v.get(positionInVector + 1);
                                int l = 0;
                                while (onew[0] != atomAfterCenterInOriginalSmiles) {
                                    Object placeholder = onew[onew.length - 1];
                                    for (int k = onew.length - 2; k > -1; k--) {
                                        onew[k + 1] = onew[k];
                                    }
                                    onew[0] = placeholder;
                                    l++;
                                    if (l > onew.length) {
                                        break;
                                    }
                                }
                            }
                            //This cares about ring openings. Here the ring closure (represendted by a figure) must be the first atom. In onew the closure is null.
                            if (getRingOpenings(atom, null).size() > 0) {
                                int l = 0;
                                while (onew[0] != null) {
                                    Object placeholder = onew[0];
                                    for (int k = 1; k < onew.length; k++) {
                                        onew[k - 1] = onew[k];
                                    }
                                    onew[onew.length - 1] = placeholder;
                                    l++;
                                    if (l > onew.length) {
                                        break;
                                    }
                                }
                            }
                            //The last in onew is a vector: This means we need to exchange the rest of the original smiles with the rest of this vector.
                            if (onew[numberOfAtoms - 1] instanceof List) {
                                for (int i = 0; i < numberOfAtoms; i++) {
                                    if (onew[i] instanceof IAtom) {
                                        List vtemp = new Vector();
                                        vtemp.add(onew[i]);
                                        for (int k = positionInVector + 1 + numberOfAtoms; k < v.size(); k++) {
                                            vtemp.add(v.get(k));
                                        }
                                        onew[i] = vtemp;
                                        for (int k = v.size() - 1; k > positionInVector + 1 + numberOfAtoms - 1; k--) {
                                            v.remove(k);
                                        }
                                        for (int k = 1; k < ((Collection) onew[numberOfAtoms - 1]).size(); k++) {
                                            v.add(((List) onew[numberOfAtoms - 1]).get(k));
                                        }
                                        onew[numberOfAtoms - 1] = ((List) onew[numberOfAtoms - 1]).get(0);
                                        break;
                                    }
                                }
                            }
                            //Put the onew objects in the original Vector
                            int k = 0;
                            for (Object onew1 : onew) {
                                if (onew1 != null) {
                                    v.set(positionInVector + 1 + k, onew1);
                                    k++;
                                }
                            }
                        }
                    }
                }
                parent = atom;
            } else {
                //Have Vector
                //LOGGER.debug("in parseChain after else");
                boolean brackets = true;
                List result = new Vector();
                addAtoms((List) o, result);
                IAtom prevAtom;

                /*
                 * Got to find last atom that was processed. This is to check the relative position of the current
                 * atom/chain with respect to its parent
                 */
                prevAtom = (IAtom) ((Vector) atomsInOrderOfSmiles).lastElement();
                int maxConnectedBondCount = 4;
                /**
                 * If the parent atom of this new chain is the very first atom
                 * in the SMILES string and this chain is placed immediately
                 * after the parent atom then the max connected bond count for
                 * the parent should be 3 instead of 4.
                 */
                if (atomsInOrderOfSmiles.indexOf(parent) == 0 && prevAtom == parent) {
                    maxConnectedBondCount = 3;
                }
                if (isRingOpening(parent, result) && container.getConnectedBondsCount(parent) < maxConnectedBondCount) {
                    brackets = false;
                }
                if (brackets) {
                    buffer.append('(');
                }
                parseChain((List) o, buffer, container, parent, chiral, doubleBondConfiguration, atomsInOrderOfSmiles, useAromaticity);
                if (brackets) {
                    buffer.append(')');
                }
            }

            positionInVector++;
            //LOGGER.debug("in parseChain after positionVector++");
        }
    }

    /**
     * Append the symbol for the bond order between <code>a1</code> and
     * <code>a2</code> to the <code>line</code>.
     *
     * @param line the StringBuffer that the bond symbol is appended to.
     * @param a1 Atom participating in the bond.
     * @param a2 Atom participating in the bond.
     * @param atomContainer the GraphAtomContainer that the SMILES string is
     * generated for.
     * @param useAromaticity true=aromaticity or sp2 will trigger lower case
     * letters, wrong=only sp2
     */
    private void parseBond(StringBuffer line, IAtom a1, IAtom a2, IAtomContainer atomContainer, boolean useAromaticity) {
        //LOGGER.debug("in parseBond()");
        if (useAromaticity && a1.getFlag(ISAROMATIC) && a2.getFlag(ISAROMATIC)) {
            return;
        }
        if (atomContainer.getBond(a1, a2) == null) {
            return;
        }
        IBond.Order type = atomContainer.getBond(a1, a2).getOrder();
        if (null != type) {
            switch (type) {
                case SINGLE:
                    break;
                case DOUBLE:
                    line.append("=");
                    break;
                case TRIPLE:
                    line.append("#");
                    break;
                // //LOGGER.debug("Unknown bond type");
                default:
                    break;
            }
        }
    }

    /**
     * Generates the SMILES string for the atom
     *
     * @param a the atom to generate the SMILES for.
     * @param buffer the string buffer that the atom is to be apended to.
     * @param container the GraphAtomContainer to analyze.
     * @param chiral is a chiral smiles wished?
     * @param parent the atom we came from.
     * @param atomsInOrderOfSmiles a vector containing the atoms in the order
     * they are in the smiles.
     * @param currentChain The chain we currently deal with.
     * @param useAromaticity true=aromaticity or sp2 will trigger lower case
     * letters, wrong=only sp2
     */
    private void parseAtom(IAtom a, StringBuffer buffer, IAtomContainer container, boolean chiral, boolean[] doubleBondConfiguration, IAtom parent, List atomsInOrderOfSmiles, List currentChain, boolean useAromaticity) {
        String symbol = a.getSymbol();
        if (a instanceof PseudoAtom) {
            symbol = "*";
        }

        boolean stereo = false;
        if (chiral) {
            stereo = isStereo(container, a);
        }
        boolean brackets = symbol.equals("B") || symbol.equals("C") || symbol.equals("N") || symbol.equals("O") || symbol.equals("P") || symbol.equals("S") || symbol.equals("F") || symbol.equals("Br") || symbol.equals("I") || symbol.equals("Cl");
        brackets = !brackets;
        //LOGGER.debug("in parseAtom()");
        //Deal with the start of a double bond configuration
        if (chiral && isStartOfDoubleBond(container, a, parent, doubleBondConfiguration)) {
            buffer.append('/');
        }

        String mass = generateMassString(a);
        brackets |= !mass.isEmpty();

        String charge = generateChargeString(a);
        brackets |= !charge.isEmpty();

        if (chiral && stereo && (isTrigonalBipyramidalOrOctahedral(container, a) != 0 || isSquarePlanar(container, a) || isTetrahedral(container, a, false) != 0 || isSquarePlanar(container, a))) {
            brackets = true;
        }
        if (brackets) {
            buffer.append('[');
        }
        buffer.append(mass);
        if ((useAromaticity && a.getFlag(ISAROMATIC))) {
            // we put in a special check for N.planar3 cases such
            // as for indole and pyrrole, which require an explicit
            // H on the nitrogen. However this only makes sense when
            // the connectivity is not 3 - so for a case such as n1ncn(c1)CC
            // the PLANAR3 N already has 3 bonds, so don't add a H for this case
            if (a.getSymbol().equals("N") && a.getHybridization() == PLANAR3 && container.getConnectedAtomsList(a).size() != 3) {
                buffer.append("[").append(a.getSymbol().toLowerCase()).append("H]");
            } else {
                buffer.append(a.getSymbol().toLowerCase());
            }
        } else {
            buffer.append(symbol);
            if (symbol.equals("*") && a.getImplicitHydrogenCount() != null && a.getImplicitHydrogenCount() > 0) {
                buffer.append("H").append(a.getImplicitHydrogenCount());
            }
        }
        if (a.getProperty(RING_CONFIG) != null && a.getProperty(RING_CONFIG).equals(UP)) {
            buffer.append('/');
        }
        if (a.getProperty(RING_CONFIG) != null && a.getProperty(RING_CONFIG).equals(DOWN)) {
            buffer.append('\\');
        }
        if (chiral && stereo && (isTrigonalBipyramidalOrOctahedral(container, a) != 0 || isSquarePlanar(container, a) || isTetrahedral(container, a, false) != 0)) {
            buffer.append('@');
        }
        if (chiral && stereo && isSquarePlanar(container, a)) {
            buffer.append("SP1");
        }
        //chiral
        //hcount
        buffer.append(charge);
        if (brackets) {
            buffer.append(']');
        }

        //LOGGER.debug("in parseAtom() after dealing with Pseudoatom or not");
        //Deal with the end of a double bond configuration
        if (chiral && isEndOfDoubleBond(container, a, parent, doubleBondConfiguration)) {
            IAtom viewFrom = null;
            for (int i = 0; i < currentChain.size(); i++) {
                if (currentChain.get(i) == parent) {
                    int k = i - 1;
                    while (k > -1) {
                        if (currentChain.get(k) instanceof IAtom) {
                            viewFrom = (IAtom) currentChain.get(k);
                            break;
                        }
                        k--;
                    }
                }
            }
            if (viewFrom == null) {
                for (int i = 0; i < atomsInOrderOfSmiles.size(); i++) {
                    if (atomsInOrderOfSmiles.get(i) == parent) {
                        viewFrom = (IAtom) atomsInOrderOfSmiles.get(i - 1);
                    }
                }
            }
            boolean afterThisAtom = false;
            IAtom viewTo = null;
            for (int i = 0; i < currentChain.size(); i++) {
                if (afterThisAtom && currentChain.get(i) instanceof IAtom) {
                    viewTo = (IAtom) currentChain.get(i);
                    break;
                }
                if (afterThisAtom && currentChain.get(i) instanceof List) {
                    viewTo = (IAtom) ((List) currentChain.get(i)).get(0);
                    break;
                }
                if (a == currentChain.get(i)) {
                    afterThisAtom = true;
                }
            }
            try {
                if (isCisTrans(viewFrom, a, parent, viewTo, container)) {
                    buffer.append('\\');
                } else {
                    buffer.append('/');
                }
            } catch (CDKException ex) {
                //If the user wants a double bond configuration, where there is none, we ignore this.
            }
        }
        List v = new Vector();
        Iterator it = getRingOpenings(a, v).iterator();
        Iterator it2 = v.iterator();
        //LOGGER.debug("in parseAtom() after checking for Ring openings");
        while (it.hasNext()) {
            Integer integer = (Integer) it.next();
            IAtom a2 = (IAtom) it2.next();
            IBond b = container.getBond(a2, a);
            IBond.Order type = b.getOrder();
            if (!(useAromaticity
                    && a.getFlag(ISAROMATIC)
                    && a2.getFlag(ISAROMATIC))) {
                if (type == DOUBLE) {
                    buffer.append("=");
                } else if (type == TRIPLE) {
                    buffer.append("#");
                }
            }
            if (integer >= 10) {
                buffer.append("%").append(integer);
            } else {
                buffer.append(integer);
            }
        }
        atomsInOrderOfSmiles.add(a);
        //LOGGER.debug("End of parseAtom()");
    }

    /**
     * Creates a string for the charge of atom <code>a</code>. If the charge is
     * 1 + is returned if it is -1 - is returned. The positive values all have +
     * in front of them.
     *
     * @return string representing the charge on <code>a</code>
     */
    private String generateChargeString(IAtom a) {
        int charge = Objects.equals(a.getFormalCharge(), UNSET) ? 0 : a.getFormalCharge();
        StringBuilder buffer = new StringBuilder(3);
        if (charge > 0) {
            //Positive
            buffer.append('+');
            if (charge > 1) {
                buffer.append(charge);
            }
        } else if (charge < 0) {
            //Negative
            if (charge == -1) {
                buffer.append('-');
            } else {
                buffer.append(charge);
            }
        }
        return buffer.toString();
    }

    /**
     * Creates a string containing the mass of the atom <code>a</code>. If the
     * mass is the same as the majour isotope an empty string is returned.
     *
     * @param a the atom to create the mass
     */
    private String generateMassString(IAtom a) {
        if (isotopeFactory == null) {
            setupIsotopeFactory(a.getBuilder());
        }

        if (a instanceof IPseudoAtom) {
            if (a.getMassNumber() != null) {
                return Integer.toString(a.getMassNumber());
            } else {
                return "";
            }
        }

        IIsotope majorIsotope = isotopeFactory.getMajorIsotope(a.getSymbol());
        if (majorIsotope == null || Objects.equals(majorIsotope.getMassNumber(), a.getMassNumber())) {
            return "";
        } else if (a.getMassNumber() == null) {
            return "";
        } else {
            return Integer.toString(a.getMassNumber());
        }
    }

    private void setupIsotopeFactory(IChemObjectBuilder builder) {
        try {
            isotopeFactory = getInstance();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * Returns the current AllRingsFinder instance
     *
     * @return the current AllRingsFinder instance
     */
    public AllRingsFinder getRingFinder() {
        return ringFinder;
    }

    /**
     * Sets the current AllRingsFinder instance Use this if you want to
     * customize the timeout for the AllRingsFinder. AllRingsFinder is stopping
     * its quest to find all rings after a default of 5 seconds.
     *
     * @see org.openscience.cdk.ringsearch.AllRingsFinder
     *
     * @param ringFinder The value to assign ringFinder.
     */
    public void setRingFinder(AllRingsFinder ringFinder) {
        this.ringFinder = ringFinder;
    }

    /**
     * Indicates whether output should be an aromatic SMILES.
     *
     * @param useAromaticityFlag if false only SP2-hybridized atoms will be
     * lower case (default), true=SP2 or aromaticity trigger lower case
     */
    public void setUseAromaticityFlag(boolean useAromaticityFlag) {
        this.useAromaticityFlag = useAromaticityFlag;
    }

    class BrokenBond {

        /**
         * The atoms which close the ring
         */
        private final IAtom a1;
        private final IAtom a2;
        /**
         * The number of the marker
         */
        private final int marker;

        /**
         * Construct a BrokenBond between <code>a1</code> and <code>a2</code>
         * with the marker <code>marker</code>.
         *
         * @param marker the ring closure marker. (Great comment!)
         */
        BrokenBond(IAtom a1, IAtom a2, int marker) {
            this.a1 = a1;
            this.a2 = a2;
            this.marker = marker;
        }

        /**
         * Getter method for a1 property
         *
         * @return The a1 value
         */
        public IAtom getA1() {
            return a1;
        }

        /**
         * Getter method for a2 property
         *
         * @return The a2 value
         */
        public IAtom getA2() {
            return a2;
        }

        /**
         * Getter method for marker property
         *
         * @return The marker value
         */
        public int getMarker() {
            return marker;
        }

        @Override
        public String toString() {
            return Integer.toString(marker);
        }

        @Override
        public boolean equals(Object o) {
            if (!(o instanceof BrokenBond)) {
                return false;
            }
            BrokenBond bond = (BrokenBond) o;
            return (a1.equals(bond.getA1()) && a2.equals(bond.getA2())) || (a1.equals(bond.getA2()) && a2.equals(bond.getA1()));
        }

        @Override
        public int hashCode() {
            int hash = 5;
            hash = 67 * hash + Objects.hashCode(this.a1);
            hash = 67 * hash + Objects.hashCode(this.a2);
            hash = 67 * hash + this.marker;
            return hash;
        }
    }
}
