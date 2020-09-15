/*
 * Copyright (c) 2018. BioInception Labs Pvt. Ltd.
 */
package org.openscience.smsd.algorithm.mcsplus;

import java.io.IOException;
import java.util.*;
import java.util.concurrent.RecursiveTask;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomType;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.isomorphism.matchers.IQueryAtom;
import org.openscience.cdk.isomorphism.matchers.IQueryBond;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import org.openscience.smsd.algorithm.matchers.DefaultMatcher;
import org.openscience.smsd.helper.LabelContainer;

/**
 * This class generates compatibility graph between query and target molecule.
 * It also marks edges in the compatibility graph as c-edges or d-edges.
 *
 *
 *
 *
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public class GenerateCompatibilityGraphFJ extends RecursiveTask<List<Result>> {

    private final static ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(GenerateCompatibilityGraphFJ.class);
    private final static boolean DEBUG = false;
    private static final int THRESHOLD = 20;
    private static final int COMPLEX_MAX_GRAPH_NODE_COUNT = 30;
    private final int startIndex;
    private final int endIndex;

    private final IAtomContainer source;
    private final IAtomContainer target;
    private final boolean shouldMatchRings;
    private final boolean shouldMatchBonds;
    private final boolean matchAtomType;

    /**
     *
     * @param startIndex
     * @param endIndex
     * @param source
     * @param target
     * @param shouldMatchBonds
     * @param shouldMatchRings
     * @param matchAtomType
     */
    public GenerateCompatibilityGraphFJ(int startIndex,
            int endIndex,
            IAtomContainer source,
            IAtomContainer target,
            boolean shouldMatchBonds,
            boolean shouldMatchRings,
            boolean matchAtomType) {
        this.endIndex = endIndex;
        this.source = source;
        this.target = target;
        this.startIndex = startIndex;
        this.shouldMatchRings = shouldMatchRings;
        this.shouldMatchBonds = shouldMatchBonds;
        this.matchAtomType = matchAtomType;
    }

    @Override

    protected List<Result> compute() {

        if (endIndex - startIndex < THRESHOLD) {
            List<Result> arrayList = new ArrayList<>();
            arrayList.add(processing(startIndex, endIndex));
            return new ArrayList<>(new HashSet<>(arrayList));//remove any duplicates
        } else {

            if (DEBUG) {
                System.out.println("Splitting workLoad startIndex: " + startIndex + ", endIndex: " + endIndex);
            }

            List<GenerateCompatibilityGraphFJ> subtasks
                    = new ArrayList<>();
            subtasks.addAll(createSubtasks());

            //Collection<CustomRecursiveTask> invokeAll = invokeAll(subtasks);
            subtasks.forEach((subtask) -> {
                subtask.fork();
            });

            List<Result> result = new ArrayList<>();
            subtasks.forEach((subtask) -> {
                result.addAll(subtask.join());
            });
            return new ArrayList<>(new HashSet<>(result));//remove any duplicates;
        }
    }

    private List<GenerateCompatibilityGraphFJ> createSubtasks() {
        List<GenerateCompatibilityGraphFJ> dividedTasks = new ArrayList<>();
        int middle = (endIndex + startIndex) / 2;

        GenerateCompatibilityGraphFJ partOne = new GenerateCompatibilityGraphFJ(startIndex, middle, source, target, shouldMatchBonds, shouldMatchRings, matchAtomType);
        GenerateCompatibilityGraphFJ partTwo = new GenerateCompatibilityGraphFJ(middle, endIndex, source, target, shouldMatchBonds, shouldMatchRings, matchAtomType);
        dividedTasks.add(partOne);
        dividedTasks.add(partTwo);

        return dividedTasks;
    }

    private Result processing(int startIndex, int endIndex) {
        Result result;
        if (DEBUG) {
            System.out.println(" GenerateCompatibilityGraphFJ ");
            System.out.println("Splitting workLoad startIndex: " + startIndex + ", endIndex: " + endIndex);
        }

        if ((!shouldMatchBonds || !matchAtomType)
                && source.getAtomCount() > COMPLEX_MAX_GRAPH_NODE_COUNT
                || target.getAtomCount() > COMPLEX_MAX_GRAPH_NODE_COUNT) {
            result = new Result();
            if (DEBUG) {
                System.out.println("CASE LARGE GRAPH");
            }
            List<Integer> compGraphNodesCZero = new ArrayList<>(); //Initialize the compGraphNodesCZero List
            compatibilityGraphNodesIfCEdgeIsZero(startIndex, endIndex, result, compGraphNodesCZero);
            compatibilityGraphCEdgeZero(result, compGraphNodesCZero);
            compGraphNodesCZero.clear();
        } else {
            result = new Result();
            if (DEBUG) {
                System.out.println("Calling Compatibility Graph Nodes ");
            }
            compatibilityGraphNodes(startIndex, endIndex, result);
            if (DEBUG) {
                System.out.println("Calling Compatibility Graph ");
            }
            compatibilityGraph(result);
            if (DEBUG) {
                System.out.println("c-edges " + result.cEdges.size());
            }
            if (DEBUG) {
                System.out.println("d-edges " + result.dEdges.size());
            }

            if (result.cEdges.isEmpty()) {
                result = new Result();
                List<Integer> compGraphNodesCZero = new ArrayList<>(); //Initialize the compGraphNodesCZero List
                compatibilityGraphNodesIfCEdgeIsZero(startIndex, endIndex, result, compGraphNodesCZero);
                compatibilityGraphCEdgeZero(result, compGraphNodesCZero);
                compGraphNodesCZero.clear();
            }
        }
        return result;
    }

    /**
     * compGraphNodesCZero is used to build up of the edges of the compatibility
     * graph
     *
     * @return
     * @throws IOException
     */
    private Integer compatibilityGraphNodesIfCEdgeIsZero(int startIndex, int endIndex, Result result, List<Integer> compGraphNodesCZero) {

        int count_nodes = 1;
        List<String> list = new ArrayList<>();
        LabelContainer labelContainer = LabelContainer.getInstance();

        for (int i = startIndex; i < endIndex; i++) {
            for (int j = 0; j < target.getAtomCount(); j++) {
                IAtom atom1 = source.getAtom(i);
                IAtom atom2 = target.getAtom(j);

                //You can also check object equal or charge, hydrogen count etc
                if ((atom1 instanceof IQueryAtom)
                        && ((IQueryAtom) atom1).matches(atom2)
                        && !list.contains(i + "_" + j)) {
                    compGraphNodesCZero.add(i);
                    compGraphNodesCZero.add(j);
                    compGraphNodesCZero.add(labelContainer.getLabelID(atom2.getSymbol())); //i.e C is label 1
                    compGraphNodesCZero.add(count_nodes);
                    result.compGraphNodes.add(i);
                    result.compGraphNodes.add(j);
                    result.compGraphNodes.add(count_nodes);
                    count_nodes += 1;
                    list.add(i + "_" + j);
                } else if (atom1.getSymbol().equalsIgnoreCase(atom2.getSymbol())
                        && !list.contains(i + "_" + j)) {
                    compGraphNodesCZero.add(i);
                    compGraphNodesCZero.add(j);
                    compGraphNodesCZero.add(labelContainer.getLabelID(atom1.getSymbol())); //i.e C is label 1
                    compGraphNodesCZero.add(count_nodes);
                    result.compGraphNodes.add(i);
                    result.compGraphNodes.add(j);
                    result.compGraphNodes.add(count_nodes);
                    count_nodes += 1;
                    list.add(i + "_" + j);
                }
            }
        }
        list.clear();
        if (DEBUG) {
            System.out.println("count_nodes " + count_nodes);
        }
        return count_nodes;
    }

    /**
     * compatibilityGraphCEdgeZero is used to build up of the edges of the
     * compatibility graph BIS
     *
     * @return
     * @throws IOException
     */
    private int compatibilityGraphCEdgeZero(Result result, List<Integer> compGraphNodesCZero) {

        int compGraphNodesCZeroListSize = compGraphNodesCZero.size();

        for (int a = 0; a < compGraphNodesCZeroListSize; a += 4) {
            int index_a = compGraphNodesCZero.get(a);
            int index_aPlus1 = compGraphNodesCZero.get(a + 1);
            for (int b = a + 4; b < compGraphNodesCZeroListSize; b += 4) {
                int index_b = compGraphNodesCZero.get(b);
                int index_bPlus1 = compGraphNodesCZero.get(b + 1);

                // if element atomCont !=jIndex and atoms on the adjacent sides of the bonds are not equal
                if ((a != b) && (index_a != index_b)
                        && (index_aPlus1 != index_bPlus1)) {

                    IBond reactantBond;
                    IBond productBond;

                    reactantBond = source.getBond(source.getAtom(index_a), source.getAtom(index_b));
                    productBond = target.getBond(target.getAtom(index_aPlus1), target.getAtom(index_bPlus1));

                    if (reactantBond != null && productBond != null) {
                        addZeroEdges(result.cEdges, result.dEdges, reactantBond, productBond, a, b);
                    } else if (reactantBond == null && productBond == null
                            && ((source.getAtomCount() < (COMPLEX_MAX_GRAPH_NODE_COUNT)
                            && target.getAtomCount() < (COMPLEX_MAX_GRAPH_NODE_COUNT))
                            || (result.dEdges.size() < result.compGraphNodes.size()))) {
                        //50 unique condition to speed up the AAM
                        Edge edge = new Edge((a / 4) + 1, (b / 4) + 1);
                        if (!result.dEdges.contains(edge)) {
                            result.dEdges.add(edge);
                        }
                    }
                }
            }
        }
        if (DEBUG) {
            //Size of C and D edges of the compatibility graph
            int cEdgesSize = result.cEdges.size();
            int dEdgesSize = result.dEdges.size();
            System.out.println("cEdgesSize " + cEdgesSize);
            System.out.println("dEdgesSize " + dEdgesSize);
        }
        return 0;
    }

    private void addZeroEdges(List<Edge> cEdges, List<Edge> dEdges,
            IBond reactantBond, IBond productBond,
            int indexI, int indexJ) {
        if (isMatchFeasible(reactantBond, productBond, shouldMatchBonds, shouldMatchRings, matchAtomType)) {
            Edge edge = new Edge((indexI / 4) + 1, (indexJ / 4) + 1);
            if (!cEdges.contains(edge)) {
                cEdges.add(edge);
            }
        } else {
            Edge edge = new Edge((indexI / 4) + 1, (indexJ / 4) + 1);
            if (!dEdges.contains(edge)) {
                dEdges.add(edge);
            }
        }
    }

    /**
     *
     * @param bondA1
     * @param bondA2
     * @param shouldMatchBonds
     * @param shouldMatchRings
     * @return
     */
    private boolean isMatchFeasible(
            IBond bondA1,
            IBond bondA2,
            boolean shouldMatchBonds,
            boolean shouldMatchRings,
            boolean matchAtomType) {

        if (bondA1 instanceof IQueryBond) {
            if (((IQueryBond) bondA1).matches(bondA2)) {
                IQueryAtom atom1 = (IQueryAtom) (bondA1.getAtom(0));
                IQueryAtom atom2 = (IQueryAtom) (bondA1.getAtom(1));
                return atom1.matches(bondA2.getAtom(0)) && atom2.matches(bondA2.getAtom(1))
                        || atom1.matches(bondA2.getAtom(1)) && atom2.matches(bondA2.getAtom(0));
            }
            return false;
        } else {
            /*
             This one also matches atom type, not just symbols
             */
            return DefaultMatcher.matches(bondA1, bondA2, shouldMatchBonds, shouldMatchRings, matchAtomType);
        }
    }

    private Map<IAtom, List<String>> labelAtomsBySymbol(int startIndex, int endIndex, IAtomContainer atomCont) {
        Map<IAtom, List<String>> label_list = new HashMap<>();

        for (int i = startIndex; i < endIndex; i++) {
            List<String> label = new ArrayList<>(7);
            for (int a = 0; a < 7; a++) {
                label.add(a, "Z9");
            }

            IAtom refAtom = atomCont.getAtom(i);
            if (refAtom == null) {
                return label_list;
            }
            /*
             * Important Step: Discriminate between source atom types
             */
            String referenceAtom;
            if (refAtom instanceof IQueryAtom) {
                referenceAtom = ((IQueryAtom) refAtom).getSymbol() == null ? "*" : ((IQueryAtom) refAtom).getSymbol();
                if (DEBUG) {
                    System.out.println("referenceAtom " + referenceAtom);
                }
            } else if (!(refAtom instanceof IQueryAtom) && this.matchAtomType) {
                referenceAtom = refAtom.getAtomTypeName() == null ? refAtom.getSymbol() : refAtom.getAtomTypeName();
            } else {
                referenceAtom = refAtom.getSymbol();
            }
            label.set(0, referenceAtom);
            List<IAtom> connAtoms = atomCont.getConnectedAtomsList(refAtom);

            int counter = 1;

            for (IAtom negAtom : connAtoms) {
                String neighbouringAtom;
                if (refAtom instanceof IQueryAtom) {
                    neighbouringAtom = ((IQueryAtom) negAtom).getSymbol() == null ? "*" : ((IQueryAtom) negAtom).getSymbol();
//                    System.out.println("neighbouringAtom " + neighbouringAtom);
                } else if (!(negAtom instanceof IQueryAtom) && this.matchAtomType) {
                    neighbouringAtom = negAtom.getAtomTypeName() == null ? negAtom.getSymbol() : negAtom.getAtomTypeName();
                } else {
                    neighbouringAtom = negAtom.getSymbol();
                }
                label.set(counter, neighbouringAtom);
                counter += 1;
            }
            if (DEBUG) {
                System.out.println("label " + label);
            }
            bubbleSort(label);
            label_list.put(refAtom, label);
        }
        return label_list;
    }

    private void bubbleSort(List<String> num) {
        int j;
        boolean flag = true;   // set flag to true to begin first pass
        String temp;   //holding variable

        while (flag) {
            flag = false;    //set flag to false awaiting a possible swap
            for (j = 0; j < (num.size() - 1); j++) {
                if (num.get(j).compareTo(num.get(j + 1)) > 0) // change to < for descending sort
                {
                    temp = num.get(j);                //swap elements
                    num.set(j, num.get(j + 1));
                    num.set(j + 1, temp);
                    flag = true;              //shows a swap occurred  
                }
            }
        }
    }

    /**
     * Generate Compatibility Graph Nodes
     *
     * @return
     * @throws IOException
     */
    private int compatibilityGraphNodes(int startIndex, int endIndex, Result result) {

        Set<Edge> edges = new HashSet<>();

        int nodeCount = 1;
        Map<IAtom, List<String>> labelAtomsBySymbolA = labelAtomsBySymbol(startIndex, endIndex, source);
        Map<IAtom, List<String>> labelAtomsBySymbolB = labelAtomsBySymbol(0, target.getAtomCount(), target);

//        for (Map.Entry<IAtom, List<String>> labelA : labelAtomsBySymbolA.entrySet()) {
//            if (DEBUG) {
//                System.out.println("labelA.getValue() " + labelA.getValue());
//            }
//            for (Map.Entry<IAtom, List<String>> labelB : labelAtomsBySymbolB.entrySet()) {
//                IAtom atom = labelA.getKey();
//                if (((atom instanceof IQueryAtom) && ((IQueryAtom) atom).matches(labelB.getKey()))
//                        || (!(atom instanceof IQueryAtom) && atom.getSymbol().equals(labelB.getKey().getSymbol()))) {
//                    if (DEBUG) {
//                        System.out.println("labelB.getValue() " + labelB.getValue());
//                    }
//                    int atomNumberI = source.indexOf(labelA.getKey());
//                    int atomNumberJ = target.indexOf(labelB.getKey());
//                    Edge e = new Edge(atomNumberI, atomNumberJ);
//                    if (!edges.contains(e)) {
//                        edges.add(e);
//                        result.compGraphNodes.add(atomNumberI);
//                        result.compGraphNodes.add(atomNumberJ);
//                        result.compGraphNodes.add(nodeCount);
//                        nodeCount += 1;
//                    }
//                }
//            }
//        }
        int i;
        IAtom[] atoms_a = new IAtom[labelAtomsBySymbolA.keySet().size()];
        IAtom[] atoms_b = new IAtom[labelAtomsBySymbolB.keySet().size()];
        i = 0;
        for (IAtom atom : labelAtomsBySymbolA.keySet()) atoms_a[i++] = atom;
        i = 0;
        for (IAtom atom : labelAtomsBySymbolB.keySet()) atoms_b[i++] = atom;
        Arrays.sort(atoms_a, new AtomComparator());
        Arrays.sort(atoms_b, new AtomComparator());
        for (IAtom atom : atoms_a) {
            if (DEBUG)
                System.out.println("AtomA ID: " + atom.getID() + ", label: " + labelAtomsBySymbolA.get(atom));
            for (IAtom atom_2 : atoms_b)
                if (((atom instanceof IQueryAtom) && ((IQueryAtom) atom).matches(atom_2)) || (!(atom instanceof IQueryAtom) && atom.getSymbol().equals(atom_2.getSymbol()))) {
                    int atom_number_I = source.indexOf(atom);
                    int atom_number_J = target.indexOf(atom_2);
                    Edge e = new Edge(atom_number_I, atom_number_J);
                    if (!edges.contains(e)) {
                        edges.add(e);
                        result.compGraphNodes.add(atom_number_I);
                        result.compGraphNodes.add(atom_number_J);
                        result.compGraphNodes.add(nodeCount);
                        nodeCount++;
                    }
                }
        }
        return 0;
    }

    /**
     * Generate Compatibility Graph Nodes Bond Insensitive
     *
     * @return
     * @throws IOException
     */
    private int compatibilityGraph(Result result) {
        int comp_graph_nodes_List_size = result.compGraphNodes.size();
        Map<String, Edge> c_edges = new TreeMap<>();
        Map<String, Edge> d_edges = new TreeMap<>();
        for (int i = 0; i < comp_graph_nodes_List_size; i += 3) {
            int i_node_0 = result.compGraphNodes.get(i);
            int i_node_1 = result.compGraphNodes.get(i + 1);
            for (int j = i + 3; j < comp_graph_nodes_List_size; j += 3) {
                int j_node_0 = result.compGraphNodes.get(j);
                int j_node_1 = result.compGraphNodes.get(j + 1);
                if (i_node_0 == j_node_0 || i_node_1 == j_node_1) continue;
                IBond reactantBond = source.getBond(source.getAtom(i_node_0), source.getAtom(j_node_0));
                IBond productBond = target.getBond(target.getAtom(i_node_1), target.getAtom(j_node_1));

                if (reactantBond != null && productBond != null)
                    addEdges(c_edges, d_edges, reactantBond, productBond, i, j);
                else if (reactantBond == null && productBond == null) {
                    Edge edge = new Edge((i / 3) + 1, (j / 3) + 1);
                    String key = edge.toString();
                    d_edges.put(key, edge);
                }

            }
        }
        result.cEdges.addAll(c_edges.values());
        result.dEdges.addAll(d_edges.values());
        return 0;
    }

    private void addEdges(List<Edge> cEdges, List<Edge> dEdges, IBond rBond, IBond pBond, int iIndex, int jIndex) {
        if (!shouldMatchBonds && !shouldMatchRings && !matchAtomType) {
            if (isRawMatch(rBond, pBond)) {
                cEdges.add(new Edge((iIndex / 3) + 1, (jIndex / 3) + 1));
            }
        } else if (isMatchFeasible(rBond, pBond, shouldMatchBonds, shouldMatchRings, matchAtomType)) {
            Edge edge = new Edge((iIndex / 3) + 1, (jIndex / 3) + 1);
            if (!cEdges.contains(edge)) {
                cEdges.add(edge);
            }
        } else {
            Edge edge = new Edge((iIndex / 3) + 1, (jIndex / 3) + 1);
            if (!dEdges.contains(edge)) {
                dEdges.add(edge);
            }
        }
    }

    private void addEdges(Map<String, Edge> c_edges, Map<String, Edge> d_edges, IBond rbond, IBond pbond, int i, int j) {
        Edge e = new Edge((i / 3) + 1, (j / 3) + 1);
        String key = e.toString();
        if (!shouldMatchBonds && !shouldMatchRings && !matchAtomType) {
            if (isRawMatch(rbond, pbond)) c_edges.put(key, e);
        }
        else if (isMatchFeasible(rbond, pbond, shouldMatchBonds, shouldMatchRings, matchAtomType)) {
            c_edges.put(key, e);
        }
        else d_edges.put(key, e);
    }

    private boolean isRawMatch(IBond reactantBond, IBond productBond) {
        String r_at_0_sym = reactantBond.getAtom(0).getSymbol();
        String r_at_1_sym = reactantBond.getAtom(1).getSymbol();
        String p_at_0_sym = productBond.getAtom(0).getSymbol();
        String p_at_1_sym = productBond.getAtom(1).getSymbol();
        IAtomType.Hybridization r_at_0_h = reactantBond.getAtom(0).getHybridization();
        IAtomType.Hybridization r_at_1_h = reactantBond.getAtom(1).getHybridization();
        IAtomType.Hybridization p_at_0_h = productBond.getAtom(0).getHybridization();
        IAtomType.Hybridization p_at_1_h = productBond.getAtom(1).getHybridization();
        if (r_at_0_sym.equals(p_at_0_sym) && r_at_1_sym.equals(p_at_1_sym))
            return r_at_0_h.equals(p_at_0_h) && r_at_1_h.equals(p_at_1_h);
        else if (r_at_1_sym.equals(p_at_0_sym) && r_at_0_sym.equals(p_at_1_sym))
            return r_at_1_h.equals(p_at_0_h) && r_at_0_h.equals(p_at_1_h);
        return false;
    }
}
