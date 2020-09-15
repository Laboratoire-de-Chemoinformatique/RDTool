/* Copyright (R) 2009-2018  Syed Asad Rahman <asad at ebi.ac.uk>
 *
 * Contact: cdk-devel@lists.sourceforge.net
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * All we ask is that proper credit is given for our work, which includes
 * - but is not limited to - adding the above copyright notice to the beginning
 * of your source code files, and to any copyright notice that you may distribute
 * with programs based on this work.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package org.openscience.smsd.algorithm.mcsplus.mcsplus3;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

/**
 * This class implements Bron-Kerbosch clique detection algorithm as it is
 * described in [F. Cazals, R. Karande: An Algorithm for reporting maximal
 * c-cliques; processedVertex.Comp. Sc. (2005); vol 349; pp. 484-490]
 *
 *
 * BronKerboschCazalsKarandeKochCliqueFinder.java
 *
 *
 *
 * @author Syed Asad Rahman <asad at ebi.ac.uk>
 */
public class MCSPlus extends Filter {

    private boolean DEBUG = false;

    /**
     * Creates a new instance of SearchCliques
     *
     *
     * @param f1
     * @param f2
     * @param shouldMatchBonds
     * @param shouldMatchRings
     * @param matchAtomType
     */
    public MCSPlus(IAtomContainer f1, IAtomContainer f2, boolean shouldMatchBonds, boolean shouldMatchRings, boolean matchAtomType) {

        super(f1, f2, shouldMatchBonds, shouldMatchRings, matchAtomType);

    }

    private List<List<Integer>> label_atoms(List<Integer> basic_atom_vector, int bond_num, List<IAtom> atoms, List<Integer[]> i_tab, List<String[]> c_tab) {

        ArrayList<List<Integer>> label_list = new ArrayList<>();

        for (int i = 0; i < basic_atom_vector.size(); i++) {

            List<Integer> label = new ArrayList<>(7);
            /*
             * Initialize the vector
             */
            for (int j = 0; j < 7; j++) {
                label.add(0);
            }

            IAtom atom1 = atoms.get(i);
            String atom1_type = matchAtomType && atom1.getAtomTypeName() != null ? atom1.getAtomTypeName() : atom1.getSymbol();
            if (!SYMBOL_VALUE.containsKey(atom1_type)) {
                int value = atom1.getAtomicNumber() == null ? atom1.hashCode() : atom1.getAtomicNumber();
                SYMBOL_VALUE.put(atom1_type, value + 1000);
            }
            label.set(0, SYMBOL_VALUE.get(atom1_type));
            int count_neighbors = 1;
            for (int j = 0; j < bond_num; j++) {
                Integer[] table_entry = i_tab.get(j);
                String[] c_table_entry = c_tab.get(j);
                int ind;
                if (basic_atom_vector.get(i).equals(table_entry[0])) ind = 0;
                else if (basic_atom_vector.get(i).equals(table_entry[1])) ind = 1;
                else continue;
                int other_ind = (ind + 1) % 2;
                IAtom atom2 = atoms.get(table_entry[other_ind] - 1);
                String atom2_type = matchAtomType && atom2.getAtomTypeName() != null ? atom2.getAtomTypeName() : c_table_entry[other_ind];
                if (!SYMBOL_VALUE.containsKey(atom2_type)) {
                    int value = atom2.getAtomicNumber() == null ? atom2.hashCode() : atom2.getAtomicNumber();
                    SYMBOL_VALUE.put(atom2_type, value);
                }
                label.set(count_neighbors, SYMBOL_VALUE.get(atom2_type));
                count_neighbors++;
            }
            List<Integer> bubbleSort = Utility.getBubbleSort(label);
            label_list.add(bubbleSort);
        }
        return label_list;
    }

    private List<Integer> reduce_atomset(int atom_num, int bond_num, List<IAtom> a_str, List<Integer[]> i_table, List<String[]> c_table) {
        // TODO: this stuff could be rewritten involving symbols information (c_table) which is not used at the moment
        List<Integer> phosphate_O_atoms = new ArrayList<>();
        List<Integer> h_atoms = new ArrayList<>();

        for (int i = 0; i < atom_num; i++) {
            if ("O".equals(a_str.get(i).getSymbol())) {
                int O_neighbor_num = 0;
                boolean P_neighbor = false;
                for (int j = 0; j < bond_num; j++) {
                    Integer[] table_entry = i_table.get(j);
                    if (i + 1 == table_entry[0]) {
                        O_neighbor_num++;
                        if ("P".equals(a_str.get(table_entry[1] - 1).getSymbol()) && table_entry[2] != 2)
                            P_neighbor = true;
                    }
                    else if (i + 1 == table_entry[1]) {
                        O_neighbor_num++;
                        if ("P".equals(a_str.get(table_entry[0] - 1).getSymbol()) && table_entry[2] != 2)
                            P_neighbor = true;
                    }
                }
                if (O_neighbor_num == 1 && P_neighbor) phosphate_O_atoms.add(i + 1);
            }
            if ("H".equals(a_str.get(i).getSymbol())) h_atoms.add(i + 1);
        }

        List<Integer> basic_atoms = new ArrayList<>();
        for (int i = 0; i < atom_num; i++)
            if (!phosphate_O_atoms.contains(i + 1) && !h_atoms.contains(i + 1)) basic_atoms.add(i + 1);
        return basic_atoms;
    }

    private int generate_compatibility_graph_nodes() {

        List<Integer> basic_atom_vec_A = reduce_atomset(num_heavy_atoms_1, num_bonds_1, atomstr1, int_table_1, char_table_1);
        List<Integer> basic_atom_vec_B = reduce_atomset(num_heavy_atoms_2, num_bonds_2, atomstr2, int_table_2, char_table_2);

        List<List<Integer>> label_list_molA = label_atoms(basic_atom_vec_A, num_bonds_1, atomstr1, int_table_1, char_table_1);
        List<List<Integer>> label_list_molB = label_atoms(basic_atom_vec_B, num_bonds_2, atomstr2, int_table_2, char_table_2);

        int i = 0;

        for (List<Integer> labelA : label_list_molA) {
            int j = 0;
            for (List<Integer> labelB : label_list_molB) {
                if (labelA.equals(labelB))
                    comp_graph_nodes.add(new Integer[]{basic_atom_vec_A.get(i), basic_atom_vec_B.get(j)});
                j++;
            }
            i++;
        }

        if (DEBUG) {
            System.out.println("comp_graph_nodes: " + comp_graph_nodes.size());
        }

        return 0;
    }

    private int generate_compatibility_graph() {
        for (int i = 0; i < comp_graph_nodes.size(); i++)
            for (int j = i + 1; j < comp_graph_nodes.size(); j++) {
                Integer[] i_entry = comp_graph_nodes.get(i);
                Integer[] j_entry = comp_graph_nodes.get(j);
                if (i_entry[0].equals(j_entry[0]) || i_entry[1].equals(j_entry[1])) continue;

                IBond bond1 = find_bond(int_table_1, ac1, i_entry[0], j_entry[0]);
                IBond bond2 = find_bond(int_table_2, ac2, i_entry[1], j_entry[1]);
                boolean mol1_pair_connected = bond1 != null;
                boolean mol2_pair_connected = bond2 != null;
                boolean connectedFlag = mol1_pair_connected && mol2_pair_connected;
                boolean disConnectedFlag = !mol1_pair_connected && !mol2_pair_connected;
                boolean matchBondFlag = connectedFlag && isMatchFeasible(bond1, bond2, shouldMatchBonds, shouldMatchRings, matchAtomType);
                if (connectedFlag && matchBondFlag)
                    c_edges.add(new Integer[]{i + 1, j + 1});
                if (disConnectedFlag || (connectedFlag && !matchBondFlag))
                    d_edges.add(new Integer[]{i + 1, j + 1});
            }

        c_edges_size = c_edges.size();
        d_edges_size = d_edges.size();

        return 0;
    }

    private IBond find_bond(List<Integer[]> int_table, IAtomContainer ac, int i_val, int j_val) {
        for (Integer[] table_entry : int_table) {
            List<Integer> entry_as_list = Arrays.asList(table_entry);
            if (entry_as_list.contains(i_val) && entry_as_list.contains(j_val)) {
                IAtom a1 = ac.getAtom(i_val - 1);
                IAtom a2 = ac.getAtom(j_val - 1);
                return ac.getBond(a1, a2);
            }
        }
        return null;
    }

    // comp_graph_nodes_C_zero is used to build up of the edges of the compatibility graph
    private int generate_compatibility_graph_nodes_if_C_edge_number_is_zero() {
//        System.out.println("Show atoms symbols");
//        for (int i = 0; i < num_heavy_atoms_1; i++) System.out.print(atomstr1.get(i).getSymbol() + " ");
//        System.out.println();
//        for (int i = 0; i < num_heavy_atoms_2; i++) System.out.print(atomstr2.get(i).getSymbol() + " ");
//        System.out.println();

//        System.out.println("Adding item to comp_graph_nodes");
        for (int i = 0; i < num_heavy_atoms_1; i++) {
            String atom1_type = atomstr1.get(i).getSymbol();
            int value = atomstr1.get(i).getAtomicNumber() == null ? atomstr1.get(i).hashCode() : atomstr1.get(i).getAtomicNumber();
            SYMBOL_VALUE.put(atom1_type, value + 1000);
            for (int j = 0; j < num_heavy_atoms_2; j++) {
                String atom2_type = atomstr2.get(j).getSymbol();
                if (atom1_type.equals(atom2_type)) {
                    Integer[] item1 = new Integer[]{i + 1, j + 1, SYMBOL_VALUE.get(atom1_type)};
                    Integer[] item2 = new Integer[]{i + 1, j + 1};
                    comp_graph_nodes_C_zero.add(item1);
                    comp_graph_nodes.add(item2);
                }
            }
        }
//        System.out.println("comp_graph_nodes #3");
//        for (Integer[] p : comp_graph_nodes) System.out.println("\t" + p[0] + ", " + p[1]);
//        System.exit(1);

        return 0;
    }

    private int generate_compatibility_graph_if_C_edge_number_is_zero() {

        for (int i = 0; i < comp_graph_nodes_C_zero.size(); i++) {
            Integer[] i_entry = comp_graph_nodes_C_zero.get(i);
            for (int j = i + 1; j < comp_graph_nodes_C_zero.size(); j++) {
                Integer[] j_entry = comp_graph_nodes_C_zero.get(j);
                if (i_entry[0].equals(j_entry[0]) || i_entry[1].equals(j_entry[1])) continue;;
                IBond b1 = find_bond(int_table_1, ac1, i_entry[0], j_entry[0]);
                IBond b2 = find_bond(int_table_2, ac2, i_entry[1], j_entry[1]);
                boolean mol1_pair_connected = b1 != null;
                boolean mol2_pair_connected = b2 != null;
                boolean connectedFlag = mol1_pair_connected && mol2_pair_connected;
                boolean disConnectedFlag = !mol1_pair_connected && !mol2_pair_connected;
                boolean matchBondFlag = connectedFlag && isMatchFeasible(b1, b2, shouldMatchBonds, shouldMatchRings, matchAtomType);
                if (connectedFlag && matchBondFlag)
                    c_edges.add(new Integer[]{i + 1, j + 1});
                if (disConnectedFlag || (connectedFlag && !matchBondFlag))
                    d_edges.add(new Integer[]{i + 1, j + 1});
            }
        }
        c_edges_size = c_edges.size();
        d_edges_size = d_edges.size();
        return 0;
    }

    //extract atom mapping from the clique vector and print it on the screen
    int extract_mapping(List<Integer> clique_vector) {

        List<Integer[]> temp_vector = new ArrayList<>();
        for (int i = 0; i < comp_graph_nodes.size(); i++)
            for (Integer cl_vec : clique_vector)
                if (cl_vec.equals(i))
                    temp_vector.add(comp_graph_nodes.get(i));
                    // TODO: need a break here?

        getFinalMappings().add(temp_vector);
        return 0;
    }

    //extract atom mapping from the clique vector and store it in vector clique_MAPPING_Local
    private List<Integer> extract_clique_MAPPING(List<Integer> clique_vector) {
        List<Integer> clique_mapping_local = new ArrayList<>();
        for (int i = 0; i < comp_graph_nodes.size(); i++) {
            Integer[] entry = comp_graph_nodes.get(i);
            for (Integer cliq : clique_vector)
                if (cliq.equals(i))
                    // TODO: break?
                    clique_mapping_local.addAll(Arrays.asList(comp_graph_nodes.get(i)));
        }
        return clique_mapping_local;
    }

//Function is called by the main program and serves as a starting point for the comparision procedure.
    public int search_cliques() {
//        boolean troubled = this.ac1.getID().equals("training_balanced_49:M00004") || this.ac2.getID().equals("training_balanced_49:M00004");

        generate_compatibility_graph_nodes();
        generate_compatibility_graph();
        if (c_edges_size == 0
                || (c_edges_size < this.ac1.getAtomCount() / 2 && c_edges_size < this.ac2.getAtomCount() / 2)
                && (this.ac1.getAtomCount() / 2 < 30 && this.ac2.getAtomCount() / 2 < 30)) {
            comp_graph_nodes.clear();
            c_edges.clear();
            d_edges.clear();
            c_edges_size = 0;
            d_edges_size = 0;
            generate_compatibility_graph_nodes_if_C_edge_number_is_zero();
            generate_compatibility_graph_if_C_edge_number_is_zero();
            comp_graph_nodes_C_zero.clear();
        }

        BKKCKCF cliqueFinder = new BKKCKCF(comp_graph_nodes, c_edges, d_edges);
        cliqueFinder.init_Algorithm();
        this.max_cliques_set = cliqueFinder.getMax_Cliques_Set();

        best_mapping_size = 0;

        int clique_number = 1;
//        System.out.println("c_edges:");
//        for (Integer[] p : c_edges) System.out.println("\t" + p[0] + ", " + p[1]);
//        System.out.println("d_edges:");
//        for (Integer[] p : d_edges) System.out.println("\t" + p[0] + ", " + p[1]);
        while (!max_cliques_set.empty()) {
            List<Integer> clique_vector = max_cliques_set.peek();
            int clique_size = clique_vector.size();
            // Is the number of mappings smaller than the number of atoms of molecule A and B?
            // In this case the clique is given to the McGregor algorithm
            if (clique_size < num_atoms_1 && clique_size < num_atoms_2) {
                try {
                    McGregor_IterationStart(clique_vector);
                } catch (Exception e) {
                    e.printStackTrace();
                }

            } else {
                extract_mapping(clique_vector);
            }
            max_cliques_set.pop();
//            break;
        }
//        System.out.println("\nsearch_cliques, show results:");
//        for (List<Integer[]> thelist : getFinalMappings()) {
//            System.out.println("Mapping:");
//            for (Integer[] x : thelist) System.out.println(x[0] + ", " + x[1]);
//            System.out.println("----------");
//        }
//        System.out.println("# of solutions before filter: " + getFinalMappings().size());

        postfilter();

//        System.out.println("# of solutions after filter: " + getFinalMappings().size());

        return 0;
    }

    private void clear() {
        this.max_cliques_set.clear();
        this.comp_graph_nodes.clear();
        this.comp_graph_nodes_C_zero.clear();
        this.char_table_1.clear();
        this.char_table_2.clear();
        this.c_edges.clear();
        this.d_edges.clear();
        this.c_edges_size = 0;
        this.d_edges_size = 0;
    }
}
