package org.openscience.smsd.algorithm.mcsplus.mcsplus3;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

import javax.print.attribute.Size2DSyntax;
import java.util.*;
import java.util.stream.Collectors;

public class McGregor extends Utility {
    protected final List<Integer[]> c_edges;
    protected final List<Integer[]> d_edges;
    protected int c_edges_size;
    protected int d_edges_size;
    private BinaryTree last, first;
    private final List<Integer[]> int_global_1;
    private final List<Integer[]> int_global_2;
    private final List<String[]> char_global_1;
    private final List<String[]> char_global_2;
    private List<Integer> MARCS;
    private final Stack<List<Integer>> BESTARCS;
    private int bestarcsleft;

    private boolean new_matrix;
    protected Stack<List<Integer>> max_cliques_set;

    protected final IAtomContainer ac1;
    protected final IAtomContainer ac2;
    protected final boolean shouldMatchBonds;
    protected final boolean shouldMatchRings;
    protected final boolean matchAtomType;

    protected int num_heavy_atoms_1;
    protected int num_heavy_atoms_2;
    protected List<IAtom> atomstr1;
    protected List<IAtom> atomstr2;
    protected int num_bonds_1;
    protected int num_bonds_2;
    protected int num_atoms_1;
    protected int num_atoms_2;
    protected List<Integer[]> int_table_1;
    protected List<Integer[]> int_table_2;
    protected List<String[]> char_table_1;
    protected List<String[]> char_table_2;
    private List<String[]> char_table_1_copy;
    private List<String[]> char_table_2_copy;

    protected final List<Integer[]> comp_graph_nodes;
    protected final List<Integer[]> comp_graph_nodes_C_zero;
    private final List<String> SignROW;

    protected int best_mapping_size;
    protected int best_clique_size;
    private final List<List<Integer[]>> final_mappings;
    protected final Map<String, Integer> SYMBOL_VALUE;
    private int iter_counter;

    protected McGregor(IAtomContainer ac1, IAtomContainer ac2, boolean shouldMatchBonds, boolean shouldMatchRings, boolean matchAtomType) {
        this.shouldMatchBonds = shouldMatchBonds;
        this.shouldMatchRings = shouldMatchRings;
        this.matchAtomType = matchAtomType;
        this.SYMBOL_VALUE = new TreeMap<>();
        MoleculeHandler handler_1 = new MoleculeHandler(ac1, shouldMatchBonds);
        this.num_heavy_atoms_1 = handler_1.get_num_H_atoms();
        this.num_bonds_1 = handler_1.getBondNumber();
        this.num_atoms_1 = handler_1.indexOf();
        this.int_table_1 = handler_1.int_table;
        this.char_table_1 = handler_1.char_table;
        this.ac1 = handler_1.getAtomContainer();
        this.atomstr1 = handler_1.getAtomString();
        MoleculeHandler handler_2 = new MoleculeHandler(ac2, shouldMatchBonds);
        this.num_heavy_atoms_2 = handler_2.get_num_H_atoms();
        this.num_bonds_2 = handler_2.getBondNumber();
        this.num_atoms_2 = handler_2.indexOf();
        this.int_table_2 = handler_2.int_table;
        this.char_table_2 = handler_2.char_table;
        this.ac2 = handler_2.getAtomContainer();
        this.atomstr2 = handler_2.getAtomString();

        this.comp_graph_nodes = new ArrayList<>();
        this.comp_graph_nodes_C_zero = new ArrayList<>();
        this.c_edges = new ArrayList<>();
        this.d_edges = new ArrayList<>();

        this.int_global_1 = new ArrayList<>();
        this.int_global_2 = new ArrayList<>();
        this.char_global_1 = new ArrayList<>();
        this.char_global_2 = new ArrayList<>();

        this.MARCS = new ArrayList<>();
        this.BESTARCS = new Stack<>();

        this.max_cliques_set = new Stack<>();
        this.final_mappings = new ArrayList<>();
        String[] characters = {"D", "G", "J", "T", "V", "W", "Y", "Z", "$", "%", "&", "*", "#", "?", "!", "~", "^", "<", ">", "=", "(", ")", "[", "]"};
        this.SignROW = Arrays.asList(characters);
        iter_counter = 0;
    }

    private List<String[]> copy_char_table(List<String[]> char_table, int num_bonds) {
        List<String[]> res = new ArrayList<>();
        for (int i = 0; i < num_bonds; i++) {
//            res.add(char_table.get(i));
            String[] entry = char_table.get(i);
            res.add(new String[]{entry[0], entry[1]});
            res.add(new String[]{"X", "X"});
        }
        return res;
    }

    private List<String[]> copy_existing_char_table(List<String[]> char_table, int num_bonds) {
        List<String[]> res = new ArrayList<>();
        for (int i = 0; i < num_bonds; i++) {
//            res.add(char_table.get(2 * i));
            String[] entry = char_table.get(2 * i);
            res.add(new String[]{entry[0], entry[1]});
            res.add(new String[]{"X", "X"});
        }
        return res;
    }

    private List<Integer> find_unmapped_atoms(List<Integer[]> mapped_atoms, int sz, int num_heavy_atoms, int atom_index) {
        assert (atom_index == 0 || atom_index == 1);
        List<Integer> res = new ArrayList<>();
        for (int i = 1; i <= num_heavy_atoms; i++) {  // TODO: why this range, not [0; num_heavy_atoms)?
            boolean is_unmapped = true;
            for (Integer[] pair : mapped_atoms)
                if (pair[atom_index].equals(i)) {
                    is_unmapped = false;
                    break;
                }
            if (is_unmapped) res.add(i);
        }
        return res;
    }

    private int search_corresponding_atom(int other_atom, int atom_index, List<Integer[]> mapped_atoms) {
        assert (atom_index == 0 || atom_index == 1);

        int corr_atom = 0;
        int other_index = atom_index == 0 ? 1 : 0;
        for (Integer[] mapped_atom : mapped_atoms)
            if (mapped_atom[atom_index] == other_atom)
                return mapped_atom[other_index];
        return corr_atom;
    }

    private List<String[]> change_char_bonds(int corr_atom, String new_sym, int nbr_bond_num, List<Integer[]> int_bond_nbrs, List<String[]> char_bond_nbrs) {
        List<String[]> res = new ArrayList<>(char_bond_nbrs);
        for (int i = 0; i < nbr_bond_num; i++) {
            if (int_bond_nbrs.get(i)[0].equals(corr_atom) && res.get(2 * i + 1)[0].equals("X")) {
                res.get(2 * i + 1)[0] = res.get(2 * i)[0];
                res.get(2 * i)[0] = new_sym;
            }
            else if (int_bond_nbrs.get(i)[1].equals(corr_atom) && res.get(2 * i + 1)[1].equals("X")) {
                res.get(2 * i + 1)[1] = res.get(2 * i)[1];
                res.get(2 * i)[1] = new_sym;
            }
        }
        return res;
    }

    protected int McGregor_IterationStart(List<Integer> clique_vector) {
        char_table_1_copy = copy_char_table(char_table_1, num_bonds_1);
        char_table_2_copy = copy_char_table(char_table_2, num_bonds_2);
        List<Integer[]> mapped_atoms = new ArrayList<>();
        for (Integer vec_val : clique_vector)
            mapped_atoms.add(comp_graph_nodes.get(vec_val - 1));
//            for (int i = 0; i < comp_graph_nodes.size(); i++)
//                if (vec_val.equals(i + 1)) {
//                    mapped_atoms.add(comp_graph_nodes.get(i));
//                    break;
//                }

        // find unmapped atoms and extract bonds for 1st molecule
        List<Integer> unmapped_atoms_mol_1 = find_unmapped_atoms(mapped_atoms, clique_vector.size(), num_heavy_atoms_1, 0);
        List<Integer[]> int_bond_nbrs_1 = new ArrayList<>();
        List<Integer[]> int_bond_set_1 = new ArrayList<>();
        List<String[]> char_bond_nbrs_1 = new ArrayList<>();
        List<String[]> char_bond_set_1 = new ArrayList<>();

        int SR_count = extract_bonds(
                0, num_bonds_1, num_bonds_2, 0, unmapped_atoms_mol_1, mapped_atoms,
                int_table_1, int_table_2, char_table_1_copy, char_table_2_copy,
                int_bond_nbrs_1, char_bond_nbrs_1, int_bond_set_1, char_bond_set_1, false
        );
        int nbr_bondnum_1 = int_bond_nbrs_1.size();
        int set_bondnum_1 = int_bond_set_1.size();

        // same for 2nd molecule
        List<Integer> unmapped_atoms_mol_2 = find_unmapped_atoms(mapped_atoms, clique_vector.size(), num_heavy_atoms_2, 1);
        List<Integer[]> int_bond_nbrs_2 = new ArrayList<>();
        List<Integer[]> int_bond_set_2 = new ArrayList<>();
        List<String[]> char_bond_nbrs_2 = new ArrayList<>();
        List<String[]> char_bond_set_2 = new ArrayList<>();

        extract_bonds(
                SR_count, num_bonds_2, nbr_bondnum_1, 1, unmapped_atoms_mol_2, mapped_atoms,
                int_table_2, int_bond_nbrs_1, char_table_2_copy, char_bond_nbrs_1,
                int_bond_nbrs_2, char_bond_nbrs_2, int_bond_set_2, char_bond_set_2, false
        );
        int nbr_bondnum_2 = int_bond_nbrs_2.size();
        int set_bondnum_2 = int_bond_set_2.size();
        Iterator(
                false, mapped_atoms.size(), mapped_atoms, nbr_bondnum_1, nbr_bondnum_2,
                int_bond_nbrs_1, int_bond_nbrs_2, char_bond_nbrs_1, char_bond_nbrs_2, set_bondnum_1, set_bondnum_2,
                int_bond_set_1, int_bond_set_2, char_bond_set_1, char_bond_set_2
        );
//        System.out.println("Showing final mappings for algorithm #3");
//        for (List<Integer[]> m : getFinalMappings()) {
//            System.out.println("Next mapping:");
//            for (Integer[] x : m) System.out.println(x[0] + ", " + x[1]);
//            System.out.println("-----");
//        }
        return 0;
    }

    private int Iterator(
            boolean mapping_check, int num_mapped_atoms, List<Integer[]> mapped_atoms, int nbr_bondnum_1, int nbr_bondnum_2,
            List<Integer[]> int_bond_nbrs_1, List<Integer[]> int_bond_nbrs_2, List<String[]> char_bond_nbrs_1, List<String[]> char_bond_nbrs_2,
            int set_bondnum_1, int set_bondnum_2, List<Integer[]> int_bond_set_1, List<Integer[]> int_bond_set_2,
            List<String[]> char_bond_set_1, List<String[]> char_bond_set_2
    ) {
        iter_counter++;
        int local_iter = iter_counter;
//        System.out.println("Entering iterator, iter = " + iter_counter);
//        boolean stop_mapping = false;
//        if (!mapping_check && nbr_bondnum_1 != 0 && nbr_bondnum_2 != 0) {
//            for (int i = 0; i < nbr_bondnum_1; i++) {
//                String[] atoms_1 = char_bond_nbrs_1.get(2 * i);
//                IBond bond_1 = get_bond(ac1, int_bond_nbrs_1.get(i));
//                for (int j = 0; j < nbr_bondnum_2; j++) {
//                    String[] atoms_2 = char_bond_nbrs_2.get(2 * j);
//                    IBond bond_2 = get_bond(ac2, int_bond_nbrs_2.get(j));
//                    if (!isMatchFeasible(bond_1, bond_2, shouldMatchBonds, shouldMatchRings, matchAtomType)) continue;
//                    if (atoms_1[0].equals(atoms_2[0]) && atoms_1[1].equals(atoms_2[1])) {
//                        stop_mapping = true;
//                        break;
//                    }
//                    else if (atoms_1[0].equals(atoms_2[1]) && atoms_1[1].equals(atoms_2[0])) {
//                        stop_mapping = true;
//                        break;
//                    }
//                }
//                if (stop_mapping) break;
//            }
//        }
//        if (stop_mapping) {
//            if (num_mapped_atoms > best_mapping_size) getFinalMappings().clear();
//            if (num_mapped_atoms >= best_mapping_size) {
//                best_mapping_size = Math.max(best_mapping_size, num_mapped_atoms);
//                getFinalMappings().add(mapped_atoms);
//            }
//            System.out.println("\tQuitting iterator, cond #1, iter = " + local_iter_num);
//            return 0;
//        }
        boolean no_map = true;
        for (int i = 0; i < nbr_bondnum_1; i++) {
            String[] atoms_1 = char_bond_nbrs_1.get(2 * i);
            IBond bond_1 = get_bond(ac1, int_bond_nbrs_1.get(i));
            for (int j = 0; j < nbr_bondnum_2; j++) {
                List<String> atoms_2 = Arrays.asList(char_bond_nbrs_2.get(2 * j));
                IBond bond_2 = get_bond(ac2, int_bond_nbrs_2.get(j));
                boolean match_feasible = isMatchFeasible(bond_1, bond_2, shouldMatchBonds, shouldMatchRings, matchAtomType);
                if (!match_feasible) continue;
                if (atoms_1[0].equals(atoms_2.get(0)) && atoms_1[1].equals(atoms_2.get(1))) {
                    no_map = false;
                    break;
                }
                else if (atoms_1[0].equals(atoms_2.get(1)) && atoms_1[1].equals(atoms_2.get(0))) {
                    no_map = false;
                    break;
                }
            }
            if (!no_map) break;
        }
        if (nbr_bondnum_1 == 0 || nbr_bondnum_2 == 0 || mapping_check || no_map) {
//            System.out.println("Alg3, iter = " + iter_counter);
//            System.out.println("best_mapping_size = " + best_mapping_size);
//            System.out.println("num_mapped_atoms = " + num_mapped_atoms);
//            System.out.println("List size = " + mapped_atoms.size() + "\n");
            if (num_mapped_atoms >= best_mapping_size) {
                if (num_mapped_atoms > best_mapping_size) {
                    int sz = getFinalMappings().size();
                    getFinalMappings().clear();
                    best_mapping_size = num_mapped_atoms;
//                    System.out.println("Cleaned mappings, size was " + sz + " (iter = " + local_iter + ")");
                }
                getFinalMappings().add(mapped_atoms);
//                System.out.println("Added mapping (iter = " + local_iter + ")");
//                for (Integer[] entry : mapped_atoms) System.out.println(entry[0] + ", " + entry[1]);
//                System.out.println("---");
            }
//            System.out.println("\tQuitting, condition #1");
            return 0;
        }

        int_global_1.clear();
        int_global_2.clear();
        char_global_1.clear();
        char_global_2.clear();

        // redefining of global vectors and variables
        int_global_1.addAll(int_bond_nbrs_1);
        int_global_2.addAll(int_bond_nbrs_2);
        char_global_1.addAll(char_bond_nbrs_1);
        char_global_2.addAll(char_bond_nbrs_2);
        this.MARCS.clear();
        this.MARCS = new ArrayList<>(nbr_bondnum_1 * nbr_bondnum_2);
        for (int i = 0; i < nbr_bondnum_1 * nbr_bondnum_2; i++)
            MARCS.add(i, 0);
        for (int i = 0; i < nbr_bondnum_1; i++) {
            String[] atoms_1 = char_bond_nbrs_1.get(2 * i);
            for (int j = 0; j < nbr_bondnum_2; j++) {
                String[] atoms_2 = char_bond_nbrs_2.get(2 * j);
                if (Arrays.equals(atoms_1, atoms_2) || Arrays.equals(atoms_1, new String[] {atoms_2[1], atoms_2[0]}))
                    MARCS.set(i * nbr_bondnum_2 + j, 1);
            }
        }
        first = last = new BinaryTree();
        last.setValue(-1);
        bestarcsleft = 0;
        start_search(nbr_bondnum_1, nbr_bondnum_2);
        Stack<List<Integer>> BESTARCS_copy = (Stack<List<Integer>>) BESTARCS.clone();
        BESTARCS.clear();

        while (!BESTARCS_copy.empty()) {
            List<Integer> MARCS_vector = BESTARCS_copy.peek();
            List<Integer[]> new_mapping = find_mcgregor_MAPPING(
                    MARCS_vector, num_mapped_atoms, mapped_atoms, nbr_bondnum_1, int_bond_nbrs_1, nbr_bondnum_2, int_bond_nbrs_2
            );
            boolean no_further_mapping = num_mapped_atoms == mapped_atoms.size();
            List<Integer> unmapped_atoms_1 = find_unmapped_atoms(new_mapping, new_mapping.size(), num_heavy_atoms_1, 0);
            List<String[]> char_set_1_copy = new ArrayList<>();
            for (String[] s : char_bond_set_1) char_set_1_copy.add(new String[]{s[0], s[1]});
            List<String[]> char_set_2_copy = copy_existing_char_table(char_bond_set_2, set_bondnum_2);
            List<Integer[]> new_int_nbrs_1 = new ArrayList<>();
            List<String[]> new_char_nbrs_1 = new ArrayList<>();
            List<Integer[]> new_int_bond_set_1 = new ArrayList<>();
            List<String[]> new_char_bond_set_1 = new ArrayList<>();
            int SR_count = extract_bonds(
                    0, set_bondnum_1, set_bondnum_2, 0, unmapped_atoms_1, new_mapping,
                    int_bond_set_1, int_bond_set_2, char_set_1_copy, char_set_2_copy, new_int_nbrs_1, new_char_nbrs_1,
                    new_int_bond_set_1, new_char_bond_set_1, false
            );
            List<Integer> unmapped_atoms_2 = find_unmapped_atoms(new_mapping, new_mapping.size(), num_heavy_atoms_2, 1);
            List<Integer[]> new_int_nbrs_2 = new ArrayList<>();
            List<String[]> new_char_nbrs_2 = new ArrayList<>();
            List<Integer[]> new_int_bond_set_2 = new ArrayList<>();
            List<String[]> new_char_bond_set_2 = new ArrayList<>();
            extract_bonds(
                    SR_count, set_bondnum_2, new_int_nbrs_1.size(), 1, unmapped_atoms_2, new_mapping,
                    int_bond_set_2, new_int_nbrs_1, char_set_2_copy, new_char_nbrs_1, new_int_nbrs_2, new_char_nbrs_2,
                    new_int_bond_set_2, new_char_bond_set_2, false
            );
            Iterator(
                    no_further_mapping, new_mapping.size(), new_mapping, new_int_nbrs_1.size(), new_int_nbrs_2.size(),
                    new_int_nbrs_1, new_int_nbrs_2, new_char_nbrs_1, new_char_nbrs_2,
                    new_char_bond_set_1.size() / 2, new_char_bond_set_2.size() / 2, new_int_bond_set_1, new_int_bond_set_2,
                    new_char_bond_set_1, new_char_bond_set_2
            );
            if (!BESTARCS_copy.isEmpty()) BESTARCS_copy.pop();
        }
//        System.out.println("\tQuitting, condition #2");
        return 0;
    }

    private static void print_str_list(List<String[]> list, String alias, boolean add_dashes) {
        if (list.isEmpty()) return;
        for (String[] x : list)
            System.out.println("\t" + alias + ": " + x[0] + ", " + x[1]);
        if (add_dashes) System.out.println("----------");
    }

    private static void print_int_list(List<Integer[]> list, String alias, boolean add_dashes) {
        if (list.isEmpty()) return;
        for (Integer[] x : list) {
            System.out.print("\t" + alias + ": ");
            for (int i = 0; i < x.length; i++)
                System.out.print(x[i] + (i == x.length - 1 ? "\n" : ", "));
        }
        if (add_dashes) System.out.println("----------");
    }

    private static boolean list_contains_simple(List<Integer[]> list, Integer[] entry) {
        // for pairs only
        for (Integer[] x : list) if (x[0].equals(entry[0]) && x[1].equals(entry[1])) return true;
        return false;
    }

    private int extract_bonds(
            int SR_count, int num_bonds, int num_bonds_2, int atom_index,
            List<Integer> unmapped_atoms, List<Integer[]> mapped_atoms,
            List<Integer[]> int_table, List<Integer[]> int_table_2,
            List<String[]> char_table_1_copy, List<String[]> char_table_2_copy,
            List<Integer[]> int_bond_nbrs, List<String[]> char_bond_nbrs,
            List<Integer[]> int_bond_set, List<String[]> char_bond_set, boolean troubled
    ) {
        for (int i = 0; i < num_bonds; i++) {
            Integer[] table_entry = int_table.get(i);
            boolean unmapped_contains_fst = unmapped_atoms.contains(table_entry[0]);
            if (!unmapped_contains_fst && !unmapped_atoms.contains(table_entry[1])) continue;
            boolean normal_bond = true;
            int unmapped_ind = unmapped_contains_fst ? 0 : 1;
            int other_ind = (unmapped_ind + 1) % 2;
            for (Integer[] mapped_pair : mapped_atoms)
                if (mapped_pair[atom_index].equals(table_entry[other_ind])) {
                    int_bond_nbrs.add(table_entry);
                    if (char_table_1_copy.get(2 * i + 1)[other_ind].equals("X")) {
                        String[] fst_char_item = new String[2];
                        String[] snd_char_item = new String[2];
                        fst_char_item[unmapped_ind] = char_table_1_copy.get(2 * i)[unmapped_ind];
                        fst_char_item[other_ind] = SignROW.get(SR_count);
                        snd_char_item[unmapped_ind] = "X";
                        snd_char_item[other_ind] = char_table_1_copy.get(2 * i)[other_ind];
                        char_bond_nbrs.add(fst_char_item);
                        char_bond_nbrs.add(snd_char_item);
                        char_table_1_copy = change_char_bonds(table_entry[other_ind], SignROW.get(SR_count), num_bonds, int_table, char_table_1_copy);
                        int corr_atom = search_corresponding_atom(table_entry[other_ind], atom_index, mapped_atoms);
                        if (char_table_2_copy != null) {
                            char_table_2_copy = change_char_bonds(corr_atom, SignROW.get(SR_count), num_bonds_2, int_table_2, char_table_2_copy);
                        }
                        SR_count++;
                    }
                    else {
                        char_bond_nbrs.add(char_table_1_copy.get(2 * i));
                        char_bond_nbrs.add(char_table_1_copy.get(2 * i + 1));
                    }
                    normal_bond = false;
//                    break; // WTF???
                }
            if (normal_bond) {
                int_bond_set.add(table_entry);
                char_bond_set.add(char_table_1_copy.get(2 * i));
                char_bond_set.add(new String[] {"X", "X"});
            }
        }
        return SR_count;
    }

    private List<Integer[]> find_mcgregor_MAPPING(
            List<Integer> MARCS_vector, int mapped_atoms_num, List<Integer[]> current_MAPPING,
            int bondnum_A, List<Integer[]> i_bonds_A, int bondnum_B, List<Integer[]> i_bonds_B) {
//        print_int_list(current_MAPPING, "current_mapping", true);
//        print_int_list(i_bonds_A, "bonds_A", true);
//        print_int_list(i_bonds_B, "bonds_B", true);
        List<Integer[]> additional_mapping = new ArrayList<>();
        for (int i = 0; i < bondnum_A; i++) {
            for (int j = 0; j < bondnum_B; j++) {
                if (MARCS_vector.get(i * bondnum_B + j) != 1) continue;
                Integer[] mol_A = i_bonds_A.get(i);
                Integer[] mol_B = i_bonds_B.get(j);
                for (int im = 0; im < mapped_atoms_num; im++) {
                    Integer[] mapping = current_MAPPING.get(im);
                    if (mapping[0] == mol_A[0] && mapping[1] == mol_B[0]) {
                        additional_mapping.add(new Integer[]{mol_A[1], mol_B[1]});
                    }
                    if (mapping[0] == mol_A[0] && mapping[1] == mol_B[1]) {
                        additional_mapping.add(new Integer[]{mol_A[1], mol_B[0]});
                    }
                    if (mapping[0] == mol_A[1] && mapping[1] == mol_B[0]) {
                        additional_mapping.add(new Integer[]{mol_A[0], mol_B[1]});
                    }
                    if (mapping[0] == mol_A[1] && mapping[1] == mol_B[1]) {
                        additional_mapping.add(new Integer[]{mol_A[0], mol_B[0]});
                    }
                }
            }
        }
//        for (int k = 0; k < MARCS_vector.size(); k++)
//            if (MARCS_vector.get(k) == 1) {
//                int j = k % bondnum_B;
//                int i = (k - j) / bondnum_B;
//                Integer[] mol_A = i_bonds_A.get(i);
//                Integer[] mol_B = i_bonds_B.get(j);
//                for (int im = 0; im < mapped_atoms_num; im++) {
//                    Integer[] mapping = current_MAPPING.get(im);
//                    if (mapping[0] == mol_A[0] && mapping[1] == mol_B[0]) {
//                        additional_mapping.add(new Integer[]{mol_A[1], mol_B[1]});
//                    }
//                    if (mapping[0] == mol_A[0] && mapping[1] == mol_B[1]) {
//                        additional_mapping.add(new Integer[]{mol_A[1], mol_B[0]});
//                    }
//                    if (mapping[0] == mol_A[1] && mapping[1] == mol_B[0]) {
//                        additional_mapping.add(new Integer[]{mol_A[0], mol_B[1]});
//                    }
//                    if (mapping[0] == mol_A[1] && mapping[1] == mol_B[1]) {
//                        additional_mapping.add(new Integer[]{mol_A[0], mol_B[0]});
//                    }
//                }
//            }
        current_MAPPING.addAll(additional_mapping);
//        System.out.println("New mapping before removal");
//        for (Integer[] p : current_MAPPING) System.out.println("\t" + p[0] + ", " + p[1]);
        return remove_recurring_mappings(current_MAPPING);
    }

    private List<Integer[]> remove_recurring_mappings(List<Integer[]> mappings) {
        List<Integer[]> res = new ArrayList<>();
        for (int i = 0; i < mappings.size(); i++) {
            boolean exist = false;
            Integer[] i_map = mappings.get(i);
            for (int j = i + 1; j < mappings.size(); j++) {
                Integer[] j_map = mappings.get(j);
                if (i_map[0].equals(j_map[0])) {
                    exist = true;
                    break;
                }
            }
            if (!exist) res.add(i_map);
        }
        return res;
    }

    private List<Integer> remove_redundant_arcs(int row, int col, List<Integer> MARCS, int nbr_bondnum_1, int nbr_bondnum_2) {
        List<Integer> res = new ArrayList<>(MARCS);
        Integer[] row_atoms = int_global_1.get(row);
        Integer[] col_atoms = int_global_2.get(col);
        for (int i = 0; i < nbr_bondnum_1; i++) {
            List<Integer> i_atoms = new ArrayList<>();
            i_atoms.add(int_global_1.get(i)[0]);
            i_atoms.add(int_global_1.get(i)[1]);
            boolean row_match = i_atoms.contains(row_atoms[0]) || i_atoms.contains(row_atoms[1]);
            for (int j = 0; j < nbr_bondnum_2; j++) {
                List<Integer> j_atoms = new ArrayList<>();
                j_atoms.add(int_global_2.get(j)[0]);
                j_atoms.add(int_global_2.get(j)[1]);
                boolean col_match = j_atoms.contains(col_atoms[0]) || j_atoms.contains(col_atoms[1]);
                if ((row_match && !col_match) || (col_match && !row_match)) res.set(i * nbr_bondnum_2 + j, 0);
            }
        }
        for (int i = 0; i < nbr_bondnum_1; i++) res.set(i * nbr_bondnum_2 + col, 0);
        for (int j = 0; j < nbr_bondnum_2; j++) res.set(row * nbr_bondnum_2 + j, 0);
        res.set(row * nbr_bondnum_2 + col, 1);
        return res;
    }

    private void start_search(int nbr_bondnum_1, int nbr_bondnum_2) {
        int i = 0, j = 0;
        for (int k = 0; k < MARCS.size(); k++)
            if (MARCS.get(k) == 1) {
                j = k % nbr_bondnum_2;
                i = (k - j) / nbr_bondnum_2;
                break;
            }
        if (MARCS.get(i * nbr_bondnum_2 + j) == 0) {
            part_search(i, j, MARCS, nbr_bondnum_1, nbr_bondnum_2);
        }
        if (MARCS.get(i * nbr_bondnum_2 + j) != 0) {
            part_search(i, j, MARCS, nbr_bondnum_1, nbr_bondnum_2);
            MARCS.set(i * nbr_bondnum_2 + j, 0);
            part_search(i, j, MARCS, nbr_bondnum_1, nbr_bondnum_2);
        }
    }

    private void part_search(int x0, int y0, List<Integer> TEMPMARCS, int nbr_bondnum_1, int nbr_bondnum_2) {
        int x = x0, y = y0;
        if (TEMPMARCS.get(x0 * nbr_bondnum_2 + y0) == 1) {
            TEMPMARCS = remove_redundant_arcs(x0, y0, TEMPMARCS, nbr_bondnum_1, nbr_bondnum_2);
            int arcsleft = (int) TEMPMARCS.stream().filter(s -> s == 1).count();
            if (arcsleft >= bestarcsleft) {
                do {
                    y++;
                    if (y == nbr_bondnum_2) {
                        y = 0;
                        x++;
                    }
                } while (x < nbr_bondnum_1 && TEMPMARCS.get(x * nbr_bondnum_2 + y) != 1);
                if (x < nbr_bondnum_1) {
                    part_search(x, y, TEMPMARCS, nbr_bondnum_1, nbr_bondnum_2);
                    TEMPMARCS.set(x * nbr_bondnum_2 + y, 0);
                    part_search(x, y, TEMPMARCS, nbr_bondnum_1, nbr_bondnum_2);
                }
                else {
                    if (arcsleft > bestarcsleft) {
                        BinaryTree.remove_tree_structure(first);
                        first = last = new BinaryTree();
                        last.setValue(-1);
                        last.equal = null;
                        last.not_equal = null;
                        BESTARCS.clear();
                    }
                    bestarcsleft = arcsleft;
                    if (check_MARCS(TEMPMARCS, nbr_bondnum_1, nbr_bondnum_2)) BESTARCS.push(TEMPMARCS);
                }
            }
        }
        else {
            do {
                y++;
                if (y == nbr_bondnum_2) {
                    y = 0;
                    x++;
                }
            } while (x < nbr_bondnum_1 && TEMPMARCS.get(x * nbr_bondnum_2 + y) != 1);
            if (x < nbr_bondnum_1) {
                part_search(x, y, TEMPMARCS, nbr_bondnum_1, nbr_bondnum_2);
                TEMPMARCS.set(x * nbr_bondnum_2 + y, 0);
                part_search(x, y, TEMPMARCS, nbr_bondnum_1, nbr_bondnum_2);
            }
            else {
                int arcsleft = (int) TEMPMARCS.stream().filter(s -> s == 1).count();
                if (arcsleft >= bestarcsleft) {
                    if (arcsleft > bestarcsleft) {
                        BinaryTree.remove_tree_structure(first);
                        first = last = new BinaryTree();
                        last.setValue(-1);
                        last.equal = null;
                        last.not_equal = null;
                        BESTARCS.clear();
                    }
                    bestarcsleft = arcsleft;
                    if (check_MARCS(TEMPMARCS, nbr_bondnum_1, nbr_bondnum_2)) BESTARCS.push(TEMPMARCS);
                }
            }
        }
    }

    private boolean check_MARCS(List<Integer> MARCS, int nbr_bondnum_1, int nbr_bondnum_2) {
        List<Integer> postnum_list = new ArrayList<>(nbr_bondnum_1 * nbr_bondnum_2);
        for (int i = 0; i < nbr_bondnum_1 * nbr_bondnum_2; i++) postnum_list.add(i, 0);
        int j = 0;
        int count = 0;
        for (int i = 0; i < postnum_list.size(); i++)
            if (MARCS.get(i) == 1) {
                postnum_list.set(j++, i);
                count++;
            }
        verify_nodes(postnum_list, first, 0, count);
        return new_matrix;
    }

    private boolean verify_nodes(List<Integer> matrix, BinaryTree curr_structure, int x, int field_length) {
        int item = matrix.get(x);
        if (item == curr_structure.getValue() && x < field_length && curr_structure.equal != null) {
            new_matrix = false;
            verify_nodes(matrix, curr_structure.equal, x + 1, field_length);
        }
        if (item != curr_structure.getValue()) {
            if (curr_structure.not_equal != null) verify_nodes(matrix, curr_structure.not_equal, x, field_length);
            else {
                curr_structure.not_equal = new BinaryTree();
                curr_structure.not_equal.setValue(item);
                BinaryTree last_one = curr_structure.not_equal;
                for (int i = 0; i + x + 1 < field_length; i++) {
                    last_one.equal = new BinaryTree();
                    last_one = last_one.equal;
                    last_one.setValue(matrix.get(i + x + 1));
                    last_one.not_equal = null;
                }
                last_one.equal = null;
                new_matrix = true;
            }
        }
       return true;
    }

    private IBond get_bond(IAtomContainer ac, Integer[] inds) {
        IAtom a1 = ac.getAtom(inds[0] - 1);
        IAtom a2 = ac.getAtom(inds[1] - 1);
        return ac.getBond(a1, a2);
    }

    public List<List<Integer[]>> getFinalMappings() {
        return final_mappings;
    }
}
