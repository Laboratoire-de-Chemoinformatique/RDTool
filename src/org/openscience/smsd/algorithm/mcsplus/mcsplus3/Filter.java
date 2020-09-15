/* Copyright (C) 2009-2018  Syed Asad Rahman <asad at ebi.ac.uk>
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

/**
 *
 * @author Syed Asad Rahman <asad at ebi.ac.uk>
 */
public class Filter extends McGregor {

    /**
     *
     * @param f1
     * @param f2
     * @param shouldMatchBonds
     * @param shouldMatchRings
     * @param matchAtomType
     */
    public Filter(IAtomContainer f1,IAtomContainer f2, boolean shouldMatchBonds, boolean shouldMatchRings, boolean matchAtomType) {
        super(f1, f2,  shouldMatchBonds,  shouldMatchRings,  matchAtomType);
    }

    private void find_co_groups(
            List<Integer> C_inds, List<Integer> carb_vec,
            List<Integer[]> int_table, List<String[]> char_table,
            int num_bonds, int num_heavy_atoms, boolean debug) {
        for (int ia = 0; ia < num_heavy_atoms; ia++) {
            int O_num = 0;
            boolean c_group = true;
            List<Integer> vector = new ArrayList<>();

            for (int ib = 0; ib < num_bonds; ib++) {
                String[] char_table_entry = char_table.get(ib);
                if (!Arrays.asList(char_table_entry).containsAll(Arrays.asList("C", "O"))) continue;
                Integer[] int_table_entry = int_table.get(ib);
                int c_ind = char_table_entry[0].equals("C") ? 0 : 1;
                int o_ind = (c_ind + 1) % 2;
                if (ia + 1 != int_table_entry[c_ind]) continue;

                if (check(int_table_entry[o_ind], num_bonds, int_table)) {
                    vector.add(int_table_entry[o_ind]);
                    vector.add(int_table_entry[2]);
                    O_num++;
                } else c_group = false;
            }
            if (O_num == 2 && c_group) {
                C_inds.add(ia + 1);
                carb_vec.add(vector.get(1) != 2 ? vector.get(2) : vector.get(0));
            }
        }
    }

    private void find_po_groups(
            List<Integer> P_inds, List<Integer> phos_vec,
            List<Integer[]> int_table, List<String[]> char_table,
            int num_bonds, int num_heavy_atoms) {
        for (int ib = 0; ib < num_bonds; ib++) {
            String[] char_table_entry = char_table.get(ib);
            if (!Arrays.asList(char_table_entry).containsAll(Arrays.asList("P", "O"))) continue;

            int P_num = 0;
            List<Integer[]> vector = new ArrayList<>();
            Integer[] int_table_entry = int_table.get(ib);
            int p_ind = char_table_entry[0].equals("P") ? 0 : 1;
            int o_ind = (p_ind + 1) % 2;
            for (int ia = 0; ia < num_heavy_atoms; ia++) {
                if (ia + 1 != int_table_entry[p_ind]) continue;
                if (check(int_table_entry[o_ind], num_bonds, int_table)) {
                    vector.add(new Integer[]{int_table_entry[o_ind], int_table_entry[2]});
                    P_num++;
                }
            }
            if (P_num > 0 && P_num < 5) {
                P_inds.add(int_table_entry[p_ind]);
                boolean double_found = false;
                for (Integer[] v : vector)
                    if (v[1] == 2) {
                        phos_vec.add(v[0]);
                        double_found = true;
                        break;
                    }
                if (!double_found) phos_vec.add(vector.get(0)[0]);
            }
        }
    }

    private void find_amino_groups(
            List<Integer> N_inds, List<Integer> amino_vec,
            List<Integer[]> int_table, List<String[]> char_table,
            int num_bonds, int num_heavy_atoms) {
        for (int ib = 0; ib < num_bonds; ib++) {
            String[] char_table_entry = char_table.get(ib);
            if (!Arrays.asList(char_table_entry).containsAll(Arrays.asList("C", "N"))) continue;

            int N_num = 0;
            List<Integer[]> vector = new ArrayList<>();
            Integer[] int_table_entry = int_table.get(ib);
            int n_ind = char_table_entry[0].equals("N") ? 0 : 1;
            int c_ind = (n_ind + 1) % 2;
            for (int ia = 0; ia < num_heavy_atoms; ia++) {
                if (ia + 1 != int_table_entry[c_ind]) continue;
                if (check(int_table_entry[n_ind], num_bonds, int_table)) {
                    vector.add(new Integer[]{int_table_entry[n_ind], int_table_entry[2]});
                    N_num++;
                }
            }
            if (N_num > 1 && N_num < 4) {
                N_inds.add(int_table_entry[c_ind]);
                boolean double_found = false;
                for (Integer[] v : vector)
                    if (v[1] == 2) {
                        amino_vec.add(v[0]);
                        double_found = true;
                        break;
                    }
                if (!double_found) amino_vec.add(vector.get(0)[0]);
            }
        }
    }

    private void find_so_groups(
            List<Integer> S_inds, List<List<Integer>> sulf_vec,
            List<Integer[]> int_table, List<String[]> char_table,
            int num_bonds, int num_heavy_atoms, List<IAtom> atoms, boolean debug) {
        for (int ia = 0; ia < num_heavy_atoms; ia++) {
            if (!atoms.get(ia).getSymbol().equals("S")) continue;
            int O_num = 0;
            List<Integer[]> vector = new ArrayList<>();
            for (int ib = 0; ib < num_bonds; ib++) {
                String[] char_table_entry = char_table.get(ib);
                Integer[] int_table_entry = int_table.get(ib);
                int o_ind = char_table_entry[0].equals("O") ? 0 : 1;
                int s_ind = (o_ind + 1) % 2;
                if (ia + 1 != int_table_entry[s_ind]) continue;
                if (check(int_table_entry[o_ind], num_bonds, int_table)) {
                    vector.add(new Integer[]{int_table_entry[o_ind], int_table_entry[2]});
                    O_num++;
                }
            }
            if (O_num > 1 && O_num < 5) S_inds.add(ia + 1);
            if (O_num == 2) {
                List<Integer> temp_vec = new ArrayList<>();
                if (vector.get(0)[1] != 2) {
                    temp_vec.add(vector.get(1)[0]);
                    temp_vec.add(vector.get(0)[0]);
                } else {
                    temp_vec.add(vector.get(0)[0]);
                    temp_vec.add(vector.get(1)[0]);
                }
                sulf_vec.add(temp_vec);
            }
            if (O_num == 3) {
                boolean no_single_bond = true;
                List<Integer> temp_vec = new ArrayList<>();
                if (vector.get(2)[1] == 1) {
                    temp_vec.add(vector.get(0)[0]);
                    temp_vec.add(vector.get(1)[0]);
                    no_single_bond = false;
                }
                if (vector.get(0)[1] == 1) {
                    temp_vec.add(vector.get(1)[0]);
                    temp_vec.add(vector.get(2)[0]);
                    temp_vec.add(vector.get(0)[0]);
                    no_single_bond = false;
                }
                if (vector.get(1)[1] == 1) {
                    temp_vec.add(vector.get(0)[0]);
                    temp_vec.add(vector.get(2)[0]);
                    temp_vec.add(vector.get(1)[0]);
                    no_single_bond = false;
                }
                if (no_single_bond)
                    for (Integer[] entry : vector)
                        temp_vec.add(entry[0]);
                sulf_vec.add(temp_vec);
            }
            if (O_num == 4) {
                List<Integer> single_bond_pos = new ArrayList<>();
                List<Integer> double_bond_pos = new ArrayList<>();
                for (Integer[] v : vector)
                    if (v[1] == 1) single_bond_pos.add(v[0]);
                    else if (v[1] == 2) double_bond_pos.add(v[0]);
                List<Integer> temp_vec = new ArrayList<>();
                temp_vec.addAll(double_bond_pos);
                temp_vec.addAll(single_bond_pos);
                sulf_vec.add(temp_vec);
            }
        }
    }

    private void find_no3_groups(
            List<Integer> N_inds, List<List<Integer>> nitro_vec,
            List<Integer[]> int_table, List<String[]> char_table,
            int num_bonds, int num_heavy_atoms, List<IAtom> atoms) {
        for (int ia = 0; ia < num_heavy_atoms; ia++) {
            if (!atoms.get(ia).getSymbol().equals("N")) continue;

            int O_num = 0;
            List<Integer[]> vector = new ArrayList<>();
            for (int ib = 0; ib < num_bonds; ib++) {
                Integer[] int_table_entry = int_table.get(ib);
                String[] char_table_entry = char_table.get(ib);
                int o_ind = char_table_entry[0].equals("O") ? 0 : 1;
                int n_ind = (o_ind + 1) % 2;
                if (ia + 1 != int_table_entry[n_ind]) continue;

                if (check(int_table_entry[o_ind], num_bonds, int_table)) {
                    vector.add(new Integer[]{int_table_entry[o_ind], int_table_entry[2]});
                    O_num++;
                }
            }
            if (O_num == 2) {
                List<Integer> temp_vec = new ArrayList<>();
                if (vector.get(0)[1] != 2) {
                    temp_vec.add(vector.get(1)[0]);
                    temp_vec.add(vector.get(0)[0]);
                }
                else {
                    temp_vec.add(vector.get(0)[0]);
                    temp_vec.add(vector.get(1)[0]);
                }
                nitro_vec.add(temp_vec);
                N_inds.add(ia + 1);
            }
            else if (O_num == 3) {
                boolean no_single_bond = true;
                List<Integer> temp_vec = new ArrayList<>();
                for (Integer[] pair : vector)
                    if (pair[1] == 2)
                        temp_vec.add(pair[0]);
                for (Integer[] pair : vector)
                    if (pair[1] == 1) {
                        temp_vec.add(pair[0]);
                        no_single_bond = false;
                    }
                if (no_single_bond) {
                    temp_vec = new ArrayList<>();
                    for (Integer[] pair : vector)
                        temp_vec.add(pair[0]);
                }
                nitro_vec.add(temp_vec);
                N_inds.add(ia + 1);
            }
        }

//        for (int ib = 0; ib < num_bonds; ib++) {
//            String[] char_table_entry = char_table.get(ib);
//            if (!Arrays.asList(char_table_entry).containsAll(Arrays.asList("N", "O"))) continue;
//
//            int O_num = 0;
//            List<Integer[]> vector = new ArrayList<>();
//            Integer[] int_table_entry = int_table.get(ib);
//            int o_ind = char_table_entry[0].equals("O") ? 0 : 1;
//            int n_ind = (o_ind + 1) % 2;
//            for (int ia = 0; ia < num_heavy_atoms; ia++) {
//                if (ia + 1 != int_table_entry[n_ind]) continue;
//                if (check(int_table_entry[o_ind], num_bonds, int_table)) {
//                    vector.add(new Integer[]{int_table_entry[o_ind], int_table_entry[2]});
//                    O_num++;
//                }
//            }
//            if (O_num == 2) {
//                List<Integer> temp_vec = new ArrayList<>();
//                if (vector.get(0)[1] != 2) {
//                    temp_vec.add(vector.get(1)[0]);
//                    temp_vec.add(vector.get(0)[0]);
//                }
//                else {
//                    temp_vec.add(vector.get(0)[0]);
//                    temp_vec.add(vector.get(1)[0]);
//                }
//                nitro_vec.add(temp_vec);
//                N_inds.add(int_table_entry[n_ind]);
//            }
//            else if (O_num == 3) {
//                boolean no_single_bond = true;
//                List<Integer> temp_vec = new ArrayList<>();
//                for (Integer[] pair : vector)
//                    if (pair[1] == 2)
//                        temp_vec.add(pair[0]);
//                for (Integer[] pair : vector)
//                    if (pair[1] == 1) {
//                        temp_vec.add(pair[0]);
//                        no_single_bond = false;
//                    }
//                if (no_single_bond) {
//                    temp_vec = new ArrayList<>();
//                    for (Integer[] pair : vector)
//                        temp_vec.add(pair[0]);
//                }
//                nitro_vec.add(temp_vec);
//                N_inds.add(int_table_entry[n_ind]);
//            }
//        }
    }

    private static Integer find_index(int item, List<Integer> items) {
        for (int i = 0; i < items.size(); i++)
            if (items.get(i).equals(item)) return i;
        return -1;
    }

    private List<List<Integer[]>> find_final_mappings(
            List<List<Integer[]>> initial_mappings,
            List<Integer> inds_1, List<Integer> inds_2, List<Integer> vectors_1, List<Integer> vectors_2) {
        if (inds_1.size() == 0 || inds_2.size() == 0)
            return new ArrayList<>(initial_mappings);
        List<List<Integer[]>> res = new ArrayList<>();
        for (List<Integer[]> final_solution : initial_mappings) {
            boolean map_correct = true;
            for (int i = 0; i < best_mapping_size && map_correct; i++) {
                Integer[] solution_entry = final_solution.get(i);
                int ind1 = find_index(solution_entry[0], inds_1);
                int ind2 = find_index(solution_entry[1], inds_2);
                if (ind1 != -1 && ind2 != -1) {
                    int fst_1 = vectors_1.get(ind1);
                    int fst_2 = vectors_2.get(ind2);
                    boolean miss_map = true;
                    for (int j = 0; j < best_mapping_size; j++) {
                        Integer[] other_entry = final_solution.get(j);
                        if (other_entry[0].equals(fst_1) && other_entry[1].equals(fst_2)) {
                            miss_map = false;
                            break;
                        }
                    }
                    map_correct = !miss_map;
                }
            }
            if (map_correct) res.add(final_solution);
        }
        return res;
    }

    private List<List<Integer[]>> find_final_mappings_2(
            List<List<Integer[]>> initial_mappings,
            List<Integer> inds_1, List<Integer> inds_2,
            List<List<Integer>> vectors_1, List<List<Integer>> vectors_2, boolean debug) {
        if (inds_1.size() == 0 || inds_2.size() == 0)
             return new ArrayList<>(initial_mappings);
        boolean no_correct_mappings = true;
        List<List<Integer[]>> res = new ArrayList<>();
        for (List<Integer[]> final_solution : initial_mappings) {
            boolean map_correct = true;
            for (int i = 0; i < best_mapping_size && map_correct; i++) {
                Integer[] solution_entry = final_solution.get(i);
                int ind1 = inds_1.indexOf(solution_entry[0]);
                int ind2 = inds_2.indexOf(solution_entry[1]);
                if (ind1 != -1 && ind2 != -1) {
                    int fst_1 = vectors_1.get(ind1).get(0);
                    int snd_1 = vectors_1.get(ind1).get(1);
                    int fst_2 = vectors_2.get(ind2).get(0);
                    int snd_2 = vectors_2.get(ind2).get(1);
                    boolean miss_map_1 = true;
                    boolean miss_map_2 = true;
                    for (int j = 0; j < best_mapping_size; j++) {
                        Integer[] other_entry = final_solution.get(j);
                        if (other_entry[0] == fst_1 && other_entry[1] == fst_2) miss_map_1 = false;
                        if (other_entry[0] == snd_1 && other_entry[1] == snd_2) miss_map_2 = false;
                        map_correct = !miss_map_1 && !miss_map_2;
                        if (map_correct) break;
                    }
                }
            }
            if (map_correct) {
//                System.out.println("Adding item, condition #1");
                res.add(final_solution);
                no_correct_mappings = false;
            }
        }
        if (no_correct_mappings) {
//            if (debug) System.out.println("No correct mapping");
            List<Integer> mol_1_Os = new ArrayList<>();
            for (int i = 0; i < inds_1.size(); i++)
                mol_1_Os.addAll(vectors_1.get(i));
            List<Integer> mol_2_Os = new ArrayList<>();
            for (int i = 0; i < inds_2.size(); i++)
                mol_2_Os.addAll(vectors_2.get(i));
            List<List<Integer[]>> temp_s_f_m = new ArrayList<>();
            for (List<Integer[]> solution : initial_mappings) {
                List<Integer[]> t_map = new ArrayList<>();
                for (Integer[] pair : solution) {
                    int ind_1 = find_index(pair[0], mol_1_Os);
                    int ind_2 = find_index(pair[1], mol_2_Os);
                    if (ind_1 == -1 || ind_2 == -1) t_map.add(pair);
                }
                temp_s_f_m.add(t_map);
            }
            List<List<Integer[]>> temp_s_f_m2 = new ArrayList<>();
            for (int i = 0; i < temp_s_f_m.size(); i++) {
                List<Integer[]> map_1 = temp_s_f_m.get(i);
                boolean unique_map = true;
                for (int j = i + 1; j < temp_s_f_m.size() && unique_map; j++) {
                    List<Integer[]> map_2 = temp_s_f_m.get(j);
                    boolean map_contained = true;
                    for (Integer[] pair : map_2) {
                        map_contained = map_1.contains(pair);
                        if (!map_contained) break;
                    }
                    unique_map = !map_contained;
                }
                if (unique_map) temp_s_f_m2.add(map_1);
            }
//            System.out.println("Showing temp_s_f_m_2");
//            for (List<Integer[]> item : temp_s_f_m2) {
//                System.out.println("Next:");
//                for (Integer[] p : item) System.out.println("\t" + p[0] + ", " + p[1]);
//                System.out.println("End next");
//            }
            for (int i = 0; i < temp_s_f_m2.size(); i++) {
                List<Integer[]> map_element = temp_s_f_m2.get(i);
//                List<Integer[]> new_element = temp_s_f_m2.get(i);
                List<Integer[]> new_element = new ArrayList<>();
                for (Integer[] pair : map_element) {
                    int ind_1 = find_index(pair[0], inds_1);
                    int ind_2 = find_index(pair[1], inds_2);
                    if (ind_1 != -1 && ind_2 != -1) {
                        List<Integer> Os_1 = vectors_1.get(ind_1);
                        List<Integer> Os_2 = vectors_2.get(ind_2);
                        int min_size = Math.min(Os_1.size(), Os_2.size());
                        for (int j = 0; j < min_size; j++)
                            new_element.add(new Integer[]{Os_1.get(j), Os_2.get(j)});
                    }
                }
                new_element.addAll(0, temp_s_f_m2.get(i));
                res.add(new_element);
            }
        }
        return res;
    }

    int postfilter() {
//        System.out.println("best_mapping_size = " + best_mapping_size);
//        System.out.println("# of mappings = " + getFinalMappings().size());
        if (best_mapping_size == 0 && best_clique_size != 0) {
            java.util.Iterator<List<Integer[]>> iter = getFinalMappings().iterator();
            List<Integer[]> vec = iter.next();
            best_mapping_size = vec.size();
        }
        if (best_mapping_size == 0 && best_clique_size == 0) {
            return 0;
        }

        // 1. Search for carboxyl groups
        List<Integer> c_inds_1 = new ArrayList<>();
        List<Integer> c_inds_2 = new ArrayList<>();
        List<Integer> carb_vec_1 = new ArrayList<>();
        List<Integer> carb_vec_2 = new ArrayList<>();
        find_co_groups(c_inds_1, carb_vec_1, int_table_1, char_table_1, num_bonds_1, num_heavy_atoms_1, true);
        find_co_groups(c_inds_2, carb_vec_2, int_table_2, char_table_2, num_bonds_2, num_heavy_atoms_2, false);
//        System.out.println("Showing c_inds_1");
//        for (Integer v : c_inds_1) System.out.print(v + " ");
//        System.out.println();
//        System.out.println("Showing c_inds_2");
//        for (Integer v : c_inds_2) System.out.print(v + " ");
//        System.out.println();
//        System.out.println("Showing carb_vec_1");
//        for (Integer v : carb_vec_1) System.out.println("\t" + v);
//        System.out.println();
//        System.out.println("Showing carb_vec_2");
//        for (Integer v : carb_vec_2) System.out.println("\t" + v);
//        System.out.println();
        List<List<Integer[]>> co_final_mappings = find_final_mappings(getFinalMappings(), c_inds_1, c_inds_2, carb_vec_1, carb_vec_2);
//        System.out.println("# of CO solutions: " + co_final_mappings.size());

        // 2. Search for phospate groups
        List<Integer> p_inds_1 = new ArrayList<>();
        List<Integer> p_inds_2 = new ArrayList<>();
        List<Integer> phos_vec_1 = new ArrayList<>();
        List<Integer> phos_vec_2 = new ArrayList<>();
        find_po_groups(p_inds_1, phos_vec_1, int_table_1, char_table_1, num_bonds_1, num_heavy_atoms_1);
        find_po_groups(p_inds_2, phos_vec_2, int_table_2, char_table_2, num_bonds_2, num_heavy_atoms_2);
        List<List<Integer[]>> po_final_mappings = find_final_mappings(co_final_mappings, p_inds_1, p_inds_2, phos_vec_1, phos_vec_2);
//        System.out.println("# of PO solutions: " + po_final_mappings.size());

        // 3. Search for amino groups
        List<Integer> n_inds_1 = new ArrayList<>();
        List<Integer> n_inds_2 = new ArrayList<>();
        List<Integer> amino_vec_1 = new ArrayList<>();
        List<Integer> amino_vec_2 = new ArrayList<>();
        find_amino_groups(n_inds_1, amino_vec_1, int_table_1, char_table_1, num_bonds_1, num_heavy_atoms_1);
        find_amino_groups(n_inds_2, amino_vec_2, int_table_2, char_table_2, num_bonds_2, num_heavy_atoms_2);
        List<List<Integer[]>> amino_final_mappings = find_final_mappings(po_final_mappings, n_inds_1, n_inds_2, amino_vec_1, amino_vec_2);
//        System.out.println("# of amino solutions: " + amino_final_mappings.size());

        // 4. Search for sulfo groups
        List<Integer> s_inds_1 = new ArrayList<>();
        List<Integer> s_inds_2 = new ArrayList<>();
        List<List<Integer>> sulf_vec_1 = new ArrayList<>();
        List<List<Integer>> sulf_vec_2 = new ArrayList<>();
        find_so_groups(s_inds_1, sulf_vec_1, int_table_1, char_table_1, num_bonds_1, num_heavy_atoms_1, atomstr1, false);
        find_so_groups(s_inds_2, sulf_vec_2, int_table_2, char_table_2, num_bonds_2, num_heavy_atoms_2, atomstr2, false);
        List<List<Integer[]>> sulfo_final_mappings = find_final_mappings_2(amino_final_mappings, s_inds_1, s_inds_2, sulf_vec_1, sulf_vec_2, false);
//        System.out.println("# of sulfo solutions: " + sulfo_final_mappings.size());

        // 5. Search for nitro groups
        List<Integer> no_inds_1 = new ArrayList<>();
        List<Integer> no_inds_2 = new ArrayList<>();
        List<List<Integer>> nitro_vec_1 = new ArrayList<>();
        List<List<Integer>> nitro_vec_2 = new ArrayList<>();
        find_no3_groups(no_inds_1, nitro_vec_1, int_table_1, char_table_1, num_bonds_1, num_heavy_atoms_1, atomstr1);
//        System.out.println("Showing no_inds_1");
//        for (Integer v : no_inds_1) System.out.print(v + " ");
//        System.out.println();
//        System.out.println("Showing nitro_vec_1");
//        for (List<Integer> item : nitro_vec_1) {
//            System.out.println("Next:");
//            System.out.print("\t");
//            for (Integer v : item) System.out.print(v + " ");
//            System.out.println();
//        }
//        System.out.println();
        find_no3_groups(no_inds_2, nitro_vec_2, int_table_2, char_table_2, num_bonds_2, num_heavy_atoms_2, atomstr2);
//        System.out.println("Showing no_inds_2");
//        for (Integer v : no_inds_2) System.out.print(v + " ");
//        System.out.println();
//        System.out.println("Showing nitro_vec_2");
//        for (List<Integer> item : nitro_vec_2) {
//            System.out.println("Next:");
//            System.out.print("\t");
//            for (Integer v : item) System.out.print(v + " ");
//            System.out.println();
//        }
//        System.out.println();
        List<List<Integer[]>> nitro_final_mappings = find_final_mappings_2(sulfo_final_mappings, no_inds_1, no_inds_2, nitro_vec_1, nitro_vec_2, true);
//        System.out.println("# of nitro solutions: " + nitro_final_mappings.size());

        //6. Searching for redundant Methyl-group mappings
        getFinalMappings().clear();
        getFinalMappings().addAll(nitro_final_mappings);

//        System.out.println("Show result");
//        for (List<Integer[]> item : nitro_final_mappings) {
//            System.out.println("Next:");
//            for (Integer[] p : item) System.out.println("\t" + p[0] + ", " + p[1]);
//            System.out.println("End next");
//        }

        return 0;
    }

    private boolean check(int atom, int num_bonds, List<Integer[]> int_table) {
        int count_nbrs = 0;
        for (int i = 0; i < num_bonds; i++) {
            Integer[] table_entry = int_table.get(i);
            if (table_entry[0] == atom || table_entry[1] == atom) {
                count_nbrs++;
                if (count_nbrs > 1) return false;
            }
        }
        return count_nbrs == 1;
    }
}
