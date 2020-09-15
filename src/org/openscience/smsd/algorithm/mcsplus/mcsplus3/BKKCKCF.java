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
import java.util.Stack;
import java.util.stream.Stream;

import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import org.xmlcml.euclid.Int;
import uk.ac.ebi.reactionblast.graphics.direct.MoleculeLabelDrawer;

/**
 *
 * @author Syed Asad Rahman, EMBL-EBI, Cambridge, UK
 * @contact asad@ebi.ac.uk
 */
public class BKKCKCF {

    private static final ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(MoleculeLabelDrawer.class);

    private final List<Integer[]> comp_graph_nodes;
    private final List<Integer[]> c_edges;
    private final List<Integer[]> d_edges;
    private final Stack<List<Integer>> max_Cliques_Set;
    private int best_clique_size;

    /*
     *T: is a set of vertices which have already been used for the
     * initialization of ENUMERATE_CLIQUES
     */
    protected final List<Integer> T;

    /**
     *
     * @param comp_graph_nodes
     * @param C_edges
     * @param D_edges
     */
    public BKKCKCF(List<Integer[]> comp_graph_nodes, List<Integer[]> C_edges, List<Integer[]> D_edges) {
        this.comp_graph_nodes = comp_graph_nodes;
        this.c_edges = C_edges;
        this.d_edges = D_edges;
        this.best_clique_size = 0;
        this.max_Cliques_Set = new Stack<>();
        this.T = new Stack<>();
    }

    /*
   
     * R: set of vertices belonging to the current clique
     
     * X: set of vertices which are not allowed to be added
     * to R, defined as X in paper
    
     * P: is a set of vertices which <b>can</b> be added to R, because they are
     * neighbours of vertex u via <i>c-edges</i>
    
     * Q: is a set of vertices which <b>cannot</b> be added to R, because they are
     * neighbours of vertex u via <i>d-edges</i>
     
     * V: stored all the vertices for the Graph G
     * V[G]: nodes of vector comp_graph_nodes are stored in V
     
     */
    int init_Algorithm() {
        List<Integer> R = new ArrayList<>();
        Stack<Integer> Q = new Stack<>();
        List<Integer> X = new ArrayList<>();
        List<Integer> N1 = new ArrayList<>();
        List<Integer> N2 = new ArrayList<>();
        Stack<Integer> P = new Stack<>();

        // nodes of vector comp_graph_nodes are stored in V
        Stack<Integer> V = new Stack<>();//Initialization of Stack V
        int V_size = comp_graph_nodes.size();
        for (int i = 1; i <= V_size; i++) V.push(i);
        V.push(0);

        int b = 0;

        while (V.get(b) != 0) { // V[b] is node u
            int central_node = V.get(b);

            P.clear();
            Q.clear();
            X.clear();
            N1.clear();
            N2.clear();
            R.clear();

            // find the neighbors of the central node from V
            N1 = find_neighbors(central_node, c_edges);
            N2 = find_neighbors(central_node, d_edges);
            for (Integer val : N1) {
                if (T.contains(val)) X.add(val);
                else P.push(val);
                int nbr_pos = find_nbr_pos(V, val, -1);
                if (nbr_pos != -1) {
                    for (int i = nbr_pos; i < V.size() - 1; i++) V.set(i, V.get(i + 1));
                    V.pop();
                    if (nbr_pos < b) b--;
                }
            }
            for (Integer val : N2) {
                Q.push(val);
                int nbr_pos = find_nbr_pos(V, val, -1);
                if (nbr_pos != -1) {
                    for (int i = nbr_pos; i < V.size() - 1; i++) V.set(i, V.get(i + 1));
                    V.pop();
                    if (nbr_pos < b) b--;
                }
            }
            P.add(0);
            R.add(central_node);
            Enumerate_Cliques(R, P, Q, X);
            T.add(central_node);
            b++;
        }

        return 0;
    }

    private static int find_nbr_pos(Stack<Integer> V, int val, int def) {
        for (int i = 0; i < V.size(); i++)
            if (V.get(i).equals(val)) return i;
        return def;
    }

    private int Enumerate_Cliques(List<Integer> R, Stack<Integer> P, Stack<Integer> Q, List<Integer> X) {
        List<Integer> N1 = new ArrayList<>();
        List<Integer> N2 = new ArrayList<>();
        Stack<Integer> P_Prime = new Stack<>();

        P_Prime.addAll(P);

        List<Integer> R_copy = new ArrayList<>();
        Stack<Integer> P_copy = new Stack<>();
        Stack<Integer> Q_copy = new Stack<>();
        List<Integer> X_copy = new ArrayList<>();

        if (P.size() == 1 && X.isEmpty()) {
            int clique_size = R.size();
            if (clique_size >= best_clique_size) {
                if (clique_size > best_clique_size) {
                    max_Cliques_Set.clear();
                    best_clique_size = clique_size;
                }
                max_Cliques_Set.push(R);
            }
            return 0;
        }
        int a = 0;
        while (P_Prime.get(a) != 0) {

            int ui = P_Prime.get(a);
            int P_size = P.size();
            int ut_node_pos = find_nbr_pos(P, ui, Integer.MAX_VALUE);
            if (ut_node_pos == Integer.MAX_VALUE) {
                LOGGER.debug("ut_node_pos = " + Integer.MAX_VALUE);
            }
            // delete P_Prime node in P
            for (int counter = ut_node_pos; counter < P_size - 1; counter++)
                P.setElementAt(P.get(counter + 1), counter);
            P.pop(); //TO DO

            R_copy.clear();
            P_copy.clear();
            Q_copy.clear();
            X_copy.clear();
            N1.clear();
            N2.clear();

            R_copy.addAll(R);
            P_copy.addAll(P);
            Q_copy.addAll(Q);
            X_copy.addAll(X);
            P_copy.pop();

            // find the neighbors of the central node from P
            N1 = find_neighbors(ui, c_edges);
            N2 = find_neighbors(ui, d_edges);

            for (Integer n1_val : N1) {
                for (Integer q_val : Q) {
                    if (!n1_val.equals(q_val)) continue;
                    if (T.contains(n1_val)) X_copy.add(n1_val);
                    else P_copy.push(n1_val);
                    int nbr_pos = find_nbr_pos(Q_copy, n1_val, Integer.MAX_VALUE);
                    for (int i = nbr_pos; i < Q_copy.size() - 1; i++)
                        Q_copy.set(i, Q_copy.get(i + 1));
                    Q_copy.pop();
                }
                int ut_set_size = P_Prime.size();
                int nbr_pos = find_nbr_pos(P_Prime, n1_val, -1);
                if (nbr_pos != -1) {
                    for (int i = nbr_pos; i < ut_set_size - 1; i++)
                        P_Prime.setElementAt(P_Prime.get(i + 1), i);
                    P_Prime.pop();
                    if (nbr_pos < a) a--;
                }
            }
            for (Integer n2_val : N2) {
                int ut_set_size = P_Prime.size();
                int nbr_pos = find_nbr_pos(P_Prime, n2_val, -1);
                if (nbr_pos != -1) {
                    for (int i = nbr_pos; i < ut_set_size - 1; i++)
                        P_Prime.setElementAt(P_Prime.get(i + 1), i);
                    P_Prime.pop();
                    if (nbr_pos < a) a--;
                }
            }

            Stack<Integer> P_copy_N_intersec = new Stack<>();
            Stack<Integer> Q_copy_N_intersec = new Stack<>();
            List<Integer> X_copy_N_intersec = new ArrayList<>();
            for (Integer nval : Stream.concat(N1.stream(), N2.stream()).toArray(Integer[]::new)) {
                if (P_copy.contains(nval)) P_copy_N_intersec.push(nval);
                if (Q_copy.contains(nval)) Q_copy_N_intersec.push(nval);
                if (X_copy.contains(nval)) X_copy_N_intersec.add(nval);
            }

            P_copy_N_intersec.push(0);
            R_copy.add(ui);
            Enumerate_Cliques(R_copy, P_copy_N_intersec, Q_copy_N_intersec, X_copy_N_intersec);
            X.add(ui);
            a++;
        }

        return 0;
    }

    private List<Integer> find_neighbors(int central_node, List<Integer[]> edges) {
        List<Integer> res = new ArrayList<>();
        for (Integer[] edge : edges)
            if (edge[0].equals(central_node)) res.add(edge[1]);
            else if (edge[1].equals(central_node)) res.add(edge[0]);
        return res;
    }

    public Stack<List<Integer>> getMax_Cliques_Set() {
        return max_Cliques_Set;
    }
}
