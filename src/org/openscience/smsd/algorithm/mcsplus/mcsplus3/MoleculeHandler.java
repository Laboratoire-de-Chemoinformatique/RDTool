package org.openscience.smsd.algorithm.mcsplus.mcsplus3;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

public class MoleculeHandler {
    private final IAtomContainer ac;
    private int num_atoms;
    private int num_bonds;
    private int num_heavy_atoms;
    private List<IAtom> atomString;
    public List<Integer[]> int_table = new LinkedList<>();
    public List<String[]> char_table = new LinkedList<>();
    private List<Integer[]> specified_int_tab = new LinkedList<>();

    private List<Integer[]> int_tab = new LinkedList<>();
    private List<String[]> char_tab = new LinkedList<>();
    private final boolean matchBonds;

    MoleculeHandler(IAtomContainer atomContainer, boolean matchBonds) {
        this.ac = atomContainer;
        this.matchBonds = matchBonds;
        this.num_bonds = atomContainer.getBondCount();
        this.num_atoms = atomContainer.getAtomCount();

        set_atom_string();
        set_connection_tables();
        set_H_number();
        boolean is_H2 = num_atoms == 2 && num_bonds == 1 && atomString.get(0).getSymbol().equals("H") && atomString.get(1).getSymbol().equals("H");
        discard_H_bonds(!is_H2);
        int_tab_specifier(num_atoms * num_bonds);
    }

    int getBondNumber() { return num_bonds; }

    public int indexOf() { return num_atoms; }

    List<IAtom> getAtomString() { return atomString; }

    public IAtomContainer getAtomContainer() { return ac; }

    private void set_atom_string() {
        atomString = new ArrayList<>();
        for (int i = 0; i < num_atoms; i++) atomString.add(ac.getAtom(i));
    }

    int get_num_H_atoms() {
        return num_heavy_atoms;
    }

    private void set_H_number() {
        num_heavy_atoms = num_atoms;
        for (IAtom atom : atomString)
            if (atom.getSymbol().equals("H")) num_heavy_atoms--;
    }

    private void set_connection_tables() {
        for (int i = 0; i < num_bonds; i++) {
            IBond bond = ac.getBond(i);
            int fst_ind = ac.indexOf(bond.getAtom(0)) + 1;
            int snd_ind = ac.indexOf(bond.getAtom(1)) + 1;
            int order = matchBonds ? bond.getOrder().numeric() : 1;
            String fst_sym = bond.getAtom(0).getSymbol();
            String snd_sym = bond.getAtom(1).getSymbol();
            int_table.add(new Integer[]{fst_ind, snd_ind, order});
            char_table.add(new String[]{fst_sym, snd_sym});
        }
    }

    private void discard_H_bonds(boolean is_no_H2) {
        int count_bonds = 0;
        if (is_no_H2)
            for (int i = 0; i < num_bonds; i++) {
                String[] entry = char_table.get(i);
                if (entry[0].equals("H") || entry[1].equals("H")) num_atoms--;
                if (!entry[0].equals("H") && !entry[1].equals("H")) {
                    char_tab.add(entry);
                    int_tab.add(int_table.get(i));
                    count_bonds++;
                }
            }
        else
            for (int i = 0; i < num_bonds; i++) {
                char_tab.add(char_table.get(i));
                int_tab.add(int_table.get(i));
                count_bonds++;
            }
        num_bonds = count_bonds;
    }

    private void int_tab_specifier(int specifier_value) {
        for (int i = 0; i < num_bonds; i++) {
            Integer[] entry = int_tab.get(i);
            specified_int_tab.add(new Integer[]{entry[0] + specifier_value, entry[1] + specifier_value, entry[2]});
        }
    }
}
