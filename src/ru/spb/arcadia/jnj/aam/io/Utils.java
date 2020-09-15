package ru.spb.arcadia.jnj.aam.io;

import org.openscience.cdk.interfaces.*;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

public class Utils {
    public static void show_mappings(String title, IReaction reaction) {
        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter("mappings.txt", true));
            String msg = "Reactions " + reaction.getID() + " mappings: ";
//            System.out.println(title);
//            System.out.println(msg);
            bw.write(msg);
            bw.newLine();
            for (IMapping m : reaction.mappings()) {
                IAtom fst = (IAtom) m.getChemObject(0), snd = (IAtom) m.getChemObject(1);
                String fst_sym = fst.getSymbol(), snd_sym = snd.getSymbol();
                String fst_id = fst.getID(), snd_id = snd.getID();
                int fst_ind = fst.getIndex(), snd_ind = snd.getIndex();
//                msg = "\tatoms: " + fst_sym + ", " + snd_sym + "; " + fst_ind + ", " + snd_ind;
                msg = "\tatoms: " + fst_sym + ", " + snd_sym + "; " + fst_id + ", " + snd_id;
                bw.write(msg);
                bw.newLine();
//                System.out.println(msg);
            }
            bw.newLine();
            bw.close();
//            System.out.println();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }

    public static void show_bonds(String title, IAtomContainer mol) {
        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter("bonds.txt", true));
            bw.write("----- " + title + " -----");
            bw.newLine();
            String msg = "Molecule " + mol.getID() + " bonds: ";
            bw.write(msg);
            bw.newLine();
            for (IBond b : mol.bonds()) {
                IBond.Order order = b.getOrder();
                IAtom fst = b.getAtom(0), snd = b.getAtom(1);
                String fst_sym = fst.getSymbol(), snd_sym = snd.getSymbol();
                String fst_id = fst.getID(), snd_id = snd.getID();
                msg = "\t" + fst_sym + " (" + fst_id + ") " + order.toString() + " " + snd_sym + " (" + snd_id + ")";
                bw.write(msg);
                bw.newLine();
            }
            bw.newLine();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }

    public static void show_atoms(IAtomContainer mol) {
        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter("atoms.txt", true));
            String msg = "Molecule " + mol.getID() + " atoms: ";
            bw.write(msg);
            bw.newLine();
            for (IAtom a : mol.atoms()) {
                msg = "\t" + a.getSymbol() + ": " + a.getID();
                bw.write(msg);
                bw.newLine();
            }
            bw.newLine();
            bw.close();
        }
        catch (IOException e) {
            e.printStackTrace();
        }
    }
}
