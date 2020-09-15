package org.openscience.smsd.algorithm.mcsplus;

import org.openscience.cdk.interfaces.IAtom;

import java.util.Comparator;

public class AtomComparator implements Comparator<IAtom> {
    @Override
    public int compare(IAtom fst, IAtom snd) {
        if (fst.getID() == null && snd.getID() == null) return 0;
        if (fst.getID() == null) return -1;
        if (snd.getID() == null) return 1;
        return Integer.parseInt(fst.getID()) - Integer.parseInt(snd.getID());
    }
}
