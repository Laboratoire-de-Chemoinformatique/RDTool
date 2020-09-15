package org.openscience.smsd.algorithm.mcsplus;

import java.util.Comparator;

public class EdgeComparator implements Comparator<Edge> {
    @Override
    public int compare(Edge fst, Edge snd) {
        return fst.getSource() == snd.getSource() ? fst.getSink() - snd.getSink() : fst.getSource() - snd.getSource();
    }
}
