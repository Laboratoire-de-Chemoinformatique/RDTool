package ru.spb.arcadia.jnj.aam;

import java.io.IOException;
import java.io.Serializable;
import java.util.HashMap;
import java.util.Map;

import static java.lang.Math.abs;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IReaction;
import static org.openscience.cdk.interfaces.IBond.Order.DOUBLE;
import static org.openscience.cdk.interfaces.IBond.Order.QUADRUPLE;
import static org.openscience.cdk.interfaces.IBond.Order.SINGLE;
import static org.openscience.cdk.interfaces.IBond.Order.TRIPLE;
import org.openscience.smsd.tools.BondEnergies;
import static org.openscience.smsd.tools.BondEnergies.getInstance;

import org.xmlcml.euclid.Int;
import uk.ac.ebi.reactionblast.fingerprints.interfaces.IFeature;
import uk.ac.ebi.reactionblast.fingerprints.interfaces.IPatternFingerprinter;
import uk.ac.ebi.reactionblast.mechanism.BondChangeCalculator;

public class MappingSolution implements Serializable {

    private static final String NEW_LINE = System.getProperty("line.separator");
    private static final long serialVersionUID = 1678787866L;

    private final String reaction_id;
    private final String algorithm_name;
    private final double bond_energy_sum;
    private final double energy_delta;
    private final int total_bond_changes;
    private final int total_fragment_changes;
    private final int total_stereo_changes;
    private final int smallest_fragment_count;
    private final int total_carbon_bond_changes;
    private final int total_changes;
    private final BondChangeCalculator bcc;
    private IReaction mapped_reaction;
    private boolean generate_2D;
    private boolean generate_3D;

    public MappingSolution(String algorithm_name, IReaction mapped_reaction, int total_fragment_changes, boolean generate_2D, boolean generate_3D) throws Exception {
        this.reaction_id = mapped_reaction.getID();
        this.algorithm_name = algorithm_name;
        this.mapped_reaction = mapped_reaction;
        this.generate_2D = generate_2D;
        this.generate_3D = generate_3D;
        this.bcc = new BondChangeCalculator(this.mapped_reaction, generate_2D, generate_3D);
//        this.bcc = new BondChangeCalculator(this.mapped_reaction);
//        this.bcc.computeBondChanges(generate_2D, generate_3D);
        this.bond_energy_sum = get_total_bond_change_energy(this.bcc.getFormedCleavedWFingerprint(), true);
        this.energy_delta = bcc.getEnergyDelta();
        this.total_bond_changes = (int) (get_total_bond_change(bcc.getFormedCleavedWFingerprint()) + get_total_bond_change(bcc.getOrderChangesWFingerprint()));
        this.total_fragment_changes = total_fragment_changes;
        this.total_stereo_changes = (int) get_total_bond_change(bcc.getStereoChangesWFingerprint());
        this.smallest_fragment_count = bcc.getTotalSmallestFragmentSize();
        this.total_carbon_bond_changes = get_total_carbon_bond_change(bcc.getFormedCleavedWFingerprint());
        this.total_changes = this.total_bond_changes + this.total_fragment_changes;
    }

    public MappingSolution(uk.ac.ebi.reactionblast.mechanism.MappingSolution solution) throws Exception {
        this.bcc = solution.getBondChangeCalculator();
        this.algorithm_name = solution.getAlgorithmID().toString();
        this.mapped_reaction = bcc.getReaction();
        this.reaction_id = this.mapped_reaction.getID();
        this.generate_2D = true;
        this.generate_3D = true;
        this.bond_energy_sum = get_total_bond_change_energy(this.bcc.getFormedCleavedWFingerprint(), true);
        this.energy_delta = this.bcc.getEnergyDelta();
        this.total_bond_changes = (int) (get_total_bond_change(this.bcc.getFormedCleavedWFingerprint()) + get_total_bond_change(this.bcc.getOrderChangesWFingerprint()));
        this.total_fragment_changes = solution.getReactor().getDelta();
        this.total_stereo_changes = (int) get_total_bond_change(this.bcc.getStereoChangesWFingerprint());
        this.smallest_fragment_count = this.bcc.getTotalSmallestFragmentSize();
        this.total_carbon_bond_changes = get_total_carbon_bond_change(this.bcc.getFormedCleavedWFingerprint());
        this.total_changes = this.total_bond_changes + this.total_fragment_changes;
    }

    public MappingSolution(String reaction_id, String algorithm_name, boolean generate_2D, boolean generate_3D, double bond_energy_sum, double energy_delta,
                           int total_bond_changes, int total_fragment_changes, int total_stereo_changes, int smallest_fragment_count,
                           int total_carbon_bond_changes, int total_changes) {
        this.reaction_id = reaction_id;
        this.algorithm_name = algorithm_name;
        this.mapped_reaction = null;
        this.generate_2D = generate_2D;
        this.generate_3D = generate_3D;
        this.bond_energy_sum = bond_energy_sum;
        this.bcc = null;
        this.energy_delta = energy_delta;
        this.total_bond_changes = total_bond_changes;
        this.total_fragment_changes = total_fragment_changes;
        this.total_stereo_changes = total_stereo_changes;
        this.smallest_fragment_count = smallest_fragment_count;
        this.total_carbon_bond_changes = total_carbon_bond_changes;
        this.total_changes = total_changes;
    }

    private static MappingSolution from_meta(String id_field, Map<String, String> meta, IReaction reaction) {
        MappingSolution solution;
        try {
            IReaction fixed_reaction = new MappingFixer(reaction).get_result();
            solution = new MappingSolution(
                    meta.get(id_field), meta.get("Algorithm"),
                    Boolean.parseBoolean(meta.get("Generate_2D")), Boolean.parseBoolean(meta.get("Generate_3D")),
                    Double.parseDouble(meta.get("Bond_energy_sum")), Double.parseDouble(meta.get("Energy_delta")),
                    Integer.parseInt(meta.get("Total_bond_changes")), Integer.parseInt(meta.getOrDefault("Total_fragment_changes", "0")),
                    Integer.parseInt(meta.get("Total_stereo_changes")), Integer.parseInt(meta.get("Smallest_fragment_count")),
                    Integer.parseInt(meta.get("Total_carbon_bond_changes")), Integer.parseInt(meta.get("Total_changes"))
            );
            solution.setReaction(fixed_reaction);
        }
        catch (Exception e) {
            e.printStackTrace();
            solution = null;
        }
        return solution;
    }

    private static MappingSolution from_mapped_reaction(IReaction mapped_reaction) {
        MappingSolution solution;
        try {
            IReaction fixed_reaction = new MappingFixer(mapped_reaction).get_result();
            solution = new MappingSolution("custom", fixed_reaction, 0, true, true);
        }
        catch (Exception e) {
            e.printStackTrace();
            solution = null;
        }
        return solution;
    }

    public static MappingSolution from_rdf_entry(String id_field, Map<String, String> meta, IReaction reaction) throws Exception {
        String[] fields = new String[]{
                id_field,
                "Algorithm",
                "Bond_energy_sum",
                "Energy_delta",
                "Total_bond_changes",
                "Total_fragment_changes",
                "Total_stereo_changes",
                "Smallest_fragment_count",
                "Total_carbon_bond_changes",
                "Total_changes"
        };
        boolean is_full_meta = true;
        for (String field : fields)
            if (!meta.containsKey(field)) {
                is_full_meta = false;
                break;
            }
        return is_full_meta ? from_meta(id_field, meta, reaction) : from_mapped_reaction(reaction);
    }

    private void setReaction(IReaction reaction) {
        // TODO: check mapping
        this.mapped_reaction = reaction;
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append(NEW_LINE);
        sb.append("-----------------------------------");
        sb.append(NEW_LINE);
        sb.append("Algorithm name is = ").append(getAlgorithmName());
        sb.append(", ");
        sb.append(NEW_LINE).append("Scores = " + "(Chaos Delta:");
        sb.append(this.getTotalFragmentChanges()).append(", Sigma: ");
        sb.append(this.getTotalBondChanges()).append(", Energy: ");
        sb.append(this.getBondEnergySum()).append(")");
        sb.append(NEW_LINE);
        sb.append("-----------------------------------");
        sb.append(NEW_LINE);
        sb.append("MappingSolution{" + "algorithm name = ").append(algorithm_name).append(", bondEnergyChange = ").append(bond_energy_sum);
        sb.append(", totalBondChanges = ").append(total_bond_changes).append(", totalFragmentChanges = ").append(total_fragment_changes);
        sb.append(", totalStereoChanges = ").append(total_stereo_changes).append(", smallestFragmentCount = ").append(smallest_fragment_count);
        sb.append(", totalChanges = ").append(total_changes).append(", generate2D = ").append(generate_2D).append(", generate3D = ").append(generate_3D).append("}");
        sb.append(NEW_LINE);
        return sb.toString();
    }

    boolean other_is_better_ignore_tfc(MappingSolution other) {
        // ignores value of "total_fragment_changes" field
        if (other.total_bond_changes == 0 && other.total_stereo_changes == 0)
            return false;
        if (other.bond_energy_sum == 0
                && other.total_bond_changes == 0
                && total_stereo_changes >= other.total_stereo_changes)
            return true; // Condition 1
        if (total_bond_changes > other.total_bond_changes
                && total_carbon_bond_changes > 0
                && total_carbon_bond_changes > other.total_carbon_bond_changes)
            return true; // Condition 2
        if (total_bond_changes > other.total_bond_changes)
            return true; // Condition 3
        if (smallest_fragment_count >= other.smallest_fragment_count
                && bond_energy_sum > other.bond_energy_sum
                && total_carbon_bond_changes >= other.total_carbon_bond_changes) {
            return true; // Condition 4
        }
        if (smallest_fragment_count > other.smallest_fragment_count)
            return true; // Condition 5
        if (smallest_fragment_count == other.smallest_fragment_count
                && bond_energy_sum > other.bond_energy_sum
                && total_carbon_bond_changes >= other.total_carbon_bond_changes)
            return true; // Condition 6
        if (bond_energy_sum > other.bond_energy_sum)
            return true; // Condition 7
        if (total_bond_changes == other.total_bond_changes)
            return true; // Condition 8
        // Conditions 9, 10 are always false at this point
        if (bond_energy_sum > other.bond_energy_sum && total_carbon_bond_changes > other.total_carbon_bond_changes)
            return true; // Condition 11
        // Conditions 12 - 15 are always false at this point
        return false;
    }

    boolean other_is_better(MappingSolution other) {
        // this is reimplementation of "isChangeFeasible" in uk.ac.ebi.reactionblast.ReactionMechanismTool
//        System.out.println("compare " + this.getAlgorithmName() + " to " + other.getAlgorithmName());

        if (other.getTotalBondChanges() == 0 && other.getTotalStereoChanges() == 0) {
            return false;
        }
        if (other.bond_energy_sum == 0
                && other.total_fragment_changes == 0
                && other.total_bond_changes == 0
                && total_stereo_changes >= other.total_stereo_changes) {
            return true; // Condition 1
        }
        if (total_bond_changes > other.total_bond_changes
                && total_carbon_bond_changes > 0
                && total_carbon_bond_changes > other.total_carbon_bond_changes
                && total_fragment_changes > other.total_fragment_changes) {
            return true; // Condition 2
        }
        if (total_bond_changes > other.total_bond_changes && total_fragment_changes > 0 && other.total_fragment_changes > 0) {
            return true; // Condition 3
        }
        if (total_fragment_changes >= other.total_fragment_changes
                && smallest_fragment_count >= other.smallest_fragment_count
                && bond_energy_sum > other.bond_energy_sum
                && total_carbon_bond_changes >= other.total_carbon_bond_changes) {
            return true; // Condition 4
        }
        if (total_fragment_changes > other.total_fragment_changes && smallest_fragment_count > other.smallest_fragment_count) {
            return true; // Condition 5
        }
        if (total_fragment_changes == other.total_fragment_changes
                && smallest_fragment_count == other.smallest_fragment_count
                && bond_energy_sum > other.bond_energy_sum
                && total_carbon_bond_changes >= other.total_carbon_bond_changes) {
            return true; // Condition 6
        }
        if (total_fragment_changes > other.total_fragment_changes && bond_energy_sum > other.bond_energy_sum) {
            return true; // Condition 7
        }
        if (total_bond_changes == other.total_bond_changes && total_fragment_changes > other.total_fragment_changes) {
            return true; // Condition 8
        }
        if (total_fragment_changes == other.total_fragment_changes
                && bond_energy_sum == other.bond_energy_sum
                && total_bond_changes > other.total_bond_changes) {
            return true; // Condition 9
        }
        if (bond_energy_sum == other.bond_energy_sum
                && total_bond_changes == other.total_bond_changes
                && total_stereo_changes > other.total_stereo_changes) {
            return true; // Condition 10
        }
        if (bond_energy_sum > other.bond_energy_sum && total_carbon_bond_changes > other.total_carbon_bond_changes) {
            return true; // Condition 11
        }
        if (total_bond_changes < other.total_bond_changes
                && bond_energy_sum < other.bond_energy_sum
                && total_carbon_bond_changes > 0
                && total_carbon_bond_changes > other.total_carbon_bond_changes
                && smallest_fragment_count > other.smallest_fragment_count) {
            return true; // Condition 12
        }
        if (total_bond_changes > other.total_bond_changes
                && total_carbon_bond_changes > other.total_carbon_bond_changes
                && smallest_fragment_count > other.smallest_fragment_count) {
            return true; // Condition 13
        }
        if (total_bond_changes == other.total_bond_changes
                && total_carbon_bond_changes == other.total_carbon_bond_changes
                && smallest_fragment_count > other.smallest_fragment_count) {
            return true; // Condition 14
        }
        if (total_bond_changes == other.total_bond_changes
                && total_carbon_bond_changes == other.total_carbon_bond_changes
                && bond_energy_sum > other.bond_energy_sum) {
            return true; // Condition 15
        }
        return false;
    }

    public String getAlgorithmName() { return algorithm_name; }

    public double getBondEnergySum() { return bond_energy_sum; }

    public int getTotalBondChanges() { return total_bond_changes; }

    public int getTotalFragmentChanges() { return total_fragment_changes; }

    public int getTotalStereoChanges() { return total_stereo_changes; }

    public int getSmallestFragmentCount() { return smallest_fragment_count; }

    public IReaction getReaction() { return mapped_reaction; }

    public String getReactionId() { return reaction_id; }

    public int getTotalChanges() { return total_changes; }

    BondChangeCalculator getBondChangeCalculator() { return bcc; }

    public boolean isGenerate3D() { return generate_3D; }

    public boolean isGenerate2D() { return generate_2D; }

    public int getTotalCarbonBondChanges() { return total_carbon_bond_changes; }

    public double getEnergyDelta() { return energy_delta; }

    private double get_total_bond_change(IPatternFingerprinter fingerprint) {
        double total = 0;
        total = fingerprint.getFeatures().stream().map((key) -> key.getWeight()).filter((val) -> (val > 0.)).map((val) -> val).reduce(total, (accumulator, _item) -> accumulator + _item);
        return total;
    }

    private int get_total_carbon_bond_change(IPatternFingerprinter fingerprint) {
        double total = 0;
        total = fingerprint.getFeatures().stream().filter((key) -> (key.getPattern().contains("C-C")
                || key.getPattern().contains("C=C")
                || key.getPattern().contains("C#C")
                || key.getPattern().contains("C%C")
                || key.getPattern().contains("C@C"))).map((key) -> key.getWeight()).filter((val) -> (val > 0.)).map((val) -> val).reduce(total, (accumulator, _item) -> accumulator + _item);
        return (int) total;
    }

    private int get_total_bond_change_energy_2(IPatternFingerprinter fingerprinter, boolean skipHydrogen) {
        int total = 0;
        Map<String, IBond.Order> bond_orders = new HashMap<>();
        bond_orders.put("-", SINGLE);
        bond_orders.put("%", SINGLE);
        bond_orders.put("@", SINGLE);
        bond_orders.put("=", DOUBLE);
        bond_orders.put("#", TRIPLE);
        bond_orders.put("$", QUADRUPLE);
        try {
            BondEnergies be =getInstance();
            for (IFeature feature : fingerprinter.getFeatures()) {
                String feature_key = feature.getPattern();
                double feature_val = feature.getWeight();
                for (String key : bond_orders.keySet()) {
                    if (feature_key.contains(key)) {
                        String[] syms = feature_key.split(key);
                        IBond.Order order = bond_orders.get(key);
                        int high_ring = (key.equals("%")) ? 1 : 0;
                        if (!skipHydrogen && !syms[0].equals("H") && !syms[1].equals("H")) {
                            int energy = be.getEnergies(syms[0], syms[1], order);
                            if (energy > 0) {
                                total += feature_val * (energy - 5.0 * high_ring);
                            }
                        }
                        break;
                    }
                }
            }
        }
        catch (CDKException e) {

        }
        return abs(total);
    }

    private int get_total_bond_change_energy(IPatternFingerprinter fingerprint, boolean skipHydrogen) {
        int total = 0;
        try {
            BondEnergies be = getInstance();
            for (IFeature feature : fingerprint.getFeatures()) {
                double val = feature.getWeight();
                String key = feature.getPattern();
                if (val > 0) {
                    if (key.contains("-") || key.contains("%") || key.contains("@")) {
                        String[] temp = null;
                        if (key.contains("-")) {
                            temp = key.split("-");
                        } else if (key.contains("%")) {
                            temp = key.split("%");
                        } else if (key.contains("@")) {
                            temp = key.split("@");
                        }
                        if (skipHydrogen && (temp[0].equals("H") || temp[1].equals("H"))) {
                            continue;
                        }
                        int energy = be.getEnergies(temp[0], temp[1], SINGLE);
                        if (energy > 0) {
                            if (key.contains("%")) {
                                total += val * (energy - 5.0);

                            } else {
                                total += val * energy;
                            }
                        }
                    } else if (key.contains("=")) {
                        String[] temp = key.split("=");
                        if (skipHydrogen && (temp[0].equals("H") || temp[1].equals("H"))) {
                            continue;
                        }
                        int energy = be.getEnergies(temp[0], temp[1], DOUBLE);
                        if (energy > 0) {
                            total += val * energy;
                        }
                    } else if (key.contains("#")) {
                        String[] temp = key.split("#");
                        int energy = be.getEnergies(temp[0], temp[1], TRIPLE);
                        if (skipHydrogen && (temp[0].equals("H") || temp[1].equals("H"))) {
                            continue;
                        }
                        if (energy > 0) {
                            total += val * energy;
                        }
                    } else if (key.contains("$")) {
                        String[] temp = key.split("$");
                        if (skipHydrogen && (temp[0].equals("H") || temp[1].equals("H"))) {
                            continue;
                        }
                        int energy = be.getEnergies(temp[0], temp[1], QUADRUPLE);
                        if (energy > 0) {
                            total += val * energy;
                        }
                    }
                }
            }
        } catch (CDKException ex) {
        }
        return abs(total);
    }
}
