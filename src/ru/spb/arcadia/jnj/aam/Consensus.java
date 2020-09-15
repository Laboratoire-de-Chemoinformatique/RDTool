package ru.spb.arcadia.jnj.aam;

import org.openscience.cdk.exception.CDKException;
import ru.spb.arcadia.jnj.aam.io.Commons;
import ru.spb.arcadia.jnj.aam.io.RDFReader;
import uk.ac.ebi.reactionblast.mapping.interfaces.IMappingAlgorithm;

import java.io.IOException;
import java.util.*;

public class Consensus {
    private MappingSolution best = null;
    private final boolean ignore_tfc;

    public Consensus(Map<String, MappingSolution> mapping_solutions, boolean ignore_tfc) throws Exception {
        this.ignore_tfc = ignore_tfc;
        // check standard mappings in specified order
        for (IMappingAlgorithm algorithm : IMappingAlgorithm.values()) {
            String algorithm_name = algorithm.name();
            if (mapping_solutions.containsKey(algorithm_name)) {
                MappingSolution other = mapping_solutions.get(algorithm_name);
                if (best == null || (ignore_tfc ? best.other_is_better_ignore_tfc(other) : best.other_is_better(other)))
                    best = other;
                mapping_solutions.remove(algorithm_name, other);
            }
        }
        // check custom mappings if any
        for (MappingSolution solution : mapping_solutions.values())
            if (best == null || (ignore_tfc ? best.other_is_better_ignore_tfc(solution) : best.other_is_better(solution)))
                best = solution;
        if (best == null)
            throw new Exception("Failed to find best solution");
    }

    private MappingSolution getConsensus() {
        return best;
    }

    static void entry_point(String[] input_files, String output_file, String id_field, boolean ignore_tfc) {
        Map<String, List<MappingSolution>> inputs = new HashMap<>();
        for (String input_file : input_files) {
            MappingSolution[] solutions;
            try {
                RDFReader reader = new RDFReader(input_file, id_field);
                solutions = reader.read_solutions();
                for (MappingSolution solution : solutions) {
                    String reaction_id = solution.getReactionId();
                    if (!inputs.containsKey(reaction_id)) {
                        List<MappingSolution> new_list = new ArrayList<>();
                        inputs.put(reaction_id, new_list);
                    }
                    inputs.get(reaction_id).add(solution);
                }
            } catch (Exception e) {
                System.out.println("Failed to read file " + input_file);
            }
        }
        ArrayList<MappingSolution> consensuses = new ArrayList<>();
        for (Map.Entry<String, List<MappingSolution>> entry : inputs.entrySet()) {
            Map<String, MappingSolution> mapping_solutions = new HashMap<>();
            for (MappingSolution solution : entry.getValue())
                mapping_solutions.put(solution.getAlgorithmName(), solution);
            try {
                Consensus consensus = new Consensus(mapping_solutions, ignore_tfc);
                consensuses.add(consensus.getConsensus());
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        }
        try {
            MappingSolution[] consensuses_to_write = new MappingSolution[consensuses.size()];
            for (int i = 0; i < consensuses.size(); i++) consensuses_to_write[i] = consensuses.get(i);
            Commons.write_rdf(output_file, "$RDFILE 1", "Reaction_ID", consensuses_to_write);
        }
        catch (IOException | CDKException e) {
            e.printStackTrace();
        }
    }
}
