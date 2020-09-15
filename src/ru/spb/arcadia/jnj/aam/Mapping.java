package ru.spb.arcadia.jnj.aam;

import java.io.File;
import java.util.*;
import java.util.concurrent.*;

import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IMapping;
import org.openscience.cdk.interfaces.IReaction;

import org.openscience.smsd.tools.ExtAtomContainerManipulator;
import ru.spb.arcadia.jnj.aam.io.Commons;
import ru.spb.arcadia.jnj.aam.io.RDFReader;
import ru.spb.arcadia.jnj.aam.io.Utils;
import ru.spb.arcadia.jnj.aam.utils.ThreadSafeCache;
import uk.ac.ebi.reactionblast.mapping.MappingThread;
import uk.ac.ebi.reactionblast.mapping.Reactor;
import uk.ac.ebi.reactionblast.mapping.interfaces.IMappingAlgorithm;
import uk.ac.ebi.reactionblast.mechanism.BondChangeCalculator;
import uk.ac.ebi.reactionblast.mechanism.ReactionMechanismTool;
import uk.ac.ebi.reactionblast.tools.StandardizeReaction;

public class Mapping {
    private int timeout;

    public Mapping(
            IReaction[] reactions, String output_dir, int timeout,
            boolean skip_min, boolean skip_max, boolean skip_mixture) throws Exception {
        this.timeout = timeout;
        if (this.timeout <= 0)
            System.out.println("No timeout is set for mapping");
        System.out.println("Will map " + reactions.length + " reactions");
        ThreadSafeCache<String, Integer> cache = ThreadSafeCache.getInstance();
        cache.put("timeout", this.timeout);
        Map<IMappingAlgorithm, String> algorithms = new HashMap<>();
        File dir = new File(output_dir);
        if (!dir.exists())
            dir.mkdirs();
        if (!skip_min)
            algorithms.put(IMappingAlgorithm.MIN, "MIN");
        if (!skip_max)
            algorithms.put(IMappingAlgorithm.MAX, "MAX");
        if (!skip_mixture)
            algorithms.put(IMappingAlgorithm.MIXTURE, "MIXTURE");
        for (IMappingAlgorithm algorithm : algorithms.keySet())
        {
            System.out.println("New mapping: start algorithm " + algorithm.toString());
            long start_alg = System.currentTimeMillis();
            Reactor[] solutions = run_method(reactions, algorithm);
            long end_alg = System.currentTimeMillis();
            write_method_results(solutions, algorithm, algorithms.get(algorithm), output_dir);
            System.out.println(algorithm + " algorithm time (s): " + 1.0 * (end_alg - start_alg) / 1000 + "\n");
        }
        cache.cleanup();
    }

    public Mapping(IReaction[] reactions, String output_dir) throws Exception {
        File dir = new File(output_dir);
        if (!dir.exists())
            dir.mkdirs();
        boolean forced_mapping = true;
        boolean generate_2D = true;
        boolean generate_3D = true;
        StandardizeReaction standardizer = new StandardizeReaction();
        List<MappingSolution> solutions = new ArrayList<>();
        Map<IMappingAlgorithm, ArrayList<MappingSolution>> all_solutions = new HashMap<>();
        all_solutions.put(IMappingAlgorithm.MAX, new ArrayList<>());
        all_solutions.put(IMappingAlgorithm.MIN, new ArrayList<>());
        all_solutions.put(IMappingAlgorithm.MIXTURE, new ArrayList<>());
//        all_solutions.put(IMappingAlgorithm.RINGS, new ArrayList<>());
        for (IReaction reaction : reactions) {
            try {
                long start = System.currentTimeMillis();
                ReactionMechanismTool rmt = new ReactionMechanismTool(reaction, forced_mapping, generate_2D, generate_3D, standardizer);
                uk.ac.ebi.reactionblast.mechanism.MappingSolution rdt_solution = rmt.getSelectedSolution();
//                ru.spb.arcadia.jnj.aam.io.Utils.show_mappings(rdt_solution.getReaction());
//                for (IAtomContainer mol : rdt_solution.getReaction().getReactants().atomContainers())
//                    ru.spb.arcadia.jnj.aam.io.Utils.show_atoms(mol);
//                for (IAtomContainer mol : rdt_solution.getReaction().getProducts().atomContainers())
//                    ru.spb.arcadia.jnj.aam.io.Utils.show_atoms(mol);
                System.out.println("Reaction ID: " + reaction.getID() + ", algorithm: " + rdt_solution.getAlgorithmID());
                solutions.add(new MappingSolution(rdt_solution));
                long end = System.currentTimeMillis();
                System.out.println("\tElapsed time (s): " + 1.0 * (end - start) / 1000);
                for (uk.ac.ebi.reactionblast.mechanism.MappingSolution sub_rdt_solution : rmt.getAllSolutions()) {
                    MappingSolution solution = new MappingSolution(sub_rdt_solution);
                    all_solutions.get(sub_rdt_solution.getAlgorithmID()).add(solution);
                }
            }
            catch (Exception e) {
                System.out.println("Failed to process reaction: " + reaction.getID());
            }
        }
        String rxn_file = output_dir + "//consensus.rdf";
        MappingSolution[] consensuses = new MappingSolution[solutions.size()];
        for (int i = 0; i < solutions.size(); i++) consensuses[i] = solutions.get(i);
        System.out.println(consensuses.length + " consensuses will be written");
        if (consensuses.length > 0) Commons.write_rdf(rxn_file, "$RDFILE 1", "Reaction_ID", consensuses);
        for (IMappingAlgorithm algorithm : all_solutions.keySet()) {
            String alg_file = output_dir + "//" + algorithm.toString() + "_solutions.rdf";
            MappingSolution[] alg_solutions = new MappingSolution[all_solutions.get(algorithm).size()];
            for (int i = 0; i < alg_solutions.length; i++) alg_solutions[i] = all_solutions.get(algorithm).get(i);
            System.out.println("\t" + alg_solutions.length + " will be written for algorithm " + algorithm);
            if (alg_solutions.length > 0) Commons.write_rdf(alg_file, "$RDFILE 1", "Reaction_ID", alg_solutions);
        }
    }

    private Reactor[] run_method(IReaction[] reactions, IMappingAlgorithm algorithm) throws Exception {
        ArrayList<Reactor> solutions = new ArrayList<>();
        for (IReaction reaction : reactions) {
            Reactor solution = run_single_mapping(reaction, algorithm);
            if (solution != null) solutions.add(solution);
            else
                System.out.println("Failed to acquire solution for reaction " + reaction.getID());
        }
        Reactor[] result = new Reactor[solutions.size()];
        for (int i = 0; i < solutions.size(); i++) result[i] = solutions.get(i);
        return result;
    }

    private Reactor run_single_mapping(IReaction reaction, IMappingAlgorithm algorithm) {
        String rid = reaction.getID();
        Reactor reactor = null;
//        if (timeout <= 0) {
//            try {
//                long start = System.currentTimeMillis();
//                reactor = new Reactor(reaction, true, algorithm);
//                long end = System.currentTimeMillis();
//                System.out.println("Reaction " + rid + " time: " + 1.0 * (end -start) / 1000);
//            } catch (Exception e) {
//                e.printStackTrace();
//            }
//        }
//        else {
//            ExecutorService es = Executors.newCachedThreadPool();
//            Callable<Reactor> task = () -> new Reactor(reaction, true, algorithm);
//            Future<Reactor> future = es.submit(task);
//            try {
//                reactor = future.get(timeout, TimeUnit.SECONDS);
//            } catch (TimeoutException e) {
//                System.out.println("TimeoutException: " + rid);
//                future.cancel(true);
//                reactor = null;
//            } catch (InterruptedException e) {
//                System.out.println("InterruptedException: " + rid);
//            } catch (ExecutionException e) {
//                System.out.println("ExecutionException: " + rid);
//            } finally {
//                es.shutdown();
//            }
//        }
        try {
            long start = System.currentTimeMillis();
            reactor = new Reactor(reaction, true, algorithm);
            long end = System.currentTimeMillis();
            System.out.println("Reaction " + rid + " time: " + 1.0 * (end - start) / 1000);
        }
        catch (Exception e) { e.printStackTrace(); }
        return reactor;
    }

    private void write_method_results(Reactor[] reactors, IMappingAlgorithm algorithm, String algo_alias, String output_dir) throws Exception {
        String alias = output_dir + "//" + algo_alias;
        ArrayList<MappingSolution> solutions = new ArrayList<>();
        for (Reactor reactor : reactors) {
            IReaction reaction = reactor.getReactionWithAtomAtomMapping();
//            ru.spb.arcadia.jnj.aam.io.Utils.show_mappings(reaction);
//            for (IAtomContainer mol : reaction.getReactants().atomContainers())
//                ru.spb.arcadia.jnj.aam.io.Utils.show_atoms(mol);
//            for (IAtomContainer mol : reaction.getProducts().atomContainers())
//                ru.spb.arcadia.jnj.aam.io.Utils.show_atoms(mol);
            int total_fragment_changes = reactor.getDelta();
            try {
                MappingSolution solution = new MappingSolution(algorithm.toString(), reaction, total_fragment_changes, true, true);
                solutions.add(solution);
            }
            catch (Exception e) {
                System.out.println("Failed to acquire solution for reaction " + reaction.getID());
            }
        }
        String rxn_file = alias + "_reactions.rdf";
        MappingSolution[] solutions_to_write = new MappingSolution[solutions.size()];
        for (int i = 0; i < solutions.size(); i++) solutions_to_write[i] = solutions.get(i);
        Commons.write_rdf(rxn_file, "$RDFILE 1", "Reaction_ID", solutions_to_write);
    }

    private static IReaction[] prepare_reactions(IReaction[] reactions) {
        StandardizeReaction standardizer = new StandardizeReaction();
        ArrayList<IReaction> clean_reactions = new ArrayList<>();
        for (IReaction reaction : reactions) {
            IReaction res = reaction;
            for (IAtomContainer a : res.getReactants().atomContainers())
                ExtAtomContainerManipulator.setNullHCountToZero(a);
            for (IAtomContainer a : res.getProducts().atomContainers())
                ExtAtomContainerManipulator.setNullHCountToZero(a);
            try {
                res = standardizer.standardize(res);
            }
            catch (Exception e) {
                System.out.println("Failed to standardize: " + reaction.getID());
            }
            String rid = res.getID();
            for (IAtomContainer a : res.getReactants().atomContainers())
                a.setID(rid + ":" + a.getID());
            for (IAtomContainer a : res.getProducts().atomContainers())
                a.setID(rid + ":" + a.getID());
            clean_reactions.add(res);
        }
        IReaction[] result = new IReaction[clean_reactions.size()];
        for (int i = 0; i < clean_reactions.size(); i++) result[i] = clean_reactions.get(i);
        return result;
    }

    static <T> T[] list_to_array(List<T> the_list, T[] the_array) {
        for (int i = 0; i < the_list.size(); i++) the_array[i] = the_list.get(i);
        return the_array;
    }

    static void entry_point(String input_file, String output_dir, String format,
                            String reaction_id, boolean skip_min, boolean skip_max,
                            boolean skip_mixture, int timeout, boolean is_old) {
        try {
            IReaction[] reactions;
            if (format.equalsIgnoreCase("RDF")) {
                RDFReader reader = new RDFReader(input_file, reaction_id);
                reactions = reader.read_reactions();
            }
            else if (format.equalsIgnoreCase("SMI")) {
                reactions = Commons.read_smiles(input_file, false);
            }
            else
                throw new Exception("Unknown format: " + format + ", exiting.");
            reactions = prepare_reactions(reactions);
            if (is_old)
                new Mapping(reactions, output_dir);
            else
                new Mapping(reactions, output_dir, timeout, skip_min, skip_max, skip_mixture);
        } catch (Exception e) {
            System.out.println(e);
            System.exit(1);
        }
    }
}
