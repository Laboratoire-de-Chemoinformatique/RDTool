package ru.spb.arcadia.jnj.aam.io;

import org.apache.commons.cli.*;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import uk.ac.ebi.reactionblast.tools.StandardizeReaction;
import uk.ac.ebi.reactionblast.tools.rxnfile.MDLV2000RXNWriter;

import ru.spb.arcadia.jnj.aam.MappingSolution;

import java.io.*;
import java.util.ArrayList;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.HashMap;
import java.util.Map;

public class Commons {
    public static IReaction[] read_smiles(String in_file, boolean standardize) throws InvalidSmilesException {
        SmilesParser smiles_parser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        StandardizeReaction standardizer = new StandardizeReaction();
        ArrayList<IReaction> reactions = new ArrayList<>();
        try {
            BufferedReader reader = new BufferedReader(new FileReader(in_file));
            int counter = 0;
            while (true) {
                String line = reader.readLine();
                if (line == null) {
                    break;
                }
                String[] strs = line.split("\\s+");
                String smiles = strs[0];
                String id = strs.length > 1 ? strs[1] : Integer.toString(counter);
                IReaction reaction = smiles_parser.parseReactionSmiles(smiles);
                if (standardize)
                    reaction = standardizer.standardize(reaction);
                reaction.setID(id);
                reactions.add(reaction);
                counter++;
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        IReaction[] result = new IReaction[reactions.size()];
        for (int i = 0; i < reactions.size(); i++) result[i] = reactions.get(i);
        return result;
    }

    public static void write_smiles(String out_file, IReaction[] reactions) throws IOException, CDKException {
        SmilesGenerator sg = new SmilesGenerator(SmiFlavor.AtomAtomMap);
        BufferedWriter bw = new BufferedWriter(new FileWriter(out_file));
        for (IReaction reaction: reactions)
        {
            String str_reaction = sg.create(reaction);
            String id_reaction = reaction.getID();
            bw.write(str_reaction + " " + id_reaction);
            bw.newLine();
        }
        bw.close();
    }

    public static void write_rdf(String out_file, String file_header, String rxn_header, MappingSolution[] solutions) throws IOException, CDKException {
        if (solutions.length == 0) {
            System.out.println("No reactions to write! Exiting");
            System.exit(1);
        }
        BufferedWriter bw = new BufferedWriter(new FileWriter(out_file));
        file_header = file_header == null ? "$RDFILE 1" : file_header;
        bw.write(file_header + "\n");
        String current_date_time = new SimpleDateFormat("MM/dd/yy HH:mm").format(new Date());
        bw.write("$DATM " + current_date_time + "\n");
        MDLV2000RXNWriter rxn_writer = new MDLV2000RXNWriter(bw);
        for (MappingSolution solution : solutions) {
            bw.write("$RFMT\n");
            rxn_writer.write(solution.getReaction());
            Map<String, Object> meta = new HashMap<>();
            meta.put(rxn_header, solution.getReactionId());
            meta.put("Algorithm", solution.getAlgorithmName());
            meta.put("Bond_energy_sum", solution.getBondEnergySum());
            meta.put("Energy_delta", solution.getEnergyDelta());
            meta.put("Total_bond_changes", solution.getTotalBondChanges());
            meta.put("Total_fragment_changes", solution.getTotalFragmentChanges());
            meta.put("Total_stereo_changes", solution.getTotalStereoChanges());
            meta.put("Smallest_fragment_count", solution.getSmallestFragmentCount());
            meta.put("Total_carbon_bond_changes", solution.getTotalCarbonBondChanges());
            meta.put("Total_changes", solution.getTotalChanges());
            for (String key : meta.keySet()) {
                bw.write("$DTYPE " + key + "\n");
                bw.write("$DATUM " + meta.get(key) + "\n");
            }
        }
        rxn_writer.close();
        bw.close();
    }

    private static MappingSolution line2mapping_solution(String line) {
        String[] strs = line.trim().split(",");
        try {
            String reaction_id = strs[0];
            String name = strs[1];
            boolean generate_2D = Boolean.parseBoolean(strs[2]);
            boolean generate_3D = Boolean.parseBoolean(strs[3]);
            double bond_energy_sum = Double.parseDouble(strs[4]);
            double energy_delta = Double.parseDouble(strs[5]);
            int total_bond_changes = Integer.parseInt(strs[6]);
            int total_fragment_changes = Integer.parseInt(strs[7]);
            int total_stereo_changes = Integer.parseInt(strs[8]);
            int smallest_fragment_count = Integer.parseInt(strs[9]);
            int total_carbon_changes = Integer.parseInt(strs[10]);
            int total_changes = Integer.parseInt(strs[11]);
            return new MappingSolution(
                    reaction_id, name, generate_2D, generate_3D, bond_energy_sum, energy_delta, total_bond_changes,
                    total_fragment_changes, total_stereo_changes, smallest_fragment_count, total_carbon_changes,
                    total_changes
            );
        }
        catch (Exception e) {
            return null;
        }
    }

    public static void main(String[] args) throws IOException, CDKException {
        Options options = new Options();
        Option inp = new Option("i", "input", true, "inputs");
        inp.setRequired(true);
        options.addOption(inp);
        Option out = new Option("o", "output", true, "outputs");
        out.setRequired(true);
        options.addOption(out);
        CommandLineParser parser = new DefaultParser();
        HelpFormatter formatter = new HelpFormatter();
        CommandLine cmd = null;
        try {
            cmd = parser.parse(options, args);
        } catch (Exception e) {
            System.out.println(e.getMessage());
            formatter.printHelp("aam-tools: EntryPoint", options);
            System.exit(1);
        }
        String input_file = cmd.getOptionValue("input");
        String output_file = cmd.getOptionValue("output");
        RDFReader reader = new RDFReader(input_file, "Reaction_ID");
        MappingSolution[] solutions = reader.read_solutions();
        IReaction[] reactions = new IReaction[solutions.length];
        for (int i = 0; i < reactions.length; i++) reactions[i] = solutions[i].getReaction();
        write_smiles(output_file, reactions);
    }
}
