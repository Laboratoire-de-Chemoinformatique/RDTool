package ru.spb.arcadia.jnj.aam;

import org.apache.commons.cli.*;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IReaction;
import uk.ac.ebi.reactionblast.tools.StandardizeReaction;

import java.io.IOException;
import java.util.ArrayList;

import ru.spb.arcadia.jnj.aam.io.Commons;

public class Standardize {
    public static void main(String[] args){
        Options options = new Options();
        Option input = new Option("i", "input_file", true, "input file path");
        input.setRequired(true);
        options.addOption(input);
        Option output = new Option("o", "output_file", true, "output file path");
        output.setRequired(true);
        options.addOption(output);

        CommandLineParser parser = new DefaultParser();
        HelpFormatter formatter = new HelpFormatter();
        CommandLine cmd = null;
        try {
            cmd = parser.parse(options, args);
        } catch (Exception e) {
            System.out.println(e.getMessage());
            formatter.printHelp("aam-tools:standardize", options);
            System.exit(1);
        }
        String input_file = cmd.getOptionValue("input_file");
        String output_file = cmd.getOptionValue("output_file");
        IReaction[] inp_reactions = null;
        try {
            inp_reactions = Commons.read_smiles(input_file, false);
        }
        catch (InvalidSmilesException e) {
            e.printStackTrace();
            System.exit(1);
        }
        ArrayList<IReaction> standardized_reactions = new ArrayList<>();
        StandardizeReaction standardizer = new StandardizeReaction();
        for (IReaction reaction : inp_reactions) {
            try {
                IReaction clean_reaction = standardizer.standardize(reaction);
                standardized_reactions.add(clean_reaction);
            }
            catch (Exception e) {
                System.out.println("Failed to standardize reaction " + reaction.getID());
                System.out.println(e);
            }
        }
        try {
            IReaction[] reactions_to_write = new IReaction[standardized_reactions.size()];
            for (int i = 0; i < standardized_reactions.size(); i++) reactions_to_write[i] = standardized_reactions.get(i);
            Commons.write_smiles(output_file, reactions_to_write);
        }
        catch (IOException | CDKException e) {
            e.printStackTrace();
            System.exit(1);
        }
    }
}
