package ru.spb.arcadia.jnj.aam;

import org.apache.commons.cli.*;

public class EntryPoint {
    public static void main(String[] args) {
        Options options = new Options();
        Option job_opt = new Option("j", "job", true, "job to run");
        job_opt.setRequired(true);
        options.addOption(job_opt);
        Option inp = new Option("i", "input", true, "inputs");
        inp.setRequired(true);
        options.addOption(inp);
        Option out = new Option("o", "output", true, "outputs");
        out.setRequired(true);
        options.addOption(out);
        Option format = new Option("fmt", "format", true, "IO format");
        format.setRequired(false);
        options.addOption(format);
        Option rdf_id_field = new Option("rdf_id", "rdf_id", true, "RDF reaction id");
        rdf_id_field.setRequired(false);
        options.addOption(rdf_id_field);
        Option no_min = new Option("min", "no_min_algo", false, "skip MIN algorithm");
        no_min.setRequired(false);
        options.addOption(no_min);
        Option no_max = new Option("max", "no_max_algo", false, "skip MAX algorithm");
        no_max.setRequired(false);
        options.addOption(no_max);
        Option no_mixture = new Option("mixture", "no_mixture_algo", false, "skip MIXTURE algorithm");
        no_mixture.setRequired(false);
        options.addOption(no_mixture);
        Option timeout = new Option("t", "timeout", true, "timeout for single reaction");
        timeout.setRequired(false);
        options.addOption(timeout);
        Option old = new Option("old", "old", false, "Use old workflow");
        old.setRequired(false);
        options.addOption(old);
        Option tfc = new Option("ignore_tfc", "ignore_total_fragment_changes", false, "Ignore metric total_fragment_changes in consensus (= 0)");
        tfc.setRequired(false);
        options.addOption(tfc);
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
        String job = cmd.getOptionValue("job");
        if (job.equalsIgnoreCase("MAPPING")) {
            String input_file = cmd.getOptionValue("input");
            String output_dir = cmd.getOptionValue("output");
            String io_format = cmd.hasOption("format") ? cmd.getOptionValue("format") : "RDF";
            String reaction_id = cmd.hasOption("rdf_id") ? cmd.getOptionValue("rdf_id") : null;
            int timeout_val = Integer.parseInt(cmd.getOptionValue("timeout", "-1"));
            boolean skip_min = cmd.hasOption("min");
            boolean skip_max = cmd.hasOption("max");
            boolean skip_mixture = cmd.hasOption("mixture");
            boolean is_old = cmd.hasOption("old");
            if (!is_old & skip_min & skip_max & skip_mixture) {
                System.out.println("No algorithm was specified, exiting.");
                System.exit(1);
            }
            if (is_old & timeout_val != -1)
                System.out.println("Timeout option is ignored when an \"-old\" flag is set.");
            Mapping.entry_point(input_file, output_dir, io_format, reaction_id, skip_min, skip_max, skip_mixture, timeout_val, is_old);
        } else if (job.equalsIgnoreCase("CONSENSUS")) {
            String[] input_files = cmd.getOptionValue("input").split(";");
            String output_file = cmd.getOptionValue("output");
            String reaction_id = cmd.hasOption("rdf_id") ? cmd.getOptionValue("rdf_id") : null;
            boolean ignore_tfc = cmd.hasOption("ignore_tfc");
            if (ignore_tfc)
                System.out.println("Will ignore metric \"total_fragment_changes\" in consensus");
            Consensus.entry_point(input_files, output_file, reaction_id, ignore_tfc);
        } else if (job.equalsIgnoreCase("RENDERING")) {
            String input_file = cmd.getOptionValue("input");
            String output_dir = cmd.getOptionValue("output");
            String reaction_id = cmd.hasOption("rdf_id") ? cmd.getOptionValue("rdf_id") : null;
            MappingRenderer.entry_point(reaction_id, input_file, output_dir);
        } else
            System.out.println("Job not implemented: " + job);
    }
}
