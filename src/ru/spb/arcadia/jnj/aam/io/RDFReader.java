package ru.spb.arcadia.jnj.aam.io;

import org.openscience.cdk.CDK;
import org.openscience.cdk.Reaction;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.*;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.smsd.tools.ExtAtomContainerManipulator;
import ru.spb.arcadia.jnj.aam.MappingSolution;
import uk.ac.ebi.reactionblast.tools.rxnfile.MDLRXNV2000Reader;

import java.io.*;
import java.util.*;
import java.util.stream.Collectors;

public class RDFReader {
    class RDFEntry {
        String id;
        Map<String, String> meta;
        IReaction reaction;

        RDFEntry(String[] reaction_lines) throws IOException {
            this.meta = get_reaction_meta(reaction_lines);
            if (!meta.containsKey(id_field)) throw new IOException("No " + id_field + " field in meta");
            this.id = meta.get(id_field);
            String data = get_reaction_data(reaction_lines);
            this.reaction = new Reaction();
            MDLRXNV2000Reader reader = new MDLRXNV2000Reader(new StringReader(data));
            IReaction reaction = new Reaction();
            try {
                reaction = reader.read(reaction);
                for (IAtomContainer a : reaction.getReactants().atomContainers())
                    ExtAtomContainerManipulator.setNullHCountToZero(a);
                for (IAtomContainer a : reaction.getProducts().atomContainers())
                    ExtAtomContainerManipulator.setNullHCountToZero(a);
                if (convert) reaction = convertRoundTripRXNSMILES(reaction);
                reaction.setID(this.id);
                this.reaction = reaction;
            }
            catch (CDKException e) {
                System.out.println("Failed to construct reaction " + this.id);
                e.printStackTrace();
                this.reaction = null;
            }
        }

        private String get_reaction_data(String[] lines) {
            String[] data_lines = Arrays.stream(lines).filter(s -> !s.contains("$DTYPE") && !s.contains("$DATUM")).toArray(String[]::new);
            return String.join("\n", data_lines);
        }

        private Map<String, String> get_reaction_meta(String[] lines) throws IOException {
            int num_lines = lines.length;
            Map<String, String> meta = new HashMap<>();
            boolean skip_line = false;
            for (int i = 0; i < num_lines; i++) {
                if (skip_line) {
                    skip_line = false;
                    continue;
                }
                String line = lines[i].trim();
                if (line.startsWith("$DTYPE")) {
                    if (i == num_lines - 1) {
                        throw new IOException("No data for line <<" + line + ">>");
                    }
                    String datum_line = lines[i + 1].trim();
                    if (!datum_line.startsWith("$DATUM")) {
                        throw new IOException("Not a datum line <<" + datum_line + ">>");
                    }
                    String key = line.split("\\s+")[1].trim();
                    String val = datum_line.split("\\s+")[1].trim();
                    meta.put(key, val);
                    skip_line = true;
                }
            }
            return meta;
        }

        private IReaction convertRoundTripRXNSMILES(IReaction r) throws CDKException {
            final SmilesGenerator sg = new SmilesGenerator(
                    SmiFlavor.AtomAtomMap
                            | SmiFlavor.UseAromaticSymbols
                            | SmiFlavor.Stereo);
            String createSmilesFromReaction = sg.create(r);
            final SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
            IReaction parseReactionSmiles = smilesParser.parseReactionSmiles(createSmilesFromReaction);
            parseReactionSmiles.setID(r.getID());
            for (int i = 0; i < r.getReactantCount(); i++) {
                parseReactionSmiles.getReactants().getAtomContainer(i).setID(r.getReactants().getAtomContainer(i).getID());
            }
            for (int i = 0; i < r.getProductCount(); i++) {
                parseReactionSmiles.getProducts().getAtomContainer(i).setID(r.getProducts().getAtomContainer(i).getID());
            }
            return parseReactionSmiles;
        }
    }

    private BufferedReader br;
    private String id_field;
    private boolean is_read;
    private boolean convert;

    public RDFReader(String filename, String id_field) throws FileNotFoundException {
        this.br = new BufferedReader(new FileReader(filename));
        this.id_field = id_field;
        this.is_read = false;
    }

    public IReaction[] read_reactions() {
        this.convert = true;
        IReaction[] reactions = null;
        try {
            reactions = Arrays.stream(read_entries()).map(e -> e.reaction).toArray(IReaction[]::new);
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        return reactions;
    }

    public MappingSolution[] read_solutions() {
        this.convert = false;
        List<MappingSolution> solutions = new ArrayList<>();
        try {
            for (RDFEntry entry : read_entries()) {
                MappingSolution solution = MappingSolution.from_rdf_entry(id_field, entry.meta, entry.reaction);
                if (solution != null)
                    solutions.add(solution);
                else
                    System.out.println("Failed to read solution, reaction id = " + entry.id);
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        MappingSolution[] result = new MappingSolution[solutions.size()];
        for (int i = 0; i < result.length; i++) result[i] = solutions.get(i);
        return result;
    }

    private RDFEntry[] read_entries() throws Exception {
        if (is_read) {
            throw new Exception("File was already read.");
        }
        ArrayList<RDFEntry> entries = new ArrayList<>();
        while (true) {
            try {
                RDFEntry entry = read_next();
                if (entry.reaction == null) continue;
                entries.add(entry);
            }
            catch (RuntimeException e) {
                break;
            }
        }
        is_read = true;
        Set<String> ids = entries.stream().map(entry -> entry.id).collect(Collectors.toSet());
        if (ids.size() != entries.size()) {
            throw new Exception("Non-unique reaction IDs found.");
        }
        RDFEntry[] result = new RDFEntry[entries.size()];
        for (int i = 0; i < result.length; i++) result[i] = entries.get(i);
        return result;
    }

    private RDFEntry read_next() {
        RDFEntry entry = null;
        try {
            ArrayList<String> lines = new ArrayList<>();
            boolean fetch_lines = false;
            while (true) {
                String line = br.readLine();
                if (line == null) break;
                if (line.trim().equals("$RXN")) fetch_lines = true;
                if (fetch_lines && (line.trim().startsWith("$RFMT"))) break;
                if (fetch_lines) lines.add(line);
            }
            if (lines.size() > 0) {
                String[] reaction_lines = new String[lines.size()];
                for (int i = 0; i < reaction_lines.length; i++)reaction_lines[i] = lines.get(i);
                entry = new RDFEntry(reaction_lines);
            }
            else
                throw new RuntimeException("No data to parse");
        }
        catch (IOException e) {
            e.printStackTrace();
        }
        return entry;
    }
}
