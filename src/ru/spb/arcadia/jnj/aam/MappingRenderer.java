package ru.spb.arcadia.jnj.aam;

import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IReaction;
import uk.ac.ebi.reactionblast.tools.ImageGenerator;
import ru.spb.arcadia.jnj.aam.io.RDFReader;

import java.io.File;

public class MappingRenderer {
    public static void entry_point(String id_field, String input_file, String output_dir) {
        try {
            RDFReader reader = new RDFReader(input_file, id_field);
            File dir = new File(output_dir);
            if (!dir.exists()) dir.mkdirs();
            MappingSolution[] solutions = reader.read_solutions();
            for (MappingSolution solution : solutions) {
                IReaction mapped_reaction = solution.getBondChangeCalculator().getReactionWithCompressUnChangedHydrogens();
                int counter = 0;
                System.out.println("Reaction ID: " + mapped_reaction.getID());
                System.out.println("Reactants info");
                for (IAtomContainer ac : mapped_reaction.getReactants().atomContainers()) {
                    ac.setID(Integer.toString(counter++));
                    System.out.println(ac.getBondCount());
                }
                System.out.println("Products info");
                for (IAtomContainer ac : mapped_reaction.getProducts().atomContainers()) {
                    ac.setID(Integer.toString(counter++));
                    System.out.println(ac.getBondCount());
                }
//                try {
//                    File file = new File(output_dir);
//                    new ImageGenerator().drawTopToBottomReactionLayout(file, mapped_reaction, mapped_reaction.getID());
//                    System.out.println("Rendered reaction " + mapped_reaction.getID());
//                }
//                catch (Exception e2) {
//                    System.out.println("Failed to render: " + mapped_reaction.getID());
//                    e2.printStackTrace();
//                }
                break;
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
}
