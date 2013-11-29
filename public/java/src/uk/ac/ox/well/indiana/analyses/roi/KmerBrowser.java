package uk.ac.ox.well.indiana.analyses.roi;

import uk.ac.ox.well.indiana.tools.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexMap;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;

public class KmerBrowser extends Module {
    @Argument(fullName="cortexGraph", shortName="cg", doc="Cortex graph")
    public CortexMap CORTEX_MAP;

    @Override
    public void execute() {
        BufferedReader br = new BufferedReader(new InputStreamReader(System.in));

        String command = "noop noop";
        do {
            System.out.print("kb> ");

            try {
                command = br.readLine();
                String[] pieces = command.split("\\s+");

                String op = pieces[0];

                if (op.contains("show")) {
                    for (int i = 1; i < pieces.length; i++) {
                        String kmer = pieces[1];


                    }
                }
            } catch (IOException e) {
                throw new RuntimeException("Error reading command from stdin.", e);
            }
        } while (!command.contains("quit"));
    }
}
