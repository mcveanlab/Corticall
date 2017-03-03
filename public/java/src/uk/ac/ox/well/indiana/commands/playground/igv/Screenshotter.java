package uk.ac.ox.well.indiana.commands.playground.igv;

import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;

import java.io.*;
import java.net.Socket;
import java.net.UnknownHostException;
import java.util.ArrayList;

/**
 * Created by kiran on 16/11/2016.
 */
public class Screenshotter extends Module {
    @Argument(fullName="locus", shortName="l", doc="Loci")
    public ArrayList<String> LOCI;

    @Override
    public void execute() {
        try {
            Socket socket = new Socket("127.0.0.1", 60151);
            PrintWriter out = new PrintWriter(socket.getOutputStream(), true);
            BufferedReader in = new BufferedReader(new InputStreamReader(socket.getInputStream()));

            out.println("goto chr1:65,827,301");
            String response = in.readLine();
            log.info("response: {}", response);

            /*
            out.println("load na12788.bam,n12788.tdf");
            String response = in.readLine();
            System.out.println(response);

            out.println("genome hg18");
            response = in.readLine();
            System.out.println(response);

            out.println("goto chr1:65,827,301");
            response = in.readLine();
            System.out.println(response);

            out.println("snapshotDirectory /screenshots");
            response = in.readLine();
            System.out.println(response);

            out.println("snapshot");
            response = in.readLine();
            System.out.println(response);
            */
        } catch (UnknownHostException e) {
            throw new IndianaException("Unknown host", e);
        } catch (IOException e) {
            throw new IndianaException("IOException", e);
        }
    }
}
