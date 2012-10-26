package uk.ac.ox.well.indiana.utils.arguments;

import com.google.common.base.Joiner;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.Parser;
import org.apache.commons.cli.UnrecognizedOptionException;

import java.util.ArrayList;
import java.util.HashMap;

public class IndianaParser extends Parser {
    @Override
    protected String[] flatten(Options options, String[] strings, boolean stopAtNonOption) {
        HashMap<Option, ArrayList<String>> arguments = new HashMap<Option, ArrayList<String>>();

        Option currentOption = null;

        for (String token : strings) {
            if (token.startsWith("-")) {
                if (options.hasOption(token)) {
                    currentOption = options.getOption(token);

                    if (!arguments.containsKey(currentOption)) {
                        arguments.put(currentOption, new ArrayList<String>());
                    }
                } else {
                    throw new RuntimeException("The option '" + token + "' is not a recognized option");
                }
            } else {
                if (currentOption != null) {
                    arguments.get(currentOption).add(token);
                }
            }
        }

        ArrayList<String> args = new ArrayList<String>();
        for (Option option : arguments.keySet()) {
            args.add("--" + option.getLongOpt());
            args.add(Joiner.on(",").join(arguments.get(option)));
        }

        return args.toArray(new String[args.size()]);
    }
}
