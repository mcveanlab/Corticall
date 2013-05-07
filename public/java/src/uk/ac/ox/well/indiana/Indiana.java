package uk.ac.ox.well.indiana;

import ch.qos.logback.classic.Level;
import ch.qos.logback.classic.Logger;
import ch.qos.logback.classic.LoggerContext;
import ch.qos.logback.classic.encoder.PatternLayoutEncoder;
import ch.qos.logback.classic.spi.ILoggingEvent;
import ch.qos.logback.core.ConsoleAppender;
import org.apache.commons.lang.time.DurationFormatUtils;
import org.slf4j.LoggerFactory;
import processing.core.PApplet;
import uk.ac.ox.well.indiana.sketches.Sketch;
import uk.ac.ox.well.indiana.tools.Tool;
import uk.ac.ox.well.indiana.utils.packageutils.IRunner;
import uk.ac.ox.well.indiana.utils.packageutils.PackageInspector;
import uk.ac.ox.well.indiana.utils.performance.PerformanceUtils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.Map;

/**
 * Main class for INDIANA.  Sets up the logger, handles help message, selects the module to run and passes through command-line arguments.
 */
public class Indiana {
    private static Logger log = configureLogger();

    /**
     * Main method for INDIANA.  First argument must be the module to run.  All other arguments are passed through to the module for processing.
     *
     * @param args  Command-line arguments
     * @throws Exception
     */
    public static void main(String[] args) throws Exception {
        log.debug("Started up");

        if (args.length == 0 || args[0].equals("-h") || args[0].equals("--help")) {
            showPrimaryHelp();
        } else if (args.length > 0) {
            Date startTime = new Date();

            String moduleName = args[0];
            String[] moduleArgs = Arrays.copyOfRange(args, 1, args.length);

            Map<String, Class<? extends IndianaModule>> modules = new PackageInspector<IndianaModule>(IndianaModule.class).getExtendingClassesMap();

            if (!modules.containsKey(moduleName)) {
                showInvalidModuleMessage(moduleName);
            } else {
                Class module = modules.get(moduleName);

                if (Tool.class.isAssignableFrom(module)) {
                    IRunner.main(module.getName(), moduleArgs);
                } else if (Sketch.class.isAssignableFrom(module)) {
                    ArrayList<String> newArgs = new ArrayList<String>();
                    newArgs.add(module.getName());
                    newArgs.addAll(Arrays.asList(moduleArgs));

                    PApplet.main(newArgs.toArray(new String[newArgs.size()]));
                }
            }

            Date elapsedTime = new Date((new Date()).getTime() - startTime.getTime());

            log.info("Performance:");
            log.info("\ttime: {}", DurationFormatUtils.formatDurationHMS(elapsedTime.getTime()));
            log.info("\t mem: {}", PerformanceUtils.getCompactMemoryUsageStats());
        }

        log.debug("Finished");
    }

    /**
     * Configure a logger that contains the log level, timestamp, module, method, and line number of the logging statement.
     *
     * @return  A fully-configured logger
     */
    private static Logger configureLogger() {
        Logger rootLogger = (Logger) LoggerFactory.getLogger(Indiana.class);

        LoggerContext loggerContext = rootLogger.getLoggerContext();
        loggerContext.reset();

        PatternLayoutEncoder encoder = new PatternLayoutEncoder();
        encoder.setContext(loggerContext);
        encoder.setPattern("%level [%date{dd/MM/yy HH:mm:ss} %class{0}.%M:%L]: %message%n");
        encoder.start();

        ConsoleAppender<ILoggingEvent> appender = new ConsoleAppender<ILoggingEvent>();
        appender.setContext(loggerContext);
        appender.setEncoder(encoder);
        appender.start();

        rootLogger.addAppender(appender);

        String logLevel = System.getProperty("indiana.loglevel");
        if (logLevel != null) {
            if      (logLevel.equalsIgnoreCase("OFF"))   { rootLogger.setLevel(Level.OFF);   }
            else if (logLevel.equalsIgnoreCase("TRACE")) { rootLogger.setLevel(Level.TRACE); }
            else if (logLevel.equalsIgnoreCase("DEBUG")) { rootLogger.setLevel(Level.DEBUG); }
            else if (logLevel.equalsIgnoreCase("INFO"))  { rootLogger.setLevel(Level.INFO);  }
            else if (logLevel.equalsIgnoreCase("WARN"))  { rootLogger.setLevel(Level.WARN);  }
            else if (logLevel.equalsIgnoreCase("ERROR")) { rootLogger.setLevel(Level.ERROR); }
            else if (logLevel.equalsIgnoreCase("ALL"))   { rootLogger.setLevel(Level.ALL);   }
        } else {
            rootLogger.setLevel(Level.INFO);
        }

        return rootLogger;
    }

    /**
     * Get the configured logger.
     *
     * @return  A fully-configured logger
     */
    public static Logger getLogger() {
        return log;
    }

    /**
     * List all of the available modules, grouped by package.
     */
    private static void showPrimaryHelp() {
        Map<String, Map<String, Class<? extends Tool>>> tools = new PackageInspector<Tool>(Tool.class).getExtendingClassTree();
        Map<String, Map<String, Class<? extends Sketch>>> sketches = new PackageInspector<Sketch>(Sketch.class).getExtendingClassTree();

        System.out.println();
        System.out.println("usage: java -jar indiana.jar [-h|--help]");
        System.out.println("                             <command> [<args>]");
        System.out.println();

        System.out.println("tools:");
        for (String p : tools.keySet()) {
            System.out.println("   " + p);

            for (String t : tools.get(p).keySet()) {
                System.out.println("      " + t);
            }
        }
        System.out.println();

        System.out.println("sketches:");
        for (String p : sketches.keySet()) {
            System.out.println("   " + p);

            for (String t : sketches.get(p).keySet()) {
                System.out.println("      " + t);
            }
        }
        System.out.println();
    }

    /**
     * Show error message when a requested module is not available.
     *
     * @param module  The name of the requested module
     */
    private static void showInvalidModuleMessage(String module) {
        System.out.println("indiana: '" + module + "' is not a valid INDIANA module. See 'java -jar indiana.jar --help'.");

        System.exit(1);
    }
}