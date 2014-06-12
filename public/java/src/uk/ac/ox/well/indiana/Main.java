package uk.ac.ox.well.indiana;

import ch.qos.logback.classic.Level;
import ch.qos.logback.classic.Logger;
import ch.qos.logback.classic.LoggerContext;
import ch.qos.logback.classic.encoder.PatternLayoutEncoder;
import ch.qos.logback.classic.spi.ILoggingEvent;
import ch.qos.logback.core.ConsoleAppender;
import org.slf4j.LoggerFactory;
import processing.core.PApplet;
import uk.ac.ox.well.indiana.commands.Command;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.commands.Sketch;
import uk.ac.ox.well.indiana.utils.arguments.Description;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.packageutils.Dispatch;
import uk.ac.ox.well.indiana.utils.packageutils.PackageInspector;

import java.io.IOException;
import java.io.InputStream;
import java.lang.management.ManagementFactory;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Map;
import java.util.Properties;

/**
 * Main class for INDIANA.  Sets up the logger, handles help message, selects the module to run and passes through command-line arguments.
 */
public class Main {
    private static Logger log = configureLogger();

    public static String progName;
    public static String progDesc;
    public static String rootPackage;
    public static String commandPackage;

    /**
     * Main method for INDIANA.  First argument must be the module to run.  All other arguments are passed through to the module for processing.
     *
     * @param args  Command-line arguments
     * @throws Exception
     */
    public static void start(String newProgName, String newProgDesc, String newRootPackage, String newCommandPackage, String[] args) throws Exception {
        progName = newProgName;
        progDesc = newProgDesc;
        rootPackage = newRootPackage;
        commandPackage = newCommandPackage;

        org.apache.log4j.PropertyConfigurator.configure("/dev/null");

        if (args.length == 0 || args[0].equals("-h") || args[0].equals("--help")) {
            showPrimaryHelp();
        } else if (args.length > 0) {
            String moduleName = args[0];
            String[] moduleArgs = Arrays.copyOfRange(args, 1, args.length);

            //Map<String, Class<? extends Command>> modules = new PackageInspector<Command>(Command.class, rootPackage).getExtendingClassesMap();
            Map<String, Class<? extends Module>> modules = new PackageInspector<Module>(Module.class, rootPackage).getExtendingClassesMap();

            if (!modules.containsKey(moduleName)) {
                showInvalidModuleMessage(moduleName);
            } else {
                Class module = modules.get(moduleName);

                if (Module.class.isAssignableFrom(module)) {
                    Dispatch.main(module.getName(), moduleArgs);
                } else if (Sketch.class.isAssignableFrom(module)) {
                    ArrayList<String> newArgs = new ArrayList<String>();
                    newArgs.add(module.getName());
                    newArgs.addAll(Arrays.asList(moduleArgs));

                    PApplet.main(newArgs.toArray(new String[newArgs.size()]));
                }
            }

        }
    }

    /**
     * Print the startup banner
     */
    private static void printBanner() {
        Properties prop = getBuildProperties();

        String version = prop.get("major.version") + "." + prop.get("minor.version") + "-" + prop.get("git.version");
        String gitDate = (String) prop.get("git.date");
        String buildDate = (String) prop.get("build.date");
        String dates = "(repo) " + gitDate + "; (build) " + buildDate;
        log.info(progName + " {}; {}", version, dates);
    }

    /**
     * Get the process id for this instance
     *
     * @return the process id
     */
    private static String getProcessID() {
        String vmName = ManagementFactory.getRuntimeMXBean().getName();

        if (vmName.contains("@")) {
            return vmName.substring(0, vmName.indexOf('@'));
        }

        return vmName;
    }

    /**
     * Configure a logger that contains the log level, timestamp, module, method, and line number of the logging statement.
     *
     * @return  A fully-configured logger
     */
    private static Logger configureLogger() {
        Logger rootLogger = (ch.qos.logback.classic.Logger) LoggerFactory.getLogger(Indiana.class);

        String logLevel = System.getProperty("loglevel");
        LoggerContext loggerContext = rootLogger.getLoggerContext();
        loggerContext.reset();

        PatternLayoutEncoder encoder = new PatternLayoutEncoder();
        encoder.setContext(loggerContext);

        if (logLevel != null && logLevel.equals("DEBUG")) {
            encoder.setPattern("%level [%date{dd/MM/yy HH:mm:ss} " + getProcessID() + " %class{0}:%L] %message%n");
        } else {
            encoder.setPattern("%.-1level [%date{yyyy-MM-dd HH:mm} " + getProcessID() + "] %message%n");
        }

        encoder.start();

        ConsoleAppender<ILoggingEvent> appender = new ConsoleAppender<ILoggingEvent>();
        appender.setContext(loggerContext);
        appender.setEncoder(encoder);
        appender.start();

        rootLogger.addAppender(appender);

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
     * Extract the build properties from the file automatically built at compile time.
     *
     * @return a populated Properties object
     */
    public static Properties getBuildProperties() {
        InputStream propStream = Indiana.class.getClassLoader().getResourceAsStream("build.properties");

        if (propStream != null) {
            try {
                Properties prop = new Properties();
                prop.load(propStream);

                return prop;
            } catch (IOException e) {
                throw new IndianaException("Unable to read build.properties file from within package", e);
            }
        }

        return null;
    }

    /**
     * List all of the available modules, grouped by package.
     */
    private static void showPrimaryHelp() {
        //Map<String, Class<? extends Command>> commands = new PackageInspector<Command>(Command.class, commandPackage).getExtendingClassesMap();
        Map<String, Class<? extends Module>> commands = new PackageInspector<Module>(Module.class, commandPackage).getExtendingClassesMap();

        int maxlength = 0;
        for (String t : commands.keySet()) {
            maxlength = (t.length() > maxlength) ? t.length() : maxlength;
        }

        Properties prop = getBuildProperties();

        System.out.println();
        System.out.println("Program: " + progName + " (" + progDesc + ")");
        System.out.println();
        System.out.format("Version: %s.%s.%s (%s)%n", prop.get("major.version"), prop.get("minor.version"), prop.get("git.version"), prop.getProperty("git.version.long"));
        System.out.println("Times:   (repo) " + prop.get("git.date") + "; (build) " + prop.get("build.date"));
        System.out.println("Contact: Kiran V Garimella <kiran@well.ox.ac.uk>");
        System.out.println();

        System.out.println("Usage:   java -jar " + progName.toLowerCase() + ".jar <command> [options]");
        System.out.println();

        System.out.print("Command:");
        int padwidth = 1;
        for (String t : commands.keySet()) {
            Description d = commands.get(t).getAnnotation(Description.class);
            String description = (d == null) ? "no description available" : d.text();

            System.out.format("%" + padwidth + "s%-" + maxlength + "s     %s%n", "", t, description);

            padwidth = 9;
        }

        System.out.println();
    }

    /**
     * Show error message when a requested module is not available.
     *
     * @param module  The name of the requested module
     */
    private static void showInvalidModuleMessage(String module) {
        System.out.println(progName.toLowerCase() + ": '" + module + "' is not a valid module. See 'java -jar " + progName.toLowerCase() + ".jar --help'.");

        System.exit(1);
    }
}
