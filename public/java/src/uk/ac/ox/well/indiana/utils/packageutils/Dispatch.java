package uk.ac.ox.well.indiana.utils.packageutils;

import com.google.common.base.Joiner;
import org.apache.commons.lang.time.DurationFormatUtils;
import uk.ac.ox.well.indiana.Main;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.performance.PerformanceUtils;

import java.io.PrintStream;
import java.lang.annotation.Annotation;
import java.lang.management.ManagementFactory;
import java.lang.management.RuntimeMXBean;
import java.lang.reflect.Field;
import java.util.*;

/**
 * Instantiates a module and passes it command-line arguments.
 */
public class Dispatch {
    /**
     * Private constructor - this class cannot be instantiated!
     */
    private Dispatch() {}

    /**
     * Main method for instantiating INDIANA modules
     *
     * @param moduleName  the name of the module
     * @param moduleArgs  the arguments for the module
     */
    @SuppressWarnings("unchecked")
    public static void main(String moduleName, String[] moduleArgs) {
        try {
            Class<? extends Module> module = (Class<? extends Module>) Class.forName(moduleName);

            Module instance = module.newInstance();
            instance.args = moduleArgs;

            Set<String> specifiedArgs = new HashSet<>();
            for (int i = 0; i < moduleArgs.length - 1; i+=2) {
                specifiedArgs.add(moduleArgs[i]);
            }

            Field[] instanceFields = instance.getClass().getDeclaredFields();
            Field[] superFields = instance.getClass().getSuperclass().getDeclaredFields();

            List<Field> fields = new ArrayList<>();
            fields.addAll(Arrays.asList(instanceFields));
            fields.addAll(Arrays.asList(superFields));

            List<String> defaultArgs = new ArrayList<>();
            for (Field field : fields) {
                for (Annotation annotation : field.getDeclaredAnnotations()) {
                    if (annotation.annotationType().equals(Argument.class)) {
                        Argument arg = (Argument) annotation;

                        if (!specifiedArgs.contains(arg.shortName()) && !specifiedArgs.contains(arg.fullName())) {
                            if (field.get(instance) != null) {
                                Object defaultValue = field.get(instance);

                                defaultArgs.add("-" + arg.shortName() + " " + defaultValue.toString());
                            }
                        }
                    }
                }
            }

            RuntimeMXBean runtimeMxBean = ManagementFactory.getRuntimeMXBean();
            List<String> arguments = runtimeMxBean.getInputArguments();
            String jvmArgs = Joiner.on(" ").join(arguments);
            String jar = System.getProperty("sun.java.command").split("\\s+")[0];
            String indianaCmd = instance.getClass().getSimpleName();
            String modArgs = Joiner.on(" ").join(moduleArgs);
            String defArgs = Joiner.on(" ").join(defaultArgs);
            String fullCmd = Joiner.on(" ").join("java", jvmArgs, "-jar", jar, indianaCmd, modArgs, defArgs);

            instance.init();

            Main.getLogger().info("{}", getBanner());
            //Main.getLogger().info("Invocation: {}", fullCmd);
            Main.getLogger().info("");

            Date startTime = new Date();

            instance.execute();

            for (Field field : fields) {
                for (Annotation annotation : field.getDeclaredAnnotations()) {
                    if (annotation.annotationType().equals(Output.class)) {
                        if (field.get(instance) instanceof PrintStream) {
                            ((PrintStream) field.get(instance)).close();
                        }
                    }
                }
            }

            Date elapsedTime = new Date((new Date()).getTime() - startTime.getTime());

            Main.getLogger().info("");
            Main.getLogger().info("Complete. (time) {}; (mem) {}",
                                  DurationFormatUtils.formatDurationHMS(elapsedTime.getTime()),
                                  PerformanceUtils.getCompactMemoryUsageStats() );
        } catch (ClassNotFoundException | IllegalAccessException | InstantiationException e) {
            e.printStackTrace();
        }
    }

    private static String getBanner() {
        Properties prop = Main.getBuildProperties();

        if (prop != null) {
            String version = prop.get("major.version") + "." + prop.get("minor.version") + "-" + prop.get("git.version");
            String gitDate = (String) prop.get("git.date");
            String buildDate = (String) prop.get("build.date");
            String dates = "(repo) " + gitDate + "; (build) " + buildDate;
            //String banner = "INDIANA " + version + "; " + dates;
            String banner = Main.progName + " " + version + "; " + dates;

            return banner;
        } else {
            return Main.progName + " unknown version; unknown build time";
        }
    }
}
