package uk.ac.ox.well.indiana;

import processing.core.PApplet;
import uk.ac.ox.well.indiana.sketches.Sketch;
import uk.ac.ox.well.indiana.tools.Tool;
import uk.ac.ox.well.indiana.utils.packageutils.IRunner;
import uk.ac.ox.well.indiana.utils.packageutils.PackageInspector;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Set;

public class IndianaMain {
    public static void main(String[] args) throws Exception {
        if (args.length == 0 || args[0].equals("-h") || args[0].equals("--help") || args[0].equals("-H")) {
            showPrimaryHelp();
        } else if (args.length > 0) {
            String moduleName = args[0];
            String[] moduleArgs = Arrays.copyOfRange(args, 1, args.length);

            HashMap<String, Class<? extends IndianaModule>> modules = new PackageInspector<IndianaModule>(IndianaModule.class).getExtendingClassesMap();

            if (!modules.containsKey(moduleName)) {
                showInvalidModuleMessage(moduleName);
            } else {
                Class module = modules.get(moduleName);

                if (Tool.class.isAssignableFrom(module)) {
                    IRunner.main(module.getName(), moduleArgs);
                } else if (Sketch.class.isAssignableFrom(module)) {
                    PApplet.main(module.getName(), moduleArgs);
                }
            }
        }
    }

    private static void showPrimaryHelp() {
        Set<Class<? extends Tool>> tools = new PackageInspector<Tool>(Tool.class).getExtendingClasses();
        Set<Class<? extends Sketch>> sketch = new PackageInspector<Sketch>(Sketch.class).getExtendingClasses();

        System.out.println("usage: java -jar indiana.jar [-h|--help|-H]");
        System.out.println("                             <command> [<args>]");
        System.out.println();

        System.out.println("tools:");
        for (Class t : tools) {
            System.out.println("   " + t.getSimpleName());
        }
        System.out.println();

        System.out.println("sketch:");
        for (Class s : sketch) {
            System.out.println("   " + s.getSimpleName());
        }
        System.out.println();

        System.exit(1);
    }

    private static void showInvalidModuleMessage(String module) {
        System.out.println("indiana: '" + module + "' is not a valid INDIANA module. See 'java -jar indiana.jar --help'.");

        System.exit(1);
    }
}