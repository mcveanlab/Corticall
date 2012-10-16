package uk.ac.ox.well.indiana.utils.packageutils;

import uk.ac.ox.well.indiana.tools.Tool;

public class IRunner {
    private IRunner() {}

    @SuppressWarnings("unchecked")
    public static void main(String moduleName, String[] moduleArgs) {
        try {
            Class<? extends Tool> module = (Class<? extends Tool>) Class.forName(moduleName);

            Tool instance = module.newInstance();
            instance.args = moduleArgs;

            instance.init();
            instance.execute();
        } catch (ClassNotFoundException e) {
            e.printStackTrace();
        } catch (InstantiationException e) {
            e.printStackTrace();
        } catch (IllegalAccessException e) {
            e.printStackTrace();
        }
    }
}
