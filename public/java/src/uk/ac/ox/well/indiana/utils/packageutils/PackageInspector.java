package uk.ac.ox.well.indiana.utils.packageutils;

import org.reflections.Reflections;
import org.reflections.scanners.SubTypesScanner;
import org.reflections.util.ClasspathHelper;
import org.reflections.util.ConfigurationBuilder;

import java.net.URL;
import java.util.HashMap;
import java.util.Set;

public class PackageInspector<ClassType> {
    private Reflections reflections;
    private Class<ClassType> classType;

    public PackageInspector(Class<ClassType> classType) {
        Set<URL> classPath = ClasspathHelper.forPackage("uk.ac.ox.well.indiana");

        ConfigurationBuilder config = new ConfigurationBuilder();
        config.setUrls(classPath);
        config.setScanners(new SubTypesScanner());

        this.reflections = new Reflections(config);
        this.classType = classType;
    }

    public Set<Class<? extends ClassType>> getExtendingClasses() {
        return reflections.getSubTypesOf(classType);
    }

    public HashMap<String, Class<? extends ClassType>> getExtendingClassesMap() {
        HashMap<String, Class<? extends ClassType>> classHashMap = new HashMap<String, Class<? extends ClassType>>();

        for (Class<? extends ClassType> c : getExtendingClasses()) {
            classHashMap.put(c.getSimpleName(), c);
        }

        return classHashMap;
    }
}
