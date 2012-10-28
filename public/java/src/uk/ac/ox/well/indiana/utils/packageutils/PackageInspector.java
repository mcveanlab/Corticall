package uk.ac.ox.well.indiana.utils.packageutils;

import org.reflections.Reflections;
import org.reflections.scanners.SubTypesScanner;
import org.reflections.util.ClasspathHelper;
import org.reflections.util.ConfigurationBuilder;

import java.lang.reflect.Modifier;
import java.net.URL;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

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
        Set<Class<? extends ClassType>> extendingClasses = reflections.getSubTypesOf(classType);
        Set<Class<? extends ClassType>> nonAbstractClasses = new HashSet<Class<? extends ClassType>>();

        for (Class<? extends ClassType> c : extendingClasses) {
            if (!Modifier.isAbstract(c.getModifiers())) {
                nonAbstractClasses.add(c);
            }
        }

        return nonAbstractClasses;
    }

    public Map<String, Class<? extends ClassType>> getExtendingClassesMap() {
        Map<String, Class<? extends ClassType>> classHashMap = new TreeMap<String, Class<? extends ClassType>>();

        for (Class<? extends ClassType> c : getExtendingClasses()) {
            classHashMap.put(c.getSimpleName(), c);
        }

        return classHashMap;
    }
}
