package uk.ac.ox.well.cortexjdk.utils.arguments;

import java.lang.annotation.ElementType;
import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;
import java.lang.annotation.Target;

@Retention(RetentionPolicy.RUNTIME)
@Target(ElementType.FIELD)
public @interface Argument {
    String shortName();
    String fullName();
    String doc();
    boolean required() default true;
}
