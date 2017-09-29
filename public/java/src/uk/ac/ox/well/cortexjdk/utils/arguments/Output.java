package uk.ac.ox.well.cortexjdk.utils.arguments;

import java.lang.annotation.ElementType;
import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;
import java.lang.annotation.Target;

@Retention(RetentionPolicy.RUNTIME)
@Target(ElementType.FIELD)
public @interface Output {
    String shortName() default "o";
    String fullName() default "output";
    String doc() default "The output file";
    //boolean required() default false;
}
