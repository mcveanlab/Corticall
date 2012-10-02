package uk.ac.ox.well.indiana.utils.arguments;

import java.lang.annotation.ElementType;
import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;
import java.lang.annotation.Target;

@Retention(RetentionPolicy.RUNTIME)
@Target(ElementType.FIELD)
public @interface Output {
    String shortName() default "o";
    String fullName() default "out";
    String doc() default "the output file";
    //boolean required() default false;
}
