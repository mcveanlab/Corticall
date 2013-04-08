package uk.ac.ox.well.indiana.utils.io.utils;

import com.google.common.base.Joiner;

import java.util.*;

public class TableWriter {
    private String primaryKey;
    private Object primaryKeyDefault;

    //          column
    private Set<String> columnNames = new LinkedHashSet<String>();

    //          column  type
    private Map<String, Class> columnTypes = new HashMap<String, Class>();

    //          column  default value
    private Map<String, Object> columnDefaultValues = new HashMap<String, Object>();

    //          primary     column  datum
    private Map<Object, Map<String, Object>> data = new HashMap<Object, Map<String, Object>>();

    public <T> void addPrimaryKey(String primaryKey, T defaultValue) {
        this.primaryKey = primaryKey;
        this.primaryKeyDefault = defaultValue;
    }

    public <T> void addColumn(String columnName, T defaultValue) {
        columnNames.add(columnName);
        columnTypes.put(columnName, defaultValue.getClass());
        columnDefaultValues.put(columnName, defaultValue);
    }

    public void set(Object primaryKey, String columnName, Object value) {
        if (!data.containsKey(primaryKey)) {
            data.put(primaryKey, new HashMap<String, Object>());
        }

        data.get(primaryKey).put(columnName, value);
    }

    public Object get(Object primaryKey, String columnName) {
        Object value = data.get(primaryKey).get(columnName);

        return (columnTypes.get(columnName).cast(value));
    }

    public String toString() {
        StringBuilder builder = new StringBuilder();

        builder.append(primaryKey + "\t" + Joiner.on("\t").join(columnNames)).append("\n");

        for (Object primaryKey : data.keySet()) {
            ArrayList<String> row = new ArrayList<String>();

            row.add(primaryKeyDefault.getClass().cast(primaryKey).toString());

            for (String columnName : columnNames) {
                Object value = (data.containsKey(primaryKey) && data.get(primaryKey).containsKey(columnName)) ? data.get(primaryKey).get(columnName) : columnDefaultValues.get(columnName);

                //row.add(columnTypes.get(columnName).cast(value).toString());
                row.add(value.toString());
            }

            builder.append(Joiner.on("\t").join(row)).append("\n");
        }

        return builder.toString();
    }
}
