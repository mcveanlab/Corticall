package uk.ac.ox.well.indiana.utils.containers;

import com.google.common.base.Joiner;
import uk.ac.ox.well.indiana.utils.io.utils.LineReader;

import java.io.File;
import java.io.PrintStream;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Scanner;

public class DataTables {
    private Map<String, DataTable> tables = new LinkedHashMap<String, DataTable>();

    public DataTables() {}

    public DataTables(File in) {
        read(in);
    }

    public void addTable(String tableName, String description, String... columnNames) {
        if (!tables.containsKey(tableName)) {
            tables.put(tableName, new DataTable(tableName, description, columnNames));
        }
    }

    public void addTable(DataTable t) {
        if (!tables.containsKey(t.getTableName())) {
            tables.put(t.getTableName(), t);
        }
    }

    public DataTable getTable(String tableName) {
        return tables.get(tableName);
    }

    public void write(PrintStream out) {
        for (String tableName : tables.keySet()) {
            DataTable t = tables.get(tableName);

            t.write(out);
        }
    }

    public void read(File in) {
        LineReader lr = new LineReader(in);

        String currentTableName = "";
        String[] header = null;

        while (lr.hasNext()) {
            String l = lr.getNextRecord();

            if (!l.startsWith("^#") && !l.isEmpty()) {
                String[] fields = l.split("\\s+");
                String tableName = fields[0];

                if (!tableName.equals(currentTableName)) {
                    header = fields;

                    DataTable t = new DataTable(tableName, "unknown");

                    for (int i = 1; i < header.length; i++) {
                        t.addColumn(header[i]);
                    }

                    addTable(t);

                    currentTableName = tableName;
                } else {
                    DataTable t = getTable(currentTableName);
                    String pk = Joiner.on(".").join(fields);

                    for (int i = 1; i < header.length; i++) {
                        if (fields.length > i) {
                            String field = fields[i];

                            if (isInteger(field)) {
                                t.set(pk, header[i], Integer.valueOf(field));
                            } else if (isLong(field)) {
                                t.set(pk, header[i], Long.valueOf(field));
                            } else if (isFloat(field)) {
                                t.set(pk, header[i], Float.valueOf(field));
                            } else if (isDouble(field)) {
                                t.set(pk, header[i], Double.valueOf(field));
                            } else if (isBoolean(field)) {
                                t.set(pk, header[i], Boolean.valueOf(field));
                            } else {
                                t.set(pk, header[i], field);
                            }
                        } else {
                            t.set(pk, header[i], null);
                        }
                    }
                }
            }
        }
    }

    private static boolean isInteger(String s) {
        try {
            Integer.parseInt(s);
        } catch (NumberFormatException e) {
            return false;
        }

        return true;
    }

    private static boolean isLong(String s) {
        try {
            Long.parseLong(s);
        } catch (NumberFormatException e) {
            return false;
        }

        return true;
    }

    private static boolean isFloat(String s) {
        try {
            Float.parseFloat(s);
        } catch (NumberFormatException e) {
            return false;
        }

        return true;
    }

    private static boolean isDouble(String s) {
        try {
            Double.parseDouble(s);
        } catch (NumberFormatException e) {
            return false;
        }

        return true;
    }

    private static boolean isBoolean(String s) {
        return s.equals("true") || s.equals("false");
    }
}
