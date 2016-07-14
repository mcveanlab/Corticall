package uk.ac.ox.well.indiana.utils.progress;

import org.slf4j.Logger;

public class ProgressMeterFactory {
    private long currentRecord = 0;
    private long updateRecord = 1000000;
    private long maxRecord = 0;
    private String header = "";
    private String message = "";
    private String indent = "  ";

    public ProgressMeterFactory currentRecord(int currentRecord) { this.currentRecord = currentRecord; return this; }
    public ProgressMeterFactory updateRecord(int updateRecord) { this.updateRecord = updateRecord; return this; }
    public ProgressMeterFactory maxRecord(int maxRecord) { this.maxRecord = maxRecord; return this; }
    public ProgressMeterFactory header(String header) { this.header = header; return this; }
    public ProgressMeterFactory message(String message) { this.message = message; return this; }
    public ProgressMeterFactory indent(String indent) { this.indent = indent; return this; }

    public ProgressMeter make(Logger log) {
        return new ProgressMeter(log, currentRecord, updateRecord, maxRecord, header, message, indent);
    }
}
