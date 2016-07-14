package uk.ac.ox.well.indiana.utils.progress;

import org.slf4j.Logger;

public class ProgressMeter {
    private Logger log;
    private long currentRecord = 0;
    private long updateRecord = 1000000;
    private long maxRecord = 0;
    private String header = "";
    private String message = "";
    private String indent = "";

    public ProgressMeter(Logger log, long currentRecord, long updateRecord, long maxRecord, String header, String message, String indent) {
        this.log = log;
        this.currentRecord = currentRecord;
        this.updateRecord = updateRecord;
        this.maxRecord = maxRecord;
        this.header = header;
        this.message = message;
        this.indent = indent;

        log.info("{}", this.header);
    }

    public void update() {
        update(message);
    }

    public void update(String newMessage) {
        if (currentRecord % updateRecord == 0) {
            if (maxRecord > 0) {
                log.info("{}{}/{} {}", indent, currentRecord, maxRecord, newMessage);
            } else {
                log.info("{}{} {} ", indent, currentRecord, newMessage);
            }
        }

        currentRecord++;
    }

    public void reset() {
        currentRecord = 0;
    }
}
