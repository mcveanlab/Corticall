package uk.ac.ox.well.cortexjdk.utils.progress;

import org.slf4j.Logger;

public class ProgressMeter {
    private Logger log;
    private long currentRecord = 1;
    private long updateRecord = 0;
    private long maxRecord = 0;
    private String header = "";
    private String message = "";
    private String indent = "";
    private int updateTime = 120;

    private long startTime = System.currentTimeMillis();

    public ProgressMeter(Logger log, long currentRecord, long updateRecord, long maxRecord, int updateTime, String header, String message, String indent) {
        this.log = log;
        this.currentRecord = currentRecord;
        this.updateRecord = updateRecord;
        this.maxRecord = maxRecord;
        this.updateTime = updateTime;
        this.header = header;
        this.message = message;
        this.indent = indent;

        log.info("{}", this.header);
    }

    public void update() {
        update(message);
    }

    public void update(String newMessage) {
        if (currentRecord % updateRecord == 0 || currentRecord == 1 || currentRecord == maxRecord || isTimeToUpdate()) {
            if (maxRecord > 0) {
                log.info("{}{}/{} ({}%) {}", indent, currentRecord, maxRecord, String.format("%.2f", 100.0*((double)currentRecord)/((double)maxRecord)), newMessage);
            } else {
                log.info("{}{} {} ", indent, currentRecord, newMessage);
            }

            startTime = System.currentTimeMillis();
        }

        currentRecord++;
    }

    public void reset() { currentRecord = 1; }

    public long pos() { return currentRecord; }

    private boolean isTimeToUpdate() {
        if (updateTime <= 0) { return false; }

        long currentTime = System.currentTimeMillis();
        return (currentTime - startTime) / 1000 >= updateTime;
    }
}
