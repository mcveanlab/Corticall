plotScree <- function(pca, n=10, ylim=1.0) {
    propOfVar = summary(pca)$importance[2,1:n];
    cumPropOfVar = summary(pca)$importance[3,1:n];

    plot(1:n, propOfVar, type="b", lwd=2, col="black", ylim=c(0, ylim), xlab="Principal component", ylab="Proportion of variance", cex=1.3, cex.lab=1.3, cex.axis=1.3, bty="n");
    points(1:n, cumPropOfVar, type="b", lwd=2, col="red");

    legend("topleft", c("Variance", "Cumulative variance"), lwd=2, col=c("black", "red"), cex=1.3);
}

theme_kiran <- function(base_size = 12, base_family = "") {
    theme_grey(base_size = base_size, base_family = base_family) %+replace%
    theme(
        axis.text        = element_text(size = rel(1.0)),
        axis.ticks       = element_line(colour = "black"),
        legend.key       = element_rect(colour = "grey80"),
        panel.background = element_rect(fill = "white", colour = NA),
        panel.border     = element_rect(fill = NA, colour = NA),
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        strip.background = element_rect(fill = "grey80", colour = "grey50")
    )
}

vcf.getSamples <- function(vcf_file) {
    q = readLines(vcf_file, n=1000);
    h = q[grep("#CHROM", q)];

    fields = unlist(strsplit(h, "\t"));
    s = fields[10:length(fields)];

    return(s);
}

vcf.read <- function(vcf_file) {
    q = read.table(vcf_file, header=FALSE, stringsAsFactors=FALSE, sep="\t");
    vcfSamples = vcf.getSamples(vcf_file);
    header = c('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', vcfSamples);
    names(q) = header;

    return(q);
}
