# -----------------------------------------------------------------------------
# 1. PAKETLER (en basta yuklenir)
# -----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
  library(scales)
})

# -----------------------------------------------------------------------------
# 2. KOMUT SATIRI ARGUMANLARI VE VERI YUKLEME
# -----------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]

# Sample adini cikti dosyasindan al
SAMPLE <- strsplit(basename(output_file), "_qc_plots_R.png")[[1]][1]
if (is.na(SAMPLE) || SAMPLE == "") {
  SAMPLE <- basename(dirname(output_file))
}

# VERIYI OKU
if (!file.exists(input_file)) {
  stop(paste("ERROR:", input_file, "not found!"))
}
df <- read.csv(input_file)

# -----------------------------------------------------------------------------
# 3. ISTATISTIK HESAPLAMALARI
# -----------------------------------------------------------------------------
total_reads <- nrow(df)
total_bases <- sum(df$length)
total_yield_mb <- total_bases / 1e6

# N50 hesapla
sorted_lengths <- sort(df$length, decreasing = TRUE)
cumsum_lengths <- cumsum(sorted_lengths)
n50_value <- sorted_lengths[which(cumsum_lengths >= total_bases / 2)[1]]

# Kalite kategorileri
q10_pct <- sum(df$mean_quality >= 10) / total_reads * 100
q20_pct <- sum(df$mean_quality >= 20) / total_reads * 100
q30_pct <- sum(df$mean_quality >= 30) / total_reads * 100

# -----------------------------------------------------------------------------
# 4. KONSOLA YAZDIR
# -----------------------------------------------------------------------------
cat("\n=== SUMMARY STATISTICS (R) ===\n")
cat(paste("GC Content   - Mean:", round(mean(df$gc_content), 2), "| Median:", round(median(df$gc_content), 2), "\n"))
cat(paste("Read Length  - Mean:", round(mean(df$length), 2), "| Median:", round(median(df$length), 2), "| N50:", n50_value, "\n"))
cat(paste("Quality Score- Mean:", round(mean(df$mean_quality), 2), "| Median:", round(median(df$mean_quality), 2), "\n"))
cat(paste("\nQ10 Above:", round(q10_pct, 1), "% | Q20 Above:", round(q20_pct, 1), "% | Q30 Above:", round(q30_pct, 1), "%\n"))
cat(paste("Total Yield:", round(total_yield_mb, 2), "Mb\n"))

# -----------------------------------------------------------------------------
# 5. GRAFIKLER - Python ile tutarli renkler
# -----------------------------------------------------------------------------

# --- GRAFIK 1: Read Length Density (Log Scale) ---
# Renk: steelblue (Python ile ayni)
p1 <- ggplot(df, aes(x = length)) +
  geom_density(fill = "steelblue", alpha = 0.7, color = "steelblue") +
  scale_x_log10(labels = comma) +
  geom_vline(xintercept = median(df$length), color = "orange", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = n50_value, color = "red", linetype = "dashed", linewidth = 1) +
  labs(title = paste0(SAMPLE, " - Read Length Density (Log Scale)"),
       x = "Length (bp) [log10]", y = "Density") +
  theme_minimal() +
  annotate("text", x = Inf, y = Inf,
           label = paste0("Mean: ", round(mean(df$length), 1),
                          "\nMedian: ", round(median(df$length), 1),
                          "\nN50: ", n50_value,
                          "\nTotal: ", total_reads),
           hjust = 1.1, vjust = 1.1, size = 3, color = "black")

# --- GRAFIK 2: Quality Score Boxplot ---
p2 <- ggplot(df, aes(x = "", y = mean_quality)) +
  geom_boxplot(fill = "coral", alpha = 0.7, color = "black", outlier.size = 0.5) +
  geom_hline(yintercept = 10, color = "darkred", linetype = "dotted", linewidth = 1) +
  geom_hline(yintercept = 20, color = "purple", linetype = "dotted", linewidth = 1) +
  geom_hline(yintercept = 30, color = "blue", linetype = "dotted", linewidth = 1) +
  labs(title = paste0(SAMPLE, " - Quality Score Distribution"),
       x = "", y = "Mean Quality Score") +
  theme_minimal() +
  annotate("text", x = Inf, y = Inf,
           label = paste0("Mean: ", round(mean(df$mean_quality), 1),
                          "\nMedian: ", round(median(df$mean_quality), 1),
                          "\nQ10+: ", round(q10_pct, 1), "%",
                          "\nQ20+: ", round(q20_pct, 1), "%",
                          "\nQ30+: ", round(q30_pct, 1), "%"),
           hjust = 1.1, vjust = 1.1, size = 3, color = "black")

# --- GRAFIK 3: GC Content Violin Plot ---
p3 <- ggplot(df, aes(x = "", y = gc_content)) +
  geom_violin(fill = "seagreen", alpha = 0.7, color = "seagreen") +
  geom_hline(yintercept = median(df$gc_content), color = "orange", linetype = "dashed", linewidth = 1) +
  geom_hline(yintercept = 52, color = "purple", linetype = "dashed", linewidth = 1) +
  labs(title = paste0(SAMPLE, " - GC Content Distribution"),
       x = "", y = "GC Content (%)") +
  theme_minimal() +
  annotate("text", x = Inf, y = Inf,
           label = paste0("Mean: ", round(mean(df$gc_content), 1), "%",
                          "\nMedian: ", round(median(df$gc_content), 1), "%"),
           hjust = 1.1, vjust = 1.1, size = 3, color = "black")

# -----------------------------------------------------------------------------
# 6. GRAFIKLERI BIRLESTIR VE KAYDET
# -----------------------------------------------------------------------------
combined <- grid.arrange(p1, p2, p3, ncol = 3)
ggsave(output_file, combined, width = 18, height = 6, dpi = 150)

cat(paste("\n[SUCCESS] R plot saved as", output_file, "\n"))