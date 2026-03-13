# PYTHON VE R YOLLARI - Conda/Environment uyumlu 
PYTHON = "python"
RSCRIPT = "Rscript"

import os
import glob
import shutil
import time

# OTOMATIK DATA BULMA - data/ klasoru yoksa olustur, FASTQ'yu bul ve tasi

# 1. data klasoru var mi kontrol et, yoksa olustur
if not os.path.exists("data"):
    os.makedirs("data")
    print("📁 'data' klasoru olusturuldu.")

# 2. Ana klasordeki .fastq dosyalarini bul ve data/ klasorune tasi
current_dir_fastq = glob.glob("*.fastq") + glob.glob("*.fastq.gz")
if current_dir_fastq and not glob.glob("data/*.fastq") and not glob.glob("data/*.fastq.gz"):
    for f in current_dir_fastq:
        shutil.move(f, os.path.join("data", f))
        print(f"📦 {f} -> data/ klasorune tasindi.")

# 3. Data klasorundeki tum fastq'lari bul
fastq_files = glob.glob("data/*.fastq") + glob.glob("data/*.fastq.gz")
fastq_names = [os.path.basename(f).replace(".fastq", "").replace(".gz", "") for f in fastq_files]

print(f"\n📁 Data klasorunde bulunan dosyalar: {fastq_names}")


# DATA SECIMI - KULLANIM SEKLINE GORE OTOMATIK veya MANUEL
#
# KULLANIM SEKILLERI:
# -------------------
# 1️⃣ TEK DATA VARSA OTOMATIK:
#    snakemake --cores 1
#    (data klasorunde tek dosya varsa onu isler)
#
# 2️⃣ TEK DATA SECMEK ICIN:
#    snakemake --cores 1 --config sample=data2
#
# 3️⃣ COKLU DATA SECMEK ICIN:
#    snakemake --cores 1 --config samples='["data1","data2","data3"]'
#
# 4️⃣ TUM DATALARI ISLEMEK ICIN:
#    snakemake --cores 1 --config samples='["hepsi"]'
#    (data klasorundeki TUM dosyalari isler)


# Komut satirindan gelen parametreler
SAMPLE = config.get("sample", None)
SAMPLES = config.get("samples", None)

# OZEL DURUM: --config samples='["hepsi"]' yazilirsa tum datalari isle
if SAMPLES == ["hepsi"]:
    SAMPLES = fastq_names
    print(f"\n TUM DATALAR islenecek: {SAMPLES}")


# SECIM MANTIGI

if not SAMPLE and not SAMPLES:
    # Hic parametre verilmemis
    if len(fastq_names) == 1:
        TARGET_SAMPLES = fastq_names
        print(f"\n✅ Tek data bulundu -> {TARGET_SAMPLES[0]} otomatik secildi.")
    elif len(fastq_names) == 0:
        raise ValueError("❌ HATA: data/ klasorunde hic .fastq dosyasi bulunamadi!")
    else:
        print(f"\n⚠️  {len(fastq_names)} data bulundu. Lutfen secim yapin:")
        print("   Tek data icin: snakemake --cores 1 --config sample=DATA_ADI")
        print("   Coklu data icin: snakemake --cores 1 --config samples='[\"data1\",\"data2\"]'")
        print("   Tum datalar icin: snakemake --cores 1 --config samples='[\"hepsi\"]'")
        raise ValueError("Secim yapilmadi!")

elif SAMPLE and not SAMPLES:
    # Tek data secilmis
    if SAMPLE in fastq_names:
        TARGET_SAMPLES = [SAMPLE]
        print(f"\n✅ Secilen data: {SAMPLE}")
    else:
        raise ValueError(f"❌ HATA: {SAMPLE} data klasorunde bulunamadi! Mevcut datalar: {fastq_names}")

elif SAMPLES and not SAMPLE:
    # Coklu data secilmis
    missing = [s for s in SAMPLES if s not in fastq_names]
    if missing:
        raise ValueError(f"❌ HATA: Su datalar bulunamadi: {missing}\n   Mevcut datalar: {fastq_names}")
    TARGET_SAMPLES = SAMPLES
    print(f"\n✅ Secilen datalar ({len(TARGET_SAMPLES)} adet): {TARGET_SAMPLES}")

else:
    raise ValueError("❌ HATA: Ayni anda hem 'sample' hem 'samples' kullanamazsin!")


# PIPELINE KURALLARI - HER DATA KENDI KLASORUNDE, DOSYA ISIMLERI SAMPLE ICERIYOR

rule all:
    input:
        expand("results/{sample}/nanoqc_report/{sample}_nanoQC.html", sample=TARGET_SAMPLES),
        expand("results/{sample}/nanoqc_report/{sample}_nanoQC.log", sample=TARGET_SAMPLES),
        expand("results/{sample}/{sample}_read_stats.csv", sample=TARGET_SAMPLES),
        expand("results/{sample}/{sample}_summary_stats.csv", sample=TARGET_SAMPLES),
        expand("results/{sample}/{sample}_report.xlsx", sample=TARGET_SAMPLES),
        expand("results/{sample}/{sample}_summary_table.png", sample=TARGET_SAMPLES),
        expand("results/{sample}/{sample}_qc_plots_python.png", sample=TARGET_SAMPLES),
        expand("results/{sample}/{sample}_qc_plots_R.png", sample=TARGET_SAMPLES)

rule nanoqc:
    input:
        fastq = "data/{sample}.fastq"
    output:
        html = "results/{sample}/nanoqc_report/{sample}_nanoQC.html",
        log = "results/{sample}/nanoqc_report/{sample}_nanoQC.log"
    run:
        # 1. NanoQC'yi calistir
        shell("nanoqc -o results/{wildcards.sample}/nanoqc_report {input.fastq}")
        
        # 2. Dosya sisteminin tamamlanmasi icin kısa bir bekleme
        time.sleep(1)
        
        # 3. Dosya isimlerini Python ile degistir (Windows uyumlu)
        import os
        import shutil
        
        base_path = f"results/{wildcards.sample}/nanoqc_report"
        
        # HTML dosyasini tasi
        old_html = os.path.join(base_path, "nanoQC.html")
        if os.path.exists(old_html):
            shutil.move(old_html, output.html)
            
        # Log dosyasini tasi (Buyuk/Kucuk harf duyarlılıgı kontroluyle)
        for log_name in ["NanoQC.log", "nanoQC.log"]:
            old_log = os.path.join(base_path, log_name)
            if os.path.exists(old_log):
                shutil.move(old_log, output.log)
                break

rule analyze_reads:
    input:
        "data/{sample}.fastq"
    output:
        read_stats = "results/{sample}/{sample}_read_stats.csv",
        summary_stats = "results/{sample}/{sample}_summary_stats.csv",
        report = "results/{sample}/{sample}_report.xlsx",
        summary_table = "results/{sample}/{sample}_summary_table.png"
    shell:
        "{PYTHON} scripts/analyze_reads.py {input} {output.read_stats} {output.summary_stats} {output.report} {output.summary_table}"

rule visualize_python:
    input:
        read_stats = "results/{sample}/{sample}_read_stats.csv"
    output:
        plot = "results/{sample}/{sample}_qc_plots_python.png"
    shell:
        "{PYTHON} scripts/visualize.py {input.read_stats} {output.plot}"

rule visualize_R:
    input:
        read_stats = "results/{sample}/{sample}_read_stats.csv"
    output:
        plot = "results/{sample}/{sample}_qc_plots_R.png"
    shell:
        "{RSCRIPT} scripts/visualize.R {input.read_stats} {output.plot}"