# FastQ File Download Reference and Notes

**Reference:**
[Qiime2 NCBI to QIIME2 Tutorial](https://github.com/stephanbitterwolf/Qiime2/blob/main/05_NCBI_to_QIIME2)

**Personal Notes: Downloading FastQ Files**

1. **Set up QIIME2 environment:**

```bash
conda activate qiime2-amplicon-2024.5
```

2. **Download FastQ files in parallel (example using SRR accession list):**

```bash
cat SRR_Acc_List.txt | parallel fastq-dump -X 10000 --split-files {}
```

3. **Organize downloaded files:**

```bash
mkdir fastq
mv *.fastq fastq/
```

4. **Verify download count:**

```bash
ls | wc -l
```

*This ensures the correct number of files were downloaded.*

5. **Notes on specific datasets:**

* **McLeod data:** Files appear merged; final length = 253 (2\*153 with overlapping reads removed).
* **Study 15:** Some files were missing; downloaded individually using:

```bash
fastq-dump --split-files -X 10000 SRR2143837
```

6. **Separating files and listing SRR numbers:**

```bash
ls *_1.fastq | sed 's/_1.fastq//' > downloaded.txt
```

*This creates a list of SRR numbers for all `_1.fastq` files.*
