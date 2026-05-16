# Coding Exercises — Detailed Plan

Per-lecture Google Colab exercises (≈ 60 minutes each) for the 27-lecture bioinformatics course. Students are EE undergrads with minimal biology background; every exercise lands one signal-processing / linear-algebra / statistics aha.

## Cross-cutting design constraints

- **Time budget:** 60 min for a motivated student who has watched the lecture.
- **Compute:** runs on free Colab CPU; ≤ 5 min total compute per exercise; no GPU required.
- **Data sourcing — synthetic or public only.** Every example file must be either (a) generated in-notebook from a deterministic seed, or (b) fetched at runtime from a stable public URL (NCBI / UniProt / Ensembl / ChEMBL / ClinVar / gnomAD / GIAB / Scanpy or squidpy built-ins / HuggingFace public models / Pfam / JASPAR / COSMIC public release). No curator-bundled private datasets ship with the repo. If a file appears to be "provided" or "precomputed", the spec must say how it's generated or fetched. Total payload ≤ 100 MB.
- **Style:** narrative cells explain *why*, code cells have `# TODO` for the student. Last cell prints / asserts the expected answer.
- **Don't duplicate artifacts:** the interactive HTMLs already show the algorithm visually; the exercise *implements* it and measures something.

## Lecture-numbering note

The actual `<h1>` titles in `lectures/lecture-NN.html` differ from the table in the project `README.md` for L05–L11. The exercise specs below use the **actual** titles (verified against the lecture HTML). The README table needs a separate refresh; outside the scope of this plan but tracked as a follow-up.

| # | Actual title in lecture HTML |
|---|---|
| 05 | Bulk RNA-seq: From Reads to Transcript Abundances |
| 06 | Differential Expression and Count Statistics |
| 07 | Single-Cell RNA-seq Fundamentals |
| 08 | Advanced Single-Cell: Trajectories, Integration, Multi-Modal |
| 09 | ChIP-seq, ATAC-seq, and Peak Calling |
| 10 | Methylation, Hi-C, and 3D Genome Organisation |
| 11 | Long Reads and the Pangenome |

## Decisions taken

- **Every one of the 27 lectures gets an exercise.** Even L14 (Data Engineering) — recast as workflow-DAG + GIAB-style benchmarking instead of a full Snakemake pipeline (too slow on Colab).
- **One notebook per lecture, hidden solution cells.** Single `lecture-NN/exercise.ipynb` with each solution cell collapsed using Colab's `cellView: "form"` / `jupyter.source_hidden` metadata. Students see the `# TODO` cell by default; if they need the answer they expand the hidden cell beneath it. No separate `_solution.ipynb`.
- **Colab link at the top of each lecture page.** Add an "Open in Colab" link to the masthead of `lectures/lecture-NN.html` (alongside or just below the existing artifact links), pointing at `colab.research.google.com/github/.../coding_exercises/lecture-NN/exercise.ipynb`. Same convention for every lecture.
- **No difficulty or extension flags surfaced to students.** The exercises are ungraded — students self-pace. The Easy / Medium / Hard tags in the per-lecture sections and in the "Difficulty distribution" table below are internal planning notes only.
- **Directory layout:** `coding_exercises/lecture-NN/` per lecture, with the notebook plus any small data files.
- **Data hosting:** ≤ 10 MB embedded in the notebook (base64 or Python-generated). > 10 MB → small public URL (E. coli reference, Scanpy datasets, gnomAD subset).
- **Standard preamble cell:** pinned `!pip install` block at top: `numpy pandas scipy scikit-learn matplotlib biopython networkx scanpy rdkit-pypi`. Heavier deps (`squidpy`, `transformers`) only on the lectures that need them.
- **Self-check:** every exercise ends with an `assert`-style cell that flags the right answer (e.g. `assert 0.40 <= recovered_ne <= 0.60`).

---

## L01 — Foundations: From Cells to Sequences to FASTQ

- **Exercise:** Parse and QC a synthetic FASTQ; fingerprint the sequencer error profile.
- **Dataset:** 10,000 reads × 150 bp generated in ~3 lines of NumPy with position-dependent Phred decay. ~2 MB in memory; nothing to download.
- **5-step outline (60 min):**
  1. (8) Generate synthetic FASTQ; write to a string.
  2. (10) Parse the 4-line records; decode Phred+33 → per-base error probability.
  3. (12) Histogram of mean quality per read; compare to a theoretical Poisson.
  4. (15) Position-wise quality curve over the 150 cycles.
  5. (15) Quality-trim at Q10; report surviving reads + new median length.
- **Aha:** FASTQ is a lossy compression of trace intensity into symbols + per-base confidence. Quality decay is a fingerprint of the instrument.
- **EE framing:** Phred = SNR in decibels (`Q = −10 log₁₀ P_err`). Trimming = preprocessing a noisy signal before correlation.
- **Avoid duplicating:** Artifacts 04 illumina-basecaller, 05 nanopore-squiggle, 06 fastq-inspector already visualise; this is the student's *implementation*.
- **Difficulty:** Easy.

## L02 — Read Alignment: From Brute Force to FM Index and Back

- **Exercise:** Build a suffix array + BWT + FM index on E. coli; do backward search for 100 simulated reads; measure scaling vs naive scan.
- **Dataset:** E. coli K-12 MG1655 reference (~4.6 MB FASTA, single download). 100 50-bp reads sampled from the reference + 0–2 random mismatches.
- **5-step outline (60 min):**
  1. (8) Load FASTA; build suffix array via naive O(N log N) sort.
  2. (8) Derive BWT + LF mapping from the SA; verify with a hand example.
  3. (12) Implement backward search: binary-search the SA range for each query.
  4. (12) Time 100 exact-match queries; compare to brute-force scan; report speedup.
  5. (20) Vary query length (30 / 50 / 100 bp); plot time vs L; overlay theoretical O(L log N).
- **Aha:** Preprocess once (O(N log N)), query forever (O(L log N)). FM index trades RAM for query speed.
- **EE framing:** Binary search on the SA = successive-approximation lookup (SAR-ADC analogy). BWT = data-dependent transform that concentrates like symbols → low entropy → fast backward search.
- **Avoid duplicating:** Artifacts 02 suffix-array and 03 bwt-builder show structure; this exercise adds the *scaling* measurement.
- **Difficulty:** Medium.

## L03 — DNA Sequence Assembly: From Reads to Genomes

- **Exercise:** Build a de Bruijn graph (k=25) from simulated reads; extract contigs; measure N50 across k.
- **Dataset:** 10 kb random ACGT genome; 500 simulated 50-bp reads at ~100× coverage; 1% error rate on a fraction. All generated in-notebook.
- **5-step outline (60 min):**
  1. (8) Count 25-mers, drop singletons.
  2. (10) Build directed graph (nodes = (k-1)-mers, edges = k-mers).
  3. (12) Greedy walk to extract linear contigs; halt at branch points.
  4. (15) Naive alignment of contigs back to the simulated genome; report N50 + coverage fraction.
  5. (15) Sweep k ∈ {15, 25, 35}; plot N50 vs k.
- **Aha:** Assembly = overlap problem recast as graph traversal. Low k → long contigs but more branches; high k → cleaner but gaps. N50 reveals the sweet spot.
- **EE framing:** k-mer counting = hash-table sketch. Contig extraction = path-finding on a DAG, analogous to maximum-likelihood decoding through a trellis.
- **Avoid duplicating:** Artifact 03 debruijn-builder lets students hand-build a small graph; this scales to 500 reads with quality metrics.
- **Difficulty:** Medium-hard.

## L04 — Variant Calling: From Aligned Reads to Called Differences

- **Exercise:** Simulate a 2 kb diploid region with a heterozygous SNV; build a pileup; compute genotype likelihoods; call the variant; sweep error rate.
- **Dataset:** Generated in-notebook. 500 simulated 100 bp reads at 50× depth; truth = het A/G at position 1000.
- **5-step outline (60 min):**
  1. (8) Generate ref + diploid ground truth + reads (with quality-based errors).
  2. (10) Build pileup; report base composition at the variant position.
  3. (12) Compute log-likelihoods for genotypes AA / AG / GG using per-read Phred.
  4. (15) Bayes-rule normalise with uniform priors → MAP genotype.
  5. (15) Sweep error rate ∈ {0, 1, 5} %; plot genotype-accuracy vs error.
- **Aha:** Calling is Bayesian MAP inference, not thresholding. The pileup is a sufficient statistic.
- **EE framing:** Argmax of posterior = signal-detection problem with unequal priors. Pileup = data-compression to a sufficient statistic. Log-likelihood = inner product of read base × log-emission.
- **Avoid duplicating:** Artifact 03 genotype-likelihoods is hand-computed examples; this scales it and quantifies robustness to noise.
- **Difficulty:** Medium.

## L05 — Bulk RNA-seq: From Reads to Transcript Abundances

- **Exercise:** Build a k-mer compatibility-class index for 3 overlapping isoforms; run EM to recover ground-truth abundances.
- **Dataset:** Three ~2 kb simulated isoforms sharing exons; 500 50-bp reads sampled with ground-truth weights (40 / 35 / 25 %).
- **5-step outline (60 min):**
  1. (8) Build k-mer → isoform index (k=25); assign each read to its compatibility class.
  2. (10) Initialise EM with uniform abundances.
  3. (15) E-step: fractional read assignment proportional to current abundances.
  4. (12) M-step: re-estimate abundances; repeat for 10 iterations.
  5. (15) Plot convergence; compare to ground truth; perturb ground truth and re-run.
- **Aha:** Ambiguous reads are not garbage — soft EM assignments recover the truth in ~5 iterations.
- **EE framing:** EM = iterative soft-decision decoding (the same idea behind belief propagation, turbo codes, Kalman smoothers).
- **Avoid duplicating:** Artifact 04 em-visualizer animates the dynamics; this exercise measures convergence and accuracy.
- **Difficulty:** Medium.

## L06 — Differential Expression and Count Statistics

- **Exercise:** Fit a negative-binomial GLM on a 500-gene × 6-sample count matrix; apply empirical-Bayes dispersion shrinkage; control FDR; produce a volcano plot.
- **Dataset:** Simulated 500 genes × 6 samples (3 vs 3). 50 genes are true positives with log₂ FC ∈ {±1, ±2}; rest are null. Counts ~ NB(μ_g · sample factor, dispersion_g) with dispersion sampled from a realistic prior. ~50 KB in memory.
- **5-step outline (60 min):**
  1. (10) Generate the count matrix; estimate library-size factors (median-of-ratios).
  2. (12) Fit per-gene NB GLM (statsmodels.GLM with NB family); extract MLE dispersions; plot dispersion vs mean expression.
  3. (12) Apply EB shrinkage: regress log-dispersion on log-mean; shrink each gene's dispersion toward the trend by a precision-weighted average.
  4. (15) Wald test for the condition coefficient; p-value per gene; BH FDR; produce a volcano plot (log₂ FC vs −log₁₀ p-adj).
  5. (11) Compare TPR + FDR at different shrinkage strengths (no shrinkage / EB / strong-prior) to show why DESeq2-style EB beats naive per-gene MLE.
- **Aha:** With 6 samples, the per-gene dispersion MLE is hopelessly noisy. EB shrinkage borrows strength across genes → calibrated p-values + actually-controlled FDR.
- **EE framing:** Per-gene dispersion estimate = noisy measurement. EB shrinkage = optimal Bayesian estimator under a learned prior — the same idea as James-Stein shrinkage, Wiener filtering, or Tikhonov regularisation.
- **Avoid duplicating:** Artifact volcano-plot explorer (if present) visualises the result; this exercise *fits* the GLM and shows why shrinkage matters.
- **Difficulty:** Medium-hard. The hardest single exercise in the first half of the course because NB GLM + dispersion estimation packs a lot.

## L07 — Single-Cell RNA-seq Fundamentals

- **Exercise:** QC → HVG → PCA → UMAP → Leiden → marker genes on a 3k-cell PBMC dataset.
- **Dataset:** Scanpy built-in `sc.datasets.pbmc3k()` (~5 MB).
- **5-step outline (60 min):**
  1. (8) Load, inspect shape, plot QC metrics (n_counts, n_genes, pct_mt).
  2. (10) Drop low-quality cells (mt > 15 %, < 500 counts, > 6000 counts).
  3. (15) Log-normalise; select 2 000 HVGs by dispersion.
  4. (18) PCA → 30 PCs → UMAP → Leiden (res 0.5); colour-by-cluster scatter.
  5. (9) Wilcoxon marker test for one cluster; heatmap of top 5 markers.
- **Aha:** Cells are sparse vectors in 20k-dim gene space; PCA + UMAP make discrete cell types fall out as clusters.
- **EE framing:** PCA = optimal rank-k approximation. UMAP = nonlinear manifold embedding. Leiden = graph partition by modularity (spectral clustering in disguise).
- **Avoid duplicating:** Artifacts already visualise the pipeline; this exercise *runs* it end-to-end.
- **Difficulty:** Easy.

## L08 — Advanced Single-Cell: Trajectories, Integration, Multi-Modal

- **Exercise:** Compute pseudotime on a myeloid trajectory; integrate two batches with Harmony; quantify integration quality.
- **Dataset:** 5k-cell Paul et al. 2015 myeloid trajectory (Scanpy built-in or GEO download ≤ 50 MB) + a 4k-cell two-batch PBMC subset.
- **5-step outline (60 min):**
  1. (8) Standard scRNA pipeline on the trajectory dataset.
  2. (12) Slingshot / Monocle3 pseudotime; UMAP coloured by pseudotime should show a smooth gradient.
  3. (12) Concat two batches; show UMAP separates by batch.
  4. (18) Run Harmony; re-UMAP; quantify integration via batch-silhouette delta.
  5. (10) Light-touch CellPhoneDB on one cluster pair — interpret as hypothesis, not ground truth.
- **Aha:** Pseudotime linearises a snapshot. Harmony aligns batches without seeing labels. L-R inference is hypothesis-generation.
- **EE framing:** Pseudotime = 1-D coordinate on a manifold. RNA velocity = state-space estimation (spliced / unspliced as two measurements of the same hidden state). Harmony = domain adaptation by iterative scaling.
- **Avoid duplicating:** Artifacts cover the algorithms; here students score integration quality with a metric.
- **Difficulty:** Medium-hard.

## L09 — ChIP-seq, ATAC-seq, and Peak Calling

- **Exercise:** Peak-call on an ATAC bedgraph with a local-Poisson test; replicate concordance; differential accessibility NB GLM; PWM motif scan.
- **Dataset:** Precomputed 1 Mb ATAC coverage bedgraph for K562 + a treated replicate (or simulated equivalent). ≤ 10 MB.
- **5-step outline (60 min):**
  1. (6) Load bedgraph; plot signal/background histogram; explain peak = local pile-up.
  2. (12) Call peaks with a local-Poisson test at q ∈ {0.01, 0.05, 0.1}; tabulate.
  3. (12) Replicate overlap: how many peaks survive at 500 bp window?
  4. (18) NB GLM on peak counts (callback to L06) — control vs treated; volcano plot.
  5. (12) PWM scan for one JASPAR motif on diff-accessible peaks; enrichment vs background.
- **Aha:** Peak calling = CFAR detection (local-noise-adaptive threshold). Differential accessibility = the same NB GLM as RNA-seq, just on chromatin.
- **EE framing:** MACS2 = local-Poisson hypothesis test (CFAR). PWM scan = matched filter (inner product of sequence ↔ PWM). NB GLM = L06 callback.
- **Avoid duplicating:** Footprint/fragment-length artifacts are QC; this exercise calls peaks + does diff-binding.
- **Difficulty:** Medium.

## L10 — Methylation, Hi-C, and 3D Genome Organisation

- **Exercise:** Beta posterior on methylation proportions + ICE-normalise a Hi-C matrix + insulation-score TADs + A/B compartments by first eigenvector.
- **Dataset:** Tiny methylation table (~100 CpGs × 5 samples) + a 40 kb-binned 50×50 Hi-C contact matrix for a 2 Mb region. Both ≤ 5 MB.
- **5-step outline (60 min):**
  1. (8) Visualise the methylation landscape; spot a CpG island.
  2. (10) Fit Beta posterior to one sample's methylation proportions.
  3. (12) ICE-normalise the contact matrix; before / after comparison.
  4. (18) Compute insulation score → detect TAD boundaries as change-points.
  5. (12) PCA on contact matrix; sign of first eigenvector = A/B compartments.
- **Aha:** Methylation is a Beta-distributed proportion. Hi-C is a covariance-like structure; TADs = blocks; A/B compartments fall out of the leading eigenvector.
- **EE framing:** Beta = conjugate Bayesian inference on proportions. ICE = iterative-bias scaling. Insulation score = edge detection (1-D filtering). Compartments = PC1 of the contact matrix.
- **Avoid duplicating:** Bisulfite-walkthrough artifact covers decoding; this exercise focuses on Beta inference + Hi-C linear algebra.
- **Difficulty:** Hard.

## L11 — Long Reads and the Pangenome

- **Exercise:** Detect a 5 kb deletion from HiFi coverage + align reads to a toy pangenome GFA via Viterbi-on-DAG + phase reads into two haplotypes.
- **Dataset:** Simulated 100 kb region + 50 HiFi reads spanning a known 5 kb deletion + a 5-node pangenome GFA. ~5 MB.
- **5-step outline (60 min):**
  1. (8) Load reads; map to linear reference; spot the SV signature in coverage.
  2. (10) Load the GFA; visualise the graph (linear backbone + variant branch).
  3. (12) Implement minimizer seeding; locate seeds across nodes.
  4. (15) Generalise the L02 Viterbi to a topologically-sorted DAG; pick the best path per read.
  5. (15) Phase reads by allele support at the deletion locus — two source sets.
- **Aha:** Pangenome = DAG over the linear reference. Long reads phase directly by spanning the SV. Viterbi generalises from a chain to a DAG.
- **EE framing:** Pangenome alignment = Viterbi on a DAG (callback to L02). Minimizer seeding = locality-sensitive hashing / sketching. Phasing = source separation.
- **Avoid duplicating:** Linear FM-index covered in L02; this is the graph generalisation.
- **Difficulty:** Hard.

## L12 — Population Genetics Fundamentals

- **Exercise:** Wright-Fisher drift simulation; estimate effective population size from observed variance; add selection.
- **Dataset:** Generated in-notebook: N = 1000, T = 200 generations, 100 replicates.
- **5-step outline (60 min):**
  1. (8) Initialise allele frequencies; binomial sampling per generation.
  2. (15) Run 100 replicates of neutral drift; plot 5 trajectories.
  3. (12) Estimate Nₑ from `var(Δp · t / (p₀(1−p₀)))`; scatter true vs estimated.
  4. (15) Add a selection coefficient s = 0.01; overlay drift trajectories.
  5. (10) Output table: (generation, mean p, var p) for neutral; selection summary.
- **Aha:** Drift variance grows as √(t / N). Weak selection is just a directional bias on top of stochastic noise.
- **EE framing:** Drift = 1-D random walk with √-time variance. Selection = deterministic trend superposed on stochastic noise. SNR ∝ Nₑ · s.
- **Avoid duplicating:** Wright-Fisher artifact simulates; this exercise *estimates* Nₑ from observed variance.
- **Difficulty:** Easy-medium.

## L13 — GWAS and Statistical Genetics

- **Exercise:** Simulate a GWAS with population stratification; produce Manhattan + QQ plots; show λ_GC distinguishes confounding from true polygenicity.
- **Dataset:** Simulated 10k individuals × 50k SNPs with 100 causal SNPs + a 0.15-correlated ancestry axis. ≤ 10 MB.
- **5-step outline (60 min):**
  1. (10) Generate genotype matrix + ancestry labels + causal effects.
  2. (12) Per-SNP logistic regression; z-scores; p-values.
  3. (15) Manhattan plot; annotate 5×10⁻⁸ threshold.
  4. (15) QQ plot + λ_GC. Re-run *without* stratification; overlay both QQs.
  5. (8) Top-10 SNP table with causal flag; λ_GC stratified vs unstratified.
- **Aha:** Stratification lifts the whole QQ tail; polygenicity only the upper tail. Shape, not just λ_GC, tells the story.
- **EE framing:** GWAS = 50 000 parallel hypothesis tests (multi-channel detection). QQ plot = inverse-CDF probability transform. PCA correction = noise whitening.
- **Avoid duplicating:** Manhattan-plot artifact replays a fixed scan; here students simulate + observe how confounding changes the plot.
- **Difficulty:** Medium.

## L14 — Data Engineering, File Formats, and Reproducibility

- **Exercise:** Workflow-DAG resumption logic + GIAB-style VCF benchmarking. *Recast from full Snakemake to a precision/recall + cache-key exercise — Snakemake on Colab is too slow.*
- **Dataset:** Two VCFs generated in-notebook with NumPy: a 500-row "predicted" set + a 300-row "truth" set, seeded for reproducibility. (Optional extension cell: download a small public GIAB v4.2.1 high-confidence VCF slice for chr20 and re-run the matching logic against real-world data.) < 1 MB.
- **5-step outline (60 min):**
  1. (8) Load both VCFs; inspect schema.
  2. (12) Implement TP/FP/FN matching by (CHROM, POS, REF, ALT); confusion matrix.
  3. (15) Precision / recall / F1; ROC curve by varying QUAL threshold.
  4. (15) Simulate a DAG of tasks (load_truth → load_pred → match → metrics); manually fail one; show which downstream tasks must re-run; explain content-addressable cache keys.
  5. (10) Markdown summary: when resume matters (long-running tasks, large intermediates).
- **Aha:** Workflow managers = content-hash caches + schedulers. Reproducibility comes from declarative dependency graphs.
- **EE framing:** VCF parsing = tokenisation. Matching = nearest-neighbour lookup. Precision/recall = detection metrics. Hashing = integrity checksum.
- **Avoid duplicating:** Workflow-DAG visualiser and GIAB benchmarker artifacts cover the concepts; this exercise *computes* the metrics + DAG logic.
- **Difficulty:** Easy-medium.

## L15 — Protein Structure Prediction in the AlphaFold Era

- **Exercise:** DCA contact prediction from a real MSA, vs ESMFold contact map (public CPU-friendly transformer) on the same sequence. Show that the learned transformer absorbs DCA + captures more.
- **Dataset:** Pfam PF00042 (globin) full MSA fetched from InterPro / Pfam REST (`https://www.ebi.ac.uk/interpro/api/entry/pfam/PF00042/?fields=alignment`); ~200 sequences × 150 columns. Reference protein: human haemoglobin α (UniProt P69905) via UniProt REST. ESMFold-150M loaded from HuggingFace `facebook/esmfold_v1` (skip if too heavy on free Colab; alternative is the precomputed AlphaFold prediction for P69905 from the public AlphaFold DB at `https://alphafold.ebi.ac.uk/files/AF-P69905-F1-predicted_aligned_error_v4.json`). All public.
- **5-step outline (60 min):**
  1. (10) Fetch the Pfam MSA + reference sequence; compute per-column entropy + conservation heatmap.
  2. (10) DCA: empirical covariance → inverse covariance → top contacts → contact heatmap.
  3. (15) Run ESMFold inference on the reference (or load the AlphaFold-DB PAE JSON for the same UniProt accession); extract pair logits / contact probabilities.
  4. (15) Correlate DCA scores with the transformer's predicted contacts; scatter plot + Pearson correlation.
  5. (10) Markdown explanation: why transformers absorb DCA + capture higher-order correlations.
- **Aha:** DCA = inverse-covariance contact estimator. Modern transformers learn a richer pair representation by iterated axial attention.
- **EE framing:** Inverse covariance = precision matrix. Axial attention = separable convolution on a 2-D matrix. DCA → ESMFold/AlphaFold = hand-crafted kernel → learned kernel.
- **Avoid duplicating:** Coevolution / Evoformer artifacts visualise; here students *compute* DCA and compare numerically against a public learned predictor.
- **Difficulty:** Medium.

## L16 — ML in Genomics: Architectures, Pitfalls, Frontiers

- **Exercise:** Pretrained-model inference + embedding probing. Load DNABERT and ESM2-35M; embed sequences; ablate motifs; show foundation models capture local features but aren't zero-shot predictors.
- **Dataset:** 10 DNA + 10 protein sequences, all from public sources fetched in-notebook — 5 human promoter sequences from Ensembl REST (`https://rest.ensembl.org/sequence/region/human/...`), 5 random/intergenic controls generated with NumPy, 5 protein sequences from UniProt REST (haemoglobin α, lysozyme, insulin, etc.), 5 disorder-prone controls from the IDP DisProt public API. Model weights pulled from public HuggingFace endpoints (`zhihan1996/DNABERT-2-117M` for DNA; `facebook/esm2_t6_8M_UR50D` for protein — both small enough for CPU). ~50 MB.
- **5-step outline (60 min):**
  1. (10) Load DNABERT; embed a 500 bp promoter; report shape + summary stats.
  2. (12) Embed 10 DNA sequences (coding / regulatory / intergenic); pairwise cosine similarity heatmap.
  3. (15) Same for 10 protein sequences with ESM2-35M.
  4. (15) Ablation probe: remove a 6 bp motif; measure embedding distance from original.
  5. (8) Markdown: "Foundation models are feature extractors, not off-the-shelf predictors."
- **Aha:** Embeddings respect local sequence features (ablation shifts distance) but zero-shot inference ≠ fine-tuned predictor.
- **EE framing:** Embedding = learned dimensionality reduction. Cosine similarity = metric in latent space. Ablation = sensitivity analysis.
- **Avoid duplicating:** DNA-LM tokeniser artifact shows the tokeniser; this exercise runs *inference* + probing.
- **Difficulty:** Medium-hard.

## L17 — Clinical Genomics, Variant Interpretation, and Ethics

- **Exercise:** Rule-based ACMG/AMP classifier + ancestry-bias demo on PRS.
- **Dataset:** 25 variants fetched at runtime from public ClinVar (`https://www.ncbi.nlm.nih.gov/clinvar/api/...` or the FTP variant_summary.txt) cross-joined with gnomAD ancestry-stratified AFs from the public gnomAD REST (`https://gnomad.broadinstitute.org/api`) and REVEL scores from the public download. The exercise filters to a curated subset of 25 known anchor variants (BRCA1, BRCA2, TP53, CFTR, etc.) in-notebook with deterministic IDs — no private TSV ships. 100 SNPs for the PRS demo sampled the same way. < 5 MB total fetch.
- **5-step outline (60 min):**
  1. (8) Load TSV; inspect ACMG evidence fields for one variant.
  2. (12) Assign ACMG codes (PVS1, PS1–4, PM1–6, …) to two anchor variants; walk the rule-combination engine.
  3. (15) Run the classifier on all 25 variants; tally classifications vs ClinVar.
  4. (15) Recompute 5 variants using EUR-only vs full-gnomAD frequencies; flag classification flips; PRS demo from 100 SNPs across EUR vs AFR backgrounds.
  5. (10) Markdown: which codes are ancestry-sensitive? Why labs must report ancestry-stratified frequencies.
- **Aha:** ACMG is auditable evidence assembly, not a black-box algorithm. Ancestry bias hides pathogenic variants in under-represented populations.
- **EE framing:** ACMG rules = Boolean logic gate network. Evidence codes = sensor inputs. PRS-portability gap = covariate shift.
- **Avoid duplicating:** ACMG classifier + ancestry-bias artifacts are interactive; here students *run* the rules + compute PRS bias numerically.
- **Difficulty:** Easy. Aha is conceptual, not algorithmic.

## L18 — Cancer Genomics: Integrated Capstone

- **Exercise:** NMF on a 96-trinucleotide signature matrix; recover SBS1/4/7a/3 from a simulated 50-sample cohort; model-order selection.
- **Dataset:** Synthetic 50 × 96 spectrum matrix, mixed from COSMIC SBS catalogue v3.4 + Poisson noise. ~10 KB.
- **5-step outline (60 min):**
  1. (8) Load matrix; barplot one sample; label dominant trinucleotide contexts.
  2. (12) Implement Lee–Seung multiplicative-update NMF with K = 4.
  3. (15) Run 50 iterations; cosine-similarity recovered W columns vs COSMIC.
  4. (15) Sweep K ∈ {2, 3, 5}; reconstruction error + over-fit at K = 5.
  5. (10) Interpret: HRD signature (SBS3) → PARP-inhibitor candidate; UV (SBS7a) → melanoma confirmation.
- **Aha:** Signatures = non-negative blind source separation. Wrong K invents spurious sources.
- **EE framing:** NMF = matrix factorisation under positivity constraint; multiplicative updates ⊂ projected gradient descent. Model-order selection = same problem as PCA-k or filter-order selection.
- **Avoid duplicating:** Signature-decomposer artifact is interactive; here students *implement* the algorithm + run model-order sweep.
- **Difficulty:** Medium.

## L19 — BLAST and Sequence Search Statistics

- **Exercise:** Seed-and-extend ungapped BLAST + Karlin-Altschul E-values; compare against published NCBI hits.
- **Dataset:** 5 query proteins + a 1 MB UniProt subset (~20k sequences).
- **5-step outline (60 min):**
  1. (10) Load DB; BLOSUM62 lookup; hand-compute one alignment.
  2. (12) Enumerate query 3-mer neighbourhoods at T = 11 vs T = 5.
  3. (12) Scan DB for seeds + X-drop ungapped extension; top-10 HSPs per query.
  4. (15) Bit score from `(λ·S − ln K)/ln 2`; E-value = `m·n·2⁻ᴮ`.
  5. (11) Run real BLASTP on one query; compare top-5 hits.
- **Aha:** BLAST = matched-filter cascade. Bit score is a log-likelihood ratio; E-value is its false-alarm rate.
- **EE framing:** Cascade detector — cheap k-mer match → expensive BLOSUM scoring → log-odds test statistic. T = sensitivity-speed knob.
- **Avoid duplicating:** Seed-neighbourhood and E-value-calculator artifacts visualise; this exercise *implements* the full pipeline + benchmarks against NCBI.
- **Difficulty:** Medium-hard.

## L20 — MSA, Phylogenetics, Comparative Genomics

- **Exercise:** Progressive MSA + NJ tree on 8 globins + dN/dS on a codon-aligned gene.
- **Dataset:** 8 globin proteins (~150 aa each) + Pfam PF00042 reference MSA + a 10-species ~300-codon gene. ≤ 5 KB.
- **5-step outline (60 min):**
  1. (10) 28 pairwise distances → NJ tree; visualise with branch lengths.
  2. (12) Progressive MSA: pairwise align, build profile, add sequence-vs-profile until 8-way.
  3. (12) SP score vs Pfam reference; conservation per column.
  4. (15) Codon-align the 10-species gene; Nei-Gojobori dN/dS; LRT vs neutrality.
  5. (11) Interpret tree topology; flag any branch with dN/dS > 1 as positive selection.
- **Aha:** MSA is NP-hard; progressive is a greedy heuristic that depends on a good guide tree. Trees summarise distance, not "truth". dN/dS is an LRT on the codon-substitution model.
- **EE framing:** MSA = constrained DP on a tree-induced order. Tree = factor graph encoding conditional independence. dN/dS = LRT (χ²-distributed null).
- **Avoid duplicating:** Progressive-MSA, NJ-vs-ML, and MSA-quality artifacts cover the algorithms; this exercise implements them + computes dN/dS.
- **Difficulty:** Hard.

## L21 — HMMs, Profile HMMs, and Gene Finding

- **Exercise:** Viterbi decoding on a 5-state gene HMM + profile-HMM scoring (HMMER-style) on a small test set.
- **Dataset:** Synthetic 5-state HMM parameters hardcoded in the notebook as a Python dict; 2 kb BRCA1 region fetched from Ensembl REST (`https://rest.ensembl.org/sequence/region/human/17:43044295..43046294`); Pfam PF00042 profile HMM downloaded from the public Pfam FTP (`https://www.ebi.ac.uk/interpro/api/entry/pfam/PF00042/?annotation=hmm`); 50 test sequences sampled from UniProt — 25 globin hits via the UniProt REST search filter `family:globin` and 25 random non-globin proteins. ≤ 200 KB.
- **5-step outline (60 min):**
  1. (10) Load HMM; describe emissions; initialise the Viterbi table.
  2. (12) Forward pass + backtrack; visualise predicted state path.
  3. (12) Compare to RefSeq exon annotation; sensitivity / specificity.
  4. (15) Score 50 sequences against the profile HMM; threshold at E = 0.01; ROC.
  5. (11) Discuss why HMMs still beat DL on sparse-data segmentation.
- **Aha:** Viterbi = MAP path on a state-space model; the most probable path ≠ sum of per-position marginals.
- **EE framing:** HMM = discrete-state-space model (Kalman analogue). Profile HMM = position-specific filter bank. Forward / Viterbi = belief propagation.
- **Avoid duplicating:** Viterbi-stepper + profile-HMM artifacts walk through; this exercise *implements* both algorithms in log-space.
- **Difficulty:** Hard.

## L22 — Network Biology and Pathway Analysis

- **Exercise:** Random-walk-with-restart on a 500-node STRING-like PPI + Louvain community detection + hypergeometric pathway enrichment.
- **Dataset:** Synthetic 500-node graph with planted communities, generated by NetworkX (~30 KB).
- **5-step outline (60 min):**
  1. (8) Build the graph; validate connectivity + degree distribution.
  2. (12) RWR from a 20-gene seed; sweep restart probability ∈ {0.3, 0.5, 0.7}.
  3. (15) Louvain; modularity Q vs ground-truth labels.
  4. (15) Hypergeometric pathway-enrichment on RWR's top 30.
  5. (10) Interpret RWR as Laplacian low-pass filter; restart probability = cutoff frequency.
- **Aha:** Network propagation = low-pass filter on the graph Laplacian. Louvain = modularity max ≈ spectral clustering.
- **EE framing:** RWR = heat diffusion on a resistor network. Louvain = k-means in Laplacian eigenspace. Hypergeometric = 2-sample binomial.
- **Avoid duplicating:** RWR-walker + spectral-cluster artifacts are toy; here students *scale to 500 nodes* and add enrichment.
- **Difficulty:** Medium.

## L23 — Metagenomics and the Microbiome

- **Exercise:** DADA2-style ASV inference + alpha / beta diversity + permutation test on a simulated 16S dataset.
- **Dataset:** 1 000 simulated 16S reads × 6 samples (3 case + 3 control), 10 ground-truth taxa, 2 % per-position error. ~150 KB.
- **5-step outline (60 min):**
  1. (8) Parse reads; length filter; drop singletons.
  2. (12) ASV inference: edit-distance clustering + replicate-based error model.
  3. (12) k-mer taxonomy assignment to a 10-reference DB; 6 × 10 abundance matrix.
  4. (15) Alpha (Shannon, Simpson, richness) + beta (Bray-Curtis, Jaccard); PCoA.
  5. (13) Paired t-test on Shannon; ANOSIM on BC matrix; interpret compositional bias.
- **Aha:** Sequencing depth is a population census. Bray-Curtis (abundance) vs Jaccard (presence/absence) tell different stories. Compositional data lives on a simplex.
- **EE framing:** ASV inference = error-correction coding on biological strings (edit distance = Hamming in discrete space). Alpha diversity = entropy. PCoA = eigen-decomposition of the distance matrix.
- **Avoid duplicating:** Diversity-calc artifact handles 5 taxa; here students go to 10 taxa + statistical tests.
- **Difficulty:** Medium-hard.

## L24 — CRISPR Functional Screens and DepMap

- **Exercise:** MAGeCK-style negative-binomial test + RRA gene-level aggregation on a simulated 10k sgRNA × 6 sample screen.
- **Dataset:** Simulated 10 000 sgRNAs × 6 samples; 50 essential genes (4-fold depletion) + 20 resistant (2-fold gain) + 430 nulls. ~100 KB.
- **5-step outline (60 min):**
  1. (8) Generate counts (NB with overdispersion + batch effect); normalise.
  2. (12) Per-sgRNA NB test; LFC + p-value; BH FDR.
  3. (14) Gene-level RRA aggregation; volcano plot of gene LFC vs −log₁₀ FDR.
  4. (14) sgRNA concordance per gene; expected vs observed false-discovery rate.
  5. (12) Interpret screen as compressed sensing — sparse signal (50 hits) from 10 000 measurements.
- **Aha:** Rank aggregation turns noisy per-sgRNA measurements into reliable gene hits. Multiple sgRNAs per gene = built-in replication.
- **EE framing:** Screen = sparse-signal recovery from noisy measurements (compressed sensing). NB test = matched-filter in the count domain. RRA = rank-based voting.
- **Avoid duplicating:** MAGeCK artifact handles 50-gene toys; here students scale to 10k sgRNAs + concordance checks.
- **Difficulty:** Medium.

## L25 — Causal Inference and Mendelian Randomisation

- **Exercise:** Two-sample MR with 100 SNP instruments — implement IVW, MR-Egger, weighted median; identify outliers; sensitivity triangulation.
- **Dataset:** Simulated 100 SNPs; 80 valid + 20 pleiotropic; true causal effect β = 0.5. ~5 KB.
- **5-step outline (60 min):**
  1. (10) Generate (β_ZX, β_ZY); visualise the Wald scatter.
  2. (12) IVW estimator + CI.
  3. (13) MR-Egger with free intercept; flag outliers (high Cook's distance).
  4. (12) Weighted median; compare all three estimators.
  5. (13) Triangulation: when do methods agree? Steiger directionality flag.
- **Aha:** MR = orthogonal regression under instrument validity. MR-Egger absorbs pleiotropy at the cost of power. Triangulating estimators reveals hidden pleiotropy.
- **EE framing:** IV = 2-stage least squares. IVW = weighted least squares with intercept fixed at 0. Egger = WLS with free intercept. Weighted median = L1 regression in rank space (robust to outliers).
- **Avoid duplicating:** MR-sensitivity artifact tries different estimators; this exercise *implements* them + adds outlier detection.
- **Difficulty:** Hard.

## L26 — Drug Discovery and Chemoinformatics

- **Exercise:** Morgan fingerprints + Tanimoto similarity matrix + hierarchical clustering + Lipinski / QED filter on 100 ChEMBL-style molecules.
- **Dataset:** 100 SMILES fetched at runtime from public ChEMBL — pulled via the ChEMBL REST API (`https://www.ebi.ac.uk/chembl/api/data/molecule.json?molecule_properties__mw_freebase__range=150,500&max_phase=4&limit=100`) so the student always gets a fresh, drug-like, MW-filtered sample. Fallback: a hard-coded list of 20 anchor drugs (aspirin, ibuprofen, caffeine, dabigatran, ezetimibe, etc.) generated by SMILES literals in the notebook if the network call fails. < 10 KB payload.
- **5-step outline (60 min):**
  1. (8) Parse SMILES with RDKit; compute Lipinski descriptors (MW, logP, HBD, HBA).
  2. (12) Morgan fingerprints (2048 bits, radius 2); pairwise Tanimoto similarity matrix.
  3. (12) Hierarchical clustering on 1 − Tanimoto; dendrogram + scaffold colouring.
  4. (15) Lipinski filter; QED score; scatter on (logP, MW) plane.
  5. (13) Interpret: Tanimoto = Jaccard on bit-set features; where Lipinski fails (dabigatran, ezetimibe).
- **Aha:** Similar fingerprints → similar properties → activity-cliff exceptions are rare but real. Lipinski is a cheap, useful filter — not a final answer.
- **EE framing:** Morgan FP = local feature extraction with bounded receptive field. Tanimoto = Jaccard on binary feature vectors (LSH-style). Hierarchical clustering = vector quantisation in feature space.
- **Avoid duplicating:** Tanimoto and Lipinski-QED artifacts are pairwise / single-molecule; here students do *pairwise on 100* + clustering.
- **Difficulty:** Easy-medium. RDKit does the heavy lifting.

## L27 — Mass-Spectrometry Proteomics + Metabolomics

- **Exercise:** Target-decoy FDR on 20 000 simulated PSMs + protein inference + label-free quantification + volcano plot.
- **Dataset:** Simulated 20 000 PSMs (10 000 target + 10 000 decoy) + 100 proteins × 6 samples count table. ~200 KB.
- **5-step outline (60 min):**
  1. (8) Parse the PSM list; explain target-decoy.
  2. (10) Sweep score threshold; compute FDR = #decoys / #targets at each threshold; ROC.
  3. (12) Protein inference: parsimony / Occam's razor on shared peptides; ≥ 2 unique-peptide filter.
  4. (15) Label-free quant: TopN intensity → median-of-ratios normalise → log₂ FC → t-test; volcano plot.
  5. (15) Discuss missing-value imputation + isotopic-label alternatives (SILAC, TMT).
- **Aha:** Decoy hits = direct empirical estimate of false positives. FDR is the best we can do without ground truth. LFQ is noisy; isotopic labelling reduces noise at higher cost.
- **EE framing:** PSM scoring = matched filter in the m/z domain. Target-decoy = null-hypothesis testing with an empirical null. FDR = q-value, the right metric when running 10⁵ parallel tests.
- **Avoid duplicating:** Target-decoy and diff-abundance artifacts cover pieces; this exercise stitches them into one pipeline.
- **Difficulty:** Medium-hard.

---

## Difficulty distribution

| Tier | Lectures |
|---|---|
| Easy | L01, L07, L17 |
| Easy-medium | L12, L14, L26 |
| Medium | L02, L04, L05, L09, L13, L15, L18, L22, L24 |
| Medium-hard | L06, L08, L16, L19, L23, L27 |
| Hard | L03, L10, L11, L20, L21, L25 |

Roughly normal-shaped — most exercises are medium / medium-hard, with a thin tail of easy + hard at the ends.

## EE-framing index

Every exercise lands one signal-processing / linear-algebra / statistics analogy. Highlights:

- **L01** Phred = SNR in dB.
- **L02** SA + binary search = SAR-ADC successive approximation.
- **L05** EM = iterative soft-decision decoding.
- **L06** EB shrinkage = James-Stein / Wiener filtering.
- **L09** Peak calling = CFAR detection.
- **L10** Hi-C compartments = PC1 of contact covariance.
- **L11** Graph alignment = Viterbi on a DAG (generalisation of L02).
- **L13** GWAS = parallel hypothesis tests; PC correction = noise whitening.
- **L18** NMF = blind source separation.
- **L19** BLAST = matched-filter cascade with LLR test statistic.
- **L21** Viterbi = MAP decoding on a state-space model.
- **L22** RWR = heat diffusion / spectral filtering.
- **L24** Screen = compressed sensing.
- **L25** MR = 2SLS / orthogonal regression.
- **L26** Tanimoto = Jaccard on LSH-style feature vectors.

## Next steps

1. **Stand up the directory tree** — one folder per lecture under `coding_exercises/lecture-NN/`.
2. **Build the notebook template** — single `exercise.ipynb` with:
   - Preamble cell (pinned `!pip install`).
   - Narrative markdown cell per step.
   - Visible `# TODO` cell per step.
   - Hidden solution cell beneath each TODO (`metadata.jupyter.source_hidden = true`).
   - Final self-check `assert` cell.
3. **Wire the "Open in Colab" link into each lecture page** — once we have a notebook URL pattern (probably the `github.com/.../coding_exercises/lecture-NN/exercise.ipynb` path through `colab.research.google.com/github/`), add the link to the masthead of every `lectures/lecture-NN.html` in a single batched edit.
4. **Start with 3 pilot lectures** — pick one Easy (L01), one Medium (L19 BLAST), one Hard (L21 HMM); shake out the template + the hidden-cells convention before batching the other 24.
5. **Refresh the README lecture table** — L05–L11 titles are stale (see the numbering note above).
