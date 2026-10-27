# HaploWISE
**H**aplotype **W**indow-based **I**nference of **S**election and **E**volution

HaploWISE is a reference-guided framework for genome-wide analyses of haplotype-specific genetic diversity and divergence. Designed for phased genomic data, HaploWISE combines local Principal Component Analysis (PCA) and discriminant analysis to classify haplotypes into alternative genetic classes using a reference population.

The tool focuses on highly differentiated genomic regions, such as polymorphic inversions and divergent haplotype blocks, where distinct haplotype classes can be reliably identified. Haplotypes are assigned to predefined classes within genomic windows, while ambiguous assignments are masked or excluded from downstream analyses.

Following classification, HaploWISE calculates population genetic statistics, including nucleotide diversity (π), Tajima's D, and absolute divergence (Dxy), within and between haplotype classes across the genome. By comparing haplotypes rather than relying exclusively on geographic population assignments, HaploWISE enables the investigation of genomic differentiation, diversity, and evolutionary processes associated with alternative haplotypes.

## Key Features

- **Phased genomic data:** Uses phased VCF files to distinguish alternative haplotypes within individuals.
- **Local PCA:** Uses window-based PCA to identify regions exhibiting distinct genetic clusters.
- **Reference-guided classification:** Applies discriminant analysis using reference individuals to assign haplotypes to predefined classes.
- **Uncertainty filtering:** Masks or excludes haplotypes that cannot be confidently classified.
- **Haplotype-based population genetics:** Calculates diversity and divergence statistics within and between haplotype classes.
- **Population-independent comparisons:** Enables comparisons of haplotype classes across geographic populations.
- **Applications to divergent genomic regions:** Particularly suited for studying polymorphic inversions and other highly differentiated genomic regions.
