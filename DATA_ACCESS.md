# Data Access

This project used summary-level genetic association data. Raw GWAS/pQTL summary-statistics files are not redistributed in the public release. Users should obtain the required inputs from the original resources under their terms of use.

## Primary Exposure

UK Biobank Pharma Proteomics Project (UKB-PPP) Olink Explore 3072 pQTL resource.

- UKB-PPP Synapse project: https://www.synapse.org/Synapse:syn51364943
- Registry of Open Data on AWS: https://registry.opendata.aws/ukbppp/
- Main publication: Sun BB, Chiou J, Traylor M, et al. Plasma proteomic associations with genetics and health in the UK Biobank. Nature. 2023. https://doi.org/10.1038/s41586-023-06592-6

Notes: UKB-PPP raw/regional pGWAS files may be access-controlled. Do not commit Synapse tokens, login files, or raw pGWAS archives to a public repository.

## Primary Outcomes

### Atrial Fibrillation

Nielsen et al. 2018 AF GWAS summary statistics.

- Download route: http://csg.sph.umich.edu/willer/public/afib2018/nielsen-thorolfsdottir-willer-NG2018-AFib-gwas-summary-statistics.tbl.gz
- Repository record: https://doi.org/10.18710/VC2PSH
- KP4CD record: https://www.kp4cd.org/node/1415
- Main publication: Nielsen JB, Thorolfsdottir RB, Fritsche LG, et al. Nature Genetics. 2018. https://doi.org/10.1038/s41588-018-0171-3

### Heart Failure

HERMES heart failure GWAS summary statistics.

- Download route: https://personal.broadinstitute.org/ryank/HERMES_Jan2019_HeartFailure_summary_data.txt.zip
- HERMES consortium: https://www.hermesconsortium.org/
- Main publication: Shah S, Henry A, Roselli C, et al. Nature Communications. 2020. https://doi.org/10.1038/s41467-019-13690-5

## Replication Outcomes

### FinnGen R12

- FinnGen R12: https://r12.finngen.fi/
- AF endpoint I9_AF: https://r12.risteys.finngen.fi/endpoints/I9_AF
- HF endpoint I9_HEARTFAIL: https://r12.risteys.finngen.fi/endpoints/I9_HEARTFAIL

### UKB/OpenGWAS

- `ukb-b-964`: https://gwas.mrcieu.ac.uk/datasets/ukb-b-964/
- `ukb-d-HEARTFAIL`: https://gwas.mrcieu.ac.uk/datasets/ukb-d-HEARTFAIL/
- OpenGWAS portal: https://gwas.mrcieu.ac.uk/

## Public Release Policy

The public repository should include:

- Custom analysis scripts.
- Derived result summaries, figures, and supplementary tables.
- Data-source links and file manifests.
- Software environment and citation metadata.

The public repository should exclude:

- Raw GWAS/pQTL files.
- VCF/GZ/TBI/ZIP/TAR archives downloaded from data providers.
- Any Synapse, OpenGWAS, GitHub, or Zenodo tokens.
- Local cache files, login files, and Python bytecode.
