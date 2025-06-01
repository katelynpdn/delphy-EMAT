# EMAT Tree Parser

This project provides a simplified Python implementation of the **EMAT tree format** introduced in Delphy, and tools to convert between EMAT and Newick representations of phylogenetic trees.

Note: Several functions directly adapted from files in delphy/core: tree.h, phylo_tree.h, and mutations.h

---

```bash
├── src
│   ├── alternativePhyloTree
│   │   └── newickemat.py
│   ├── phyloTree.py            # Core script for EMAT trees and conversions
│   └── phyloTreeTest.py        # Test phyloTree.py
└── tests
    └── exampleTest
        ├── ematTest.fasta      # Sample FASTA file, referenced in phyloTreeTest.py
        └── ematTest.txt        # Delphy EMAT representation of ematTest.fasta (from running gdb on Delphy)

```

# Dependencies

-DendroPy

-TreeSwift
