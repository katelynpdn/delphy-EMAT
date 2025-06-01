# EMAT Tree Parser

This project provides a simplified Python implementation of the **EMAT tree format** introduced in Delphy, and tools to convert between EMAT and Newick representations of phylogenetic trees.

Note: Several functions directly adapted from files in delphy/core: tree.h, phylo_tree.h, and mutations.h

---

src/
├── phyloTree.py # Core script for EMAT trees and conversions

tests/
├── phyloTreeTest.py # Test for EMAT to Newick, and Newick to EMAT conversion
└── exampleTest/
├──── ematTest.fasta # Sample FASTA file for testing
├──── ematTest.txt # Actual delphy EMAT representation of ematTest.fasta, from running gdb on Delphy

# Dependencies

-DendroPy
-treeswift
