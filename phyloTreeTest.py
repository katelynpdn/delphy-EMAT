import phyloTree as tree

def main():
    myTree = tree.PhyloTree("ACGT")
    myTree.insertNode(tree.PhyloNode("X", -1, [1, 2]))
    myTree.insertNode(tree.PhyloNode("Y", 0))
    myTree.insertNode(tree.PhyloNode("Z", 0))
    myTree.printTree()
    print(myTree.outputNewickTree())

if __name__ == "__main__":
    main()