import phyloTree as tree

def main():
    myTree = tree.PhyloTree("TCA")
    myTree.root = 6

    # Node 0
    myTree.insertNode(tree.PhyloNode(name="sample1|2019-12-31", parent=4, children=[], tMin=-1, tMax=-1, t=-1))
    myTree.at(0).addMutation(tree.Mutation(fromSeqLetter='A', site=2, toSeqLetter='C', t=-47.998690527716469))
    myTree.at(0).addMutation(tree.Mutation(fromSeqLetter='C', site=1, toSeqLetter='G', t=-5.1611070803280157))
    # Node 1
    myTree.insertNode(tree.PhyloNode(name = "sample2|2019-11-28", parent=4, children=[], tMin = -34, tMax = -34, t = -34))
    myTree.at(1).addMutation(tree.Mutation(fromSeqLetter='C', site=1, toSeqLetter='G', t=-36.913161628586991))
    # Node 2
    myTree.insertNode(tree.PhyloNode(name = "sample3|2019-12-28", parent=5, children=[], tMin = -4, tMax = -4, t = -4))
    # Node 3
    myTree.insertNode(tree.PhyloNode(name = "sample4|2019-12-30", parent=6, children=[], tMin = -2, tMax = -2, t = -2))
    myTree.at(3).addMutation(tree.Mutation(fromSeqLetter='A', site=2, toSeqLetter='G', t=-39.840408597973386))
    # Node 4
    myTree.insertNode(tree.PhyloNode(name = "", parent=5, children=[0, 1], tMin = -3.40282347e+38, tMax = 3.40282347e+38, t = -48.467903755611424))
    # Node 5
    myTree.insertNode(tree.PhyloNode(name = "", parent=6, children=[2, 4], tMin = -3.40282347e+38, tMax = 3.40282347e+38, t = -49.147026742181303))
    # Ndoe 6
    myTree.insertNode(tree.PhyloNode(name = "", parent=-1, children=[3, 5], tMin = -3.40282347e+38, tMax = 3.40282347e+38, t = -50.29117326075027))
    
    # myTree.printTree()
    print(myTree.returnNewickTree())
    myTree.writeNewickFile("newickEx.nwk")
    # print(myTree.drawTree())

    print("-----")


    myTreeCopy = tree.PhyloTree("TCA", "newickEx.nwk")
    print(myTreeCopy.returnNewickTree())
    # myTreeCopy.printTree()

if __name__ == "__main__":
    main()
