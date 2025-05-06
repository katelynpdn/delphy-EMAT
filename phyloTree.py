# Directly adapted from delphy/core/tree.h, phylo_tree.h, and mutations.h

from treeswift import Tree as swiftTree, Node as swiftNode

K_NO_NODE = -1  # Represents blank reference to node
class Node:
   def __init__(self, parent=K_NO_NODE, children=[]):
      self.parent = parent   # Parent index (int)
      self.children = children        # List of child indices (int)
   def __str__(self):
      return f"Node: parent={self.parent}, children={self.children}"
   def is_inner_node(self) -> bool:
      return self.children
   def is_tip(self) -> bool:
      return not self.children

class Tree:
   def __init__(self, nodes):
      self.root = K_NO_NODE   # Parent index (int)
      self.nodes = nodes      # List of Nodes
   def __init__(self, numNodes):
      self.root = K_NO_NODE   
      self.nodes = []
      for i in range(numNodes):
         self.nodes.append(Node())
   def insertNode(self, node) -> int:
      newNodeIndex = len(self.nodes)
      self.nodes.append(node)
      if (newNodeIndex) == 0:
         self.root = newNodeIndex
      return newNodeIndex
   def at(self, i) -> Node:
      return self.nodes[i]
   def size(self):
      return len(self.nodes)


# Mutation specifies specific time and place
class Mutation:
   def __init__(self, fromSeqLetter='', site=-1, toSeqLetter='', t=0.0):
      self.fromSeqLetter = fromSeqLetter
      self.site = site
      self.toSeqLetter = toSeqLetter
      self.t = t

# MissationMap "encodes which sites are missing in sequences downstream of a point X on a phylo tree." (delphy)
# Also encodes state at those sites at X.
#   a. Interval set of gaps: vector of [start,end) pairs for the ranges of missing site indices
#   b. Map containing delta between states of missing sites at X and in reference sequence.
class MissationMap:
   def __init__(self):
      self.intervals = []       # List of tuples (start,end)
      self.from_states = {}     # Map: Key: Site_index (int) Value: Real_seq_letter

class PhyloNode(Node):
   def __init__(self, name="", parent=K_NO_NODE, children=[], tMin=42.0, tMax=-42.0, t=0.0):
      super().__init__(parent, children)
      self.name = name
      self.tMin = tMin
      self.tMax = tMax
      self.t = t
      self.mutations = []
      self.missations = []

class PhyloTreeLoc:
   def __init__(self, branch=K_NO_NODE, t=0.0):
      self.branch = branch        # Branch index - referred to by endpoint nodes
      self.t = t

class PhyloTree(Tree):     
   def __init__(self, refSequence=[]):
      self.refSequence = refSequence    # List of chars (A,C,G, or T)
      super().__init__(0)
   def printTree(self):
      print(f"Tree: refSequence={self.refSequence}, root={self.root}, nodes={self.nodes}")
   
   # Traversal helper function
   def _traversal(self, curNode, curSwiftNode):
      for child in curNode.children:
         childNode = self.at(child)
         newSwiftNode = swiftNode(childNode.name)
         curSwiftNode.add_child(newSwiftNode)
         curSwiftNode = newSwiftNode
         self._traversal(childNode, curSwiftNode)
   def convertTreeSwift(self) -> swiftTree:
      tree = swiftTree()
      tree.root = swiftNode(self.at(self.root).name)
      if (self.size() != 0):
         self._traversal(self.at(self.root), tree.root)
      return tree
   def outputNewickTree(self):
      tree = self.convertTreeSwift()
      # for node in tree.traverse_preorder():
      #    print(node)
      return tree.newick()

   # Print in newick format using post-order traversal
   # def output_newick_tree(self):
   #    if (self.size == 0):
   #       return
      
   #    work_stack = []
   #    work_stack.append(self.root)
   #    while (work_stack):
   #       curNode = work_stack.pop()
   #       children = curNode.children
         
   #       if (self.at(curNode).is_inner_node()) {
   #          *os_ << "(";
   #       }