# Directly adapted from delphy/core/tree.h, phylo_tree.h, and mutations.h

from treeswift import Tree as swiftTree, Node as swiftNode
import dendropy

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
   def __init__(self, nodes=[]):
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
      return newNodeIndex
   def at(self, i) -> Node:
      return self.nodes[i]
   def size(self):
      return len(self.nodes)


# Mutation specifies specific time and place
class Mutation:
   def __init__(self, fromSeqLetter='', site=-1, toSeqLetter='', t=0.0):
      self.fromSeqLetter = fromSeqLetter
      self.site = site     #Site index (int)
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
      self.tMin = float(tMin)
      self.tMax = float(tMax)
      self.t = float(t)
      self.mutations = []
      self.missations = []
   def addMutation(self, mutation):
      self.mutations.append(mutation)
   def addMissation(self, missation):
      self.missations.append(missation)

class PhyloTreeLoc:
   def __init__(self, branch=K_NO_NODE, t=0.0):
      self.branch = branch        # Branch index - referred to by endpoint nodes
      self.t = t

# Class to encode a simplified EMAT format
# Note: Can only read from newickFile if it has comments in the form [&mutations="A1A;1.00"][&node="1,1,1"]
class PhyloTree(Tree):     
   def __init__(self, refSequence=None, newickFile=None):
      self.refSequence = []

      if isinstance(refSequence, str): # Given string input "ACG"
         self.refSequence = list(refSequence)
      elif isinstance(refSequence, list): # Given list of characters
         self.refSequence = refSequence
      elif refSequence is not None:
         raise ValueError("refSequence has to be a string or list of characters")

      super().__init__(0)

      # Read from newickFile if given
      if newickFile:
         with open(newickFile, "r") as f:
               self.readNewickString(f.read())
   
   def printTree(self):
      print(f"Tree: refSequence={self.refSequence}, root={self.root}, nodes={[i.name for i in self.nodes]}")
   def getEdgeLength(self, nodeIndex): #Get (time - parent time) at specified node, from days to years
      node = self.at(nodeIndex)
      length = (node.t - self.at(node.parent).t)/365
      return length if (length >= 0) else 0.0
   
   # Traversal helper function
   def _phyloTraversal(self, curNode, curSwiftNode):
      for child in curNode.children:
         childNode = self.at(child)
         newSwiftNode = swiftNode(childNode.name, self.getEdgeLength(child))
         curSwiftNode.add_child(newSwiftNode)
         self._phyloTraversal(childNode, newSwiftNode)

   def convertTreeSwift(self) -> swiftTree:
      tree = swiftTree()
      tree.root = swiftNode(self.at(self.root).name, self.getEdgeLength(self.root))
      if (self.size() != 0):
         self._phyloTraversal(self.at(self.root), tree.root)
      return tree
   
   def returnNexusTree(self):
      tree = self.convertTreeSwift()
      return tree.nexus()
   
   def drawTree(self):
      tree = self.convertTreeSwift()
      return tree.draw()

   #------ Functions for converting to Newick tree ------
   # Return string of all mutations for a specific node
   # String representation following the format [From Site To, Time - ParentTime]
   def mutationString(self, nodeIndex):
      node = self.at(nodeIndex)
      if (node.mutations):
         mutationList = []
         for mutation in node.mutations:
            tParent = 0
            if (node.parent != K_NO_NODE):
               tParent = self.at(node.parent).t
            else:
               tParent = self.at(self.root).t
            mutationList.append(f"{mutation.fromSeqLetter}{mutation.site}{mutation.toSeqLetter};{mutation.t - tParent}")
         return "[&mutations=\"" + ";".join(mutationList) + "\"]"
      return ""

   # Helper traversal function for returnNewickTree
   def _newickHelper(self, curNodeIndex):
      if (self.size == 0):
         return ""

      curNode = self.at(curNodeIndex)
      children = curNode.children
      nodeString = f"[&node=\"{curNode.tMin};{curNode.tMax};{curNode.t}\"]"

      if (curNode.is_inner_node()): # Internal node
         childList = []
         for i, childNodeIndex in enumerate(children):
            childList.append(self._newickHelper(childNodeIndex))
         childrenString = ",".join(childList)
         return f"({childrenString}){curNode.name}:{self.getEdgeLength(curNodeIndex)}{nodeString}{self.mutationString(curNodeIndex)}"
      else: # Leaf node
         return f"{curNode.name}:{self.getEdgeLength(curNodeIndex)}{nodeString}{self.mutationString(curNodeIndex)}"
      
   # Return tree in Newick format (with node times and mutation data as comments)
   # Note: Some tree visualizers may not accept this comment format, such as FigTree. Others such as ETE Tool Kit work fine.
   def returnNewickTree(self):
      return self._newickHelper(self.root) + ";"
   
   def writeNewickFile(self, filename):
      with open(filename, "w") as f:
         f.write(self.returnNewickTree())


   #------ Functions for reading Newick format into EMAT ------
   def _readNewickHelper(self, dendroNode, parentIndex = K_NO_NODE):

      # Access comments on each node
      if (dendroNode.annotations):
         mutationList = []
         # Extract node, mutations data from comments
         for comment in dendroNode.annotations:
            if (comment.name == "mutations"):
               commentList = comment.value.split(";")
               for i in range(0, len(commentList), 2):
                  mutationTime = float(commentList[i+1])
                  if (parentIndex != K_NO_NODE):
                     mutationTime = mutationTime + self.at(parentIndex).t  # Add commented time to parent time to get node time
                  mutationList.append(Mutation(commentList[i][0], int(commentList[i][1:-1]), commentList[i][-1], t=mutationTime))
            elif (comment.name == "node"):
               commentList = comment.value.split(";")
               tMin = float(commentList[0])
               tMax = float(commentList[1])
               t = float(commentList[2])

      # Add node to tree
      nodeName = dendroNode.taxon.label if dendroNode.taxon != None else ""
      curNodeIndex = len(self.nodes)
      numChildren = len(dendroNode.child_nodes())
      childIndices = list(range(curNodeIndex + 1, curNodeIndex + numChildren + 1))
      curNode = PhyloNode(nodeName, parentIndex, childIndices, tMin=tMin, tMax=tMax, t=t)
      self.insertNode(curNode)
      # Add each mutation to node
      for mutation in mutationList:
         curNode.addMutation(mutation)

      for child in dendroNode.child_node_iter():
         self._readNewickHelper(child, curNodeIndex)

   # Overwrites current tree data with passed in newick string
   def readNewickString(self, treeStr):
      self.root = 0
      dendroTree = dendropy.Tree.get(data=treeStr, schema='newick', preserve_underscores=True, suppress_internal_node_taxa=False)
      self._readNewickHelper(dendroTree.seed_node)
