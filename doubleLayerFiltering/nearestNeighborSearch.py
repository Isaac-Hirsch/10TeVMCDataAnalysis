import numpy as np

#binary tree object
class biTree(object):
    def __init__(self, value):
        self.value=value
        self.left=None
        self.right=None

#nodes that will make up a binary tree
class node(object):
    def __init__(self,value):
        self.value=value
        self.left=None
        self.right=None

#end of binary tree will be a list of theta and phis
class leaf(object):
    def __init__(self):
        #leafs will be itialized empty
        self.value=[]
    def add(self,value):
        #Adding new items to a list
        self.value.append(value)

#specific binary tree built for taking all of the hits and splitting them up into small subregions
class biTreethetaPhi(biTree):
    def __init__(self,depth):
        super().__init__(0)
        self.depth=depth
        self.fill(self, 0, np.pi, -np.pi, np.pi, depth)
    
    def fill(self,curNode,thetaMin,thetaMax,phiMin,phiMax,depth):

        if depth!=0:

            if (depth//2)==0:
                curNode.left=node(thetaMin+(thetaMin+thetaMax)/4)
                curNode.right=node(thetaMin+(thetaMin+thetaMax)*3/4)

                self.fill(curNode.left, thetaMin, thetaMin+(thetaMin+thetaMax)/2,phiMin,phiMax, depth-1)
                self.fill(curNode.right, thetaMin+(thetaMin+thetaMax)/2, thetaMax, phiMin,phiMax, depth-1)

            else:
                curNode.left=node(phiMin+(phiMin+phiMax)/4)
                curNode.right=node(phiMin+(phiMin+phiMax)*3/4)

                self.fill(curNode.left, thetaMin, thetaMax, phiMin, phiMin+(thetaMin+thetaMax)/2, depth-1)
                self.fill(curNode.right, thetaMin, thetaMax, phiMin+(thetaMin+thetaMax)/2, phiMax, depth-1)

        else:
            curNode.left=leaf()
            curNode.right=leaf()
    
    #insert a tuple of (theta, phi) into the binary tree
    def insert(self,value):
        assert type(value)==tuple, "tried to inser a value that was not a tuple"
        curNode=self

        for i in self.depth:
            #theta seach
            if i//2==0:
                if value[0] > curNode.value:
                    curNode=curNode.right
                else:
                    curNode=curNode.left

            #phi search
            else:
                if value[1] > curNode.value:
                    curNode=curNode.right
                else:
                    curNode=curNode.left
            
            #appending theta onto the left node of the leaf
            curNode.left.append(value[0])
            #appending phi onto the right node of the leaf
            curNode.right.append(value[1])

tree=biTreethetaPhi(10)