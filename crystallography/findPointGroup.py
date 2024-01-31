import numpy as np


#this script finds all the point groups of a crystal
# reference:  	arXiv:1808.01590

W1=np.array([[1,0,0],
             [0,1,0],
             [0,0,1]],dtype=float)

W1bar=-W1

W2=np.array([[-1,0,0],
             [0,-1,0],
             [0,0,1]],dtype=float)

W2bar=-W2

W3=np.array([[0,-1,0],
             [1,-1,0],
             [0,0,1]],dtype=float)

W3bar=-W3

W4=np.array([[0,-1,0],
             [1,0,0],
             [0,0,1]],dtype=float)

W4bar=-W4

W6=np.array([[1,-1,0],
             [1,0,0],
             [0,0,1]])

W6bar=-W6


WMats={"1":W1,"-1":W1bar,"2":W2,"-2":W2bar,"3":W3,"-3":W3bar,"4":W4,"-4":W4bar,"6":W6,"-6":W6bar}

def pointGroupChecker(convParamsAndInfo):
    """

    :param convParamsAndInfo: information about the conventional cell, with angles in degrees
    :return: the rotation operators that leave the metric tensor invariant
    """
    a = convParamsAndInfo["Basis a"]
    b = convParamsAndInfo["Basis b"]
    c = convParamsAndInfo["Basis c"]

    B=np.array([a,b,c]).T
    G=B.T@B

