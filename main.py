import sys
import numpy as np
import matplotlib.pyplot as plt
import os
import pathlib
from datetime import datetime

plt.close("all")
# Main
from cd.SymGroup import GetSpaceGroupPrimitive
from cd.NbrAtom  import FindNeighbor
from cd.SymAtom import FindAtomSymmetry
from cd.HopRel   import FindRelation
from cd.HopVal   import GetHoppingValue
from cd.HmtReal  import GetHamiltonianReal
from cd.KPoint   import GetKPoint
from cd.HmtK     import GetHamiltonianK
from cd.Eig      import GetEigenSolution
# Read & Write
from rw.ReadTBIN import ReadInput
from rw.ReadHop  import ReadHopping
from rw.ReadRel  import ReadRelation
from rw.ReadHmt  import ReadHamiltonianReal
from rw.ReadKpt  import ReadKPoint
from rw.WriteRel import WriteRelation
from rw.WriteHmt import WriteHamiltonianReal
# Plot
from pl.PlotEB   import PlotEnergyBand
from pl.PlotAt   import PlotAtoms
from pl.PlotHop  import PlotHoppingTerm
# from pl.PlotBZ  import PlotBrillouinZone

#This is the demo program for Auh-17-2023
#1. It reads config file containing information of the crystal,
# the reader may give either a primitive cell or a conventional cell,
# the program will compute primitive from conventional, or conventional from primitive.
#2. When the

# material = 'data/ABO3/primitive_TBIN_ABO3.txt'
material = 'data/Graphene/primitive_TBIN_Graphene.txt'
# material = 'data/h-BN/primitive_TBIN_h-BN.txt'
# material = 'data/NaCl/primitive_TBIN_NaCl.txt'
# material = 'data/Si/primitive_TBIN_Si.txt'
lenParams=len(sys.argv)
if lenParams<2:
    sys.argv.append(material)
else:
    sys.argv[1] = material

lenParams=len(sys.argv)
if lenParams!=2:
    raise RuntimeError("Wrong number of arguments.")

inConfigName=sys.argv[1]
#check if given path exist
if not os.path.exists(inConfigName):
    raise NotADirectoryError("path does not exist.")
#check if given file exists
if not os.path.isfile(inConfigName):
    raise FileNotFoundError("Not a file.")
# print(inConfigName)
# inConfigFolder=pathlib.Path(inConfigName).parent
# print(inConfigFolder)
#read info
ParaIn=ReadInput(inConfigName)
ParaSym = GetSpaceGroupPrimitive(ParaIn)
ParaIn["origin Bilbao"]=ParaSym["origin Bilbao"]
ParaNbr    = FindNeighbor(ParaIn)
Name=ParaIn["Name"]
# PlotAtoms(ParaIn,ParaNbr,Name)
tFindingRelationStart=datetime.now()
# ParaSymAt  = FindAtomSymmetry(ParaIn,ParaSym,ParaNbr)
# ParaRel    = FindRelation(ParaIn,ParaSym,ParaNbr,ParaSymAt)
ParaRel    = ReadRelation(ParaIn["Folder"]+ "/HopRel.txt")
# PlotHoppingTerm(ParaIn,ParaNbr,ParaRel,Name,[5,6])
tFindingRelationEnd=datetime.now()
print("Finding symmetry relations: ",tFindingRelationEnd-tFindingRelationStart)
# WriteRelation(ParaIn,ParaRel)
# From hopping relations to Hamiltonian in real space
HopValIn   = ReadHopping(ParaIn["Folder"]+"/HopValIN_"+ParaIn["Name"]+".txt")
ParaRel    = ReadRelation(ParaIn["Folder"]+ "/HopRel.txt")
HopValClas = GetHoppingValue(HopValIn, ParaRel)
ParaHmtR   = GetHamiltonianReal(ParaRel,HopValClas)
WriteHamiltonianReal(ParaIn, ParaHmtR)
# From Hamiltonian in real space to Hamiltonian in k space
ParaHmtR   = ReadHamiltonianReal(ParaIn)
ParaKptIn = ReadKPoint(ParaIn["Folder"]+"/KptIN_" + Name + ".txt")
ParaKpt = GetKPoint(ParaKptIn)
ParaHmtK   = GetHamiltonianK(ParaHmtR,ParaKpt)
tEigStart=datetime.now()
ParaEig    = GetEigenSolution(ParaHmtK)
ParaEigPlt = PlotEnergyBand(ParaKpt,ParaEig,ParaIn["Folder"])
tEigEnd=datetime.now()
print(Name+ " enerygy band: ",tEigEnd-tEigStart)