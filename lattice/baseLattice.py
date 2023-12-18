import os
import numpy as np


class baseLattice:
    # base class for a lattice
    def __init__(self):
        # parameters in base lattice
        self.ParaIn = {}
        self.ParaSym = {}
        self.ParaNbr = {}
        self.ParaSymAt = {}
        self.ParaRel = {}
        self.HopValClas = []

    def atomName2Num(self, atm):
        """

        :param atm: atom name
        :return: number of atoms in base cell
        """
        ind = self.ParaIn["AtomName"].index(atm)
        return self.ParaIn["AtomNumber"][ind]

    def checkSupercellInfoSanity(self):
        """
        Check if the information given in the supercellXXX parts is valid
        :return:
        """

        if len(self.ParaIn["supercell"]) == 0:
            return

        if not "supercellSize" in self.ParaIn["supercell"]:
            raise ValueError("Supercell size not specified.")
        if (not "supercellVacancy" in self.ParaIn["supercell"]) \
                and (not "supercellSubstitution" in self.ParaIn["supercell"]) \
                and (not "supercellInterstitial" in self.ParaIn["supercell"]):
            raise ValueError(
                "At least give one of the values: supercellVacancy, supercellSubstitution, supercellInterstitial.")
        # check if supercellSize are positive integers
        # supercellSize is already int by the means of reading
        for num in self.ParaIn["supercell"]["supercellSize"]:
            if num <= 0:
                raise ValueError("supercell size should be >0.")

        # check if supercellVacancy is valid
        if "supercellVacancy" in self.ParaIn["supercell"]:
            rowSet = set()  # to ckeck duplicated rows
            # positionAtmCountDict={}#to check vacancy numbers
            for row in self.ParaIn["supercell"]["supercellVacancy"]:
                n0n1n2, j, atm = row
                n0, n1, n2 = n0n1n2
                rowSet.add(tuple([n0, n1, n2, j, atm]))
                # posAtmKey=tuple([n0,n1,n2,atm])
                # if not posAtmKey in positionAtmCountDict:
                #     positionAtmCountDict[posAtmKey]=1
                # else:
                #     positionAtmCountDict[posAtmKey]+=1

                # n0,n1,n2 is the base cell with the vacancy
                if (n0 < 0 or n0 >= self.ParaIn["supercell"]["supercellSize"][0]) \
                        or (n1 < 0 or n1 >= self.ParaIn["supercell"]["supercellSize"][1]) \
                        or (n2 < 0 or n2 >= self.ParaIn["supercell"]["supercellSize"][2]):
                    raise ValueError("invalid vacancy position.")
                numOfAtm = self.atomName2Num(atm)
                # j-th atom to be replaced by vacancy, among all the same kinds of atoms
                if j < 0:
                    raise ValueError("Value of vacancy should be positive.")
                if j >= numOfAtm:
                    raise ValueError(
                        "Value of vacancy of atom " + atm + " should be smaller than " + str(numOfAtm) + ".")
            # check duplicated rows
            if len(rowSet) != len(self.ParaIn["supercell"]["supercellVacancy"]):
                raise ValueError("Please remove duplicated rows in supercellVacancy.")
            # check if total vacancy numbers of an atom is <= total number of this kind of atom
            # for key in positionAtmCountDict:
            #     atmIn=key[3]
            #     totNum=self.atomName2Num(atmIn)
            #     atmNumIn=positionAtmCountDict[key]
            #     if atmNumIn>totNum:
            #         raise ValueError("The vacancy number of "+atmIn+" is greater than the number of atom "+atmIn+".")

        # check if supercellSubstitution is valid
        if "supercellSubstitution" in self.ParaIn["supercell"]:
            rowSet = set()  # to ckeck duplicated rows
            for row in self.ParaIn["supercell"]["supercellSubstitution"]:
                n0n1n2, j, atmOld, atmNewAndOrbs = row
                n0, n1, n2 = n0n1n2
                atmNew = atmNewAndOrbs[0]
                # TODO: for now we make the restriction that the substuted atom is different from the new atom
                if atmNew==atmOld:
                    raise ValueError("Substitution must be a different atom from "+atmOld+".")
                orbsNew = atmNewAndOrbs[1:]
                # if orbitals are not provided
                if len(orbsNew) == 0:
                    raise ValueError("Please give orbitals of the substitution from " + atmOld + " to " + atmNew + ".")
                rowSet.add(tuple([n0, n1, n2, j, atmOld]))
                if len(orbsNew) == 0:
                    raise ValueError("Please give the orbitals of the new atom.")

                # n0,n1,n2 is the base cell with the vacancy
                if (n0 < 0 or n0 >= self.ParaIn["supercell"]["supercellSize"][0]) \
                        or (n1 < 0 or n1 >= self.ParaIn["supercell"]["supercellSize"][1]) \
                        or (n2 < 0 or n2 >= self.ParaIn["supercell"]["supercellSize"][2]):
                    raise ValueError("invalid substitution position.")
                atmOldNum = self.atomName2Num(atmOld)  # will raise exception if atmOld is not valid
                if j < 0:
                    raise ValueError("The 3rd argument in supercellSubstitution is negative.")
                if j >= atmOldNum:
                    raise ValueError(
                        "The 3rd argument in supercellSubstitution should be smaller than " + str(atmOldNum) + ".")
            # check duplicated rows
            if len(rowSet) != len(self.ParaIn["supercell"]["supercellSubstitution"]):
                raise ValueError("Please remove duplicated rows in supercellSubstitution.")

        # check if supercellInterstitial is valid
        if "supercellInterstitial" in self.ParaIn["supercell"]:
            rowList = []  # to ckeck duplicated rows
            eps = 1e-6
            for row in self.ParaIn["supercell"]["supercellInterstitial"]:
                n0n1n2, s0s1s2, atm, orbs = row
                if len(orbs) == 0:
                    raise ValueError("Please give orbitals for interstitial atom " + atm + ".")
                n0, n1, n2 = n0n1n2
                s0, s1, s2 = s0s1s2
                rowList.append([n0, n1, n2, s0, s1, s2])
                # n0,n1,n2 is the base cell with the vacancy
                if (n0 < 0 or n0 >= self.ParaIn["supercell"]["supercellSize"][0]) \
                        or (n1 < 0 or n1 >= self.ParaIn["supercell"]["supercellSize"][1]) \
                        or (n2 < 0 or n2 >= self.ParaIn["supercell"]["supercellSize"][2]):
                    raise ValueError("invalid interstitial position.")
                # check if fractional coordinates are valid:
                if s0 < 0 or s0 > 1:
                    raise ValueError("Value of s0 is not valid.")
                if s1 < 0 or s1 > 1:
                    raise ValueError("Value of s1 is not valid.")
                if s2 < 0 or s2 > 1:
                    raise ValueError("Value of s2 is not valid.")
                if len(orbs) == 0:
                    raise ValueError("Please give the orbitals of interstitial atom " + atm + ".")

            # sort by each element of row
            rowList = np.array(rowList)
            lengthTmp = len(rowList[0])
            for j in range(0, lengthTmp):
                rowList = sorted(rowList, key=lambda row: row[j])
            lenthTmp1 = len(rowList)
            for k in range(0, lenthTmp1 - 1):
                diff = np.linalg.norm(rowList[k] - rowList[k + 1], ord=2)
                if diff < eps:
                    raise ValueError("Please remove duplicated positions in supercellInterstitial.")

        # check if supercellNbr is > 0
        if "supercellNbr" in self.ParaIn["supercell"]:
            if self.ParaIn["supercell"]["supercellNbr"] <= 0:
                raise ValueError("supercellNbr should be greater than 0.")

        return

    def array2Text(self,arr,arrName,fileName):
        """
        write array arr to text
        :param arr: array
        :param arrName: array's name
        :param fileName: output file's name
        :return:
        """
        if os.path.exists(fileName):
            append_write = 'a+'  # append if already exists
        else:
            append_write = 'w+'  # make a new file if not
        fptr=open(fileName,append_write)
        fptr.write("\n")
        fptr.write(arrName+":\n")
        # fptr.write("[")
        count=0
        for row in arr:
            elemCount=0
            # fptr.write("[")
            for v in row:
                if elemCount<len(row)-1:
                    fptr.write(str(v)+" ")
                    elemCount+=1
                else:
                    fptr.write(str(v)+"")

            if count<len(arr)-1:
                fptr.write("\n")
            count+=1
        # fptr.write("]")
        fptr.close()

    def vec2Text(self,vec,vecName,fileName):
        """
        write vector vec to text
        :param vec: vector
        :param vecName: name of vector
        :param fileName: output file
        :return:
        """
        if os.path.exists(fileName):
            append_write = 'a+'  # append if already exists
        else:
            append_write = 'w+'  # make a new file if not
        fptr=open(fileName,append_write)
        fptr.write("\n")
        fptr.write(vecName + ":\n")
        # fptr.write("[\n")
        count = 0
        for v in vec:
            if count<len(vec)-1:
                fptr.write(str(v)+" ")
            else:
                fptr.write(str(v)+"")
            count+=1
        fptr.close()