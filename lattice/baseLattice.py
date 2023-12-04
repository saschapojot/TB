# import numpy as np





class lattice:
    #base class for a lattice
    def __init__(self):
        self.ParaIn={}
        self.ParaSym={}
        self.ParaNbr={}
        self.ParaSymAt={}
        self.ParaRel={}
        self.HopValClas=[]

    def atomName2Num(self,atm):
        """

        :param atm: atom name
        :return: number of atoms
        """
        ind= self.ParaIn["AtomName"].index(atm)
        return self.ParaIn["AtomNumber"][ind]
    def checkSupercellInfoSanity(self):
        if len(self.ParaIn["supercell"])==0:
            return

        if not "supercellSize" in self.ParaIn["supercell"]:
            raise ValueError("Supercell size not specified.")
        if (not "supercellVacancy" in self.ParaIn["supercell"])\
            and (not "supercellSubstitution" in self.ParaIn["supercell"])\
            and (not "supercellInterstitial" in self.ParaIn["supercell"]):
            raise ValueError("At least give one of the values: supercellVacancy, supercellSubstitution, supercellInterstitial.")
        # check if supercellSize are positive integers
        # supercellSize is already int by the means of reading
        for num in self.ParaIn["supercell"]["supercellSize"]:
            if num<=0:
                raise ValueError("supercell size should be >0.")

        #check if supercellVacancy is within the range of supercellSize
        if "supercellVacancy" in self.ParaIn["supercell"]:
            for row in self.ParaIn["supercell"]["supercellVacancy"]:
                n0n1n2,j,atm=row
                n0,n1,n2=n0n1n2
                #n0,n1,n2 is the base cell with the vacancy
                if (n0<0 or n0>=self.ParaIn["supercell"]["supercellSize"][0])\
                    or (n1<0 or n1 >=self.ParaIn["supercell"]["supercellSize"][1])\
                    or (n2<0 or n2 >=self.ParaIn["supercell"]["supercellSize"][2]):
                    raise ValueError("invalid vacancy position.")
                numOfAtm=self.atomName2Num(atm)
                #j-th atom to be replaced by vacancy, among all the same kinds of atoms
                if j<0:
                    raise ValueError("Value of vacancy should be positive.")
                if j>=numOfAtm:
                    raise ValueError("Value of vacancy of atom "+atm+" should be smaller than "+str(numOfAtm)+".")



        #check if supercellSubstitution is valid
        if "supercellSubstitution" in self.ParaIn["supercell"]:
            for row in self.ParaIn["supercell"]["supercellSubstitution"]:
                n0n1n2,j,atmOld,atmNewAndOrbs=row
                n0,n1,n2=n0n1n2
                atmNew=atmNewAndOrbs[0]
                orbsNew=atmNewAndOrbs[1:]
                if len(orbsNew)==0:
                    raise ValueError("Please give the orbitals of the new atom.")

                # n0,n1,n2 is the base cell with the vacancy
                if (n0 < 0 or n0 >= self.ParaIn["supercell"]["supercellSize"][0]) \
                        or (n1 < 0 or n1 >= self.ParaIn["supercell"]["supercellSize"][1]) \
                        or (n2 < 0 or n2 >= self.ParaIn["supercell"]["supercellSize"][2]):
                    raise ValueError("invalid substitution position.")
                atmOldNum=self.atomName2Num(atmOld)#will raise exception if atmOld is not valid
                if j<0:
                    raise ValueError("The 3rd argument in supercellSubstitution is negative.")
                if j>= atmOldNum:
                    raise ValueError("The 3rd argument in supercellSubstitution should be smaller than "+str(atmOldNum)+".")



        #check if supercellInterstitial is valid
        if "supercellInterstitial" in self.ParaIn["supercell"]:
            for row in self.ParaIn["supercell"]["supercellInterstitial"]:
                n0n1n2,s0s1s2,atm,orbs=row
                n0,n1,n2=n0n1n2
                s0,s1,s2=s0s1s2
                # n0,n1,n2 is the base cell with the vacancy
                if (n0 < 0 or n0 >= self.ParaIn["supercell"]["supercellSize"][0]) \
                        or (n1 < 0 or n1 >= self.ParaIn["supercell"]["supercellSize"][1]) \
                        or (n2 < 0 or n2 >= self.ParaIn["supercell"]["supercellSize"][2]):
                    raise ValueError("invalid interstitial position.")
                #check if fractional coordinates are valid:
                if s0<0 or s0>1:
                    raise ValueError("Value of s0 is not valid.")
                if s1<0 or s1>1:
                    raise ValueError("Value of s1 is not valid.")
                if s2<0 or s2>1:
                    raise ValueError("Value of s2 is not valid.")
                if len(orbs)==0:
                    raise ValueError("Please give the orbitals of interstitial atom "+atm+".")


        return