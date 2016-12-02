import numpy as np
import random
import matplotlib.pyplot as plt


'''
averages for E.coli:
translation: 15aa/s
transcription: 55nt/s
protein length: 267 aa
mRNA length: 801 nt
time of life of mRNA: 300 s   (RNA.p_lambda)
time of life of protein:  7200 s   (Protein.p_lambda)

X.prod_timer
801:55 = ~14 -> (for DNA, how many seconds it takes to produce one mRNA)
267:15= ~18 -> (for mRNA, how many seconds it takes to produce one protein)

'''


class Molecule(object):
    def __init__(self, simulation):
        self.simulation = simulation

    @staticmethod
    def get_life_time(for_whom, p_lambda):  # get time of life from Poisson distr, for lambda being av time of life
        k = np.random.poisson(p_lambda)
        return k


class DNA(Molecule):
    prod_timer = 14   # change accordingly to data

    def __init__(self, simulation):
        Molecule.__init__(self, simulation)
        self.ticks = 0

    def transcription(self):
        self.ticks += 1
        if self.ticks == DNA.prod_timer:
            new_rna = RNA(self.simulation)
            self.simulation.list_rna.append(new_rna)
            self.ticks = 0


class RNA(Molecule):
    prod_timer = 18  # change accordingly to data
    p_lambda = 300   # change accordingly to data

    def __init__(self, simulation):
        Molecule.__init__(self, simulation)
        self.ticks = 0
        self.life_time = self.get_life_time(self, RNA.p_lambda)

    def go_RNA(self):
        self.life_time -= 1
        if self.life_time >= 0:
            self.translation()
        else:
            self.dead()

    def dead(self):
        self.simulation.list_rna.remove(self)

    def translation(self):
        self.ticks += 1
        if self.ticks == RNA.prod_timer:
            new_p = Protein(self.simulation)
            self.simulation.list_prot.append(new_p)
            self.ticks = 0


class Protein(Molecule):
    p_lambda = 7200   # change accordingly to data

    def __init__(self, simulation):
        Molecule.__init__(self, simulation)
        self.life_time = self.get_life_time(self, Protein.p_lambda)

    def go_PROT(self):
        self.life_time -= 1
        if self.life_time >= 0:
            pass
        else:
            self.dead()

    def dead(self):
        self.simulation.list_prot.remove(self)


class Simulation(object):
    def __init__(self, t_sim):
        self.t_sim = t_sim
        d = DNA(self)
        self.list_dna, self.list_rna, self.list_prot, = [d], [], []  # list with objects
        self.A, self.B, self.C, self.T = [1], [0], [0], [0]  # list with amounts of molecules
        self.loop()

    def loop(self):
        for t in range(self.t_sim):
            for rna in self.list_rna:
                rna.go_RNA()
            for dna in self.list_dna:
                dna.transcription()
            for prot in self.list_prot:
                prot.go_PROT()
            self.A.append(len(self.list_dna))
            self.B.append(len(self.list_rna))
            self.C.append(len(self.list_prot))
            self.T.append(len(self.T))
        self.plot()

    def plot(self):
        a1, a2, a3, a4 = str(self.T[-1]), str(self.A[-1]), str(self.B[-1]), str(self.C[-1])
        a = " in " + a1 + " second of simulation there were:\n DNA molecules: " + a2 + ", RNA molecules: " + a3 + ", protein molecules: " + a4
        b = " max RNA: " + str(max(self.B)) + " in " + str(np.argmax(self.B)) + " sec"
        c = " max proteins: " + str(max(self.C)) + " in " + str(np.argmax(self.C)) + " sec"
        print a
        print b
        print c

        plt.figure(1)
        plt.plot(self.T, self.A, color='red')
        plt.plot(self.T, self.B, color='blue')
        plt.plot(self.T, self.C, color='yellow')
        plt.xlabel('Time')
        plt.ylabel('Number of molecules')
        plt.legend(['dna', 'rna', 'protein'],loc='upper left')

        plt.figure(2)
        plt.subplot(311)
        plt.plot(self.T, self.A, color='red')
        plt.legend(['dna'],loc='upper left')
        plt.xlabel('Time')
        plt.ylabel('Number of molecules')

        plt.subplot(312)
        plt.plot(self.T, self.B, color='blue')
        plt.legend(['rna'],loc='upper left')
        plt.xlabel('Time')
        plt.ylabel('Number of molecules')

        plt.subplot(313)
        plt.plot(self.T, self.C, color='yellow')
        plt.legend(['protein'],loc='upper left')
        plt.xlabel('Time')
        plt.ylabel('Number of molecules')

        plt.show()


if __name__ == "__main__":
    time_of_sim = 10800             # here put the time of simulation
    s = Simulation(time_of_sim)
