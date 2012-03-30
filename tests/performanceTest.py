"""Tests that check performance of algorithm on randomly permuted sets of sequences derived from a median.
"""
import random

class Chromosome:
    def __init__(self):
        """Creates an empty chromosome, a lot like a list, but with a few quirks
        """
        self.chromosome = []
    
    def getRandomBreakpoint(self):
        """Creates a split chromosome
        """
        a = Chromosome()
        b = Chromosome()
        i = random.randint(len(chromosome()))
        a.chromosome = self.chromosome[:i]
        b.chromosome = self.chromosome[i:]
        return a, b
    
    def reverse(self):
        self.chromosome.reverse()
        self.chromosome = [ -element for element in self.chromosome ]
        
    def append(self, element):
        assert int(element) == element
        self.chromosome.append(element)
        
    def fuse(self, _3End):
        fusedChromosome = Chromosome()
        fusedChromosome.chromosome = self.chromosome[:] + _3End.chromosome[:]
        return fusedChromosome
    
    def __str__(self):
        """Make chromosome string
        """
        return " ".join([ str(element) in self.chromosome ])

class Genome:
    def __init__(self, elementNumber=1, chromosomeNumber=1):
        """Creates a genome with elementNumber of elements and chromosomeNumber of chromosomes
        """
        assert chromosomeNumber > 0
        self.chromosomes = [ chromosome() for i in xrange(chromosomeNumber) ]
        for element in xrange(elementNumber):
            random.choice(self.chromosomes).append(element)
        
    def clone(self):
        """Clone the genome
        """
        g = Genome()
        g.chromosomes = [ chromosome[:] for chromosome in self.chromosomes ]
        return g
        
    def permuteByDCJ(self):
        """Apply a DCJ op randomly to the genome
        """
        chr1 = self.getRandomChromosome()
        chr2 = self.getRandomChromosome()
        a, b = chr1.getRandomBreakpoint()
        c, d = chr2.getRandomBreakpoint()
        if random.random() > 0.5:
            c, d = d.reverse(), c.reverse()
        self.replaceChromosome(chr1, a.fuse(d))
        self.replaceChromosome(chr2, c.fuse(b))
        
    def permuteByInversion(self):
        """Apply a inversion op randomly to the genome
        """
        chr1 = self.getRandomChromosome()
        a, b = chr1.getRandomBreakpoint()
        b, c = b.getRandomBreakpoint()
        self.replaceChromosome(chr1, a.fuse(b.reverse()).fuse(c))
        
    def permuteByTranslocation(self, invertProb=0.5):
        """Apply a translocation op randomly to the genome
        """
        chr1 = self.getRandomChromosome()
        chr2 = self.getRandomChromosome()
        a, b = chr1.getRandomBreakpoint()
        c, d = chr2.getRandomBreakpoint()
        d, e = d.getRandomBreakpoint()
        if random.random() > invertProb: #Invert translocated with random prob
            d = d.reverse()
        self.replaceChromosome(chr1, a.fuse(d).fuse(b))
        self.replaceChromosome(chr2, c.fuse(e)) 
    
    def getRandomChromosome(self):
        """Get random chromosome in genome
        """
        chr = random.choice(self.chromosomes)
        return chr
    
    def replaceChromosome(self, previousChr, newChr):
        """Replace chromosome present in the genome
        """
        assert previousChr in self.chromosomes
        self.chromosomes[self.chromosomes.index(previousChr)] = newChr
        
    def __str__(self):
        """Make genome string
        """
        return "&".join([ str(chromosome) for chromosome in self.chromosomes ])
    
    def getDirectAdjacencies(self):
        adjacencies = set()
        for chromosome in self.chromosomes:
            chrString = chromosome.getElements()
            for i in xrange(0, len(chrString)-1):
                adjacencies.add(Adjacency(chrString[i], chrString[i+1]))
        return adjacencies
    
    def getTransitiveAdjacencies(self):
        adjacencies = set()
        for chromosome in self.chromosomes:
            chrString = chromosome.getElements()
            for i in xrange(0, len(chrString)):
                for j in xrange(i+1, len(chrString)):
                    adjacencies.add(Adjacency(chrString[i], chrString[j]))
        return adjacencies
    
    def getOutOfOrderDistance(self, otherGenome):
        a = self.getTransitiveAdjacencies()
        b = otherGenome.getTransitiveAdjacencies()
        return len(a.intersection(b)) / len(a.union(b))
    
    def getDcjDistance(self, otherGenome):
        pass

    