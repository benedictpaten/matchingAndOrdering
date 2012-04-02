"""Tests that check performance of algorithm on randomly permuted sets of sequences derived from a median.
"""
import random

class Chromosome:
    def __init__(self):
        """Creates an empty chromosome, which is a list of signed integers, and which is immutable.
        """
        self.chromosome = []
    
    def getRandomBreakpoint(self):
        """Creates a split chromosome
        """
        a = Chromosome()
        b = Chromosome()
        i = random.randint(0, len(self))
        a.chromosome = self.chromosome[:i]
        b.chromosome = self.chromosome[i:]
        return a, b
    
    def getRandomSegment(self):
        """Creates a split chromosome
        """
        a, b = self.getRandomBreakpoint()
        b, c = b.getRandomBreakpoint()
        return a, b, c
    
    def getReverse(self):
        """Returns reverse
        """
        c = Chromosome()
        c.chromosome = [ -element for element in self.chromosome ]
        c.chromosome.reverse()
        return c
        
    def append(self, element):
        """Adds an element to the chromosome
        """
        assert int(element) == element
        assert element != 0
        self.chromosome.append(element)
        
    def fuse(self, _3End):
        """Joins two chromosomes together
        """
        fusedChromosome = Chromosome()
        fusedChromosome.chromosome = self.chromosome[:] + _3End.chromosome[:]
        return fusedChromosome
    
    def clone(self):
        """Copies the chromosome
        """
        c = Chromosome()
        c.chromosome = self.chromosome[:]
        return c
    
    def getOrderedElements(self):
        """Returns an ordered and oriented list of the elements
        """
        return self.chromosome[:]
    
    def __len__(self):
        return len(self.chromosome)
    
    def __str__(self):
        """Make chromosome string
        """
        return " ".join([ str(element) for element in self.chromosome ])

def orderedPair(element1, element2, distance):
    """Ordered pair tuple
    """
    if element1 < element2:
        return (element1, element2, distance)
    return (element2, element1, distance)

class Genome:
    def __init__(self, elementNumber=1, chromosomeNumber=1):
        """Creates a genome with elementNumber of elements and chromosomeNumber of chromosomes
        """
        assert chromosomeNumber >= 0
        assert elementNumber >= chromosomeNumber
        self.chromosomes = []
        for i in xrange(chromosomeNumber):
            chromosome = Chromosome()
            chromosome.append(i+1)
            self.chromosomes.append(chromosome)
        for element in xrange(chromosomeNumber+1, elementNumber+1):
            random.choice(self.chromosomes).append(element)
    
    def getChromosomeNumber(self):
        return len(self.chromosomes)
        
    def getElements(self):
        """Gets set of elements (ints) in genome
        """
        elements = set()
        for chromosome in self.chromosomes:
            for element in chromosome.getOrderedElements():
                elements.add(abs(element))
        return elements
        
    def clone(self):
        """Clone the genome
        """
        g = Genome()
        g.chromosomes = [ chromosome.clone() for chromosome in self.chromosomes ]
        return g
    
    def permuteByInversion(self):
        """Apply a inversion op randomly to the genome
        """
        chr1 = self.getRandomChromosome()
        a, b, c = chr1.getRandomSegment()
        self.replaceChromosome(chr1, a.fuse(b.getReverse()).fuse(c))
        
    def permuteByDcj(self):
        """Apply a DCJ op randomly to the genome
        """
        chr1 = self.getRandomChromosome()
        chr2 = self.getRandomChromosome()
        if chr1 == chr2:
            a, b, c = chr1.getRandomSegment()
            if random.random() > 0.5:
                self.replaceChromosome(chr1, a.fuse(b.getReverse()).fuse(c))
            else:
                self.replaceChromosome(chr1, a.fuse(c))
                self.addChromosome(b)
        else:
            a, b = chr1.getRandomBreakpoint()
            c, d = chr2.getRandomBreakpoint()
            if random.random() > 0.5:
                c, d = d.getReverse(), c.getReverse()
            self.replaceChromosome(chr1, a.fuse(d))
            self.replaceChromosome(chr2, c.fuse(b))
        
    def permuteByTranslocation(self, invertProb=0.5):
        """Apply a translocation op randomly to the genome
        """
        chr1 = self.getRandomChromosome()
        chr2 = self.getRandomChromosome()
        if chr1 == chr2:
            a, b, c = chr1.getRandomSegment()
            if random.random() > invertProb: #Invert translocated with random prob
                b = b.getReverse()
            if random.random() > 0.5:
                a1, a2 = a.getRandomBreakpoint()
                self.replaceChromosome(chr1, a1.fuse(b).fuse(a2).fuse(c))
            else:
                c1, c2 = c.getRandomBreakpoint()
                self.replaceChromosome(chr1, a.fuse(c1).fuse(b).fuse(c2))
        else:
            a, b, c = chr1.getRandomSegment()
            d, e = chr2.getRandomBreakpoint()
            if random.random() > invertProb: #Invert translocated with random prob
                b = b.getReverse()
            self.replaceChromosome(chr1, a.fuse(c))
            self.replaceChromosome(chr2, d.fuse(b).fuse(e)) 
        
    def __str__(self):
        """Make genome string
        """
        return " & ".join([ str(chromosome) for chromosome in self.chromosomes ])
    
    def getOutOfOrderDistance(self, otherGenome):
        """Computes fraction of transitive adjacencies shared by two.
        """
        def stripDistances(transitiveAdjacencies):
            return set([ (element1, element2) for element1, element2, distance in transitiveAdjacencies ])
        a = stripDistances(self.getTransitiveAdjacencies())
        b = stripDistances(otherGenome.getTransitiveAdjacencies())
        return 1.0 - len(a.intersection(b)) / float(len(a.union(b)))
    
    def getCircularDcjDistance(self, otherGenome):
        """Computes DCJ Distance
        """
        greyAdjacencies = self.getCircularAdjacencies()
        blackAdjacencies = otherGenome.getCircularAdjacencies()
        assert len(greyAdjacencies) == len(blackAdjacencies)
        assert len(greyAdjacencies) == 2*len(self.getElements())
        cycles = 0
        def fn(adjacencies1, adjacencies2, element):
            if element not in adjacencies1:
                return
            otherElement = adjacencies1[element]
            adjacencies1.pop(element)
            fn(adjacencies2, adjacencies1, otherElement)
        for element in list(greyAdjacencies.keys()):
            if element in greyAdjacencies:
                cycles += 1
                fn(greyAdjacencies, blackAdjacencies, element)
        assert len(greyAdjacencies.values()) == 0
        assert len(blackAdjacencies.values()) == 0
        assert cycles % 2 == 0 #Cycles get traversed in both directions
        return len(self.getElements()) - cycles/2
    
    #Remaining functions are private
    
    def addChromosome(self, chromosome):
        if len(chromosome) > 0:
            self.chromosomes.append(chromosome)
    
    def getRandomChromosome(self):
        """Get random chromosome in genome
        """
        chr = random.choice(self.chromosomes)
        return chr
    
    def replaceChromosome(self, previousChr, newChr):
        """Replace chromosome present in the genome
        """
        assert previousChr in self.chromosomes
        if len(newChr) == 0: #Do not add a zero length chromosome
            self.chromosomes.remove(previousChr)
        else:
            self.chromosomes[self.chromosomes.index(previousChr)] = newChr
    
    def getTransitiveAdjacencies(self):
        """Get set of transitive adjacencies
        """
        adjacencies = set()
        for chromosome in self.chromosomes:
            chrString = chromosome.getOrderedElements()
            for i in xrange(0, len(chrString)):
                for j in xrange(i+1, len(chrString)):
                    assert chrString[i] != chrString[j]
                    oP = orderedPair(chrString[i], chrString[j], j - i)
                    assert oP not in adjacencies
                    adjacencies.add(oP)
                    #For symmetry do in reverse orientation
                    oP = orderedPair(-chrString[j], -chrString[i], j - i)
                    assert oP not in adjacencies
                    adjacencies.add(oP)
        return adjacencies
    
    def getCircularAdjacencies(self):
        """Gets set of direct 'abutting' adjacencies as hash, treating chromosomes as circular
        """
        adjacencies = {}
        for chromosome in self.chromosomes:
            chrString = chromosome.getOrderedElements()
            assert len(chrString) > 0
            for i in xrange(len(chrString)-1):
                leftElement = -chrString[i]
                rightElement = chrString[i+1]
                assert abs(leftElement) != abs(rightElement)
                assert leftElement not in adjacencies
                assert rightElement not in adjacencies
                adjacencies[leftElement] = rightElement
                adjacencies[rightElement] = leftElement
            #Add circular adjacency
            leftElement = chrString[0]
            rightElement = -chrString[-1]
            assert leftElement not in adjacencies
            assert rightElement not in adjacencies
            adjacencies[leftElement] = rightElement
            adjacencies[rightElement] = leftElement
        return adjacencies
    
class MedianHistory:
    """Represents a simulation of a median genome and its children
    """
    def __init__(self, medianGenome, leafGenomeNumber):
        self.medianGenome = medianGenome
        self.leafGenomes = [ medianGenome.clone() for i in xrange(leafGenomeNumber) ]
        
    def permuteLeafGenomes(self, operationNumber=1, doInversion=False, doDcj=False, doTranslocation=False):
        """Permutes the child genomes using the given operator types, each for
        operationNumber of times. 
        """
        for leafGenome in self.leafGenomes:
            for i in xrange(operationNumber):
                if doInversion:
                    leafGenome.permuteByInversion()
                if doDcj:
                    leafGenome.permuteByDcj()
                if doTranslocation:
                    leafGenome.permuteByTranslocation()
                    
    def getMedianGenome(self):
        return self.medianGenome
    
    def getLeafGenomes(self):
        return self.leafGenomes[:]
        
    def getLeafGenomeString(self):
        """Gets string in GRIMM format representing the leaf genomes
        """
        return "".join([ '>\n%s\n' % str(leafGenome) for leafGenome in self.leafGenomes ])
    
    def getMedianDcjDistance(self, genome):
        """Gets the median DCJ Distance
        """
        return sum([ leafGenome.getCircularDcjDistance(genome) for leafGenome in self.leafGenomes ])/len(self.leafGenomes)
    
    def getMedianOutOfOrderDistance(self, genome):
        """Gets the median out of order distance
        """
        return sum([ leafGenome.getOutOfOrderDistance(genome) for leafGenome in self.leafGenomes ])/len(self.leafGenomes)

