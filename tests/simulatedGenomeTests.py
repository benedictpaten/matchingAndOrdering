#!/usr/bin/env python

#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
import unittest
import sys
import os

from sonLib.bioio import getLogLevelString
from sonLib.bioio import parseSuiteTestOptions
from sonLib.bioio import system, popenCatch

import matchingAndOrdering.tests.simulatedGenome

from matchingAndOrdering.tests.simulatedGenome import Chromosome, Genome, MedianHistory

class TestCase(unittest.TestCase):
    def setUp(self):
        self.tempFile = os.path.join(os.getcwd(), "simulatedGenomeTempFile.txt")
        self.elementNumbers=(100,)
        self.chromosomeNumbers=(1,) # 2, 5)
        self.leafGenomeNumbers=(3, 5, 10)
        self.operationNumber = (1, 10, 100) # 1000)
        self.operationType = ((True, False, False),) # (False, True, False), (False, False, True))
        self.replicates = 5
    
    def testReferenceAndGrimmAlgorithms(self):
        """Iterates through a list of simulation variants and prints results
        """
        headerLine = "\t".join(("elementNumber", "chromosomeNumber", "leafGenomeNumber", 
                                 "operationNumber",
                                 "doInversion", "doDcj", "doTranslocation", "replicate", 
                                 "medianDCJDistance", "medianOutOfOrderDistance", 
                                 "medianDCJDistanceForReferenceAlgorithm",
                                 "medianOutOfOrderDistanceForReferenceAlgorithm", 
                                 "medianDCJDistanceForGrimm", "medianOutOfOrderDistanceForGrimm", "medianGenome", "medianGenomeForReferenceAlgorithm", "medianGenomeForGrimm"))
        if getLogLevelString() in  ("DEBUG", "INFO" ):
            print headerLine
        for elementNumber in self.elementNumbers:
            for chromosomeNumber in self.chromosomeNumbers:
                for leafGenomeNumber in self.leafGenomeNumbers:
                    for operationNumber in self.operationNumber:
                        for doInversion, doDcj, doTranslocation in self.operationType:
                            for replicate in xrange(self.replicates):
                                medianHistory = MedianHistory(Genome(elementNumber=elementNumber, chromosomeNumber=chromosomeNumber), leafGenomeNumber=leafGenomeNumber)
                                medianHistory.permuteLeafGenomes(operationNumber=operationNumber, doInversion=doInversion, doDcj=doDcj, doTranslocation=doTranslocation)
                                medianDCJDistance = medianHistory.getMedianDcjDistance(medianHistory.getMedianGenome())
                                medianOutOfOrderDistance = medianHistory.getMedianOutOfOrderDistance(medianHistory.getMedianGenome())
                                #Now run reference problem algorithm   
                                referenceProblemMedianGenome = runReferenceMedianProblemTest(medianHistory)      
                                medianDCJDistanceForReferenceAlgorithm = medianHistory.getMedianDcjDistance(referenceProblemMedianGenome)
                                medianOutOfOrderDistanceForReferenceAlgorithm = medianHistory.getMedianOutOfOrderDistance(referenceProblemMedianGenome)
                                #Now run GRIMM
                                grimmProblemMedianGenome = runGrimmMedianProblemTest(medianHistory)      
                                medianDCJDistanceForGrimm = medianHistory.getMedianDcjDistance(grimmProblemMedianGenome)
                                medianOutOfOrderDistanceForGrimm = medianHistory.getMedianOutOfOrderDistance(grimmProblemMedianGenome)
                                #Now print a report line
                                line = "\t".join([ str(i) for i in 
                                (elementNumber, chromosomeNumber, leafGenomeNumber, 
                                 operationNumber,
                                 doInversion, doDcj, doTranslocation, replicate, 
                                 medianDCJDistance, medianOutOfOrderDistance, 
                                 medianDCJDistanceForReferenceAlgorithm,
                                 medianOutOfOrderDistanceForReferenceAlgorithm, 
                                 medianDCJDistanceForGrimm, medianOutOfOrderDistanceForGrimm,
                                 "'%s'" % str(medianHistory.getMedianGenome()),
                                 "'%s'" % str(referenceProblemMedianGenome),
                                 "'%s'" % str(grimmProblemMedianGenome)) ])
                                #Print line
                                if getLogLevelString() in ("DEBUG", "INFO"):
                                    print line
    
    def testChromosome(self):
        """Test basic functions of a chromosome.
        """
        c = Chromosome()
        #Append method
        c.append(1)
        c.append(2)
        #Test string generator
        self.assertEquals(c.getOrderedElements(), [ 1, 2 ])
        self.assertEqual(str(c), "1 2")
        
        #Test get reverse
        d = c.getReverse()
        self.assertEqual(str(c), "1 2")
        self.assertEqual(str(d), "-2 -1")
        
        #Test fuse
        e = c.fuse(d)
        self.assertEqual(str(c), "1 2")
        self.assertEqual(str(d), "-2 -1")
        self.assertEqual(str(e), "1 2 -2 -1")
        
        #Test fuse of null chromosome
        e = e.fuse(Chromosome())
        self.assertEquals(e.getOrderedElements(), [ 1, 2, -2, -1 ])
        self.assertEqual(str(e), "1 2 -2 -1")
        
        #Test random breakpoint
        for i in xrange(10):
            f, g = e.getRandomBreakpoint()
            self.assertEqual(str(f.fuse(g)), "1 2 -2 -1")
        
        self.assertEqual(str(d), str(d.clone()))
        
    def testGenome(self):
        """Test basic functions of a genome
        """
        d = Genome(elementNumber=10, chromosomeNumber=3)
        self.assertEqual(d.getChromosomeNumber(), 3)
        self.assertEqual(d.getElements(), set(range(1, 11)))
        
        #Test clone
        self.assertEqual(str(d), str(d.clone()))
        
        #Test inversions
        for i in xrange(100):
            e = d.clone()
            self.assertEquals(d.getOutOfOrderDistance(e), 0)
            self.assertEquals(d.getCircularDcjDistance(e), 0)
            e.permuteByInversion()
            self.assertEquals(e.getChromosomeNumber(), 3)
            self.assertEquals(e.getElements(), d.getElements())
            self.assertTrue(d.getOutOfOrderDistance(e) >= 0)
            self.assertTrue(d.getCircularDcjDistance(e) in [ 0, 1 ])
            
        #Test dcj
        for i in xrange(100):
            e = d.clone()
            self.assertEquals(d.getOutOfOrderDistance(e), 0)
            self.assertEquals(d.getCircularDcjDistance(e), 0)
            e.permuteByDcj()
            self.assertEquals(e.getElements(), d.getElements())
            self.assertTrue(d.getOutOfOrderDistance(e) >= 0)
            self.assertTrue(d.getCircularDcjDistance(e) in [ 0, 1, 2 ])
        
        #Test translocations
        for i in xrange(100):
            e = d.clone()
            self.assertEquals(d.getOutOfOrderDistance(e), 0)
            self.assertEquals(d.getCircularDcjDistance(e), 0)
            e.permuteByTranslocation()
            self.assertEquals(e.getElements(), d.getElements())
            self.assertTrue(d.getOutOfOrderDistance(e) >= 0)
            self.assertTrue(d.getCircularDcjDistance(e) in [ 0, 1, 2 ])

def runGrimmMedianProblemTest(medianHistory):
    """Runs Grimm.
    """
    
    system
    pass

def runReferenceMedianProblemTest(medianHistory):
    """Runs the reference problem for a given median history
    """
    #Make adjacencies
    stubNumber = 2
    nodeNumber = len(medianHistory.getMedianGenome().getElements()) * 2 + stubNumber;
    weights = {}
    def weightFn(distance):
        assert distance >= 1
        return 1.0/distance #Hack for now
    for genome in medianHistory.getLeafGenomes():
        for node1, node2, distance in genome.getTransitiveAdjacencies():
            if (node1, node2) in weights:
                weights[(node1, node2)] += weightFn(distance)
            else:
                weights[(node1, node2)] = weightFn(distance)
    def translateLeftSideOfElementToNode(element):
        assert element != 0
        if element < 0:        
            return abs(element) * 2
        return element * 2 + 1
    def translateLeftNodeToElement(node):
        assert node >= stubNumber
        assert node < nodeNumber
        element = node / 2 
        if (node % 2) == 0:
            element *= -1
        return element
    #Now print out the 
    input = "%i\t%i\t%i\t%s" % (nodeNumber, stubNumber, len(weights.keys()), "\t".join([ "%i\t%i\t%f" % (translateLeftSideOfElementToNode(-node1), translateLeftSideOfElementToNode(node2), weights[(node1, node2)]) for (node1, node2) in weights.keys()]))
    #Command
    command = os.path.join(os.path.split(os.path.abspath(matchingAndOrdering.tests.simulatedGenome.__file__))[0], "testBin", "referenceMedianProblemTest")
    output = popenCatch(command, input)
    medianChromosome = Chromosome()
    for adjacency in output.split():
        medianChromosome.append(translateLeftNodeToElement(int(adjacency)))
    medianGenome = Genome(chromosomeNumber=0, elementNumber=0)
    medianGenome.addChromosome(medianChromosome)
    assert medianGenome.getElements() == medianHistory.getMedianGenome().getElements()
    return medianGenome

def main():
    parseSuiteTestOptions()
    sys.argv = sys.argv[:1]
    unittest.main()
        
if __name__ == '__main__':
    main()