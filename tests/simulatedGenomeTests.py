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

from matchingAndOrdering.tests.simulatedGenome import Chromosome, Genome, MedianHistory, weightFn

class TestCase(unittest.TestCase):
    def setUp(self):
        self.elementNumbers=(250,)
        self.chromosomeNumbers=(1,) # 2, 5)
        self.leafGenomeNumbers=(5,) # 5, 10)
        self.operationNumber = (10,) #(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50) #(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100, 150, 200, 250)
        self.operationType = ((False, True, False, False, True),) #(True, False, False),) #(False, True, False), (False, False, True))
        self.replicates = 3
        self.greedyIterations = (100,)
        self.theta = (0.1,) #0.2,0.5,0.9)
    
    def testReferenceAndAsMedianAlgorithms(self):
        """Iterates through a list of simulation variants and prints results
        """
        headerLine = "\t".join(("elementNumber", "chromosomeNumber", "leafGenomeNumber", 
                                 "operationNumber",
                                 "totalOperationNumber",
                                 "doInversion", "doShortInversion", "doDcj", "doTranslocation", "doShortTranslocation", 
                                 "greedyIterations",
                                 "theta",
                                 "replicate", 
                                 "medianDCJDistance", "medianOutOfOrderDistance", 
                                 "weightedMedianOutOfOrderDistance", 
                                 "medianDCJDistanceForReferenceAlgorithm",
                                 "medianOutOfOrderDistanceForReferenceAlgorithm", 
                                 "weightedMedianOutOfOrderDistanceForReferenceAlgorithm",
                                 "dCJDistanceForReferenceAlgorithmFromMedian",
                                 "outOfOrderDistanceForReferenceAlgorithmFromMedian",
                                 "weightedOutOfOrderDistanceForReferenceAlgorithmFromMedian",
                                 "medianDCJDistanceForAsMedian", 
                                 "medianOutOfOrderDistanceForAsMedian", 
                                 "weightedMedianOutOfOrderDistanceForAsMedian", 
                                 "dCJDistanceForAsMedianFromMedian",
                                 "outOfOrderDistanceForAsMedianFromMedian",
                                 "weightedOutOfOrderDistanceForAsMedianFromMedian",
                                 "medianGenomeForReferenceAlgorithm", 
                                 "medianGenomeForAsMedian"))
        if getLogLevelString() in  ("DEBUG", "INFO" ):
            print headerLine
        for elementNumber in self.elementNumbers:
            for chromosomeNumber in self.chromosomeNumbers:
                for leafGenomeNumber in self.leafGenomeNumbers:
                    for operationNumber in self.operationNumber:
                        for doInversion, doShortInversion, doDcj, doTranslocation, doShortTranslocation in self.operationType:
                            for greedyIterations in self.greedyIterations:
                                for theta in self.theta:
                                    for replicate in xrange(self.replicates):
                                        medianHistory = MedianHistory(Genome(elementNumber=elementNumber, chromosomeNumber=chromosomeNumber), leafGenomeNumber=leafGenomeNumber)
                                        medianHistory.permuteLeafGenomes(operationNumber=operationNumber, doInversion=doInversion, doDcj=doDcj, doTranslocation=doTranslocation,
                                                                         doShortInversion=doShortInversion, doShortTranslocation=doShortTranslocation)
                                        medianDCJDistance = medianHistory.getMedianDcjDistance(medianHistory.getMedianGenome())
                                        medianOutOfOrderDistance = medianHistory.getMedianOutOfOrderDistance(medianHistory.getMedianGenome())
                                        weightedMedianOutOfOrderDistance = medianHistory.getWeightedMedianOutOfOrderDistance(medianHistory.getMedianGenome(), theta=theta)
                                        #Now run reference problem algorithm   
                                        referenceProblemMedianGenome = runReferenceMedianProblemTest(medianHistory, greedyIterations, theta)      
                                        medianDCJDistanceForReferenceAlgorithm = medianHistory.getMedianDcjDistance(referenceProblemMedianGenome)
                                        medianOutOfOrderDistanceForReferenceAlgorithm = medianHistory.getMedianOutOfOrderDistance(referenceProblemMedianGenome)
                                        weightedMedianOutOfOrderDistanceForReferenceAlgorithm = medianHistory.getWeightedMedianOutOfOrderDistance(referenceProblemMedianGenome, theta=theta)
                                        dCJDistanceForReferenceAlgorithmFromMedian = medianHistory.getMedianGenome().getCircularDcjDistance(referenceProblemMedianGenome)
                                        outOfOrderDistanceForReferenceAlgorithmFromMedian = medianHistory.getMedianGenome().getOutOfOrderDistance(referenceProblemMedianGenome)
                                        weightedOutOfOrderDistanceForReferenceAlgorithmFromMedian = medianHistory.getMedianGenome().getWeightedOutOfOrderDistance(referenceProblemMedianGenome, theta=theta)
                                        totalOperationNumber = operationNumber * len([  i for i in (doInversion, doShortInversion, doDcj, doTranslocation, doShortTranslocation) if i == True ])
                                        #Biomedian comparison turned off
                                        if False and leafGenomeNumber == 3 and doDcj == False and float(totalOperationNumber) / elementNumber <= 0.5:
                                            asMedianProblemMedianGenome = runAsMedianMedianProblemTest(medianHistory)
                                            medianDCJDistanceForAsMedian = medianHistory.getMedianDcjDistance(asMedianProblemMedianGenome)
                                            medianOutOfOrderDistanceForAsMedian = medianHistory.getMedianOutOfOrderDistance(asMedianProblemMedianGenome)
                                            weightedMedianOutOfOrderDistanceForAsMedian = medianHistory.getWeightedMedianOutOfOrderDistance(asMedianProblemMedianGenome, theta=theta)
                                            dCJDistanceForAsMedianFromMedian = medianHistory.getMedianGenome().getCircularDcjDistance(asMedianProblemMedianGenome)
                                            outOfOrderDistanceForAsMedianFromMedian = medianHistory.getMedianGenome().getOutOfOrderDistance(asMedianProblemMedianGenome)
                                            weightedOutOfOrderDistanceForAsMedianFromMedian = medianHistory.getMedianGenome().getWeightedOutOfOrderDistance(asMedianProblemMedianGenome, theta=theta)
                                        else:
                                            asMedianProblemMedianGenome = "n/a"
                                            medianDCJDistanceForAsMedian = "n/a"
                                            medianOutOfOrderDistanceForAsMedian = "n/a"
                                            weightedMedianOutOfOrderDistanceForAsMedian = "n/a"
                                            dCJDistanceForAsMedianFromMedian = "n/a"
                                            outOfOrderDistanceForAsMedianFromMedian = "n/a"
                                            weightedOutOfOrderDistanceForAsMedianFromMedian = "n/a"
                                        #Now prepare line to print
                                        line = "\t".join([ str(i) for i in 
                                        (elementNumber, chromosomeNumber, leafGenomeNumber, 
                                         operationNumber,
                                         totalOperationNumber,
                                         doInversion, doShortInversion, doDcj, doTranslocation, doShortTranslocation,
                                         greedyIterations,
                                         theta,
                                         replicate, 
                                         medianDCJDistance, medianOutOfOrderDistance, 
                                         weightedMedianOutOfOrderDistance,
                                         medianDCJDistanceForReferenceAlgorithm,
                                         medianOutOfOrderDistanceForReferenceAlgorithm, 
                                         weightedMedianOutOfOrderDistanceForReferenceAlgorithm, 
                                         dCJDistanceForReferenceAlgorithmFromMedian,
                                         outOfOrderDistanceForReferenceAlgorithmFromMedian,
                                         weightedOutOfOrderDistanceForReferenceAlgorithmFromMedian,
                                         medianDCJDistanceForAsMedian, 
                                         medianOutOfOrderDistanceForAsMedian,
                                         weightedMedianOutOfOrderDistanceForAsMedian,
                                         dCJDistanceForAsMedianFromMedian,
                                         outOfOrderDistanceForAsMedianFromMedian,
                                         weightedOutOfOrderDistanceForAsMedianFromMedian,
                                         "'%s'" % str(referenceProblemMedianGenome),
                                         "'%s'" % str(asMedianProblemMedianGenome)) ])
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

def runAsMedianMedianProblemTest(medianHistory):
    """Runs AsMedian, requires to be installed. I got it from:
    https://sites.google.com/site/andrewweixu/Home/software/asmedian
    """
    #Dump to disk
    tempFile = os.path.join(os.getcwd(), "simulatedGenomeTempFile.txt")
    fileHandle = open(tempFile, 'w')
    fileHandle.write(medianHistory.getLeafGenomeString())
    fileHandle.close()
    #-cp /Users/benedictpaten/Desktop/ASMedian-1.0 
    popenCatch("java BIOMedian %s" % tempFile)
    os.remove(tempFile)
    #Parse in
    fileHandle = open(tempFile + ".rst", 'r')
    input = fileHandle.readlines()
    fileHandle.close()
    os.remove(tempFile + ".rst")
    asMedianMedianGenome = Genome(chromosomeNumber=0, elementNumber=0)
    for line in input[1:]:
        if line[0] == '#':
            break
        asMedianChromosome = Chromosome()
        for element in line.split()[1:]:
            asMedianChromosome.append(int(element))
        asMedianMedianGenome.addChromosome(asMedianChromosome)
    return asMedianMedianGenome

def runReferenceMedianProblemTest(medianHistory, greedyIterations,theta):
    """Runs the reference problem for a given median history
    """
    #Make adjacencies
    stubNumber = 2
    nodeNumber = len(medianHistory.getMedianGenome().getElements()) * 2 + stubNumber;
    weights = {}
    for genome in medianHistory.getLeafGenomes():
        for node1, node2, distance in genome.getTransitiveAdjacencies():
            if (node1, node2) in weights:
                weights[(node1, node2)] += weightFn(distance, theta)
            else:
                weights[(node1, node2)] = weightFn(distance, theta)
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
    input = "%i\t%i\t%i\t%i\t%s" % (greedyIterations, nodeNumber, stubNumber, len(weights.keys()), "\t".join([ "%i\t%i\t%f" % (translateLeftSideOfElementToNode(-node1), translateLeftSideOfElementToNode(node2), weights[(node1, node2)]) for (node1, node2) in weights.keys()]))
    #Command
    command = os.path.join(os.path.split(os.path.abspath(matchingAndOrdering.tests.simulatedGenome.__file__))[0], "testBin", "referenceMedianProblemTest2")
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