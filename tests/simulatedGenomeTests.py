#!/usr/bin/env python

#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
import unittest
import sys
import os

from sonLib.bioio import getLogLevelString
from sonLib.bioio import parseSuiteTestOptions
from sonLib.bioio import system

from matchingAndOrdering.tests.simulatedGenome import Chromosome, Genome, MedianHistory

class TestCase(unittest.TestCase):
    def setUp(self):
        self.tempFile = os.path.join(os.getcwd(), "simulatedGenomeTempFile.txt")
        self.elementNumbers=(10, 100)
        self.chromosomeNumbers=(1, 2, 5)
        self.leafGenomeNumbers=(2, 3, 5)
        self.operationNumber = (1, 10, 100, 1000)
        self.operationType = ((True, False, False), (False, True, False), (False, False, True))
        self.replicates = 1
    
    def testReferenceAndGrimmAlgorithms(self):
        """Iterates through a list of simulation variants and prints results
        """
        headerLine = "\t".join(("elementNumber", "chromosomeNumber", "leafGenomeNumber", 
                                 "operationNumber",
                                 "doInversion", "doDcj", "doTranslocation", "replicate", 
                                 "medianDCJDistance", "medianOutOfOrderDistance", 
                                 "medianDCJDistanceForReferenceAlgorithm",
                                 "medianOutOfOrderDistanceForReferenceAlgorithm", 
                                 "medianDCJDistanceForGrimm", "medianOutOfOrderDistanceForGrimm"))
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
                                #Dump to disk
                                fileHandle = open(self.tempFile, 'w')
                                fileHandle.write(medianHistory.getLeafGenomeString())
                                fileHandle.close()
                                #Now run reference problem algorithm         
                                medianDCJDistanceForReferenceAlgorithm = 0
                                medianOutOfOrderDistanceForReferenceAlgorithm = 0
                                #Now run GRIMM
                                medianDCJDistanceForGrimm = 0
                                medianOutOfOrderDistanceForGrimm = 0
                                #Now cleanup
                                os.remove(self.tempFile)
                                
                                #Now print a report line
                                line = "\t".join([ str(i) for i in 
                                (elementNumber, chromosomeNumber, leafGenomeNumber, 
                                 operationNumber,
                                 doInversion, doDcj, doTranslocation, replicate, 
                                 medianDCJDistance, medianOutOfOrderDistance, 
                                 medianDCJDistanceForReferenceAlgorithm,
                                 medianOutOfOrderDistanceForReferenceAlgorithm, 
                                 medianDCJDistanceForGrimm, medianOutOfOrderDistanceForGrimm) ])
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

def main():
    parseSuiteTestOptions()
    sys.argv = sys.argv[:1]
    unittest.main()
        
if __name__ == '__main__':
    main()