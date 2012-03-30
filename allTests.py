#!/usr/bin/env python

#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
import unittest

from sonLib.bioio import parseSuiteTestOptions

from matchingAndOrdering.externalTools.matchGraph.matchGraphTest import TestCase as matchGraphTest
from matchingAndOrdering.externalTools.blossom.blossomTest import TestCase as blossomTest
from matchingAndOrdering.tests.allTests import TestCase as cTests

def allSuites(): 
    return unittest.TestSuite((unittest.makeSuite(matchGraphTest, 'test'),
                                   unittest.makeSuite(blossomTest, 'test'),
                                   unittest.makeSuite(cTests, 'test')))

def main():
    parseSuiteTestOptions()
    
    suite = allSuites()
    runner = unittest.TextTestRunner()
    i = runner.run(suite)
    return len(i.failures) + len(i.errors)
        
if __name__ == '__main__':
    import sys
    sys.exit(main())
