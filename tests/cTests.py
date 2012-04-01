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

class TestCase(unittest.TestCase):
    def testCuTest(self):
        system("matchingAndOrderingTests %s" % getLogLevelString())
    
def main():
    parseSuiteTestOptions()
    sys.argv = sys.argv[:1]
    unittest.main()
        
if __name__ == '__main__':
    main()
