'''
Unit-tests for the lasso module.
'''
import unittest
from lasso.lasso import ClassifiableCountDataset, COLNAME_SEQUENCE,\
    COLNAME_TARGET

TEST_MATRIX = './matrix.csv'

class TestClassifiableCountDataset(unittest.TestCase):

    def setUp(self):
        ''' 
        Create a bare-bones object, feed a known and valid file. 
        '''
        self.cdf = ClassifiableCountDataset()
        self.cdf.parse(TEST_MATRIX)
        
    def test__data_not_none(self):
        ''' 
        Assert the parsed data-structure is not None.
        '''
        self.assertIsNotNone(self.cdf.get_data())
    
    def test_has_sequence_feature(self):
        ''' 
        Assert user-provided file contains a sequence vector.
        '''
        self.assertTrue(COLNAME_SEQUENCE in self.cdf.get_data())
    
    def test_has_target_feature(self):
        ''' 
        Assert that the user-provided file has a target vector.
        '''
        self.assertTrue(COLNAME_TARGET in self.cdf.get_data())
    
    def test_is_sanitized(self):
        ''' 
        Assert that the user-provided matrix file is sanitized.
        '''
        self.assertTrue(self.cdf.is_sanitized())
        
    def test_target_is_binary(self):
        ''' 
        Assert the target vector only contains 0 and 1
        '''
        set_target = set(self.cdf.get_data()[COLNAME_TARGET])
        self.assertEqual(set_target, set([0, 1]))

if __name__ == "__main__":
    unittest.main()