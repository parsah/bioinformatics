'''
Unit-tests for the lasso module.
'''
import unittest
from lasso.lasso import BinaryClassificationFile, ClassifiableCountDataset, \
    COLNAME_SEQUENCE, COLNAME_TARGET

TEST_MATRIX = './matrix.csv'

class TestClassifiableCountDataset(unittest.TestCase):
    
    def setUp(self):
        self.f = BinaryClassificationFile(f = TEST_MATRIX)
        self.f.parse()
        self.cdf = ClassifiableCountDataset()
        self.cdf.set_data(self.f.get_data())
    
    def test_control_is_classifiable(self):
        ''' 
        Assert the control data-set is classifiable.
        '''
        return isinstance(self.cdf.get_control(), ClassifiableCountDataset)
    
    def test_query_is_classifiable(self):
        ''' 
        Assert the query data-set is classifiable.
        '''
        return isinstance(self.cdf.get_query(), ClassifiableCountDataset)
    
    def test_control_data_is_valid(self):
        ''' 
        Assert the control data-set is not null.
        '''
        self.assertIsNotNone(self.cdf.get_control().data())
        
    def test_query_data_is_valid(self):
        '''
        Assert the query data-set is not null.
        '''
        self.assertIsNotNone(self.cdf.get_query().data())
    
    def test_query_references_one(self):
        ''' 
        Assert that query observations reference target value of 1.
        '''
        target = set(self.cdf.get_query().data()[COLNAME_TARGET])
        self.assertTrue(len(target) == 1 and 1 in target)

    def test_control_references_zero(self):
        ''' 
        Assert that control observations reference target value of 0.
        '''
        target = set(self.cdf.get_control().data()[COLNAME_TARGET])
        self.assertTrue(len(target) == 1 and 0 in target)        
    
    def test_initial_query(self):
        ''' 
        Assert that a data-set without data is initially none.
        '''
        cdf = ClassifiableCountDataset()
        self.assertIsNone(cdf.data())

class TestBinaryClassificationFile(unittest.TestCase):

    def setUp(self):
        ''' 
        Create a bare-bones object, feed a known and valid file. 
        '''
        self.f = BinaryClassificationFile(f = TEST_MATRIX)
        self.f.parse()
    
    def test_unparsed_yields_none(self):
        ''' 
        Assert that an unparsed data-set yields no data.
        '''
        f = BinaryClassificationFile(f = None)
        self.assertIsNone(f.get_data())
        
    def test_unparsed_is_unsanitized(self):
        ''' 
        Assert that an unparsed data-set is not sanitized.
        '''
        f = BinaryClassificationFile(f = None)
        self.assertFalse(f.is_sanitized())
    
    def test__data_not_none(self):
        ''' 
        Assert the parsed data-structure is not None.
        '''
        self.assertIsNotNone(self.f.get_data())
    
    def test_has_sequence_feature(self):
        ''' 
        Assert user-provided file contains a sequence vector.
        '''
        self.assertTrue(COLNAME_SEQUENCE in self.f.get_data())
    
    def test_has_target_feature(self):
        ''' 
        Assert that the user-provided file has a target vector.
        '''
        self.assertTrue(COLNAME_TARGET in self.f.get_data())
    
    def test_is_sanitized(self):
        ''' 
        Assert that the user-provided matrix file is sanitized.
        '''
        self.assertTrue(self.f.is_sanitized())
    
    def test_target_is_binary(self):
        ''' 
        Assert the target vector only contains 0 and 1
        '''
        set_target = set(self.f.get_data()[COLNAME_TARGET])
        self.assertEqual(set_target, set([0, 1]))

if __name__ == "__main__":
    unittest.main()