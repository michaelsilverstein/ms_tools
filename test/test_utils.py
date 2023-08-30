from unittest import TestCase
from ms_tools.utils import check_len_n, check_n_sheets

class TestUtils(TestCase):
    
    def test_check_len_n(self):
        class TestClass:
            def __init__(self, x):
                self.x = x
        
        t = TestClass(['a', 'b'])
        
        self.assertTrue(check_len_n(t, 'x', 2))
        
        with self.assertRaises(ValueError) as context:
            check_len_n(t, 'x', 3)
        self.assertEqual('"x" must be of length 3', str(context.exception))
    
    def test_check_n_sheets(self):
        filepath = 'test/test_data/DE.Microresp1.xlsx'
        
        self.assertTrue(2, check_n_sheets(filepath, 2))
        
        wrong_n = 1
        with self.assertRaises(ValueError) as context:
            check_n_sheets(filepath, wrong_n)
        self.assertEqual(f'Excel file "{filepath}" must contain {wrong_n} sheets.', str(context.exception))