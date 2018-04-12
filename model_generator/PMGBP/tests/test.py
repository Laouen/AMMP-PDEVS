import unittest
from SBMLParser import SBMLParser


class TestSBMLParser(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        super(TestSBMLParser, self).__init__(*args, **kwargs)
        self.sbml_parser = SBMLParser()

    def test_parse_gene_association(self):
        self.assertEqual(self.sbml_parser.parse_gene_association(['(b1)']), ['b1'])
        self.assertEqual(self.sbml_parser.parse_gene_association(['(b1 or b2)']), ['b1', 'b2'])
        self.assertEqual(self.sbml_parser.parse_gene_association(['(b1 or b2 or b3)']), ['b1', 'b2', 'b3'])
        self.assertEqual(self.sbml_parser.parse_gene_association(['(b1 and b2)']), ['b1-b2'])
        self.assertEqual(self.sbml_parser.parse_gene_association(['(b1 and b2 and b3)']), ['b1-b2-b3'])
        self.assertEqual(self.sbml_parser.parse_gene_association(['(b1 and b2 or b3)']), ['b1-b2', 'b3'])
        self.assertEqual(self.sbml_parser.parse_gene_association(['(b1 and (b2 or b3))']), ['b1-b2', 'b1-b3'])
        self.assertEqual(self.sbml_parser.parse_gene_association(['(b1 or b2 and b3)']), ['b1-b3', 'b2-b3'])
        self.assertEqual(self.sbml_parser.parse_gene_association(['(b1 or (b2 and b3))']), ['b1', 'b2-b3'])
        self.assertEqual(self.sbml_parser.parse_gene_association(['(b1 and b2 and b3 and b4)']), ['b1-b2-b3-b4'])
        self.assertEqual(self.sbml_parser.parse_gene_association(['(b1 or b2 and b3 and b4)']), ['b1-b3-b4', 'b2-b3-b4'])
        self.assertEqual(self.sbml_parser.parse_gene_association(['(b1 and b2 or b3 and b4)']), ['b1-b2-b4', 'b3-b4'])
        self.assertEqual(self.sbml_parser.parse_gene_association(['(b1 and b2 and b3 or b4)']), ['b1-b2-b3', 'b4'])
        self.assertEqual(self.sbml_parser.parse_gene_association(['(b1 or b2 and b3 or b4)']), ['b1-b3', 'b2-b3', 'b4'])
        self.assertEqual(self.sbml_parser.parse_gene_association(['(b1 and b2 or b3 or b4)']), ['b1-b2', 'b3', 'b4'])
        self.assertEqual(self.sbml_parser.parse_gene_association(['(b1 or b2 or b3 or b4)']), ['b1', 'b2', 'b3', 'b4'])


if __name__ == '__main__':
    unittest.main()
