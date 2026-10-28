import unittest
from postprocessing import extrapolate


class RichardsonChecks(unittest.TestCase):
    def family(self, axial=(44,88,176), power=2):
        return [dict(nt=nt,nz=nz,value=1.25+.8/nt**power+3/nz**power)
                for nt,nz in zip((128,256,512),axial)]

    def test_two_direction_error_recovers_order_and_limit(self):
        value,order=extrapolate(self.family(),'value')
        self.assertAlmostEqual(order,2,places=9)
        self.assertAlmostEqual(value,1.25,places=13)

    def test_original_rounded_family_rejects(self):
        with self.assertRaisesRegex(ValueError,'both directions'):
            extrapolate(self.family((44,90,180)),'value')

    def test_second_pair_must_also_double(self):
        with self.assertRaisesRegex(ValueError,'both directions'):
            extrapolate(self.family((38,76,154)),'value')

    def test_sorting_and_three_finest_levels(self):
        family=self.family()
        value,order=extrapolate([family[2],dict(nt=64,nz=22,value=9.),family[0],family[1]],'value')
        self.assertAlmostEqual(value,1.25,places=13)
        self.assertAlmostEqual(order,2,places=9)

    def test_nonmonotonic_and_constant_sequences_reject(self):
        for values in [(1.,2.,1.5),(1.,1.,1.)]:
            with self.subTest(values=values),self.assertRaises(ValueError):
                rows=self.family()
                for row,value in zip(rows,values):
                    row['value']=value
                extrapolate(rows,'value')

    def test_insufficient_levels_reject(self):
        with self.assertRaises(ValueError):
            extrapolate(self.family()[:2],'value')

    def test_nonfinite_rejects(self):
        rows=self.family();rows[1]['value']=float('nan')
        with self.assertRaises(ValueError):
            extrapolate(rows,'value')


if __name__=='__main__':
    unittest.main()
