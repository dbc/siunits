import unittest
import hypothesis
import siunits as si
from siunits import Unit, BaseUnit, DerivedUnit, \
    Dimension, Dn, sqrt, root


class TestDimensionBasic(unittest.TestCase):
    def test_dimensionless(self):
        d = Dimension()
        self.assertEqual(d._exponents, [0, 0, 0, 0, 0, 0, 0, 0])
        s = repr(d)
        self.assertEqual(s, 'Dimension()')
        s = str(d)
        self.assertEqual(s, '')

    def test_base_units(self):
        # meter
        d = Dimension(meter=1)
        self.assertEqual(d._exponents, [1, 0, 0, 0, 0, 0, 0, 0])
        s = repr(d)
        self.assertEqual(s, 'Dimension(meter=1)')
        s = str(d)
        self.assertEqual(s, 'm')
        # kilogram
        d = Dimension(kilogram=1)
        self.assertEqual(d._exponents, [0, 1, 0, 0, 0, 0, 0, 0])
        s = repr(d)
        self.assertEqual(s, 'Dimension(kilogram=1)')
        s = str(d)
        self.assertEqual(s, 'kg')
        # second
        d = Dimension(second=1)
        self.assertEqual(d._exponents, [0, 0, 1, 0, 0, 0, 0, 0])
        s = repr(d)
        self.assertEqual(s, 'Dimension(second=1)')
        s = str(d)
        self.assertEqual(s, 's')
        # ampere
        d = Dimension(ampere=1)
        self.assertEqual(d._exponents, [0, 0, 0, 1, 0, 0, 0, 0])
        s = repr(d)
        self.assertEqual(s, 'Dimension(ampere=1)')
        s = str(d)
        self.assertEqual(s, 'A')
        # candela
        d = Dimension(candela=1)
        self.assertEqual(d._exponents, [0, 0, 0, 0, 1, 0, 0, 0])
        s = repr(d)
        self.assertEqual(s, 'Dimension(candela=1)')
        s = str(d)
        self.assertEqual(s, 'cd')
        # kelvin
        d = Dimension(kelvin=1)
        self.assertEqual(d._exponents, [0, 0, 0, 0, 0, 1, 0, 0])
        s = repr(d)
        self.assertEqual(s, 'Dimension(kelvin=1)')
        s = str(d)
        self.assertEqual(s, 'K')
        # mole
        d = Dimension(mole=1)
        self.assertEqual(d._exponents, [0, 0, 0, 0, 0, 0, 1, 0])
        s = repr(d)
        self.assertEqual(s, 'Dimension(mole=1)')
        s = str(d)
        self.assertEqual(s, 'mol')
        # radian
        d = Dimension(radian=1)
        self.assertEqual(d._exponents, [0, 0, 0, 0, 0, 0, 0, 1])
        s = repr(d)
        self.assertEqual(s, 'Dimension(radian=1)')
        s = str(d)
        self.assertEqual(s, 'rad')


    def test_dimension_construct01(self):
        # some ways to make m^2
        d1 = Dimension([2, 0, 0, 0, 0, 0, 0, 0])
        s = str(d1)
        self.assertEqual(s, 'm^2')
        d2 = Dimension('m^2')
        self.assertEqual(str(d1), str(d2))
        d3 = Dimension(meter=2)
        self.assertEqual(str(d1), str(d3))
        d4 = Dimension(d2)
        self.assertEqual(str(d1), str(d4))
        # some 1/sec a.k.a. Hertz
        d1 = Dimension(second=-1)
        s = str(d1)
        self.assertEqual(s, 'Hz')
        # Multiple units kg*m/sec^2 a.k.a. Newton
        d1 = Dimension(kilogram=1, meter=1, second=-2)
        s = str(d1)
        self.assertEqual(s, 'N')
        s = repr(d1)
        self.assertEqual(s, 'Dimension(meter=1, kilogram=1, second=-2)')
        d2 = Dimension('kgm/s^2')
        self.assertEqual(repr(d1), repr(d2))

class TestDimensionArith(unittest.TestCase):
    def setUp(self):
        self.dm = Dimension('m')
        self.dkg = Dimension('kg')
        self.ds = Dimension('s')
        self.dA = Dimension('A')
        self.dK = Dimension('K')
        self.dmol = Dimension('mol')
        self.dcd = Dimension('cd')
        self.drad = Dimension('rad')

    def test_4bang_basic(self):
        actual = self.dm + self.dm
        self.assertEqual(Dimension([1,0,0,0,0,0,0,0]), actual)
        actual = self.dm - self.dm
        self.assertEqual(Dimension([1,0,0,0,0,0,0,0]), actual)
        dm2 = self.dm * self.dm
        self.assertEqual(Dimension([2,0,0,0,0,0,0,0]), dm2)
        actual = self.dm / self.dm
        self.assertEqual(Dimension(), actual)
        actual = dm2 // self.dm
        self.assertEqual(self.dm, actual)
        actual = self.dm / self.ds
        self.assertEqual(Dimension([1,0,-1,0,0,0,0,0]), actual)

    def test_pow01(self):
        actual = pow(self.dkg, 2)
        self.assertEqual(Dimension([0,2,0,0,0,0,0,0]), actual)
        actual = pow(self.dkg, 3)
        self.assertEqual(Dimension([0,3,0,0,0,0,0,0]), actual)
        with self.assertRaises(TypeError):
            actual = pow(self.dkg, 0.5)
        with self.assertRaises(TypeError):
            actual = pow(self.dkg, -1)

    def test_root01(self):
        kg4 = pow(self.dkg, 4)
        actual = kg4.root()
        self.assertEqual(Dimension('kg^2'), actual)
        a6 = pow(self.dA, 6)
        actual = a6.root(3)
        self.assertEqual(Dimension('A^2'), actual)
        with self.assertRaises(ValueError):
            actual = kg4.root(-2)
        with self.assertRaises(ValueError):
            actual = a6.root(4)

    def test_inplace01(self):
        t = self.dK * self.dmol
        self.assertEqual(Dimension(kelvin=1, mole=1), t)
        t *= Dimension(self.dcd)
        self.assertEqual(Dimension(kelvin=1, mole=1, candela=1), t)

    def test_compare01(self):
        actual = self.dkg == Dimension(kilogram=1)
        self.assertTrue(actual)
        actual = self.dkg != Dimension(kilogram=1)
        self.assertFalse(actual)
        actual = self.dkg != Dimension(kelvin=1)
        self.assertTrue(actual)
        actual = self.dkg == Dimension(mole=1)
        self.assertFalse(actual)
        with self.assertRaises(TypeError):
            actual = self.dkg < self.dA
        with self.assertRaises(TypeError):
            actual = self.dkg > self.dkg
        with self.assertRaises(TypeError):
            actual = self.dmol <= self.dcd
        with self.assertRaises(TypeError):
            actual = self.drad < self.ds


class TestDerivedFactoring(unittest.TestCase):
    def setUp(self):
        self.uC = Unit.of['C']
        self.uN = Unit.of['N']
        self.uF = Unit.of['F']
        self.uWb = Unit.of['Wb']
        self.uS = Unit.of['S']
        self.ucd = Unit.of['cd']
        self.uW = Unit.of['W']

    def test_occurs_in01(self):
        #coulomb = Dimension(second=1, ampere=1)
        coulomb = Unit.of['C'].dimension
        self.assertEqual(self.uC.occurs_in(coulomb._exponents), 1)
        self.assertEqual(
            self.uC.occurs_in((coulomb * coulomb)._exponents), 2)
        newton = Unit.of['N'].dimension
        self.assertEqual(
            self.uC.occurs_in(newton._exponents), 0)
        nm = Unit.of['N'].dimension * Unit.of['m'].dimension
        self.assertEqual(
            self.uN.occurs_in(nm._exponents), 1)
        with self.assertRaises(AttributeError):
            self.assertEqual(
                self.um.occurs_in(nm._exponents), 0)
        nm2 = nm * nm
        self.assertEqual(
            self.uN.occurs_in(nm2._exponents), 2)
        self.assertEqual(
            self.uN.occurs_in([1,1,4,0,0,0,0]), 0)
        self.assertEqual(
            self.uN.occurs_in([1,1,-4,0,0,0,0]), 1)

    def test_extract_from01(self):
        foo = self.uC.dimension * self.uS.dimension
        expect = [-2, -1, 4, 3, 0, 0, 0, 0]
        self.assertEqual(foo._exponents, expect)
        self.assertEqual(
            self.uS.occurs_in(foo._exponents), 1)
        self.assertEqual(
            self.uC.occurs_in(foo._exponents), 3)
        residue = self.uS.extract_from(expect)
        # expect should not be modified.
        self.assertEqual(expect, [-2, -1, 4, 3, 0, 0, 0, 0])
        # residue should be a Coulomb
        self.assertEqual(residue, self.uC.dimension._exponents)

    def test_extract_from03(self):
        eut = [-2, -1, 4, 3, 1, 0, 0, 0]
        residue = Unit.of['S'].extract_from(eut)
        self.assertEqual(residue, [0, 0, 1, 1, 1, 0, 0, 0])

    def test_factor01(self):
        foo = self.uW.dimension * self.uN.dimension
        factored = DerivedUnit.factor(foo._exponents)
        self.assertEqual(factored, [(self.uW, 1), (self.uN, 1)])


if __name__ == '__main__':
    unittest.main()