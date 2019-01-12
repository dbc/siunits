"""siunit - A module to support dimensioned arithmetic using the SI system.

Dimensioned numbers are instances of siunit.Dn().  Arithmetic between Dn's
with incompatible units raises TypeError.  Arithmetic between Dn's with
compatible units produces a result with appropriate units.  For example, ::

  >>> m=Dn('3kg')
  >>> a=Dn('2m/s^2')
  >>> f=m*a
  >>> print(f)
  6.0N
"""

# MIT License

# Copyright (c) 2019 David B. Curtis

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

from numbers import Number
from math import sqrt as math_sqrt
import re
try:
    from numpy import array
except ImportError:
    pass


class Unit:
    """Base class for base unit and derived unit definitions.  Not directly
    instantiated.

    :param name: Full name of the unit.
    :param abbreviation: The official unit abbreviation.  Used for display
        and parsing.
    :param quantifies: A string (no white space allowed) that describes the
        quantity measured by the unit.
    """

    of = {}
    """Dictionary of all defined units.  Each Unit instance will be inserted
    under three keys: name, abbreviation, and quantifies."""

    _all = set()
    """Set of all defined units."""

    _measureable = dict((quantifies, unit_index)
        for (unit_index, quantifies) in enumerate(
            ['length', 'mass', 'time', 'electric_current', 'luminousity',
                'temperature', 'amount_of_substance', 'angle']))
    """Defines the entities quantified by BaseUnits. Also defines the order
    in which they appear in Dimension() exponent vectors."""

    def __init__(self, name, abbreviation, quantifies, display_order):
        self.name = name
        self.abbreviation = abbreviation
        self.quantifies = quantifies
        self.display_order = display_order
        self.of[self._irredundant(name)] = self
        self.of[self._irredundant(abbreviation)] = self
        self.of[self._irredundant(quantifies)] = self
        self._all.add(self)

    def _irredundant(self, ident):
        """Validates identifier as previously unused."""
        if ident in self.of:
            raise ValueError(' '.join([ident, 'already defined.']))
        else:
            return ident

    def __repr__(self):
        return ''.join(
            [self.__class__.__name__, '(', self.reprvals(), ')'])


class BaseUnit(Unit):
    """Defines an SI base unit.

    :param name: Full name of the unit.
    :param abbreviation: The official unit abbreviation.  Used for display
        and parsing.
    :param quantifies: A string (no white space allowed) that describes the
        quantity measured by the unit.
    :param display_order: A integer that controls the print-out order
        for this unit for the __str__ method..  Lower display_order units
        print first.
    """
    named = {}
    index_order = [None] * len(Unit._measureable)

    def __init__(self, name, abbreviation, quantifies, display_order):
        unit_index = self._measureable[quantifies]
        super().__init__(name, abbreviation, quantifies, display_order)
        self.unit_index = unit_index
        self.named[name] = self
        self.index_order[Unit._measureable[self.quantifies]] = self
        self._dimension = Dimension(**{self.name:1})

    def reprvals(self):
        return ', '.join([repr(x) for x in [self.name, self.abbreviation,
            self.quantifies, self.display_order]])

    @property
    def dimension(self):
        """Returns an instance of Dimension() representing this unit."""
        return self._dimension

    def gloss(self):
        """Glossary text for this unit."""
        return [self.name, self.abbreviation, self.quantifies, 'SI base unit']


class DerivedUnit(Unit):
    """Defines an SI derived unit.

    :param name: Full name of the unit.
    :param abbreviation: The official unit abbreviation.  Used for display
        and parsing.
    :param quantifies: A string (no white space allowed) that describes the
        quantity measured by the unit.
    :param display_order: A integer that controls the print-out order
        for this unit.  Lower display_order units print first.
    :param dimension: An instance of Dimension() that defines the derivation
        of this unit.
    """

    with_basis = {}
    """Dictionary of all DerivedUnit instances, indexed by a Dimension
    vector in tuple() form."""

    _factor_order = []
    """List of DerivedUnit instances in the prefered order for factoring
    out from complex dimension exponent vectors."""

    def __init__(self, name, abbreviation, quantifies, display_order,
            dimension):
        super().__init__(name, abbreviation, quantifies, display_order)
        self.dimension = Dimension(dimension)
        self.with_basis[tuple(self.dimension._exponents)] = self

    def reprvals(self):
        return ', '.join([repr(x) for x in [self.name, self.abbreviation,
            self.quantifies, self.display_order, self.dimension]])

    @classmethod
    def factor_order(cls):
        return cls._factor_order

    @classmethod
    def set_factor_order(cls, derived_units):
        """Sets the order for factoring out instances of DerivedUnit
        from a complex Dimension exponent vector.

        :param derived_units: A list of DerivedUnit instances, in order
            of decreasing preference.
        """
        l = []
        for du in derived_units:
            if du in l:
                raise ValueError(' '.join([du.name,
                    'duplicated in factor_order.']))
            if not isinstance(du, cls):
                raise ValueError(' '.join([repr(du),
                    'is not an instance of',
                    cls.__name__]))
            l.append(du)
        cls._factor_order = l

    def occurs_in(self, other_exponents):
        """Count the numer of times this dimension occurs in some unit
        exponent vector.

        :param other_exponents: An exponent vector.
        :returns: A count of the number of occurances.
        """
        # This is not smart enough to tease appart units where
        # one derived unit contributes a positive exponent and
        # another contributes a negative exponent.
        occurances = []
        for self_exp, other_exp in zip(
                self.dimension._exponents, other_exponents):
            if self_exp == 0:
                continue  # Not relevant.
            elif self_exp * other_exp < 0:
                return 0  # Different signs.
            elif abs(self_exp) <= abs(other_exp):
                occurances.append(other_exp // self_exp)
            else:
                return 0  # Too many.
        return min(occurances)

    def extract_from(self, exponents, count=None):
        """Extracts *count* occurances of *self's* exponent vector from
        the vector *exponents*.

        :param exponents: An exponent vector to be reduced.
        :param count: Optional.  Number of occurances to extract. Default=1.
        :returns: An updated (reduced) exponent vector.
        """
        count = 1 if count is None else int(count)
        return [y - (count * x)
            for x, y in zip(self.dimension._exponents, exponents)]

    @classmethod
    def factor(cls, exponents):
        """Factors a list of exponents into list of (dimension, count) tuples.
        This prepares an exponent list for display in maximally-factored
        form.

        :param exponents: An exponent vector.
        :returns: List of tuples of the form: (Dimension, exponent)
        """
        factored = []
        residue = exponents
        for derived_unit in cls.factor_order():
            occurances = derived_unit.occurs_in(residue)
            if occurances:
                #print ('du:', derived_unit.name)
                factored.append((derived_unit, occurances))
                #print('r before:', residue)
                residue = derived_unit.extract_from(residue, occurances)
                #print('r after:', residue)
        dims = [(d, residue[d.unit_index]) for d in BaseUnit.index_order
            if residue[d.unit_index] != 0]
        factored.extend(dims)
        return factored

    def gloss(self):
        """Glossary text for this unit."""
        return [self.name, self.abbreviation, self.quantifies,
            self.dimension.basis()]


class Dimension:
    """The computed dimensions of a number.  The dimensions are carried
    internally as a vector of exponents in a canonical order.  The vector
    contains one integer per measureable, representing the exponent on a
    unit, in the order defined by Unit._measureable.
    For example: ::

      m*kg is  [1, 1,  0, 0, 0, 0, 0, 0]
      m^2  is  [2, 0,  0, 0, 0, 0, 0, 0]
      m/s^2 is [1, 0, -2, 0, 0, 0, 0, 0]

    :param unit_spec: Can be one of:

        - Dimension() instance, which constructs a copy.
        - A string, which will be parsed, potentially raising an error.
        - An interable of integer-ish things, which will be interpreted as an
          exponent vector.
        - Omitted.

    :param kwargs: Any BaseUnit.name can be a kwarg parameter.  If unit_spec
        is omitted, a Dimension is constructed from kwargs.  If neither
        unit_spec nor any kwargs are supplied, a "dimensionless" Dimension()
        is constructed.
    """
    # Note: _parser_re is created on first attempt to parse a unit string.
    # If any new units are added after the regular expression is cached,
    # the cached pattern will be invalid. YAGN invalidation?
    _parser_re = None
    "Cached re match pattern use to parse unit strings."

    _preferred = dict()
    """Dictionary of display strings keyed by exponent vector.
    Used to trap out preferred, ie: "natural" display strings for
    the __str__ function.
    """

    def __init__(self, unit_spec=None, **kwargs):
        # If unit_spec is an instance of Dimension, construct copy.
        if isinstance(unit_spec, Dimension):
            self._exponents = unit_spec._exponents
            return
        # If unit_spec is a string, then try to parse it.
        if isinstance(unit_spec, str):
            self._exponents = self.parse(unit_spec)._exponents
            return
        # If unit_spec is not none, expect  a list of integers.
        if unit_spec is not None:
            if len(unit_spec) != len(Unit._measureable):
                raise ValueError(' '.join(['Length of exponent list is',
                    str(len(unit_spec)),
                        'but', str(len(Unit._measureable)), 'required.']))
            else:
                self._exponents = [int(x) for x in unit_spec]
            return
        # Create an empty exponent vector to fill with keywords params, or if
        # there are none, it will default to being a dimensionless quantity.
        self._exponents = [0] * len(Unit._measureable)
        for unit in kwargs:
            try:
                u = BaseUnit.named[unit]
            except IndexError:
                raise ValueError(' '.join([unit, 'is not a BaseUnit name.']))
            else:
                self._exponents[u.unit_index] = int(kwargs[unit])

    def __repr__(self):
        params = ', '.join(['='.join((nm,str(val)))
            for nm, val in [(u.name, self._exponents[u.unit_index])
                for u in BaseUnit.index_order] if val != 0])
        return ''.join([self.__class__.__name__, '(', params, ')'])

    @classmethod
    def prefer(cls, s, delete=None):
        """Add/update a preferred display string.

        :param s: A unit string.  The exponent vector that it represents
            will always be display as *s*.
        :param delete: Optional.  If truthy, *s* is deleted from the
            preferences dictionary.
        """
        dim = cls.parse(s)  # May raise ValueError.  Let caller handle.
        v = tuple(dim._exponents)
        if bool(delete):
            del cls._preferred[v]
        else:
            cls._preferred[v] = s

    def __str__(self):
        # See if there is a preferred display override.
        u = self._preferred.get(tuple(self._exponents), None)
        if u:
            return u
        # Look for an exact match to a derived unit for an easy win.
        u = DerivedUnit.with_basis.get(tuple(self._exponents), None)
        if u:
            return u.abbreviation
        # Convert to list of (Dimension, count) tuples, factoring out derived
        # units.
        dims = DerivedUnit.factor(self._exponents)
        # Prep for display.
        numerator = sorted([(d, v) for d, v in dims if v > 0],
            key=lambda tup: tup[0].display_order)
        numerator = ''.join([
            '^'.join([d.abbreviation, str(v)]) if v > 1 else d.abbreviation
                for d, v in numerator])
        denominator = sorted([(d, v) for d, v in dims if v < 0],
            key=lambda tup: tup[0].display_order)
        denominator = ''.join([
            '^'.join([d.abbreviation, str(abs(v))]) if v < -1 else d.abbreviation
                for d, v in denominator])
        if numerator and denominator:
            return '/'.join([numerator, denominator])
        elif denominator:
            return '/'.join(['1', denominator])
        else:
            return numerator

    def basis(self):
        """Return a human-readabble description of the basis for
        this dimension.
        """
        bases = [(d, self._exponents[d.unit_index])
            for d in BaseUnit.index_order
                if self._exponents[d.unit_index] != 0]
        return ' '.join(
            [''.join([d.name, '^' + str(exp) if exp != 1 else ''])
                for d, exp in bases if exp != 0])

    @classmethod
    def parse(cls, s):
        """Parse *s* as unit abbreviations and exponents, potentially in
        simple fraction notation.

        :param s: Unit specification as a string. No white space.
        :returns: Dimension, or raises ValueError.
        """
        # This parser is too ugly to live long.  But it seems to work.  The
        # better strategy is probably to replace this kludgery with a proper
        # recursive-descent parser.  But, for now.... meh.
        patt = cls._parser_pattern()
        # See if this is in fraction notation.
        f = s.split('/')
        if len(f) > 2:
            raise ValueError('Only simple fractions of units accepted.')
        numerator_str = f.pop(0)
        denominator_str = f.pop(0) if f else None
        numl = []
        while numerator_str:
            m = patt.match(numerator_str)
            if not m:
                raise ValueError(' '.join(['Syntax error in:', numerator_str]))
            numl.append(m.group(0))
            numerator_str = numerator_str[m.end(0):]
        numl2 = []
        while numl:
            d = numl.pop(0)
            if numl and numl[0][0] == '^':
                n = numl.pop(0)
                n = int(n[1:])
            else:
                n = 1
            numl2.append((d, n))
        # Start construction with a dimensionless Dimension()
        rslt = Dimension()
        for abbr, exp in numl2:
            d = Unit.of[abbr].dimension
            if exp < 0:
                d = Dimension()/d
                exp = -exp
            d = pow(d, exp)
            rslt *= d
        if denominator_str:
            denom = cls.parse(denominator_str)
            rslt = rslt / denom
        return rslt

    @classmethod
    def _parser_pattern(cls):
        "Construct the re pattern from abbreviations for defined units."
        if cls._parser_re is None:
            # re needs to have longest match strings first, so sort
            # abbreviations by decreasing length.
            tokens = sorted(set([u.abbreviation for u in Unit.of.values()]),
                key=lambda s:len(s), reverse=True)
            # In addition to unit abbreviations, we need to look for unit
            # exponents of the form:
            # up-caret, optional minus-sign, numbers,
            tokens.append(r'\^-?\d+')  # *** NO CARRIER
            pattern = '|'.join(tokens)
            cls._parser_re = re.compile(''.join(['(', pattern, ')']))
        return cls._parser_re

    # adding/subtracting mismatched units is invalid.
    def __add__(self, other):
        if self == other:
            return self
        raise TypeError('Adding mismatched units.')

    def __sub__(self, other):
        if self == other:
            return self
        raise TypeError('Subtracting mismatched units.')

    # multiply units by adding exponents.
    def __mul__(self, other):
        return self.__class__(
            [x + y for x, y in zip(self._exponents, other._exponents)])

    def __matmul__(self, other):
        return self.__mul__(other)

    def _div(self, other):
        return self.__class__(
            [x - y for x, y in zip(self._exponents, other._exponents)])

    def __truediv__(self, other):
        return self._div(other)

    def __floordiv__(self, other):
        return self._div(other)

    # Comparision operations only make sense for values of same units.
    # The __eq__ and __ne__ operators are defined so that they can be
    # used to check for equivalence of units.  The other comparison
    # operators raise.
    def __eq__(self, other):
        return min([self_exp == other_exp for self_exp, other_exp
                in zip(self._exponents, other._exponents)])

    def __ne__(self, other):
        return not self == other

    _invalid_comparison_msg = 'Dimensions can only be compared for equality.'
    def __lt__(self, other):
        raise TypeError(self._invalid_comparison_msg)

    def __gt__(self, other):
        raise TypeError(self._invalid_comparison_msg)

    def __le__(self, other):
        raise TypeError(self._invalid_comparison_msg)

    def __ge__(self, other):
        raise TypeError(self._invalid_comparison_msg)

    # For pow(), other must be a positive iteger.  To raise a Dimension
    # to a power, multiply the current exponent by other.
    def __pow__(self, other):
        o = int(other)
        if o < 0 or o != other:
            raise TypeError(
                'pow() only supported for positive integer exponentiation.')
        exp = [e * o for e in self._exponents]
        return self.__class__(exp)

    def __imul__(self, other):
        self._exponents = [x + y
            for x, y in zip(self._exponents, other._exponents)]
        return self

    def root(self, n=None):
        """Implement n-th root for dimensions.  All units must be evenly
        divisible by n.

        :param n: N-th root to take.  Default==2.
        """
        n = 2 if n is None else int(n)
        if n < 0:
            raise ValueError('Can only take positive roots of dimensions.')
        fail = max([e % n for e in self._exponents])
        if fail:
            raise ValueError(' '.join(['root(n) requires'
               'all dimension exponents to be divisble by n.']))
        return self.__class__([e // n for e in self._exponents])


# Define the base units for the SI system.
base_units = [
    BaseUnit('meter',    'm',   'length',              20),
    BaseUnit('kilogram', 'kg',  'mass',                10),
    BaseUnit('second',   's',   'time',                30),
    BaseUnit('ampere',   'A',   'electric_current',    40),
    BaseUnit('kelvin',   'K',   'temperature',         50),
    BaseUnit('mole',     'mol', 'amount_of_substance', 60),
    BaseUnit('candela',  'cd',  'luminousity',         70),
    BaseUnit('radian',   'rad', 'angle',               80),
]
"""The canonical list of SI base units."""

# Import-time self-check: Make sure that all base units are accounted for.
assert(len([x for x in BaseUnit.index_order if x is None]) == 0)


# Define some derived units for the SI system.
derived_units = [
    DerivedUnit('coulomb', 'C', 'charge', 90, Dimension(second=1, ampere=1)),
    DerivedUnit('hertz', 'Hz', 'frequency', 100, Dimension(second=-1)),
    DerivedUnit('newton', 'N', 'force', 110,
        Dimension(kilogram=1, meter=1, second=-2)),
    DerivedUnit('pascal', 'Pa', 'pressure', 130,
        Dimension(kilogram=1, meter=-1, second=-2)),
    DerivedUnit('joule', 'J', 'energy', 130,
        Dimension(kilogram=1, meter=2, second=-2)),
    DerivedUnit('watt', 'W', 'power', 140,
        Dimension(kilogram=1, meter=2, second=-3)),
    DerivedUnit('volt', 'V', 'electromotive_force', 150,
        Dimension(kilogram=1, meter=2, second=-3, ampere=-1)),
    DerivedUnit('farad', 'F', 'capacitance', 160,
        Dimension(kilogram=-1, meter=-2, second=4, ampere=2)),
    DerivedUnit('ohm', 'Ohm', 'resistance', 170,
        Dimension(kilogram=1, meter=2, second=-3, ampere=-2)),
    DerivedUnit('siemens', 'S', 'conductance', 180,
        Dimension(kilogram=-1, meter=-2, second=3, ampere=2)),
    DerivedUnit('weber', 'Wb', 'magnetic_flux', 190,
        Dimension(kilogram=1, meter=2, second=-2, ampere=-1)),
    DerivedUnit('tesla', 'T', 'magnetic_flux_density', 200,
        Dimension(kilogram=1, second=-2, ampere=-1)),
    DerivedUnit('henry', 'H', 'inductance', 210,
        Dimension(kilogram=1, meter=2, second=-2, ampere=-2)),
    DerivedUnit('lux', 'lx', 'illuminance', 220,
        Dimension(meter=-2, candela=1)),
    DerivedUnit('katal', 'kat', 'catalytic_activity', 230,
        Dimension(mole=1, second=-1)),
]
"""SI derived units."""

# Set up the order to factor out derived units.
DerivedUnit.set_factor_order([
    Unit.of['W'],
    Unit.of['N'],
])

# Set up display traps for things that come out strangely using the canonical
# conversion in __str__().
#Dimension.prefer('Ohm/m^2')



class Dn:
    """Dimensioned number.  The numeric parameter accepts any type, and passes
    through arithmetic operations.  This allows smooth operation with complex
    numbers, numpy arrays, etc.

    :param n: A number-ish thing .  Can be a numeric value, a Numpy array, or
        if a string an attempt will be made to parse it as a number with
        units.  If the string is *only* units, it constructs a dimensioned
        value of 1.
    :param units: Omitted when *n* is a string, as the units are taken from
        parsing the string.  Otherwise, if *n* is not a string and *units*
        is omitted, a dimensionless instance of Dimension() is used. Otherwise
        *units* will be passed to the Dimension() constructor.
    """
    _mismatched_units_message = 'Comparing values with different units.'
    _unsupported_op_message = 'operation unsupported for dimensioned numbers.'

    def __init__(self, n, units=None):
        if isinstance(n, str):
            n, units = self.parse(n)
        self.n = n
        # TODO: So.... do I want to eliminate the redunant copy constructor here?
        self.units = Dimension() if units is None else Dimension(units)

    def __repr__(self):
        params = ', '.join([repr(x) for x in [self.n, self.units]])
        return ''.join([self.__class__.__name__, '(', params, ')'])

    def __str__(self):
        return ''.join([str(x) for x in [self.n, self.units]])

    @classmethod
    def parse(cls, s):
        """Parse a string into tuple of number and Dimension instance. The
        string must parse completely, or else this will raise.  Parse makes
        some assumptions about type: If the string starts with a '(', the
        expression is passed to the complex() constructor.  If the string
        starts with a '[', the expression is passed to the numpy.array()
        constructor.  Otherwise, it is passed to the float() constructor.

        :param s: A string to be parsed.
        :returns: A tuple of (Number, Dimension)
        """
        if s[0] in '([':
            "split point is balancing close bracket"
            p, n = 1, 1
            while n:
                if s[p] in '([':
                    n += 1
                elif s[p] in ')]':
                    n -= 1
                p += 1
        else:
            "split point is first char not part of number."
            p = 0
            while s[p] in '+-0123456789.eE':
                p += 1
        num, dim = s[0:p], s[p:]
        if not len(num):
            number = 1
        elif num[0] == '[':
            number = array(num)
        elif num[0] == '(':
            number = complex(num)
        else:
            number = float(num)
        dimension = Dimension.parse(dim)
        return number, dimension

    def _maybe_promote(self, n):
        """If n is dimensionless numeric, turn it into a dimensionless Dn."""
        return self.__class__(n, Dimension()) if isinstance(n, Number) else n

    def __add__(self, other):
        other = self._maybe_promote(other)
        return self.__class__(self.n + other.n, self.units + other.units)

    def __sub__(self, other):
        other = self._maybe_promote(other)
        return self.__class__(self.n - other.n, self.units - other.units)

    def __mul__(self, other):
        other = self._maybe_promote(other)
        return self.__class__(self.n * other.n, self.units * other.units)

    def __matmul__(self, other):
        other = self._maybe_promote(other)
        return self.__class__(self.n @ other.n, self.units @ other.units)

    def __truediv__(self, other):
        other = self._maybe_promote(other)
        return self.__class__(self.n / other.n, self.units / other.units)

    def __floordiv__(self, other):
        other = self._maybe_promote(other)
        return self.__class__(self.n // other.n, self.units // other.units)

    def __mod__(self, other):
        raise TypeError(' '.join(['mod', self._unsupported_op_message]))

    def __divmod__(self, other):
        raise TypeError(' '.join(['divmod', self._unsupported_op_message]))

    def __pow__(self, other, modulo=None):
        if modulo is not None:
            raise TypeError(' '.join(['pow() with modulo',
                self._unsupported_op_message]))
        o = int(other)
        return self.__class__(self.n ** o, self.units ** o)

    # Comparisions only make sense for Dn()'s with same units.
    def __eq__(self, other):
        other = self._maybe_promote(other)
        if self.units == other.units:
            return self.n == other.n
        raise TypeError(self._mismatched_units_message)

    def __ne__(self, other):
        other = self._maybe_promote(other)
        if self.units == other.units:
            return self.n != other.n
        raise TypeError(self._mismatched_units_message)

    def __lt__(self, other):
        other = self._maybe_promote(other)
        if self.units == other.units:
            return self.n < other.n
        raise TypeError(self._mismatched_units_message)

    def __gt__(self, other):
        other = self._maybe_promote(other)
        if self.units == other.units:
            return self.n > other.n
        raise TypeError(self._mismatched_units_message)

    def __le__(self, other):
        other = self._maybe_promote(other)
        if self.units == other.units:
            return self.n <= other.n
        raise TypeError(self._mismatched_units_message)

    def __ge__(self, other):
        other = self._maybe_promote(other)
        if self.units == other.units:
            return self.n >= other.n
        raise TypeError(self._mismatched_units_message)

    def __radd__(self, other):
        other = self._maybe_promote(other)
        return self.__class__(other.n + self.n, other.units + self.units)

    def __rsub__(self, other):
        other = self._maybe_promote(other)
        return self.__class__(other.n - self.n, other.units - self.units)

    def __rmul__(self, other):
        other = self._maybe_promote(other)
        return self.__class__(other.n * self.n, other.units * self.units)

    def __rmatmul__(self, other):
        other = self._maybe_promote(other)
        return self.__class__(other.n @ self.n, other.units @ self.units)

    def __rtruediv__(self, other):
        other = self._maybe_promote(other)
        return self.__class__(other.n / self.n, other.units / self.units)

    def __rfloordiv__(self, other):
        other = self._maybe_promote(other)
        return self.__class__(other.n // self.n, other.units // self.units)

    def __neg__(self):
        return self.__class__(-self.n, self.units)

    def __pos__(self):
        # FIXME: Is it really necessary to call a copy constructor here???
        # Verify and simplify if possible.
        return self.__class__(self.n, self.units)

    def __abs__(self):
        return self.__class__(abs(self.n), self.units)

    def __iadd__(self, other):
        other = self._maybe_promote(other)
        self.units += other.units
        self.n += other.n
        return self

    def __isub__(self, other):
        other = self._maybe_promote(other)
        self.units -= other.units
        self.n -= other.n
        return self

    def __imul__(self, other):
        other = self._maybe_promote(other)
        self.units *= other.units
        self.n *= other.n
        return self

    def __imatmul__(self, other):
        other = self._maybe_promote(other)
        self.units @= other.units
        self.n @= other.n
        return self

    def __itruediv__(self, other):
        other = self._maybe_promote(other)
        self.units /= other.units
        self.n /= other.n
        return self

    def __ifloordiv__(self, other):
        other = self._maybe_promote(other)
        self.units //= other.units
        self.n //= other.n
        return self

    def sqrt(self):
        """Implement sqrt() for dimensioned numbers."""
        return self.__class__(math_sqrt(self.n), self.units.root())

    def root(self, n):
        """Implement n-th root for dimensioned numbers."""
        n = int(n)
        return self.__class__(pow(self.n, float(1.0/n)), self.units.root(n))

# Exported functions
def prefer(dimension_str):
    """Set an override string to display a complex dimension as
    explicitly specified.
    :param dimension_str: The dimension as a string.
    """
    Dimension.prefer(dimension_str)

def sqrt(dimensioned_number):
    return dimensioned_number.sqrt()

def root(dimensioned_number, n):
    return dimensioned_number.root(n)

def unit_of(s):
    return Unit.of[s]

def set_factor_order(derived_unit_list):
    DerivedUnit.set_factor_order(derived_unit_list)

def unit_definitions():
    return Unit._all

__all__ = ['Dn', 'sqrt', 'root', 'unit_of', 'set_factor_order',
    'unit_definitions']