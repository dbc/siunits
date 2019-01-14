siunits primary interface
=========================

The Dn class is sufficient for most work with siunits.

.. automodule:: siunits
   :members: Dn

siunits also defines some convenience functions to help set
preferred display formats and to simplify access to other
functionality.

.. currentmodule:: siunits
.. autofunction:: prefer
.. autofunction:: sqrt
.. autofunction:: root
.. autofunction:: unit_of
.. autofunction:: set_factor_order
.. autofunction:: unit_definitions

siunits Getting Started Guide
=============================

The main interface to the siunits module is the ``Dn()``, 
"dimensioned number" class. 
A dimensioned number supports normal arithmetic, but in 
addition to a numeric value, the dimenions are carried along
and updated when performing arithmetic.
In cases where the arithmetic does not make sense for numbers
with the given units, an exception is raised.

For example: ::

  from siunits import Dn

  m = Dn(2, 'kg')  # m is a mass
  a = Dn(3, 'm/s^2')  # a is an acceleration
  f = m*a  # f will get units of force: kgm/sec^2, a.k.a. Newtons (N)
  print (f)

The printed result will is: ::

  6N

But if you try arithetic with units that don't make sense: ::

  cd = Dn(4, 'cd')  # some candelas
  print(m+cd)  # What happens if you add meters and candelas?

This raises an exception: ::
 
  TypeError: Adding mismatched units.

The ``Dn()`` constructor is forgiving.
In the canonical form, it takes a number-like thing, and a dimension.
siunits takes advantage of Python "duck typing", so the number-like
thing is just an instance of anything that supports arithmetic.
It can be an integer, a float, a complex, a Numpy array, or even
a user-defined type that implements the special arithmetic 
"dunder" methods (__add__, __mul__, etc).

The dimension is a string that is a units expression, or an
instance of ``siunits.Dimension()``. 
The known units are: 

- A ampere
- C coulomb
- F farad
- H henry
- Hz hertz
- J joule
- K kelvin
- N newton
- Ohm ohm
- Pa pascal
- S siemens
- T tesla
- V volt
- W watt
- Wb weber
- cd candela
- kat katal
- kg kilogram
- lx lux
- m meter
- mol mole
- rad radian
- s second

The ``Dn()`` constructor can also be called in other forms:

- Dn('200K')  -- pass it a string, and it will be parsed as a number
  and some units. This is often convenient, but is not as general as
  passing the number seperately.
- Dn('kg')  -- pass in just a units string. The number 1 will be
  dimensioned and returned.  It can sometimes be more readable to
  write out something as: 4*Dn('kg'), which is completely equivalent
  to writing Dn(4, 'kg').

The form ``Dn('number-units')`` string form is not quite as general as the 
(value, dimension) form of the constructor because the parser makes
certain assumptions about the number contained in the string it
is parsing.  If the string begins with a digit, it calls ``float()``
on the numeric part to construct a floating point dimensioned number.
If the string begins with a paren, then it calls ``complex()`` on
the number part of the string.  If the string begins with a 
bracket, it calls ``numpy.array()`` on the string.  So: ::

  - Dn('200K') is equivalent to: Dn(float(200), 'K')
  - Dn('(3+2j)A') is equivalent to: Dn(complex(3+2j), 'A')
  - Dn('[[1,2],[3,4]]W') is equivalent to: 
    Dn(numpy.array([[1,2],[3,4]], 'W')

So you can not create dimensioned integers or dimensioned versions
of user-defined types with the single-string form of the constructor.

If you need to create a dimension-less dimensioned number, simply 
call the ``Dn()`` constructor with only a number: ::

  tau = Dn(6.28)
  >>> repr(tau)
  'Dn(6.28, Dimension())'


The ``repr()`` for ``Dn()'s`` spells out the dimension in canonical
form in SI base units. ::

  >>> repr(f)
  'Dn(6, Dimension(meter=1, kilogram=1, second=-2))'

Controlling Display Formatting
==============================

Dimensioned numbers make a valiant attempt to unwind the units that
result from a chain of arithmetic into a presentable string.
Calling str() on an instance of Dn() will return a string that
has the number concantenated with a string spelling out the units.
A lot of times, it will make good sense.
But sometimes the output will look pretty strange.
It is correct, but just looks odd.
For example, Hertz is the unit for 1/second. 
It is perfectly correct to say: meter-Hertz, but that
looks pretty odd since we are accustomed to seeing meter/second.

Even though the default constructor for unit string does a decent
job, sometimes you may want to control things yourself.
There are two ways to do that:

- You can explictly set a preferred display format using prefer().
- You can control the order that derived units are factored out
  of a complex unit using set_factor_order().

The prefer() function is very straight-forward.
Let's say you don't like the default display the canonical unit group FIXME,
which is FIXME, but would rather see FIXME. 
Simply say ``prefer(FIXME)``. 
This creates a dictionary entry mapping the unit group to FIXME.
If any unit group hits in the preferences dictionary, you get
that exact string back.

For any complex unit group that does not hit in the preferences
dictionary, the next step is to look for an exact match among
the derived units.  This is not under user control.
It handles the simple cases of derived units like Newton (N) and Weber (Wb), 
and their simple powers.

If neither of those lookups return a result, an attempt is made
to factor the unit group.
This is where things can get awkward, as there may be multiple
valid factorings, or there may be a possible factoring that is
beyond the current factoring capability of siunits.

For cases where the unit group doesn't factor at all, the best
way to get a nice display string is to use ``prefer()``.

Sometimes, the factoring for display can be improved by factoring
out derived units in some other order than the default order.
This is where ``set_factor_order()`` comes into play.
Call ``set_factor_order()`` with a list of unit abbreviations 
(or unit names) in the order in which they should be tried.
This will override the default factor order with yours.
Example: ::

  TBW

The default factor order is: TBW

