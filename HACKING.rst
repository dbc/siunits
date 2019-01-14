Hacking Bootstrap Guide
=======================

The full hacking guide is in doc/source/hacking.rst, and can be built
using Sphinx by going to doc/ and typing ``make html``.

If your system does not have an up-to-date Sphinx, create a 
virtual env in the usual way and add Sphinx. 
For example: ::

  cd ~
  mkdir pyenv 
  cd pyenv
  python3.6 -m venv siunits
  source siunits/bin/activate
  pip install sphinx
  
Then you should be able to build the documentation and get the
full hacking guide.

All development work to date has been done on Python 3.6.
PEP 465 for __matmul__ came in with Python 3.5. 
As siunits stands now, it will not work with Python earlier
than 3.5 because it uses the __matmul__ operator.

