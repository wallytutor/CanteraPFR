Development
===========

General Practices
-----------------

Since this project depends on external libraries, it is important to keep a
standard for updates. First, the target programming language for the package is
Python. That said, applications are expected to be conceived in that language.
Nonetheless, the main interface of the library used to develop these models is
provided in C++. This implies that development workflow should follow the path:

1. Conceive new models in C++:
    a. Header files shall be placed under *include/CanteraPFR*.
    b. Source files shall be placed under *src*.
    c. Test interface for C++ is provided in *test/main*.
    d. Edit *Makefile* to include test binary.
2. Provide Cython interface:
    a. Add C++ interfaces to *CanteraPFR.pdx*.
    b. Implement Python callable interface to new models in *CanteraPFR.pyx*.
    c. Add test to new model to a script under *test*.
3. Build and test project with new model:
    a. Run `make` from main directory.
    b. Execute the newly created test from its directory.
    c. Run `make` from *doc* directory for documentation.

Functionalities to implement
----------------------------

1. Models are not taking cross-section variations into account.
2. Generalize model by incorporating source (surface) terms.
3. Implement test for turbulent limit in momentum loss equation.
4. Provide analytical Jacobian for all models.
