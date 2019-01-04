# -*- coding: utf-8 -*-

from ctypes import addressof


cdef class Closure:
    """ Just-in-time compiler for Python functions.

    Allows a Python function to be called from C/C++ library.
    Based on https://stackoverflow.com/questions/51044122.

    TODO
    ----
        This class is incoherent, once the interface of `_inner_fun` is fixed
        and does not follow what `ftype` requires. This must be wrapped.

    Parameters
    ----------
    python_fun : function
        Python function to call from C/C++.
    ftype : ctypes.CFUNCTYPE
        Interface to provide the function.
    """

    cdef object python_fun
    cdef object jit_wrap

    def __cinit__(self, python_fun, ftype):
        self.python_fun = python_fun
        self.jit_wrap = ftype(lambda *args: self.python_fun(*args))

    cdef func_t get_fun_ptr(self):
        return (<func_t *><size_t>addressof(self.jit_wrap))[0]
