from os import path

from setuptools.command.develop import develop
from numpy.distutils.core import Extension, setup

compile_args_quad = {"libraries": ["quadminos"]}
compile_args_double = {"libraries": ["minos"]}

if path.isfile("libquadminos.a"):
    compile_args_quad["library_dirs"] = [path.abspath(".")]
else:
    raise Exception('Missing libquadminos.a')

if path.isfile("libminos.a"):
    compile_args_double["library_dirs"] = [path.abspath(".")]
else:
    raise Exception('Missing libminos.a')

ext_modules = [
    Extension(name="qminos.qwarmLP",
              sources=["qminos/src/lp/qwarmLP.f90"], **compile_args_quad),
    Extension(name="qminos.warmLP",
              sources=["qminos/src/lp/warmLP.f90"], **compile_args_double),
    Extension(name="qminos.qvaryME",
              sources=["qminos/src/fva/qvaryME.f90"], **compile_args_quad),
    Extension(name="qminos.qsolveME",
              sources=["qminos/src/nlp/qsolveME.f90",
                       "qminos/src/nlp/qmatrixA.f90"], **compile_args_quad),
    Extension(name="qminos.qquadprog",
              sources=["qminos/src/qp/qquadprog.f90",
                       "qminos/src/qp/qmatrixQ.f90"], **compile_args_quad),
    Extension(name="qminos.quadprog",
              sources=["qminos/src/qp/quadprog.f90",
                       "qminos/src/qp/matrixQ.f90"], **compile_args_double),
    Extension(name="qminos.qnlclp",
              sources=["qminos/src/nlclp/qnlclp.f90",
                       "qminos/src/nlclp/qmatrixA.f90"], **compile_args_quad),
    Extension(name="qminos.qnlcqp",
              sources=["qminos/src/nlcqp/qnlcqp.f90",
                       "qminos/src/nlcqp/qmatrixA.f90",
                       "qminos/src/nlcqp/qmatrixQ.f90" ], **compile_args_quad),
    ]

setup(
    name="qminos",
    ext_modules=ext_modules,
    cmdclass={"develop": develop}
    )
