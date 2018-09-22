# Fundamentals of Engineering Numerical Analysis (Python)
Examples from Fundamentals of Engineering Numerical Analysis, 2nd Edition, by Parviz Moin (ISBN: 9780521711234). Solutions from http://numerics.stanford.edu/ta/index.html, converted to python. Source code solutions for examples 6.14 and 6.15 is available from https://web.stanford.edu/group/fpc/ta/ch6/. Requirements:
* [Numpy](http://www.numpy.org/) - NumPy: array processing for numbers, strings, records, and objects.
* [Scipy](https://www.scipy.org/) - SciPy: Scientific Library for Python
* [Matplotlib](http://matplotlib.org/) - Python plotting package

Note:
The original mesh for example 6.16 is unavailable, so one was generated based on the given domain geometry using the open source finite-element meshing software [Gmsh](http://www.gmsh.info). This geometry is contained in the file ch6ex16.geo, with the mesh nodes and elements stored in files ch6ex16_37.npz (37 elements) and ch6ex16_491.npz (491 elements). These approximate the meshes used in the original example problem, and ch6ex16_37.pdf depicts the 37 element mesh as is done in figure 6.22 in the book.