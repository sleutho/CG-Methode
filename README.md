GC-Methode
==========
Fortran Implementierung der [Methode der konjungierten Gradienten](http://en.wikipedia.org/wiki/Conjugate_gradient_method) zur L�sung von linearen Gleichgungssystemen mit [symetrischer positiv definiter Matrix](http://en.wikipedia.org/wiki/Positive-definite_matrix).

##Usage:
```
dim (INTEGER n,Dimension)
kmax (INTEGER kmax,max Iterationen)
eps (REAL eps, gew�nschte Genauigkeit)
Matrix	(REAL a(i,j))
b (REAL b(i))
x0 (REAL x0(i))
```


