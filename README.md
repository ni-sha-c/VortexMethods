Vortex Methods
=============

<h2> Warm up : Norbury Problem. </h2>

<h3>To Execute</h3>

* The driver module - driver.py
* The parameter alpha and initial seed for X can be set here.



<h3> Function Synopsis </h3>

* driver module actually solves equation 2.19* of the paper.
* The matrix populated by <a href="http://www.codecogs.com/eqnedit.php?latex=\frac{\partial&space;P_\alpha(X^1,t)}{\partial&space;X_j}" target="_blank"><img src="http://latex.codecogs.com/gif.latex?\frac{\partial&space;P_\alpha(X^1,t)}{\partial&space;X_j}" title="\frac{\partial P_\alpha(X^1,t)}{\partial X_j}" /></a> s  is formed here.
I haven't evaluated the matrix another time because the paper says keeping it as it is, across iterations, gave good results.

* The <a href="http://www.codecogs.com/eqnedit.php?latex=P_\alpha" target="_blank"><img src="http://latex.codecogs.com/gif.latex?P_\alpha" title="P_\alpha" /></a> evaluation is 
in the module P. In this module, there is:


    in_integral 
    

that finds the inner integral at a given <a href="http://www.codecogs.com/eqnedit.php?latex=\hat{t}" target="_blank"><img src="http://latex.codecogs.com/gif.latex?\hat{t}" title="\hat{t}" /></a>
. The inner integral uses Gauss quadrature and the integrand is evaluated in  G.py

For lack of a better idea, I used Simpson's rule manually in 
 
    get_integral
    



<h3> Where refactoring is needed </h3>

At first glance, I could make the outer integral more efficient. 
It looks to me that other inefficiencies can be attributed to me being a Python noob.











