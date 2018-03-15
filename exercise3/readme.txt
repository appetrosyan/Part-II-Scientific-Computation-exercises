-*- mode: org-mode

* Exercise 3B: Helmholtz Coils and magnets

** Structure of the Solution

   The entire solution is contained in the file [[file:helmholtz.py][Helmholtz.py]], which
   generates plots inside the folder =figures/=. Every function is
   type-annotated, and provided with a doc-string, wherever an
   explanation might be needed. 

   The solution is rather general, as it can cope with many different
   wire geometries, so long as a generator of the =StraightWire=
   segments is provided. 

   *Note that for reasons of clarity, the notation followed here
   assumes that the coils are oriented with their normals in the z
   direction, not x as in the practical handout*

   Briefly, it was done to make it easier to make sense of the
   orientation of the coils. A clockwise current in the x-y plane
   produces a positive field in the z direction. 
   


** Single-Coil on axis:

   The plot =single_coil_on_axis-values.pdf= shows agreement with
   analytical value. One can clearly see that the discrepancy is
   small, and that the numerical method reproduces the theoretical
   results to reasonable accuracy. 

   To quantify exactly how small, consider the plot
   =single_coil_on_axis-residuals.pdf=. Here one can notice that the
   relative discrepancy is a value $O(10^{-14})$. It has a large
   constant factor, and fluctuates slightly point-to-point. 

   By re-running the simulation for different values of resolution we
   can find that the optimal resolution is $O(2^6)$. 
   
   Note, that Counter-intuitively increasing the number of subdivisions
   does not always, reduce the relative error. 
  
   This is due to the fact that for extreme resolutions, both the
   error in length of the chords (the sides of a regular polygon inscribed into the circle) and the directions of the normals are dominated by
   roundoff errors. Hence by increasing the number of iterations we
   can decrease the accuracy of our computation. 

** Single coil in the y-z plane:
   
   The plot of the magnetic field, is given in figure:
   =1_coils_yz_section=, and shows some of the expected features of
   the magnetic field of a single coiled wire. Namely, the field is
   symmetric about the $y=0$ plane, with our choice of orientation of
   the current (i.e. flowing clockwise in the x-y plane), the field  has the
   right sign. 

   Moreover, we can see a visual representation of the coil
   itself. There are two vortexes in the B field, which correspond to
   the places where the coil intersects the y-z plane. 

   Also, due to the symmetry we have imposed on the construction of
   straight-wire segments, we see that the field has only the z
   component on axis. 

** Helmholtz coils:

   The plot =helmholtz_coils_on_axis=, shows how the field varies in
   presence of two coils. The theoretical values are marked on with
   black crosses, and agree to within $O(10^-{13})$, which agrees with
   our prediction from previous investigation.

   Then we move on to plotting the field within a cylinder of $10
   \times 10 cm$, =2_coils_yz_section=. The field is demonstrably
   homogeneous, and the maximum value of the deviation was evaluated
   to $0.00419\ T$, i.e. $4%$ of the value at the peak. 

   If we were to expand the Magnetic field at the origin, using
   $\delta z = 0.1\ m$   quadratic terms of the Taylor expansion would
   be of the order of a  few per-cent, and thus we would expect for
   the field to be homogeneous to that extent. 

   It is also worthwhile to vary the separation of the coils, as then
   the reason for this homogeneity is apparent. 

   Upon moving the coils closer, the field peaks more sharply, while
   when the coils are further than one radius apart, there's a local
   minimum precisely in the middle. 

   
** Multiple coils:

   As soon as we increase the number of coils, we can see that the
   peak of the value in the median peak for the field at the centre,
   is increased. 

   This is the expected behaviour, as asymptotically, this collection
   of coils is a solenoid, where the field strength in the middle is
   proportional to the density of wiring on the periphery. 
   
** A note on optimisations: 

   The code is multi-process, to provide as much speed as
   possible. Interestingly, parallelising the =field= function in at
   when iterating over the wire segments actually slowed the code down, while in other cases it provided a near-linear scaling with the number of threads. I suspect that this is due to the task consuming too much memory, and having to resort to the swap file, based on system monitoring on my
   machine. 

   Compromises have been made between generality and performance. Code
   that runs better is less abstract, since it by definition has to do fewer checks. This sadly means that some of the more elegant solutions such as using a single function to evaluate the field of any wire, or allowing for Circular wires to have normals other than along the $z$ axis, had to be scrapped in favour of those that run quicker. 

   An additional modest speedup could be gained, by converting the
   most time consuming function into a =numpy.ufunc=. However, in that
   case the =CircularWire= class would need to be rewritten in C,
   which sacrifices much generality and readability in favour of too
   little a gain. 
   

   

   
