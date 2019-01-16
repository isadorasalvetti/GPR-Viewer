
# Exercises

**1) Curvatures**

In display curvatures, display GAUSSIAN and MEAN CURVATURE. Implemented for closed and bordered meshes.

**2) Iteractive Smoothing**

Implemented LAPLACIAN and BILAPLACIAN with uniform and cotangent weights. Next to the weight selection is the iterations selectior. Cotangent weights will not work on bordered meshes and create artifacts after a few iterations.

**3) Global Smoothing**

Implemented on uniform and cotangent weights (uses the weight selector of the previous session). The method of selecting fixed vertices is determined by the selection box and the ammount of vertices to be fixed is determined by the slider next to it (blue = fixed, between 1 and 99%).

**4) Magnify frequency details**

There are three smoothing frequencies. The three numbers represent which multiplier will be used on each.

**5) Discrete harmonic map**

"Texture" will compute the parametric map (if already not computed) and apply a procedural texture to the mesh. "Display Parametrization" will compute the parametric map (if already not computed) and display the map using the vertices of the mesh (this replaces the current mesh).

**Other functions:**

Reset Colors = Remove color and texture from the mesh.

Reload = Reverts the mesh to its original state (as loaded).
