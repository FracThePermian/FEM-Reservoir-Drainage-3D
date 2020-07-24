# FEM-Reservoir-Drainage-3D

Compares analytical/numerical results for the drainage of a single well. Modify the following variables to your liking.

It is written in Matlab, is self-contained, and no external dependencys. 

### Assumptions
  * Incompressible Fluid
  * Homogeneous Reservoir
  * Neumann boundary (Reservoir side)
  * Dirichlet boundary (Well side)
### Modify Parameters
 * Reservoir Depth, Area, and Volume (ft ; ft^2 ; ft^3)
 * Reservoir Pressure (psi)
 * Time and/or Iterations (days)
 * Porosity (phi)
 * Permeability (k)
 * Compressibility (ct)
 * Volume Flow (Q)
 * Damage (hk)

## Method and Materials
* The Finite Element Method (FEM) is an essential numerical tool for solving boundary value problems of PDE's. 

**Snippet of module that generates sparsity matrix**
```matlab
for i = 1:N  		%Generate sparsity matrix
    if i+NX <= N 
        T_diag = T_frac(i,i+NX,Ay,mu,Bw,dy);
        T(i,i+NX) = T(i,i+NX)-T_diag;
        T(i+NX,i) = T(i+NX,i)-T_diag;
    end
    
    % This is for diagonals toward the inside of sparsity
    if (mod(i,NX) ~= 0) && (i+1 <= N)  
        T_diag = T_frac(i,i+1,Ax,mu,Bw,dx);
        T(i,i+1) = T(i,i+1)-T_diag;
        T(i+1,i) = T(i+1,i)-T_diag;
    end
    T(i,i) = abs(sum(T(i,:)));%Sum up every row on every iteration
    if BC(i) ~= 0 %if and only if the boundary condition exists
        T(i,i) = T(i,i)+BC(i)*2*T_frac(i,i,Ax,mu,Bw,dx);
    end
```


## Results

![Pressure Visualization](https://github.com/FracThePermian/FEM-Reservoir-Drainage-3D/blob/master/Graphical-Output/3D_Pressure_Visualization1.PNG "Pressure Gain/Reduction Profile")

![Pressure Contours](https://github.com/FracThePermian/FEM-Reservoir-Drainage-3D/blob/master/Graphical-Output/Reservoir_Contour2.PNG "Pressure Contours")

![Pressure vs. Radius](https://github.com/FracThePermian/FEM-Reservoir-Drainage-3D/blob/master/Graphical-Output/PressureVRadius3.PNG "psi v. ft.")



### References
  * Hughes, Thomas J. R.The Finite Element Method: Linear Static and Dynamic Finite Element Analysis. Englewood Cliffs, NJ: Prentice-Hall, 1987. Print.
  * Kolukula, Siva. MatlabMesh Postprocessing. Net.

### License
[MIT](https://opensource.org/licenses/MIT "MIT")
