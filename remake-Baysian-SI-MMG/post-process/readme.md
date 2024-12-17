## Compilimg `.f90` files
```bash
FC=gnu95 f2py -m mmg_fort -c --f90flags='-O3'  subroutine/mmg_esso_osaka_verctor_input.f90 
```
macの時は以下の方が成功する。
```bash
FC=gfortran f2py -m mmg_fort -c --f90flags='-O3'  subroutine/mmg_esso_osaka_verctor_input.f90 
```