# OpenMPM3D
Our laboratory has developed [MPM3D](http://comdyn.hy.tsinghua.edu.cn/english/mpm3d) since 2004. And we provided a simplified version named [MPM3D-F90](https://github.com/xzhang66/MPM3D-F90) in 2016. It is inspiring that lots of scholars star and fork our code. We still recommend MPM beginners to read the Fortran version codes and the following two books to gain an overall comprehension about Material Point Method(MPM).

1. X Zhang, Z Chen, Y Liu. [**The Material Point Method** - A Continuum-Based Particle Method for Extreme Loading Cases](http://store.elsevier.com/The-Material-Point-Method/Xiong-Zhang/isbn-9780124077164/). Academic Press, 2016; 

2. X Zhang, YP Lian, Y Liu, X Zhou. [Material Point Method (in Chinese)](http://comdyn.hy.tsinghua.edu.cn/downloads/mpm-book), Beijing: Tsinghua University Press, 2013.

With the development of C++ language, C++ program can run as quick as(or even quicker than) the corresponding Fortran program. Besides, C++ program based on OOP is in favor of incremental development and version control. We highly recommend fellow scholars to develop your own MPM solvers and carry on corresponding research based on this repository [OpenMPM3D]() of C++ version.

This repository is a 3-dimensional Material Point Method(MPM) code of C-plus-plus language. It extends the [MPM3D-F90](https://github.com/xzhang66/MPM3D-F90) with the following characteristics:
- Input file in XML format
- Various failure model
- String function parsing support. (+, -, *, /, sin(), asin(), cos(), acos(), tan(), atan(), exp(), log(), log10(), abs(), sqrt(), cosh(), sinh(), tanh()... You can expand the function parsing library by yourself.)
- Stagger Grid Material Point Method(SGMPM) solver
- Finite Volume Method(FVM) solver
- anything will be added in the future...

Still, this repository is a simplified version and framework of [MPM3D](http://comdyn.hy.tsinghua.edu.cn/english/mpm3d). 
Please cite our books and related publications in your publication if **MPM3D-F90**, **OpenMPM3D** or **MPM3D** is used in your work. Here is the List of [our publications in Material Point Method](http://comdyn.hy.tsinghua.edu.cn/103-achievements/mpm3d-en/553-mpm-publications).

## Author
Dr. Ruichen Ni

>School of Aerospace Engineering, Tsinghua University, Beijing 100084, China.

>E-mail: thurcni@163.com

Prof. Xiong Zhang(*Corresponding author)

>School of Aerospace Engineering, Tsinghua University, Beijing 100084, China.

>E-mail: xzhang@tsinghua.edu.cn

(P.S. If you have any question about our code, don't hesitate to contact us. Also, we are glad to recieve any bug-report to help us refining the code.)

## Contributor
- Develop your own MPM-corresponding solvers. (You can follow the procedure of SGMPM solver and FVM solver. Leave the base class unchanged as possible as you can.)
- Merge your codes and several example input datafiles into this repository to make your work easy for other scholars to follow.

### Contributor List
If contribute to this repository, please add your information here as the following template.
#### Solver name

> Author name

> Your contact information, such as E-mail address

> Corresponding article citation

#### MPM(Material Point Method)

> Ruichen Ni(Code-writter), Xiong Zhang*

> E-mail: xzhang@tsinghua.edu.cn

> Xiong Zhang, Zhen Chen, and Yan Liu. The Material Point Method: A Continuum-Based Particle Method for Extreme Loading Cases. Academic Press, USA, 2016.


## Folder Structure
### Data
&emsp;Example input datafile (.xml).

### build
&emsp;Folder left for CMake compiling.

&emsp;(P.S. You can compile the whole project by following the instructions in "HowToBuild.md".)

### src
&emsp;Source codes.