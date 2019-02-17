This is a guide to import the code on board

1, Import the example project STK3700_emode
2, Replace the main file with main_board.c 
3, Add the folders "kmalloc", "matrix", "misc" and "qp_solvers" in the include path of the project under project->property->C/C++ Build->Settings->GNU ARM C Compiler->Include
4, Add the source file in the above folders under same path
5, Add math library(m) in project->property->C/C++ Build->Settings->GNU ARM C Linker->Library->Libraries(-I) 
6, Build and download the project
