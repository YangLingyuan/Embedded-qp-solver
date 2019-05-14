This qp_solver is to solve the following qp problem
	minimize     (1/2)*x'*P*x + q'*x + r
    	subject to   lb <= x <= ub

*******************************************************
How to use
*******************************************************

Use config / config_board for matrix size, test settings and to specify which tests to run.
To try on computer, best use
	make clean && make && ./main
To test against a pyhton based reference implementation (qpsolvers) python (u can specify version in config)
and the opt library are required. The relevant validation code is found at test/qp_ref.py and test/test.c.

In case u cannot use the exported version of the project for the board,
this is a guide to import the code manually for use on the board

1, Import the example project STK3700_emode
2, Delete main.c
3, Add source from repo (one can drag and drop and choose to use links instead of copying)
4, Add the folders "kmalloc", "matrix", "misc", "test" and "qp_solvers" to the include path of the project
   under project->property->C/C++ Build->Settings->GNU ARM C Compiler->Include
5, Add math (just add the letter m) in project->property->C/C++ Build->Settings->GNU ARM C Linker->Library->Libraries(-I) 
6, Build
