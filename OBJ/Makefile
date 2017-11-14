obj=BUILDINS.o INTERPRETER.o  WAREHOUSE.o SOLVERM.o FILEM.o INPUTM.o GAUSSM.o OPERATIONM.o SORTM.o PKEAST.o GENERATEM.o  OUTPUTM.o DISCOVERY.o 

FEM_MFS.out : $(obj)
	gfortran -g -fcheck=all -Wall $(obj) -o FEM_MFS.out

%.o : %.F95
	gfortran -c $^ -o $@


.PHONY : clean
clean:
	-rm FEM_MFS.out $(obj)
		
