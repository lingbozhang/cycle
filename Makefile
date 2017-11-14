obj_name=INTERPRETER.o  WAREHOUSE.o FILEM.o INPUTM.o GAUSSM.o OPERATIONM.o SORTM.o PKEAST.o GENERATEM.o SOLVERM.o OUTPUTM.o BUILDINS.o DISCOVERY.o 
obj=$(addprefix ./OBJ/,$(obj_name))

FEM_MFS.out : $(obj)
	gfortran -g -fcheck=all -Wall $(obj) -o FEM_MFS.out

%.o : %.F95
	gfortran -c $^ -o $@


.PHONY : clean
clean:
	-rm FEM_MFS.out $(obj)
		
