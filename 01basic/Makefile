# PATH=$PATH:/opt/llvm/Release+Asserts/bin 

# ssh bestname@eecs583a.eecs.umich.edu

# /opt/llvm-source/: contains the source code of llvm and clang
# /opt/llvm/: contains the regular build of llvm and clang
# /opt/llvm-debug/: contains the debug build of llvm and clang (for gdb)

# $ clang -emit-llvm test.c -c -o test.bc
# $ llc test.bc -o test.s
# $ g++ -o test test.s
# $ ./test

# cd /opt/llvm-source/projects/sample
# cp -vr * /host/bestname/mypass

# shopt -s globstar
# rename -n 's/sample/mypass/' **
# rename 's/sample/mypass/' **
# grep -R sample *

# rm -rf lampp

############# Important:
# To profile an application (e.g., test.c):

### fun folders: 
# /opt/llvm/Release+Asserts/bin
# /opt/llvm-source/include/llvm

# #simple
# fname=simple
# type=c
# input=in.txt

fname=test
type=cpp
input=in.txt

# #wc
# fname=wc
# type=c
# input=cccp.c

#compress
# fname=compress
# type=c
# input=compress.in

output1=$(fname).opcstats
output2=$(fname).ldstats

## usage: make profile -> make all
#### run mypass
all:
	opt -f -load mypass/Release+Asserts/lib/mypass.so -profile-loader \
		-profile-info-file=llvmprof.out -opcstats $(fname).bc > /dev/null 
		# 2> $(output1)

# opt -load mypass/Release+Asserts/lib/mypass.so -profile-loader \
# 	-profile-info-file=llvmprof.out -opcstats $(fname).bc
# mem:
# 	opt -f -load mypass/Release+Asserts/lib/mypass.so -profile-loader \
# 		-profile-info-file=llvmprof.out -ldstats $(fname).bc > /dev/null 2> out2.txt
# 	cat out2.txt

profile:
	rm -f llvmprof.out
	clang -emit-llvm $(fname).$(type) -c -o $(fname).bc
	opt -insert-edge-profiling $(fname).bc -o $(fname).profile.bc
	llc $(fname).profile.bc -o $(fname).profile.s
	g++ -o $(fname).profile $(fname).profile.s /opt/llvm/Release+Asserts/lib/libprofile_rt.so
	./$(fname).profile $(input) > out1.txt

# fname=$1


# lamp:
# 	# tested , oK
# 	rm -f llvmprof.out # Otherwise your profile runs are added together
# 	clang -emit-llvm -o $(fname).bc -c $(fname).$(type) || \
# 		{ echo "Failed to emit llvm bc"; exit 1; }
# 	opt -loop-simplify < $(fname).bc > $(fname).ls.bc || \
# 		{ echo "Failed to opt loop simplify"; exit 1; } 
# 	opt -insert-edge-profiling $(fname).ls.bc -o $(fname).profile.ls.bc
# 	llc $(fname).profile.ls.bc -o $(fname).profile.ls.s
# 	g++ -o $(fname).profile $(fname).profile.ls.s /opt/llvm/Release+Asserts/lib/libprofile_rt.so
# 	./$(fname).profile $(input)

# lamp2:
# 	opt -load mypass/Release+Asserts/lib/mypass.so -lamp-insts -insert-lamp-profiling \
# 	-insert-lamp-loop-profiling -insert-lamp-init < $(fname).ls.bc > $(fname).lamp.bc || \
# 		{ echo "Failed to opt load"; exit 1; }
# 	llc < $(fname).lamp.bc > $(fname).lamp.s || \
# 		{ echo "Failed to llc"; exit 1; }
# 	g++ -o $(fname).lamp.exe $(fname).lamp.s mypass/tools/lamp-profiler/lamp_hooks.o || \
# 		{ echo "Failed to g++"; exit 1; } 
# 	./$(fname).lamp.exe $(input)
# 	echo "Done generating lamp profile"

# lamp3:
# 	opt -load mypass/Release+Asserts/lib/mypass.so -lamp-inst-cnt -lamp-map-loop \
# 		-lamp-load-profile -profile-loader -profile-info-file=llvmprof.out -ldstats \
# 		< $(fname).ls.bc > /dev/null 2> $(output2) || \
# 		{ echo "Failed to opt-load"; exit 1; }	

clean:
	rm -f test test.bc test.p* test.s test.l* llvmprof.out test.opcstats

# submit:
# 	rm -f bestname_hwk1/src/*
# 	rm -f bestname_hwk1/results/*
# 	cp *opcstats bestname_hwk1/results/
# 	# cp *ldstats bestname_hwk1/results/
# 	cp mypass/lib/mypass/mypass.cpp bestname_hwk1/src/
# 	cp mypass/lib/mypass/memprofile.cpp bestname_hwk1/src/
# 	mv bestname_hwk1/src/mypass.cpp bestname_hwk1/src/operation_statistics.cpp 
# 	mv bestname_hwk1/src/memprofile.cpp bestname_hwk1/src/memory_profile.cpp 

# cands:
# 	tar -cvzf bestname_hw1.tgz bestname_hw1
# 	scp bestname_hw1.tgz bestname@eecs583a.eecs.umich.edu:/hw1_submissions/





