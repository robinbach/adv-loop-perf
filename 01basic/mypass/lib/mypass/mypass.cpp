/*
 * File: mypass.c
 *
 * Description:
 *  This is a mypass source file for a library.  It helps to demonstrate
 *  how to setup a project that uses the LLVM build system, header files,
 *  and libraries.
 */

#include <stdio.h>
#include <stdlib.h>

/* LLVM Header File
#include "llvm/Support/DataTypes.h"
*/

/* Header file global to this project */
#include "mypass.h"

#include "llvm/Pass.h"
#include "llvm/IR/Function.h"
#include "llvm/Analysis/ProfileInfo.h" 
#include "llvm/Analysis/LoopInfo.h"

#include "llvm/Transforms/Scalar.h"
#include "llvm/ADT/Statistic.h"
#include "llvm/Analysis/AliasAnalysis.h"
#include "llvm/Analysis/AliasSetTracker.h"
#include "llvm/Analysis/ConstantFolding.h"
#include "llvm/Analysis/Dominators.h"
#include "llvm/Analysis/LoopInfo.h"
#include "llvm/Analysis/LoopPass.h"
#include "llvm/Analysis/ValueTracking.h"
#include "llvm/IR/Constants.h"
#include "llvm/IR/DataLayout.h"
#include "llvm/IR/DerivedTypes.h"
#include "llvm/IR/Instructions.h"
#include "llvm/IR/IntrinsicInst.h"
#include "llvm/IR/LLVMContext.h"
#include "llvm/IR/Metadata.h"
#include "llvm/Support/CFG.h"
#include "llvm/Support/CommandLine.h"
#include "llvm/Support/Debug.h"
#include "llvm/Support/raw_ostream.h"
#include "llvm/Target/TargetLibraryInfo.h"
#include "llvm/Transforms/Utils/Local.h"
#include "llvm/Transforms/Utils/SSAUpdater.h"
#include <algorithm>

// SLICM includes
#include "llvm/Transforms/Utils/BasicBlockUtils.h"
#include "LAMP/LAMPLoadProfile.h"
#include "llvm/Analysis/ProfileInfo.h"

// 583 hw2 additional includes
#include <vector>
#include <tr1/unordered_map>
#include <tr1/unordered_set>
#include <map>


using namespace llvm;
using namespace std;

namespace {
	struct mypass: public LoopPass { 
		static char ID; // Pass identification, replacement for typeid
		ProfileInfo * PI;
		LoopInfo * LI;       // Current LoopInfo
		mypass() : LoopPass(ID) { }


		tr1::unordered_set<Instruction *> visitedloops;

		virtual bool runOnLoop(Loop *L, LPPassManager &LPM) {
			LI = &getAnalysis<LoopInfo>();
			PI = &getAnalysis<ProfileInfo>();

			// 1. Branch: br, switch, indirectbr:
			int branchCount = 0;

				// bias count
			int biasCount = 0;
			int unbiasCount = 0;

			// 2. Integer ALU counts:
			int intAddCount = 0;
			int intSubCount = 0;
			int intMulCount = 0;

			int intUDivCount = 0;
			int intSDivCount = 0;
			int intURemCount = 0;

			int intShlCount = 0;
			int intLShrCount = 0;
			int intAShrCount = 0;

			int intAndCount = 0;
			int intOrCount = 0;
			int intXorCount = 0;

			int intICmpCount = 0;
			int intSRemCount = 0;

			// 3. Floating-point ALU: fadd, fsub, fmul, fdiv, frem, fcmp:
			int flALUCount = 0;

			// 4. Memory operation: alloca, load, store, getelementptr, fence, cmpxchg, atomicrmw:
			int memOpCount = 0;

			// 5. other operations
			int otherOps = 0;

			// for (Function::iterator b = F.begin(); b != F.end(); b++) {
			// for (Loop::block_iterator b = F.begin(); b != F.end(); b++) {

			bool inserted = false;
			int loopiterations = 0;
			Instruction * head;
			for (Loop::block_iterator BB = L->block_begin(), BBE = L->block_end();
       BB != BBE; ++BB){
				// errs() << "Func\n";
				int blockbranchCount = 0;
				Instruction * bbhead = &( *(*BB)->begin() );
				if(visitedloops.count(bbhead)){
					errs() << "Inner loop of " << head;
					errs() << " with head id: " << bbhead << " detected\n";
					// return true;
				}
				else if(!inserted){
					inserted = true;
					head = bbhead;
					visitedloops.insert(bbhead);
					loopiterations = PI->getExecutionCount(*(BB+1)); 
				}
					

				// for (BasicBlock::iterator i = b->begin(); i != b->end(); i++) { 
				for (BasicBlock::iterator i = (*BB)->begin(), E = (*BB)->end();
         i != E; ++i){
					// errs() << "Block\n";
					
					int exeCount = PI->getExecutionCount(*BB) > 0 ? PI->getExecutionCount(*BB) : 0;
					switch (i->getOpcode()) {
						// 1. Branch:
						case Instruction::Br: case Instruction::Switch:
						case Instruction::IndirectBr:
							branchCount += exeCount;
							blockbranchCount = exeCount;
							break;

						// 2. Integer ALU: add, sub, mul, udiv, sdiv, 
						// urem, shl, lshr, ashr, and, or, xor,
						// icmp, srem
						case Instruction::Add:
							intAddCount += exeCount;
							break;
						case Instruction::Sub:
							intSubCount += exeCount;
							break;
						case Instruction::Mul:
							intMulCount += exeCount;
							break;
						case Instruction::UDiv:
							intUDivCount += exeCount;
							break;							
						case Instruction::SDiv:
							intSDivCount += exeCount;
							break;							
						case Instruction::URem:
							intURemCount += exeCount;
							break;							
						case Instruction::Shl:
							intShlCount += exeCount;
							break;							
						case Instruction::LShr:
							intLShrCount += exeCount;
							break;							
						case Instruction::AShr:
							intAShrCount += exeCount;
							break;													
						case Instruction::And:
							intAndCount += exeCount;
							break;							
						case Instruction::Or:
							intOrCount += exeCount;
							break;							
						case Instruction::Xor:
							intXorCount += exeCount;
							break;							
						case Instruction::ICmp:
							intICmpCount += exeCount;
							break;							
						case Instruction::SRem:
							intSRemCount += exeCount;
							break;

						// 3. Floating-point ALU: fadd, fsub, fmul, fdiv, frem, fcmp:
						case Instruction::FAdd: case Instruction::FSub:
						case Instruction::FMul: case Instruction::FDiv:
						case Instruction::FRem: case Instruction::FCmp:
							flALUCount += exeCount;
							break;	

						// 4. Memory operation: alloca, load, store, getelementptr, fence, cmpxchg, atomicrmw:
						case Instruction::Alloca: case Instruction::Load:
						case Instruction::Store: case Instruction::GetElementPtr:
						case Instruction::AtomicCmpXchg: case Instruction::AtomicRMW:
							memOpCount += exeCount;
							break;	

						default:
							otherOps += exeCount;
						break;


					} // -- end switch --


				} // -- end for --

				// check if bias
				if(blockbranchCount > 0){
					unbiasCount += blockbranchCount;
					for (Loop::block_iterator BB2 = L->block_begin(), BBE2 = L->block_end();
       BB2 != BBE2; ++BB2){
						ProfileInfo::Edge branchEdge = PI->getEdge(*BB, *BB2);
						int weight = PI->getEdgeWeight(branchEdge);
						double ratio = weight * 1.0/blockbranchCount;
						if(ratio > 0.8){
							biasCount += blockbranchCount;
							unbiasCount -= blockbranchCount;
							break;
						}
					}
				}
				
			}
			// 1. branch operation sum:

			// 2. int ALU operations sum:
			int intALUSum = intAddCount + intSubCount + intMulCount + intUDivCount 
				+ intSDivCount + intURemCount + intShlCount + intLShrCount 
				+ intAShrCount + intAndCount + intOrCount + intXorCount 
				+ intICmpCount + intSRemCount;

			// 5 total sum:
			int opTotalSum = branchCount + intALUSum + flALUCount + memOpCount + otherOps;
			// Output section
			//FuncName,DynOpCount,%IALU,%FALU,%MEM,%Biased-Branch,%Unbiased-Branch,%Oth
			// errs().write_escaped(F.getName()) << ",";
			errs() << "Loop" << head << ",";
			errs() << opTotalSum << ",";
			errs() << loopiterations << ",";
			// errs() << 1.0f * intALUSum/opTotalSum << ",";
			// errs() << 1.0f * flALUCount/opTotalSum << ",";
			// errs() << 1.0f * memOpCount/opTotalSum << ",";
			// errs() << 1.0f * biasCount/opTotalSum << ",";
			// errs() << 1.0f * unbiasCount/opTotalSum << ",";
			// errs() << 1.0f * otherOps/opTotalSum << "\n";
			errs() << format("%f,%f,%f,%f,%f,",
				1.0f * intALUSum/opTotalSum, 1.0f * flALUCount/opTotalSum, 
				1.0f * memOpCount/opTotalSum, 1.0f * biasCount/opTotalSum, 
				1.0f * unbiasCount/opTotalSum);
			errs() << format("%f\n", 1.0f * otherOps/opTotalSum);
			// part 1: basic profiler end.



			return true; 
		}

		void getAnalysisUsage(AnalysisUsage &AU) const {
			AU.addRequired<ProfileInfo>();
			AU.addRequired<LoopInfo>();
		} 
	};
}

char mypass::ID = 0;
static RegisterPass<mypass> X("opcstats", "loop operation count pass", false, false);

