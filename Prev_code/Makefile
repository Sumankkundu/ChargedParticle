LIB = -L${ROOTSYS}/lib -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl -rdynamic -L/home/tanmay/Software/clhep -lCLHEP

anal_merge_root:
#	g++ -g -I${ROOTSYS}/include -I${CLHEP_PARAM:_PATH}/include -o anal_merge_root anal_merge_root.C $(LIB)
	g++ -g -I${ROOTSYS}/include -I${CLHEP_PARAM:_PATH}/include -o anal_merge_root_v1 anal_merge_root_v1.C $(LIB)
clean:
#	anal_merge_root	
	anal_merge_root_v1
