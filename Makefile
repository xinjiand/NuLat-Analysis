NuLatRoot: NuLatRoot.cpp
	g++ NuLatRoot.cpp -L/home/xinjian/program/root/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -pthread -lm -ldl -rdynamic -pthread -std=c++11 -m64 -I/home/xinjian/program/root/include -o NuLatRoot.exe


