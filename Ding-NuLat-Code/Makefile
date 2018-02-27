
objects= anaMuon.exe NuLatAB-Plot.exe NuLatMuon-Plot.exe NuLat-Primary.exe NuLatAB.exe NuLatAB-PSD.exe anaposAB.exe NuLatABNoRoot.exe NuLatWindow.exe shell-generation.exe Map.exe PSD-Shell.exe event-list.exe NuLatNewMuon.exe

analysis: $(objects)

anaMuon.exe: anaMuon.cpp anafunction.h
	g++ anaMuon.cpp -L/home/daq/Program/root/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread  -pthread -lm -ldl -rdynamic -pthread -std=c++11 -m64 -I/home/daq/Program/root/include -o anaMuon.exe

NuLatAB-Plot.exe: NuLatAB-Plot.cpp anafunction.h
	g++ NuLatAB-Plot.cpp -L/home/daq/Program/root/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread  -pthread -lm -ldl -rdynamic -pthread -std=c++11 -m64 -I/home/daq/Program/root/include -o NuLatAB-Plot.exe

NuLatMuon-Plot.exe: NuLatMuon-Plot.cpp anafunction.h
	g++ NuLatMuon-Plot.cpp -L/home/daq/Program/root/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread  -pthread -lm -ldl -rdynamic -pthread -std=c++11 -m64 -I/home/daq/Program/root/include -o NuLatMuon-Plot.exe

NuLat-Primary.exe: NuLat-Primary.cpp anafunction.h
	g++ NuLat-Primary.cpp -L/home/daq/Program/root/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread  -pthread -lm -ldl -rdynamic -pthread -std=c++11 -m64 -I/home/daq/Program/root/include -o NuLat-Primary.exe

NuLatAB.exe: NuLatAB.cpp anafunction.h
	g++ NuLatAB.cpp -L/home/daq/Program/root/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl -rdynamic -pthread -std=c++11 -m64 -I/home/daq/Program/root/include -o NuLatAB.exe

NuLatAB-PSD.exe: NuLatAB-PSD.cpp anafunction.h
	g++ NuLatAB-PSD.cpp -L/home/daq/Program/root/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread  -pthread -lm -ldl -rdynamic -pthread -std=c++11 -m64 -I/home/daq/Program/root/include -o NuLatAB-PSD.exe

anaposAB.exe: anaposAB.cpp anafunction.h
	g++ anaposAB.cpp -L/home/daq/Program/root/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread  -pthread -lm -ldl -rdynamic -pthread -std=c++11 -m64 -I/home/daq/Program/root/include -o anaposAB.exe

#NuLatMuon.exe: NuLatMuon.cpp anafunction.h
#	g++ NuLatMuon.cpp -L/home/daq/Program/root/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread  -pthread -lm -ldl -rdynamic -pthread -std=c++11 -m64 -I/home/daq/Program/root/include -o NuLatMuon.exe

NuLatNewMuon.exe: NuLat_newMuon.cpp anafunction.h
	g++ NuLat_newMuon.cpp -L/home/daq/Program/root/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread  -pthread -lm -ldl -rdynamic -pthread -std=c++11 -m64 -I/home/daq/Program/root/include -o NuLatNewMuon.exe


NuLatABNoRoot.exe: NuLatAB-NoRoot.cpp anafunction.h
	g++ NuLatAB-NoRoot.cpp -L/home/daq/Program/root/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread  -pthread -lm -ldl -rdynamic -pthread -std=c++11 -m64 -I/home/daq/Program/root/include -o NuLatABNoRoot.exe

NuLatWindow.exe: NuLat-Window.cpp anafunction.h
	g++ NuLatAB-NoRoot.cpp -L/home/daq/Program/root/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread  -pthread -lm -ldl -rdynamic -pthread -std=c++11 -m64 -I/home/daq/Program/root/include -o NuLatWindow.exe

shell-generation.exe: shell.cpp
	g++ shell.cpp -o shell-generation.exe

Map.exe: NuLatMap.cpp
	g++ NuLatMap.cpp -o Map.exe

PSD-Shell.exe: Shell-PSD.cpp
	g++ Shell-PSD.cpp -o PSD-Shell.exe

event-list.exe: Event-All-list.cpp
	g++ Event-All-list.cpp -o event-list.exe


       
