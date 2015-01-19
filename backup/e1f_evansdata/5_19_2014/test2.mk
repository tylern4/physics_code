ProjectName            :=test
ConfigurationName      :=Release
WorkspacePath          := "/home/tylern/DATA/code/current"
ProjectPath            := "/home/tylern/DATA/code/current"
IntermediateDirectory  :=./Release
OutDir                 := $(IntermediateDirectory)
LinkerName             :=g++
SharedObjectLinkerName :=g++ -shared -fPIC
ObjectSuffix           :=.o
DependSuffix           :=.o.d
PreprocessSuffix       :=.o.i
DebugSwitch            :=-gstab
IncludeSwitch          :=-I
LibrarySwitch          :=-l
OutputSwitch           :=-o 
LibraryPathSwitch      :=-L
PreprocessorSwitch     :=-D
SourceSwitch           :=-c 
OutputFile             :=$(IntermediateDirectory)/$(ProjectName)
ObjectSwitch           :=-o 
PreprocessOnlySwitch   :=-E 
ObjectsFileList        :="/home/tylern/DATA/code/current/test.txt"
MakeDirCommand         :=mkdir -p
LinkOptions            :=  -fPIC -m64 
IncludePath            :=  $(IncludeSwitch). $(IncludeSwitch). $(IncludeSwitch)$(ROOTSYS)/include 
Libs                   := $(LibrarySwitch)Gui $(LibrarySwitch)Core $(LibrarySwitch)Cint $(LibrarySwitch)RIO $(LibrarySwitch)Net $(LibrarySwitch)Hist $(LibrarySwitch)Graf $(LibrarySwitch)Graf3d $(LibrarySwitch)Gpad $(LibrarySwitch)Tree $(LibrarySwitch)Rint $(LibrarySwitch)Postscript $(LibrarySwitch)Matrix $(LibrarySwitch)Physics $(LibrarySwitch)MathCore $(LibrarySwitch)Thread $(LibrarySwitch)m $(LibrarySwitch)dl 
ArLibs                 :=  "Gui" "Core" "Cint" "RIO" "Net" "Hist" "Graf" "Graf3d" "Gpad" "Tree" "Rint" "Postscript" "Matrix" "Physics" "MathCore" "Thread" "m" "dl" 
LibPath                := $(LibraryPathSwitch). $(LibraryPathSwitch)$(ROOTSYS)/lib 

##
## Common variables
## AR, CXX, CC, CXXFLAGS and CFLAGS can be overriden using an environment variables
##
AR       := ar rcus
CXX      := gcc
CC       := gcc
CXXFLAGS :=  -fPIC -O2 -Wall -fopenmp $(Preprocessors) 
CFLAGS   :=  -fPIC -O2 -Wall -fopenmp $(Preprocessors) 


##
## User defined environment variables
##
ROOTSYS:=/etc/root
LD_LIBRARY_PATH:=$(ROOTSYS)/lib
Objects=$(IntermediateDirectory)/test$(ObjectSuffix) \

##
## test Build Targets 
##
.PHONY: all clean PreBuild PrePreBuild PostBuild
all: $(OutputFile)

$(OutputFile): $(IntermediateDirectory)/.d $(Objects) 
	@$(MakeDirCommand) $(@D)
	@echo "" > $(IntermediateDirectory)/.d
	@echo $(Objects) > $(ObjectsFileList)
	$(LinkerName) $(OutputSwitch)$(OutputFile) @$(ObjectsFileList) $(LibPath) $(Libs) $(LinkOptions)

$(IntermediateDirectory)/.d:
	@test -d ./Release || $(MakeDirCommand) ./Release

PreBuild:


##
## Objects
##
$(IntermediateDirectory)/test$(ObjectSuffix): test.cpp $(IntermediateDirectory)/test$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/tylern/DATA/code/current/test.cpp" $(CXXFLAGS) $(ObjectSwitch) $(IntermediateDirectory)/test$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/test$(DependSuffix): test.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/test$(ObjectSuffix) -MF$(IntermediateDirectory)/test$(DependSuffix) -MM "/home/tylern/DATA/code/current/test.cpp"

$(IntermediateDirectory)/test$(PreprocessSuffix): test.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/test$(PreprocessSuffix) "/home/tylern/DATA/code/current/test.cpp"


-include $(IntermediateDirectory)/*$(DependSuffix)
##
## Clean
##
clean:
	$(RM) $(IntermediateDirectory)/test$(ObjectSuffix)
	$(RM) $(IntermediateDirectory)/test$(DependSuffix)
	$(RM) $(IntermediateDirectory)/test$(PreprocessSuffix)
	$(RM) $(OutputFile)