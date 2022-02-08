# Набор основных флагов.
CXXFLAGS :=-O3 -MMD -MP -g -std=c++11 -fopenmp -Wpedantic -Wall -Werror -Wextra -pedantic-errors -fPIC
# Флаги для линковки.
LINKER := -fopenmp -shared -fPIC

CXX = g++

# Находим все исходники.
CPP := $(shell find src -type f -name "*.cpp")

# Объектники и зависимости для swim.
OBJ := $(CPP:.cpp=.o)
DEP := $(OBJ:%.o=%.d)

# Все объектники и зависимости для clean.
ALLOBJ := $(shell find src -type f -name "*.o")
ALLDEP := $(shell find src -type f -name "*.d")

# Запустить SWIG
$(shell swig -python -o src/swim_wrap.cpp src/swim.i)


# Бинарник.
OUTPUT_SWIM = ../_swim.so

INCLUDE = -Isrc -Ideps/json/single_include $(shell pkg-config --cflags --libs python3)

all : $(OUTPUT_SWIM)
	$(shell cp src/swim.py ..)

$(OUTPUT_SWIM): $(OBJ)
	if [ ! -d bin ]; then \
	mkdir -p bin;\
	fi
	$(CXX) $(LINKER) $^ -o $@

.DELETE_ON_ERROR:

%.o: %.cpp
	@touch $@
	$(CXX) $(CXXFLAGS) $(ADDFLAGS) $(INCLUDE) -c $< -o $@

-include $(DEP)

.PHONY : clean
clean:
	rm -rf $(ALLOBJ) $(ALLDEP)
