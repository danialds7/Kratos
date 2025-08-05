#!/bin/bash

# Kratos C++ Performance Profiling Script
# This script provides multiple profiling options for Kratos applications

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="/media/dani/Work/Kratos/Simulations/Two_Fluid_Sed/New/ConDiff_Sing_Sed_Depos"
PYTHON_SCRIPT="MainKratos_Modif.py"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${GREEN}=== Kratos C++ Performance Profiling Setup ===${NC}"

# Check if project directory exists
if [ ! -d "$PROJECT_DIR" ]; then
    echo -e "${RED}Error: Project directory not found: $PROJECT_DIR${NC}"
    exit 1
fi

# Check if Python script exists
if [ ! -f "$PROJECT_DIR/$PYTHON_SCRIPT" ]; then
    echo -e "${RED}Error: Python script not found: $PROJECT_DIR/$PYTHON_SCRIPT${NC}"
    exit 1
fi

# Change to project directory
cd "$PROJECT_DIR"

# Create results directory
mkdir -p profiling_results
RESULTS_DIR="$PROJECT_DIR/profiling_results"

echo -e "${YELLOW}Choose profiling method:${NC}"
echo "1. Google CPU Profiler (gperftools)"
echo "2. Valgrind Callgrind (detailed but slow)"
echo "3. Valgrind Massif (memory profiling)"
echo "4. py-spy (Python-focused profiling)"
echo "5. perf (Linux perf profiling)"
echo "6. All methods (comprehensive analysis)"

read -p "Enter your choice (1-6): " choice

case $choice in
    1)
        echo -e "${GREEN}Running Google CPU Profiler...${NC}"
        export CPUPROFILE="$RESULTS_DIR/kratos_cpu.prof"
        export CPUPROFILE_FREQUENCY=1000
        export LD_PRELOAD="/usr/lib/x86_64-linux-gnu/libprofiler.so.0"
        python3 "$PYTHON_SCRIPT" 2>&1 | tee "$RESULTS_DIR/gperftools_output.log"
        unset LD_PRELOAD
        
        if [ -f "$RESULTS_DIR/kratos_cpu.prof" ]; then
            echo -e "${GREEN}Converting profile to text format...${NC}"
            google-pprof --text /usr/bin/python3 "$RESULTS_DIR/kratos_cpu.prof" > "$RESULTS_DIR/gperftools_analysis.txt"
            
            echo -e "${GREEN}Generating call graph...${NC}"
            google-pprof --svg /usr/bin/python3 "$RESULTS_DIR/kratos_cpu.prof" > "$RESULTS_DIR/gperftools_callgraph.svg"
            
            echo -e "${GREEN}Top functions by CPU usage:${NC}"
            google-pprof --top20 /usr/bin/python3 "$RESULTS_DIR/kratos_cpu.prof"
        fi
        ;;
    2)
        echo -e "${GREEN}Running Valgrind Callgrind...${NC}"
        valgrind --tool=callgrind --callgrind-out-file="$RESULTS_DIR/callgrind.out" \
                 --collect-jumps=yes --collect-systime=yes \
                 python3 "$PYTHON_SCRIPT" 2>&1 | tee "$RESULTS_DIR/valgrind_callgrind.log"
        
        if [ -f "$RESULTS_DIR/callgrind.out" ]; then
            echo -e "${GREEN}Converting callgrind to dot format...${NC}"
            gprof2dot --format=callgrind --output="$RESULTS_DIR/callgrind.dot" "$RESULTS_DIR/callgrind.out"
            dot -Tsvg "$RESULTS_DIR/callgrind.dot" -o "$RESULTS_DIR/callgrind_graph.svg"
            
            echo -e "${GREEN}Callgrind summary:${NC}"
            callgrind_annotate --auto=yes "$RESULTS_DIR/callgrind.out" | head -50
        fi
        ;;
    3)
        echo -e "${GREEN}Running Valgrind Massif (Memory Profiling)...${NC}"
        valgrind --tool=massif --massif-out-file="$RESULTS_DIR/massif.out" \
                 --time-unit=ms --detailed-freq=1 \
                 python3 "$PYTHON_SCRIPT" 2>&1 | tee "$RESULTS_DIR/valgrind_massif.log"
        
        if [ -f "$RESULTS_DIR/massif.out" ]; then
            echo -e "${GREEN}Memory usage summary:${NC}"
            ms_print "$RESULTS_DIR/massif.out" > "$RESULTS_DIR/massif_analysis.txt"
            head -100 "$RESULTS_DIR/massif_analysis.txt"
        fi
        ;;
    4)
        echo -e "${GREEN}Running py-spy profiling...${NC}"
        timeout 300 py-spy record -o "$RESULTS_DIR/pyspy_profile.svg" -- python3 "$PYTHON_SCRIPT" &
        PYSPY_PID=$!
        
        # Start the Python script in background
        python3 "$PYTHON_SCRIPT" 2>&1 | tee "$RESULTS_DIR/pyspy_output.log" &
        PYTHON_PID=$!
        
        # Wait for either to complete
        wait $PYSPY_PID $PYTHON_PID
        
        echo -e "${GREEN}py-spy profiling completed. Check $RESULTS_DIR/pyspy_profile.svg${NC}"
        ;;
    5)
        echo -e "${GREEN}Running perf profiling...${NC}"
        # Check if perf is available
        if command -v perf &> /dev/null; then
            perf record -g -o "$RESULTS_DIR/perf.data" python3 "$PYTHON_SCRIPT" 2>&1 | tee "$RESULTS_DIR/perf_output.log"
            
            echo -e "${GREEN}Generating perf report...${NC}"
            perf report -i "$RESULTS_DIR/perf.data" --stdio > "$RESULTS_DIR/perf_report.txt"
            
            echo -e "${GREEN}Top 20 functions by CPU usage:${NC}"
            head -50 "$RESULTS_DIR/perf_report.txt"
        else
            echo -e "${RED}perf not available. Installing...${NC}"
            sudo apt install -y linux-tools-$(uname -r) || echo -e "${RED}Could not install perf${NC}"
        fi
        ;;
    6)
        echo -e "${GREEN}Running comprehensive profiling (all methods)...${NC}"
        
        # Create timestamp for this run
        TIMESTAMP=$(date +%Y%m%d_%H%M%S)
        RUN_DIR="$RESULTS_DIR/comprehensive_$TIMESTAMP"
        mkdir -p "$RUN_DIR"
        
        echo -e "${YELLOW}1/4 Running Google CPU Profiler...${NC}"
        export CPUPROFILE="$RUN_DIR/kratos_cpu.prof"
        export CPUPROFILE_FREQUENCY=1000
        export LD_PRELOAD="/usr/lib/x86_64-linux-gnu/libprofiler.so.0"
        timeout 600 python3 "$PYTHON_SCRIPT" 2>&1 | tee "$RUN_DIR/gperftools_output.log"
        unset LD_PRELOAD
        
        if [ -f "$RUN_DIR/kratos_cpu.prof" ]; then
            google-pprof --text /usr/bin/python3 "$RUN_DIR/kratos_cpu.prof" > "$RUN_DIR/gperftools_analysis.txt"
            google-pprof --svg /usr/bin/python3 "$RUN_DIR/kratos_cpu.prof" > "$RUN_DIR/gperftools_callgraph.svg"
        fi
        
        echo -e "${YELLOW}2/4 Running py-spy profiling...${NC}"
        timeout 300 py-spy record -o "$RUN_DIR/pyspy_profile.svg" -- python3 "$PYTHON_SCRIPT" &
        
        echo -e "${YELLOW}3/4 Running Valgrind Callgrind (limited time)...${NC}"
        timeout 300 valgrind --tool=callgrind --callgrind-out-file="$RUN_DIR/callgrind.out" \
                             --collect-jumps=yes \
                             python3 "$PYTHON_SCRIPT" 2>&1 | tee "$RUN_DIR/valgrind_callgrind.log"
        
        if [ -f "$RUN_DIR/callgrind.out" ]; then
            gprof2dot --format=callgrind --output="$RUN_DIR/callgrind.dot" "$RUN_DIR/callgrind.out"
            dot -Tsvg "$RUN_DIR/callgrind.dot" -o "$RUN_DIR/callgrind_graph.svg"
        fi
        
        echo -e "${YELLOW}4/4 Running Memory profiling...${NC}"
        timeout 180 valgrind --tool=massif --massif-out-file="$RUN_DIR/massif.out" \
                             --time-unit=ms \
                             python3 "$PYTHON_SCRIPT" 2>&1 | tee "$RUN_DIR/valgrind_massif.log"
        
        if [ -f "$RUN_DIR/massif.out" ]; then
            ms_print "$RUN_DIR/massif.out" > "$RUN_DIR/massif_analysis.txt"
        fi
        
        echo -e "${GREEN}Comprehensive profiling completed in: $RUN_DIR${NC}"
        ;;
    *)
        echo -e "${RED}Invalid choice. Exiting.${NC}"
        exit 1
        ;;
esac

echo -e "${GREEN}=== Profiling Results ===${NC}"
echo "Results saved in: $RESULTS_DIR"
echo ""
echo -e "${YELLOW}Next steps for CUDA optimization:${NC}"
echo "1. Check the generated reports for CPU-intensive functions"
echo "2. Look for loops and mathematical operations that can be parallelized"
echo "3. Focus on functions in hydraulic_fluid_auxiliary_utilities.cpp"
echo "4. Check memory access patterns for GPU optimization potential"
echo ""
echo -e "${GREEN}Key files to analyze:${NC}"
echo "- Text reports: *_analysis.txt"
echo "- Call graphs: *_callgraph.svg, *_graph.svg"
echo "- Profile data: *.prof, *.out"
