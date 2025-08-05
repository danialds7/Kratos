#!/bin/bash

# Kratos C++ Function-Level Profiling Script
# Focus on identifying bottlenecks in hydraulic_fluid_auxiliary_utilities.cpp

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="/media/dani/Work/Kratos/Simulations/Two_Fluid_Sed/New/ConDiff_Sing_Sed_Depos"
PYTHON_SCRIPT="MainKratos_Modif.py"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${GREEN}=== Kratos C++ Function-Level Profiling ===${NC}"

# Check if project directory exists
if [ ! -d "$PROJECT_DIR" ]; then
    echo -e "${RED}Error: Project directory not found: $PROJECT_DIR${NC}"
    exit 1
fi

# Change to project directory
cd "$PROJECT_DIR"

# Create results directory
mkdir -p profiling_results/focused_cpp
RESULTS_DIR="$PROJECT_DIR/profiling_results/focused_cpp"

echo -e "${BLUE}This script will profile your Kratos application with focus on C++ bottlenecks${NC}"
echo -e "${BLUE}Specifically targeting functions in hydraulic_fluid_auxiliary_utilities.cpp${NC}"
echo ""

# Start comprehensive profiling
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
RUN_DIR="$RESULTS_DIR/run_$TIMESTAMP"
mkdir -p "$RUN_DIR"

echo -e "${YELLOW}=== Starting CPU Profiling with Google Perftools ===${NC}"
export CPUPROFILE="$RUN_DIR/kratos_cpu.prof"
export CPUPROFILE_FREQUENCY=1000  # Sample every 1000 cycles
export LD_PRELOAD="/usr/lib/x86_64-linux-gnu/libprofiler.so.0"

# Run the profiling
echo -e "${GREEN}Running Kratos application with profiling...${NC}"
python3 "$PYTHON_SCRIPT" 2>&1 | tee "$RUN_DIR/execution_log.txt"

# Cleanup
unset LD_PRELOAD

if [ -f "$RUN_DIR/kratos_cpu.prof" ]; then
    echo -e "${GREEN}=== Analyzing CPU Profile ===${NC}"
    
    # Generate text report
    echo -e "${YELLOW}Generating text analysis...${NC}"
    google-pprof --text /usr/bin/python3 "$RUN_DIR/kratos_cpu.prof" > "$RUN_DIR/cpu_analysis_full.txt"
    
    # Generate focused analysis on Kratos functions
    echo -e "${YELLOW}Generating focused Kratos analysis...${NC}"
    google-pprof --text /usr/bin/python3 "$RUN_DIR/kratos_cpu.prof" | grep -E "(Kratos|hydraulic|fluid|auxiliary)" > "$RUN_DIR/cpu_analysis_kratos.txt"
    
    # Generate call graph
    echo -e "${YELLOW}Generating call graph...${NC}"
    google-pprof --svg /usr/bin/python3 "$RUN_DIR/kratos_cpu.prof" > "$RUN_DIR/cpu_callgraph.svg"
    
    # Generate top functions report
    echo -e "${YELLOW}Generating top functions report...${NC}"
    google-pprof --top50 /usr/bin/python3 "$RUN_DIR/kratos_cpu.prof" > "$RUN_DIR/top_functions.txt"
    
    # Generate cumulative report
    echo -e "${YELLOW}Generating cumulative analysis...${NC}"
    google-pprof --cum /usr/bin/python3 "$RUN_DIR/kratos_cpu.prof" > "$RUN_DIR/cumulative_analysis.txt"
    
    echo -e "${GREEN}=== CPU Profile Analysis Complete ===${NC}"
    echo ""
    echo -e "${BLUE}Top 10 CPU-intensive functions:${NC}"
    google-pprof --top10 /usr/bin/python3 "$RUN_DIR/kratos_cpu.prof"
    echo ""
    
    echo -e "${BLUE}Kratos-specific functions:${NC}"
    google-pprof --text /usr/bin/python3 "$RUN_DIR/kratos_cpu.prof" | grep -E "(Kratos|hydraulic|fluid|auxiliary)" | head -10
    
else
    echo -e "${RED}CPU profile not generated. Check if profiling is working correctly.${NC}"
fi

echo -e "${GREEN}=== Starting Memory Profiling ===${NC}"
valgrind --tool=massif \
         --massif-out-file="$RUN_DIR/massif.out" \
         --time-unit=ms \
         --detailed-freq=1 \
         --max-snapshots=100 \
         python3 "$PYTHON_SCRIPT" 2>&1 | tee "$RUN_DIR/memory_log.txt"

if [ -f "$RUN_DIR/massif.out" ]; then
    echo -e "${GREEN}=== Analyzing Memory Profile ===${NC}"
    ms_print "$RUN_DIR/massif.out" > "$RUN_DIR/memory_analysis.txt"
    
    echo -e "${BLUE}Memory usage summary:${NC}"
    head -30 "$RUN_DIR/memory_analysis.txt"
fi

echo -e "${GREEN}=== Creating Analysis Report ===${NC}"

# Create comprehensive analysis report
cat > "$RUN_DIR/ANALYSIS_REPORT.md" << EOF
# Kratos C++ Performance Analysis Report

Generated: $(date)
Application: $PYTHON_SCRIPT

## Summary

This report analyzes the performance of the Kratos hydraulic fluid simulation application
with focus on identifying C++ bottlenecks suitable for CUDA optimization.

## Files Generated

- \`cpu_analysis_full.txt\` - Complete CPU profiling results
- \`cpu_analysis_kratos.txt\` - Kratos-specific function analysis
- \`cpu_callgraph.svg\` - Visual call graph (open in browser)
- \`top_functions.txt\` - Top 50 CPU-intensive functions
- \`cumulative_analysis.txt\` - Cumulative CPU usage analysis
- \`memory_analysis.txt\` - Memory usage profiling
- \`execution_log.txt\` - Application execution log
- \`memory_log.txt\` - Memory profiling log

## Key Areas for CUDA Optimization

### 1. CPU-Intensive Functions
Check \`top_functions.txt\` for functions with high CPU usage.
Look for:
- Mathematical computations
- Loop-heavy operations
- Matrix operations
- Iterative algorithms

### 2. Kratos-Specific Bottlenecks
Check \`cpu_analysis_kratos.txt\` for:
- HydraulicFluidAuxiliaryUtilities functions
- ProcessElementErosion/ProcessElementDeposition
- CalculateErosionRate
- ConnectNewConditions

### 3. Memory Patterns
Check \`memory_analysis.txt\` for:
- Memory allocation patterns
- Large data structures
- Potential memory access optimization

## Next Steps

1. **Identify GPU-suitable functions**: Look for functions with:
   - High CPU usage percentage
   - Mathematical operations
   - Loop structures
   - Independent data processing

2. **Analyze function characteristics**:
   - Data parallelism potential
   - Memory access patterns
   - Algorithm complexity

3. **Prioritize for CUDA conversion**:
   - Functions taking >5% CPU time
   - Loop-heavy operations
   - Mathematical calculations

## Files to Examine

1. Open \`cpu_callgraph.svg\` in browser for visual analysis
2. Review \`top_functions.txt\` for highest CPU usage
3. Check \`cpu_analysis_kratos.txt\` for Kratos-specific bottlenecks
4. Analyze \`memory_analysis.txt\` for memory optimization opportunities

EOF

echo -e "${GREEN}=== Analysis Complete ===${NC}"
echo ""
echo -e "${BLUE}Results saved in: $RUN_DIR${NC}"
echo ""
echo -e "${YELLOW}Quick Analysis Summary:${NC}"
echo "========================="

if [ -f "$RUN_DIR/cpu_analysis_full.txt" ]; then
    echo -e "${BLUE}Top 5 CPU-intensive functions:${NC}"
    head -15 "$RUN_DIR/top_functions.txt" | tail -10
    echo ""
fi

if [ -f "$RUN_DIR/cpu_analysis_kratos.txt" ]; then
    echo -e "${BLUE}Kratos-specific functions found:${NC}"
    cat "$RUN_DIR/cpu_analysis_kratos.txt" | head -5
    echo ""
fi

echo -e "${GREEN}=== Recommendations for CUDA Optimization ===${NC}"
echo "1. Focus on functions with >5% CPU usage"
echo "2. Look for mathematical operations and loops"
echo "3. Check hydraulic_fluid_auxiliary_utilities functions"
echo "4. Examine memory access patterns in top functions"
echo "5. Open cpu_callgraph.svg to visualize call relationships"
echo ""
echo -e "${YELLOW}View detailed report: $RUN_DIR/ANALYSIS_REPORT.md${NC}"
