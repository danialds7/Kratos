#!/bin/bash

# Quick Kratos Profiling Launcher
# Easy-to-use script for profiling Kratos applications

echo "=== Kratos Performance Profiling Launcher ==="
echo ""
echo "This will profile your Kratos application and generate CUDA optimization recommendations."
echo ""

# Check if running from the correct directory
if [ ! -f "MainKratos_Modif.py" ]; then
    echo "Error: MainKratos_Modif.py not found in current directory"
    echo "Please run this script from: /media/dani/Work/Kratos/Simulations/Two_Fluid_Sed/New/ConDiff_Sing_Sed_Depos"
    exit 1
fi

# Run the focused profiling
echo "Starting focused C++ profiling..."
/opt/Kratos/profile_kratos_focused.sh

# Get the latest results directory
RESULTS_BASE="profiling_results/focused_cpp"
LATEST_DIR=$(ls -t "$RESULTS_BASE" | grep "run_" | head -1)

if [ -n "$LATEST_DIR" ]; then
    echo ""
    echo "=== Generating CUDA Optimization Analysis ==="
    python3 /opt/Kratos/analyze_performance.py "$RESULTS_BASE"
    
    echo ""
    echo "=== Results Summary ==="
    echo "Profile data: $RESULTS_BASE/$LATEST_DIR"
    echo "View the following files for detailed analysis:"
    echo "- CUDA_OPTIMIZATION_REPORT.md - CUDA optimization recommendations"
    echo "- ANALYSIS_REPORT.md - General performance analysis"
    echo "- cpu_callgraph.svg - Visual call graph (open in browser)"
    echo "- top_functions.txt - Top CPU-intensive functions"
    echo ""
    echo "Next steps:"
    echo "1. Review CUDA_OPTIMIZATION_REPORT.md for optimization targets"
    echo "2. Focus on functions with high CPU usage percentage"
    echo "3. Implement GPU kernels for mathematical operations"
    echo "4. Profile again after CUDA optimization"
else
    echo "Error: No profiling results found"
fi
