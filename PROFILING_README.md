# Kratos C++ Performance Profiling Setup

This setup provides comprehensive performance profiling for your Kratos hydraulic fluid simulation application, with a focus on identifying C++ bottlenecks suitable for CUDA optimization.

## Quick Start

1. **Navigate to your project directory**:
   ```bash
   cd /media/dani/Work/Kratos/Simulations/Two_Fluid_Sed/New/ConDiff_Sing_Sed_Depos
   ```

2. **Run quick profiling**:
   ```bash
   /opt/Kratos/quick_profile.sh
   ```

This will automatically:
- Profile your application with Google Perftools
- Generate analysis reports
- Identify CUDA optimization candidates
- Create visual call graphs

## Available Profiling Scripts

### 1. `quick_profile.sh` - Recommended for most users
- Runs focused C++ profiling
- Generates CUDA optimization recommendations
- Creates comprehensive analysis reports

### 2. `profile_kratos_focused.sh` - Detailed C++ analysis
- CPU profiling with Google Perftools
- Memory profiling with Valgrind Massif
- Focused on Kratos-specific functions

### 3. `profile_kratos.sh` - Comprehensive profiling
- Multiple profiling methods (gperftools, valgrind, py-spy)
- Choose specific profiling tool or run all
- More detailed but takes longer

### 4. `analyze_performance.py` - Post-processing analysis
- Analyzes profiling results
- Identifies CUDA optimization candidates
- Generates optimization reports

## Understanding the Results

### Key Files Generated

1. **`CUDA_OPTIMIZATION_REPORT.md`** - Main recommendations for CUDA optimization
2. **`ANALYSIS_REPORT.md`** - General performance analysis
3. **`cpu_callgraph.svg`** - Visual call graph (open in browser)
4. **`top_functions.txt`** - Top CPU-intensive functions
5. **`cpu_analysis_kratos.txt`** - Kratos-specific function analysis

### What to Look For

#### High-Priority CUDA Candidates
- Functions with >5% CPU usage
- Mathematical operations (ProcessElementErosion, CalculateErosionRate)
- Loop-heavy operations
- Element/condition processing

#### Optimization Targets in Your Code
- `ProcessElementErosion` - Element erosion processing
- `ProcessElementDeposition` - Element deposition processing
- `CalculateErosionRate` - Erosion rate calculations
- `ConnectNewConditions` - Interface connectivity

## VSCode Integration

Your VSCode is configured with profiling launch configurations:

1. **"Python: MainKratos with CPU Profiling"** - Runs with CPU profiling enabled
2. **"Python: MainKratos with Valgrind"** - Runs with Valgrind profiling
3. **"C++: Attach"** - Attach debugger to running process

## Profiling Methods Available

### 1. Google CPU Profiler (gperftools)
- **Best for**: CPU usage analysis, function-level profiling
- **Output**: Call graphs, function timing, CPU percentages
- **Use when**: You need to identify CPU bottlenecks

### 2. Valgrind Callgrind
- **Best for**: Detailed function call analysis
- **Output**: Call graphs, instruction counts
- **Use when**: You need detailed call relationships

### 3. Valgrind Massif
- **Best for**: Memory usage analysis
- **Output**: Memory allocation patterns, heap usage
- **Use when**: You need to optimize memory usage

### 4. py-spy
- **Best for**: Python-specific profiling
- **Output**: Python call stacks, Python function timing
- **Use when**: You suspect Python bottlenecks

## CUDA Optimization Workflow

1. **Profile your application**:
   ```bash
   /opt/Kratos/quick_profile.sh
   ```

2. **Review results**:
   - Check `CUDA_OPTIMIZATION_REPORT.md`
   - Identify high-scoring functions
   - Focus on mathematical operations

3. **Implement CUDA kernels**:
   - Start with highest CPU usage functions
   - Focus on `hydraulic_fluid_auxiliary_utilities.cpp`
   - Implement GPU versions of mathematical operations

4. **Benchmark and validate**:
   - Compare CPU vs GPU performance
   - Measure speedup and memory usage
   - Profile again after optimization

## Example Analysis Results

```
Top 5 CUDA candidates:
1. ProcessElementErosion: 23.5% CPU (Score: 18)
2. CalculateErosionRate: 15.2% CPU (Score: 15)
3. ConnectNewConditions: 12.8% CPU (Score: 13)
4. ProcessElementDeposition: 8.7% CPU (Score: 11)
5. Mathematical operations: 6.3% CPU (Score: 9)
```

## Troubleshooting

### Common Issues

1. **No profile data generated**:
   - Check that gperftools is installed
   - Verify LD_PRELOAD is set correctly
   - Ensure profiling flags are in cmake configuration

2. **Profiling too slow**:
   - Use `quick_profile.sh` instead of comprehensive profiling
   - Reduce simulation time for profiling runs
   - Focus on specific functions

3. **Cannot open SVG files**:
   - Open SVG files in web browser
   - Install SVG viewer extension in VSCode
   - Convert to PNG: `convert file.svg file.png`

### Requirements

- Google Perftools (libgoogle-perftools-dev)
- Valgrind
- Python 3 with gprof2dot
- Graphviz (for generating graphs)

## Next Steps

1. **Immediate actions**:
   - Run profiling on your current application
   - Review the generated reports
   - Identify top 3-5 optimization targets

2. **CUDA implementation**:
   - Start with `ProcessElementErosion` function
   - Implement GPU kernel for erosion calculations
   - Add memory transfer optimizations

3. **Performance validation**:
   - Benchmark before/after CUDA implementation
   - Measure speedup and memory usage
   - Profile optimized version

## Contact

For questions about this profiling setup or CUDA optimization strategies, refer to the generated analysis reports or check the Kratos documentation.
