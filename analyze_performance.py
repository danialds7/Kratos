#!/usr/bin/env python3
"""
Kratos Performance Analysis Helper
Analyzes profiling results and identifies CUDA optimization candidates
"""

import os
import sys
import re
import json
from pathlib import Path

def analyze_gperftools_profile(profile_text_file):
    """Analyze Google Perftools profile output"""
    if not os.path.exists(profile_text_file):
        return None
    
    with open(profile_text_file, 'r') as f:
        content = f.read()
    
    # Extract function performance data
    functions = []
    lines = content.split('\n')
    
    for line in lines[6:]:  # Skip header
        if line.strip() and not line.startswith('ROUTINE'):
            parts = line.split()
            if len(parts) >= 4:
                try:
                    cpu_time = float(parts[1])
                    cpu_percent = float(parts[2])
                    cum_time = float(parts[3])
                    cum_percent = float(parts[4])
                    function_name = ' '.join(parts[5:])
                    
                    functions.append({
                        'name': function_name,
                        'cpu_time': cpu_time,
                        'cpu_percent': cpu_percent,
                        'cumulative_time': cum_time,
                        'cumulative_percent': cum_percent
                    })
                except (ValueError, IndexError):
                    continue
    
    return functions

def identify_cuda_candidates(functions):
    """Identify functions that are good candidates for CUDA optimization"""
    candidates = []
    
    # Criteria for CUDA optimization
    cuda_keywords = [
        'ProcessElementErosion', 'ProcessElementDeposition', 'CalculateErosionRate',
        'ConnectNewConditions', 'matrix', 'vector', 'calculate', 'compute',
        'loop', 'iterate', 'hydraulic', 'fluid', 'auxiliary', 'utilities'
    ]
    
    for func in functions:
        score = 0
        reasons = []
        
        # High CPU usage
        if func['cpu_percent'] > 5.0:
            score += 10
            reasons.append(f"High CPU usage: {func['cpu_percent']:.1f}%")
        
        # Mathematical operations
        if any(keyword in func['name'].lower() for keyword in ['calculate', 'compute', 'process']):
            score += 5
            reasons.append("Mathematical operations")
        
        # Kratos-specific functions
        if 'kratos' in func['name'].lower() or 'hydraulic' in func['name'].lower():
            score += 3
            reasons.append("Kratos-specific function")
        
        # Loop-heavy operations
        if any(keyword in func['name'].lower() for keyword in ['for', 'while', 'iterate', 'loop']):
            score += 4
            reasons.append("Loop-heavy operations")
        
        # Element processing
        if any(keyword in func['name'].lower() for keyword in ['element', 'condition', 'node']):
            score += 3
            reasons.append("Element/condition processing")
        
        if score > 5:
            candidates.append({
                'function': func,
                'score': score,
                'reasons': reasons
            })
    
    return sorted(candidates, key=lambda x: x['score'], reverse=True)

def generate_cuda_optimization_report(results_dir):
    """Generate a comprehensive CUDA optimization report"""
    
    # Find the most recent profiling results
    result_dirs = [d for d in os.listdir(results_dir) if d.startswith('run_')]
    if not result_dirs:
        print("No profiling results found")
        return
    
    latest_dir = os.path.join(results_dir, sorted(result_dirs)[-1])
    
    # Analyze CPU profile
    cpu_analysis_file = os.path.join(latest_dir, 'cpu_analysis_full.txt')
    functions = analyze_gperftools_profile(cpu_analysis_file)
    
    if not functions:
        print("Could not analyze CPU profile")
        return
    
    # Identify CUDA candidates
    cuda_candidates = identify_cuda_candidates(functions)
    
    # Generate report
    report_file = os.path.join(latest_dir, 'CUDA_OPTIMIZATION_REPORT.md')
    
    with open(report_file, 'w') as f:
        f.write("# CUDA Optimization Report\n\n")
        f.write(f"Generated from: {latest_dir}\n")
        f.write(f"Analysis Date: {os.path.basename(latest_dir)}\n\n")
        
        f.write("## Executive Summary\n\n")
        f.write(f"- Total functions analyzed: {len(functions)}\n")
        f.write(f"- CUDA optimization candidates: {len(cuda_candidates)}\n")
        f.write(f"- Top CPU-intensive functions: {len([f for f in functions if f['cpu_percent'] > 5.0])}\n\n")
        
        f.write("## Top CUDA Optimization Candidates\n\n")
        for i, candidate in enumerate(cuda_candidates[:10], 1):
            func = candidate['function']
            f.write(f"### {i}. {func['name']}\n")
            f.write(f"- **CPU Usage**: {func['cpu_percent']:.1f}%\n")
            f.write(f"- **Cumulative Time**: {func['cumulative_percent']:.1f}%\n")
            f.write(f"- **Optimization Score**: {candidate['score']}\n")
            f.write(f"- **Reasons**: {', '.join(candidate['reasons'])}\n\n")
        
        f.write("## Implementation Priority\n\n")
        f.write("### High Priority (Immediate CUDA conversion)\n")
        high_priority = [c for c in cuda_candidates if c['score'] > 15]
        for candidate in high_priority:
            f.write(f"- {candidate['function']['name']} (Score: {candidate['score']})\n")
        
        f.write("\n### Medium Priority (Secondary CUDA conversion)\n")
        medium_priority = [c for c in cuda_candidates if 10 <= c['score'] <= 15]
        for candidate in medium_priority:
            f.write(f"- {candidate['function']['name']} (Score: {candidate['score']})\n")
        
        f.write("\n### Low Priority (Future consideration)\n")
        low_priority = [c for c in cuda_candidates if c['score'] < 10]
        for candidate in low_priority[:5]:  # Show top 5
            f.write(f"- {candidate['function']['name']} (Score: {candidate['score']})\n")
        
        f.write("\n## Detailed Analysis\n\n")
        f.write("### Function Performance Breakdown\n\n")
        for func in functions[:20]:  # Top 20 functions
            f.write(f"- **{func['name']}**: {func['cpu_percent']:.1f}% CPU, {func['cumulative_percent']:.1f}% Cumulative\n")
        
        f.write("\n## Next Steps\n\n")
        f.write("1. **Profile-Guided Optimization**:\n")
        f.write("   - Focus on functions with >5% CPU usage\n")
        f.write("   - Analyze loop structures and memory access patterns\n\n")
        
        f.write("2. **CUDA Implementation Strategy**:\n")
        f.write("   - Start with highest-scoring functions\n")
        f.write("   - Implement GPU kernels for mathematical operations\n")
        f.write("   - Optimize memory transfers between CPU and GPU\n\n")
        
        f.write("3. **Performance Validation**:\n")
        f.write("   - Benchmark CPU vs GPU implementations\n")
        f.write("   - Measure speedup and memory usage\n")
        f.write("   - Profile after CUDA optimization\n\n")
        
        f.write("## Code Locations\n\n")
        f.write("Key files to examine for CUDA optimization:\n")
        f.write("- `hydraulic_fluid_auxiliary_utilities.cpp`\n")
        f.write("- `ProcessElementErosion` and `ProcessElementDeposition` functions\n")
        f.write("- `CalculateErosionRate` function\n")
        f.write("- `ConnectNewConditions` function\n")
    
    print(f"CUDA optimization report generated: {report_file}")
    
    # Print summary to console
    print("\n=== CUDA Optimization Summary ===")
    print(f"Top 5 CUDA candidates:")
    for i, candidate in enumerate(cuda_candidates[:5], 1):
        func = candidate['function']
        print(f"{i}. {func['name']}: {func['cpu_percent']:.1f}% CPU (Score: {candidate['score']})")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 analyze_performance.py <results_directory>")
        sys.exit(1)
    
    results_dir = sys.argv[1]
    if not os.path.exists(results_dir):
        print(f"Results directory not found: {results_dir}")
        sys.exit(1)
    
    generate_cuda_optimization_report(results_dir)
