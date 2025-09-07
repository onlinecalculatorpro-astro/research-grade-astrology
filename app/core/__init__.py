"""
Core research-grade astrological computation modules.

This package contains all the research-grade modules with observatory-level
precision and complete uncertainty propagation.
"""

from .timescales import TimeScales, UncertaintyBounds, build_timescales_research
from .system_integration import ResearchGradeSystem, SystemConfig, run_full_system_integration

__all__ = [
    "TimeScales",
    "UncertaintyBounds", 
    "build_timescales_research",
    "ResearchGradeSystem",
    "SystemConfig",
    "run_full_system_integration",
]
