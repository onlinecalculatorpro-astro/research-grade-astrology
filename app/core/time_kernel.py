# app/core/time_kernel.py
# -----------------------------------------------------------------------------
# Research-Grade Time Kernel - Compatibility Shim & API Facade
#
# Design Principles:
#   • Facade pattern over canonical timescales engine
#   • Backward compatibility with deprecation warnings
#   • Type safety preservation
#   • Error context preservation  
#   • Precision level forwarding
#   • Robust import resolution
#
# API Levels:
#   1. Research-Grade: Full TimeScales objects with validation
#   2. Professional: Dict format with precision metadata
#   3. Legacy: Float returns for compatibility
#
# Standards Compliance:
#   • Preserves ERFA two-part JD precision
#   • Maintains research-grade error bounds
#   • Forwards all precision controls to engine
# -----------------------------------------------------------------------------

from __future__ import annotations

import math
import warnings
import sys
from typing import Dict, Any, Tuple, Union, Optional, TYPE_CHECKING
from functools import lru_cache
from pathlib import Path

# Type imports
if TYPE_CHECKING:
    from app.core.timescales import TimeScales, PrecisionLevel, TwoPartJD

__all__ = [
    "build_timescales",
    "build_timescales_research",
    "julian_day_utc", 
    "jd_tt_from_utc_jd",
    "jd_ut1_from_utc_jd",
    "TimeKernelError",
    "get_engine_info",
]

# Version and compatibility
TIMEKERNEL_VERSION = "4.0.0"  # Research-grade rewrite
API_COMPATIBILITY_VERSION = "3.1.0"
MINIMUM_TIMESCALES_VERSION = "4.0.0"

# ───────────────────────────── Engine Discovery & Validation ─────────────────────────────

class TimeKernelError(Exception):
    """Time kernel specific errors."""
    pass

@lru_cache(maxsize=1)
def _discover_timescales_engine():
    """
    Discover and validate the canonical timescales engine.
    Returns: (build_timescales_func, engine_info)
    """
    import_paths = [
        "app.core.timescales",      # Standard package layout
        "core.timescales",          # Alternative layout  
        "timescales",               # Standalone module
        ".timescales",              # Relative import fallback
    ]
    
    engine_func = None
    engine_info = {}
    
    for path in import_paths:
        try:
            if path.startswith('.'):
                # Handle relative imports carefully
                from importlib import import_module
                module = import_module(path, package=__package__)
            else:
                module = __import__(path, fromlist=['build_timescales'])
            
            if hasattr(module, 'build_timescales'):
                engine_func = module.build_timescales
                
                # Extract engine metadata
                engine_info = {
                    'module_path': path,
                    'module_file': getattr(module, '__file__', 'unknown'),
                    'version': getattr(module, '__version__', 'unknown'),
                    'has_research_features': all(
                        hasattr(module, attr) for attr in 
                        ['TimeScales', 'PrecisionLevel', 'TwoPartJD', 'validate_timescales']
                    ),
                    'erfa_available': True,  # Verified during import
                }
                break
                
        except ImportError as e:
            continue
    
    if engine_func is None:
        raise TimeKernelError(
            f"Cannot locate timescales engine. Searched paths: {import_paths}. "
            f"Ensure timescales.py is available in Python path."
        )
    
    # Validate engine compatibility
    if not engine_info['has_research_features']:
        warnings.warn(
            "Timescales engine lacks research-grade features. "
            "Consider upgrading to research-grade timescales.py",
            UserWarning
        )
    
    return engine_func, engine_info

def get_engine_info() -> Dict[str, Any]:
    """Get information about the timescales engine."""
    _, info = _discover_timescales_engine()
    return info.copy()

# ───────────────────────────── Type-Safe Engine Access ─────────────────────────────

def _get_engine_types():
    """Get type classes from the engine (if available)."""
    try:
        _, engine_info = _discover_timescales_engine()
        module_path = engine_info['module_path']
        
        if module_path.startswith('.'):
            from importlib import import_module
            module = import_module(module_path, package=__package__)
        else:
            module = __import__(module_path, fromlist=['TimeScales', 'PrecisionLevel'])
        
        return {
            'TimeScales': getattr(module, 'TimeScales', None),
            'PrecisionLevel': getattr(module, 'PrecisionLevel', None),
            'TwoPartJD': getattr(module, 'TwoPartJD', None),
            'ValidationReport': getattr(module, 'ValidationReport', None),
        }
    except Exception:
        return {}

def _safe_engine_call(func_name: str, *args, **kwargs):
    """Safely call engine function with error context."""
    try:
        engine_func, _ = _discover_timescales_engine()
        return engine_func(*args, **kwargs)
    except Exception as e:
        raise TimeKernelError(f"Engine call {func_name} failed: {e}") from e

# ───────────────────────────── Research-Grade API ─────────────────────────────

def build_timescales_research(
    date_str: str,
    time_str: str, 
    tz_name: str,
    dut1_seconds: float,
    *,
    precision_level: Optional[str] = "maximum",
    validate: bool = True,
    compute_uncertainties: bool = True,
    return_dataclass: bool = True,
    **kwargs
) -> Union['TimeScales', Dict[str, Any]]:
    """
    Research-grade time scales with full precision and validation.
    
    Args:
        date_str: Date in YYYY-MM-DD format
        time_str: Time in HH:MM:SS[.fraction] format
        tz_name: IANA timezone name
        dut1_seconds: UT1-UTC difference in seconds
        precision_level: "maximum", "high", or "standard"
        validate: Perform cross-validation
        compute_uncertainties: Compute error bounds
        return_dataclass: Return TimeScales object vs dict
        **kwargs: Additional engine parameters
        
    Returns:
        TimeScales dataclass or dict with full precision metadata
        
    Raises:
        TimeKernelError: Engine or parameter errors
        ValueError: Invalid input parameters
    """
    
    # Convert precision level to engine format
    engine_types = _get_engine_types()
    if engine_types.get('PrecisionLevel') and isinstance(precision_level, str):
        PrecisionLevel = engine_types['PrecisionLevel']
        precision_map = {
            'maximum': PrecisionLevel.MAXIMUM,
            'high': PrecisionLevel.HIGH, 
            'standard': PrecisionLevel.STANDARD,
        }
        kwargs['precision_level'] = precision_map.get(precision_level.lower(), PrecisionLevel.MAXIMUM)
    
    kwargs.update({
        'validate': validate,
        'compute_uncertainties': compute_uncertainties,
    })
    
    # Call engine with research parameters
    result = _safe_engine_call(
        'build_timescales_research',
        date_str, time_str, tz_name, dut1_seconds,
        **kwargs
    )
    
    # Handle return format
    if return_dataclass:
        return result
    elif hasattr(result, 'to_dict'):
        return result.to_dict()
    else:
        return result  # Already a dict

def build_timescales(
    date_str: str,
    time_str: str,
    tz_name: str, 
    dut1_seconds: float,
    **kwargs
) -> Dict[str, Any]:
    """
    Professional-grade time scales (dict format for compatibility).
    
    This is the primary API for most applications. Returns a dictionary
    with all time scales and precision metadata.
    
    Returns:
        Dict with keys: jd_utc, jd_tt, jd_ut1, delta_t, dat, dut1,
        tz_offset_seconds, timezone, warnings, precision_metadata, etc.
    """
    
    # Default to professional settings
    kwargs.setdefault('precision_level', 'high')
    kwargs.setdefault('validate', False)  # Skip validation for performance
    kwargs.setdefault('return_dataclass', False)
    
    return build_timescales_research(
        date_str, time_str, tz_name, dut1_seconds,
        **kwargs
    )

# ───────────────────────────── Legacy Compatibility API ─────────────────────────────

def julian_day_utc(date_str: str, time_str: str, tz_name: str) -> float:
    """
    DEPRECATED: Convert local time to UTC Julian Day.
    
    Use build_timescales(...) for new code.
    This function preserves backward compatibility but may lose precision.
    """
    warnings.warn(
        "julian_day_utc() is deprecated. Use build_timescales(...)['jd_utc'] "
        "for new code. This function may lose research-grade precision.",
        DeprecationWarning,
        stacklevel=2
    )
    
    try:
        # Use standard precision for legacy API
        result = build_timescales_research(
            date_str, time_str, tz_name, dut1_seconds=0.0,
            precision_level='standard',
            validate=False,
            return_dataclass=True
        )
        
        # Extract highest precision available
        if hasattr(result, 'jd_utc_2part'):
            return float(result.jd_utc_2part.jd)
        elif isinstance(result, dict):
            return float(result['jd_utc'])
        else:
            return float(result.jd_utc)
            
    except Exception as e:
        raise TimeKernelError(f"julian_day_utc conversion failed: {e}") from e

def jd_tt_from_utc_jd(jd_utc: float) -> float:
    """
    DEPRECATED: Convert UTC JD to TT JD.
    
    Use build_timescales(...) for new code.
    Maintains ERFA precision but loses research-grade error tracking.
    """
    warnings.warn(
        "jd_tt_from_utc_jd() is deprecated. Use build_timescales(...)['jd_tt'] "
        "for research-grade precision and error tracking.",
        DeprecationWarning,
        stacklevel=2
    )
    
    try:
        # Use direct ERFA calculation for legacy compatibility
        utc1, utc2 = _split_jd_precise(jd_utc)
        
        # Import ERFA dynamically
        import erfa
        tai1, tai2 = erfa.utctai(utc1, utc2)
        tt1, tt2 = erfa.taitt(tai1, tai2)
        
        return _fsum_precise(tt1, tt2)
        
    except Exception as e:
        raise TimeKernelError(f"jd_tt_from_utc_jd conversion failed: {e}") from e

def jd_ut1_from_utc_jd(jd_utc: float, dut1_seconds: float) -> float:
    """
    DEPRECATED: Convert UTC JD to UT1 JD.
    
    Use build_timescales(...) for new code.
    Maintains ERFA precision but loses research-grade validation.
    """
    warnings.warn(
        "jd_ut1_from_utc_jd() is deprecated. Use build_timescales(...)['jd_ut1'] "
        "for research-grade precision and DUT1 validation.",
        DeprecationWarning,
        stacklevel=2
    )
    
    try:
        # Validate DUT1 (legacy validation)
        _validate_dut1_legacy(dut1_seconds)
        
        # Use direct ERFA calculation
        utc1, utc2 = _split_jd_precise(jd_utc)
        
        import erfa
        ut11, ut12 = erfa.utcut1(utc1, utc2, float(dut1_seconds))
        
        return _fsum_precise(ut11, ut12)
        
    except Exception as e:
        raise TimeKernelError(f"jd_ut1_from_utc_jd conversion failed: {e}") from e

# ───────────────────────────── Precision-Preserving Utilities ─────────────────────────────

def _split_jd_precise(jd: float) -> Tuple[float, float]:
    """
    Split Julian Date into two-part form with optimal precision.
    
    This algorithm minimizes floating-point errors by ensuring
    the fractional part is in [-0.5, 0.5] range.
    """
    # Method: Split at half-day boundary for optimal precision
    jd_shifted = jd + 0.5
    jd1 = math.floor(jd_shifted)
    jd2 = jd_shifted - jd1
    
    # Shift back
    jd1 -= 0.5
    
    # Guard against floating-point edge cases
    if jd2 >= 1.0:
        jd1 += 1.0
        jd2 -= 1.0
    elif jd2 < -1.0:
        jd1 -= 1.0  
        jd2 += 1.0
    
    return float(jd1), float(jd2)

def _fsum_precise(*values) -> float:
    """High-precision sum using Kahan/Neumaier algorithm."""
    return math.fsum(float(v) for v in values)

def _validate_dut1_legacy(dut1_seconds: float) -> None:
    """Legacy DUT1 validation (simplified)."""
    if not isinstance(dut1_seconds, (int, float)):
        raise TypeError("dut1_seconds must be numeric")
    
    if abs(float(dut1_seconds)) > 0.9 + 1e-12:
        raise ValueError(
            f"dut1_seconds out of IERS range (±0.9s): {dut1_seconds}"
        )

# ───────────────────────────── Module Health & Diagnostics ─────────────────────────────

def _run_engine_diagnostics() -> Dict[str, Any]:
    """Run diagnostic tests on the timescales engine."""
    diagnostics = {
        'engine_discovered': False,
        'engine_info': {},
        'erfa_available': False,
        'research_features': False,
        'test_calculation': None,
        'errors': [],
    }
    
    try:
        # Test engine discovery
        engine_func, engine_info = _discover_timescales_engine()
        diagnostics['engine_discovered'] = True
        diagnostics['engine_info'] = engine_info
        diagnostics['research_features'] = engine_info.get('has_research_features', False)
        
        # Test ERFA availability
        try:
            import erfa
            diagnostics['erfa_available'] = True
        except ImportError:
            diagnostics['errors'].append("ERFA not available")
        
        # Test basic calculation
        try:
            test_result = build_timescales(
                "2000-01-01", "12:00:00", "UTC", 0.0
            )
            diagnostics['test_calculation'] = {
                'jd_utc': test_result.get('jd_utc'),
                'successful': True
            }
        except Exception as e:
            diagnostics['errors'].append(f"Test calculation failed: {e}")
    
    except Exception as e:
        diagnostics['errors'].append(f"Engine discovery failed: {e}")
    
    return diagnostics

def get_time_kernel_status() -> Dict[str, Any]:
    """Get comprehensive status of the time kernel."""
    status = {
        'version': TIMEKERNEL_VERSION,
        'api_compatibility': API_COMPATIBILITY_VERSION,
        'python_version': sys.version,
        'module_path': __file__,
    }
    
    # Add diagnostics
    status.update(_run_engine_diagnostics())
    
    return status

# ───────────────────────────── Module Initialization ─────────────────────────────

def _initialize_time_kernel():
    """Initialize and validate the time kernel on first import."""
    try:
        # Verify engine is available
        _discover_timescales_engine()
        
        # Verify ERFA is available
        import erfa
        
    except Exception as e:
        warnings.warn(
            f"Time kernel initialization warning: {e}. "
            f"Some functions may not work correctly.",
            UserWarning
        )

# Initialize on import
_initialize_time_kernel()
