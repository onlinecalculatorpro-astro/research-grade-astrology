# test_imports.py
import sys
import os

# Add app to path
sys.path.insert(0, '.')

print("Testing Research-Grade Module Imports...")
print("=" * 50)

try:
    from app.core.timescales import TimeScales, UncertaintyBounds
    print("✓ timescales.py imported successfully")
except Exception as e:
    print(f"✗ timescales.py import failed: {e}")

try:
    from app.core.time_kernel import build_timescales_research
    print("✓ time_kernel.py imported successfully")
except Exception as e:
    print(f"✗ time_kernel.py import failed: {e}")

try:
    from app.core.leapseconds import enhance_timescales_with_leap_precision
    print("✓ leapseconds.py imported successfully")
except Exception as e:
    print(f"✗ leapseconds.py import failed: {e}")

try:
    from app.core.ephemeris_adapter import ResearchEphemerisAdapter
    print("✓ ephemeris_adapter.py imported successfully")
except Exception as e:
    print(f"✗ ephemeris_adapter.py import failed: {e}")

try:
    from app.core.astronomy import compute_chart_research_grade
    print("✓ astronomy.py imported successfully")
except Exception as e:
    print(f"✗ astronomy.py import failed: {e}")

try:
    from app.core.system_integration import run_full_system_integration
    print("✓ system_integration.py imported successfully")
except Exception as e:
    print(f"✗ system_integration.py import failed: {e}")

print("\n" + "=" * 50)
print("Import test complete!")

# Test basic functionality
print("\nTesting Basic Functionality...")
try:
    from app.core.timescales import build_timescales_research
    time_scales = build_timescales_research("2023-12-21", "12:00:00", "UTC", 0.0)
    print(f"✓ TimeScales creation successful: JD_TT = {time_scales.jd_tt}")
except Exception as e:
    print(f"✗ TimeScales creation failed: {e}")
