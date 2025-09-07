# API Documentation

## Research-Grade Core Modules

### TimeScales

```python
from app.core.timescales import TimeScales, build_timescales_research

# Create research-grade time scales
time_scales = build_timescales_research(
    date_str="2023-12-21",
    time_str="12:00:00", 
    tz_name="UTC",
    dut1_seconds=0.0
)
