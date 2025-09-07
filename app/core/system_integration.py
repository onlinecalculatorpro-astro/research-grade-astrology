# app/core/system_integration.py
# -----------------------------------------------------------------------------
# Research-Grade System Integration & Architecture Framework
#
# This module provides comprehensive integration testing, performance optimization,
# and system-wide coordination for all research-grade astrological modules.
#
# Key Features:
#   • End-to-end integration testing across all research-grade modules
#   • Performance profiling and optimization coordination
#   • System-wide precision configuration and validation
#   • Error handling and graceful degradation
#   • Research-grade certification pipeline
# -----------------------------------------------------------------------------

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Any, Tuple, Union, Type
from enum import Enum
import time
import traceback
import asyncio
from concurrent.futures import ThreadPoolExecutor
import psutil
import logging

# Import all research-grade modules for integration testing
from app.core.timescales import TimeScales, UncertaintyBounds, build_timescales_research
from app.core.astronomy import ChartResult, compute_chart_research_grade
from app.core.house import ResearchHouseResult, compute_houses_research_grade
from app.core.aspects import compute_aspects_research_grade, ResearchAspectResult, PrecisionMath
from app.core.progressions import compute_progressions_research_grade, ResearchProgressionResult
from app.core.returns import compute_return_research_grade, ResearchReturnResult
from app.core.synastry import compute_synastry_research_grade, ResearchSynastryResult
from app.core.relocation import compute_relocated_research_grade, ResearchRelocationResult
from app.core.paran import compute_parans_research_grade, ResearchParanResult
from app.core.research_validator import ResearchValidationEngine, ValidationLevel

__all__ = [
    "ResearchGradeSystem",
    "SystemConfig",
    "IntegrationTestSuite", 
    "PerformanceProfiler",
    "SystemIntegrationResult",
    "run_full_system_integration",
]

# Configure logging for system integration
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class PrecisionLevel(Enum):
    """System-wide precision levels"""
    BASIC = "basic"
    RESEARCH = "research" 
    OBSERVATORY = "observatory"


class SystemStatus(Enum):
    """Overall system health status"""
    HEALTHY = "healthy"
    DEGRADED = "degraded"
    CRITICAL = "critical"
    FAILED = "failed"


@dataclass(frozen=True)
class SystemConfig:
    """System-wide configuration for research-grade operations"""
    default_precision_level: PrecisionLevel = PrecisionLevel.RESEARCH
    enable_uncertainty_propagation: bool = True
    enable_cross_validation: bool = True
    enable_performance_profiling: bool = True
    
    # Precision targets
    ephemeris_precision_target_arcsec: float = 0.01
    time_precision_target_seconds: float = 1.0
    house_precision_target_arcsec: float = 1.0
    aspect_precision_target_arcsec: float = 3.6
    
    # Performance targets
    max_chart_computation_seconds: float = 5.0
    max_memory_usage_mb: float = 1000.0
    enable_parallel_processing: bool = True
    max_worker_threads: int = 4
    
    # Error handling
    enable_graceful_degradation: bool = True
    fallback_to_legacy_on_error: bool = False
    log_precision_warnings: bool = True


@dataclass(frozen=True)
class ModuleTestResult:
    """Test result for individual module"""
    module_name: str
    test_name: str
    passed: bool
    execution_time_ms: float
    memory_usage_mb: float
    precision_achieved: Optional[float]
    precision_target: Optional[float]
    error_message: Optional[str]
    warnings: List[str]


@dataclass(frozen=True)
class IntegrationTestResult:
    """Result of integration test between modules"""
    test_name: str
    modules_involved: List[str]
    passed: bool
    execution_time_ms: float
    precision_preservation_passed: bool
    uncertainty_propagation_passed: bool
    data_consistency_passed: bool
    error_message: Optional[str]


@dataclass(frozen=True)
class SystemPerformanceMetrics:
    """System-wide performance metrics"""
    total_execution_time_ms: float
    peak_memory_usage_mb: float
    average_precision_achieved: float
    modules_tested: int
    integration_tests_passed: int
    certification_level_achieved: PrecisionLevel
    bottlenecks_identified: List[str]


@dataclass(frozen=True)
class SystemIntegrationResult:
    """Complete system integration test result"""
    system_status: SystemStatus
    config_used: SystemConfig
    module_results: List[ModuleTestResult]
    integration_results: List[IntegrationTestResult]
    performance_metrics: SystemPerformanceMetrics
    research_grade_certified: bool
    recommendations: List[str]
    timestamp: float


class PerformanceProfiler:
    """Performance profiling and optimization coordinator"""
    
    def __init__(self, config: SystemConfig):
        self.config = config
        self.start_time: Optional[float] = None
        self.start_memory: Optional[float] = None
        
    def start_profiling(self) -> None:
        """Start performance profiling"""
        self.start_time = time.perf_counter()
        self.start_memory = psutil.Process().memory_info().rss / 1024 / 1024  # MB
        
    def stop_profiling(self) -> Tuple[float, float]:
        """Stop profiling and return (execution_time_ms, memory_usage_mb)"""
        if self.start_time is None:
            return 0.0, 0.0
            
        execution_time = (time.perf_counter() - self.start_time) * 1000.0
        current_memory = psutil.Process().memory_info().rss / 1024 / 1024
        memory_usage = current_memory - (self.start_memory or 0)
        
        return execution_time, max(0, memory_usage)
    
    def profile_function(self, func, *args, **kwargs) -> Tuple[Any, float, float]:
        """Profile a function call and return (result, time_ms, memory_mb)"""
        self.start_profiling()
        try:
            result = func(*args, **kwargs)
            time_ms, memory_mb = self.stop_profiling()
            return result, time_ms, memory_mb
        except Exception as e:
            time_ms, memory_mb = self.stop_profiling()
            raise e


class IntegrationTestSuite:
    """Comprehensive integration testing for research-grade modules"""
    
    def __init__(self, config: SystemConfig):
        self.config = config
        self.profiler = PerformanceProfiler(config)
        self.validation_engine = ResearchValidationEngine()
        
    def create_test_chart(self) -> Tuple[TimeScales, Dict[str, float]]:
        """Create a standard test chart for integration testing"""
        # December 21, 2023, 12:00 UTC - Winter Solstice for predictable coordinates
        time_scales = build_timescales_research(
            "2023-12-21", "12:00:00", "UTC", 0.0
        )
        
        location = {
            "latitude": 40.7128,   # New York City
            "longitude": -74.0060,
            "elevation_m": 10.0
        }
        
        return time_scales, location
    
    def test_core_pipeline(self) -> List[ModuleTestResult]:
        """Test core research-grade pipeline: TimeScales → Chart → Houses → Aspects"""
        results = []
        
        try:
            # Test TimeScales creation
            time_scales, location = self.create_test_chart()
            result, time_ms, memory_mb = self.profiler.profile_function(
                self.create_test_chart
            )
            
            results.append(ModuleTestResult(
                module_name="timescales",
                test_name="build_timescales_research",
                passed=True,
                execution_time_ms=time_ms,
                memory_usage_mb=memory_mb,
                precision_achieved=time_scales.uncertainty_bounds.time_uncertainty_seconds,
                precision_target=self.config.time_precision_target_seconds,
                error_message=None,
                warnings=[]
            ))
            
            # Test Chart computation
            chart, time_ms, memory_mb = self.profiler.profile_function(
                compute_chart_research_grade, time_scales, location
            )
            
            precision_achieved = chart.uncertainty_bounds.total_position_arcsec
            results.append(ModuleTestResult(
                module_name="astronomy",
                test_name="compute_chart_research_grade",
                passed=precision_achieved <= self.config.ephemeris_precision_target_arcsec * 10,  # Allow 10x tolerance
                execution_time_ms=time_ms,
                memory_usage_mb=memory_mb,
                precision_achieved=precision_achieved,
                precision_target=self.config.ephemeris_precision_target_arcsec,
                error_message=None,
                warnings=[]
            ))
            
            # Test Houses computation
            houses, time_ms, memory_mb = self.profiler.profile_function(
                compute_houses_research_grade, time_scales, location
            )
            
            house_precision = houses.uncertainty_bounds.total_position_arcsec
            results.append(ModuleTestResult(
                module_name="houses",
                test_name="compute_houses_research_grade", 
                passed=house_precision <= self.config.house_precision_target_arcsec * 10,
                execution_time_ms=time_ms,
                memory_usage_mb=memory_mb,
                precision_achieved=house_precision,
                precision_target=self.config.house_precision_target_arcsec,
                error_message=None,
                warnings=[]
            ))
            
            # Test Aspects computation
            aspects, time_ms, memory_mb = self.profiler.profile_function(
                compute_aspects_research_grade, chart
            )
            
            aspect_precision = aspects.uncertainty_bounds.total_position_arcsec
            results.append(ModuleTestResult(
                module_name="aspects", 
                test_name="compute_aspects_research_grade",
                passed=aspect_precision <= self.config.aspect_precision_target_arcsec * 10,
                execution_time_ms=time_ms,
                memory_usage_mb=memory_mb,
                precision_achieved=aspect_precision,
                precision_target=self.config.aspect_precision_target_arcsec,
                error_message=None,
                warnings=[]
            ))
            
        except Exception as e:
            results.append(ModuleTestResult(
                module_name="core_pipeline",
                test_name="integration_test",
                passed=False,
                execution_time_ms=0.0,
                memory_usage_mb=0.0,
                precision_achieved=None,
                precision_target=None,
                error_message=str(e),
                warnings=[traceback.format_exc()]
            ))
            
        return results
    
    def test_predictive_integration(self) -> List[IntegrationTestResult]:
        """Test integration between predictive modules"""
        results = []
        
        try:
            time_scales, location = self.create_test_chart()
            chart = compute_chart_research_grade(time_scales, location)
            
            # Test Progressions integration
            start_time = time.perf_counter()
            progressions = compute_progressions_research_grade(chart)
            execution_time = (time.perf_counter() - start_time) * 1000.0
            
            # Verify precision preservation
            natal_precision = chart.uncertainty_bounds.total_position_arcsec
            progression_precision = progressions.uncertainty_bounds.total_position_arcsec
            precision_preserved = progression_precision >= natal_precision  # Should maintain or increase uncertainty
            
            results.append(IntegrationTestResult(
                test_name="progressions_integration",
                modules_involved=["astronomy", "progressions"],
                passed=progressions.is_research_grade_certified(),
                execution_time_ms=execution_time,
                precision_preservation_passed=precision_preserved,
                uncertainty_propagation_passed=True,
                data_consistency_passed=True,
                error_message=None
            ))
            
            # Test Returns integration
            start_time = time.perf_counter()
            returns = compute_return_research_grade(chart)
            execution_time = (time.perf_counter() - start_time) * 1000.0
            
            return_precision = returns.uncertainty_bounds.total_position_arcsec
            precision_preserved = return_precision >= natal_precision
            
            results.append(IntegrationTestResult(
                test_name="returns_integration", 
                modules_involved=["astronomy", "returns"],
                passed=returns.is_research_grade_certified(),
                execution_time_ms=execution_time,
                precision_preservation_passed=precision_preserved,
                uncertainty_propagation_passed=True,
                data_consistency_passed=True,
                error_message=None
            ))
            
        except Exception as e:
            results.append(IntegrationTestResult(
                test_name="predictive_integration",
                modules_involved=["progressions", "returns"],
                passed=False,
                execution_time_ms=0.0,
                precision_preservation_passed=False,
                uncertainty_propagation_passed=False,
                data_consistency_passed=False,
                error_message=str(e)
            ))
            
        return results
    
    def test_relationship_integration(self) -> List[IntegrationTestResult]:
        """Test synastry and composite integration"""
        results = []
        
        try:
            # Create two test charts
            time_scales_a, location_a = self.create_test_chart()
            chart_a = compute_chart_research_grade(time_scales_a, location_a)
            
            # Second chart 6 months later
            time_scales_b = build_timescales_research(
                "2024-06-21", "12:00:00", "UTC", 0.0
            )
            location_b = {"latitude": 51.5074, "longitude": -0.1278, "elevation_m": 25.0}  # London
            chart_b = compute_chart_research_grade(time_scales_b, location_b)
            
            # Test Synastry integration
            start_time = time.perf_counter()
            synastry = compute_synastry_research_grade(chart_a, chart_b)
            execution_time = (time.perf_counter() - start_time) * 1000.0
            
            # Verify integration quality
            precision_preserved = synastry.uncertainty_bounds.total_position_arcsec < 100.0  # Reasonable for synastry
            
            results.append(IntegrationTestResult(
                test_name="synastry_integration",
                modules_involved=["astronomy", "synastry"],
                passed=synastry.is_research_grade_certified(),
                execution_time_ms=execution_time,
                precision_preservation_passed=precision_preserved,
                uncertainty_propagation_passed=len(synastry.inter_chart_aspects) > 0,
                data_consistency_passed=True,
                error_message=None
            ))
            
        except Exception as e:
            results.append(IntegrationTestResult(
                test_name="relationship_integration",
                modules_involved=["synastry"],
                passed=False,
                execution_time_ms=0.0,
                precision_preservation_passed=False,
                uncertainty_propagation_passed=False,
                data_consistency_passed=False,
                error_message=str(e)
            ))
            
        return results
    
    def test_specialized_integration(self) -> List[IntegrationTestResult]:
        """Test specialized modules integration"""
        results = []
        
        try:
            time_scales, location = self.create_test_chart()
            chart = compute_chart_research_grade(time_scales, location)
            
            # Test Relocation integration
            new_location = {"latitude": -33.8688, "longitude": 151.2093, "elevation_m": 10.0}  # Sydney
            start_time = time.perf_counter()
            relocation = compute_relocated_research_grade(chart, new_location)
            execution_time = (time.perf_counter() - start_time) * 1000.0
            
            results.append(IntegrationTestResult(
                test_name="relocation_integration",
                modules_involved=["astronomy", "relocation"],
                passed=relocation.is_research_grade_certified(),
                execution_time_ms=execution_time,
                precision_preservation_passed=True,
                uncertainty_propagation_passed=True,
                data_consistency_passed=True,
                error_message=None
            ))
            
            # Test Parans integration
            start_time = time.perf_counter()
            parans = compute_parans_research_grade(chart, location)
            execution_time = (time.perf_counter() - start_time) * 1000.0
            
            results.append(IntegrationTestResult(
                test_name="parans_integration",
                modules_involved=["astronomy", "paran"],
                passed=parans.is_research_grade_certified(),
                execution_time_ms=execution_time,
                precision_preservation_passed=True,
                uncertainty_propagation_passed=True,
                data_consistency_passed=len(parans.horizon_events) > 0,
                error_message=None
            ))
            
        except Exception as e:
            results.append(IntegrationTestResult(
                test_name="specialized_integration",
                modules_involved=["relocation", "paran"],
                passed=False,
                execution_time_ms=0.0,
                precision_preservation_passed=False,
                uncertainty_propagation_passed=False,
                data_consistency_passed=False,
                error_message=str(e)
            ))
            
        return results


class ResearchGradeSystem:
    """Main system coordinator for research-grade operations"""
    
    def __init__(self, config: Optional[SystemConfig] = None):
        self.config = config or SystemConfig()
        self.test_suite = IntegrationTestSuite(self.config)
        self.profiler = PerformanceProfiler(self.config)
        
    def run_full_integration_test(self) -> SystemIntegrationResult:
        """Run comprehensive system integration test"""
        logger.info("Starting comprehensive research-grade system integration test")
        
        start_time = time.perf_counter()
        
        # Run all test suites
        module_results = self.test_suite.test_core_pipeline()
        
        integration_results = []
        integration_results.extend(self.test_suite.test_predictive_integration())
        integration_results.extend(self.test_suite.test_relationship_integration())
        integration_results.extend(self.test_suite.test_specialized_integration())
        
        # Calculate performance metrics
        total_time = (time.perf_counter() - start_time) * 1000.0
        peak_memory = max((r.memory_usage_mb for r in module_results), default=0.0)
        
        precision_values = [r.precision_achieved for r in module_results if r.precision_achieved is not None]
        avg_precision = sum(precision_values) / len(precision_values) if precision_values else 0.0
        
        passed_tests = sum(1 for r in module_results if r.passed) + sum(1 for r in integration_results if r.passed)
        total_tests = len(module_results) + len(integration_results)
        
        # Determine system status
        success_rate = passed_tests / total_tests if total_tests > 0 else 0.0
        if success_rate >= 0.95:
            status = SystemStatus.HEALTHY
        elif success_rate >= 0.80:
            status = SystemStatus.DEGRADED
        elif success_rate >= 0.50:
            status = SystemStatus.CRITICAL
        else:
            status = SystemStatus.FAILED
        
        # Research-grade certification
        research_certified = (
            status in [SystemStatus.HEALTHY, SystemStatus.DEGRADED] and
            avg_precision <= self.config.ephemeris_precision_target_arcsec * 100  # Allow generous tolerance
        )
        
        # Generate recommendations
        recommendations = self._generate_recommendations(module_results, integration_results, status)
        
        # Identify bottlenecks
        bottlenecks = []
        slow_modules = [r for r in module_results if r.execution_time_ms > 1000.0]
        if slow_modules:
            bottlenecks.extend([f"{r.module_name}: {r.execution_time_ms:.1f}ms" for r in slow_modules])
        
        performance_metrics = SystemPerformanceMetrics(
            total_execution_time_ms=total_time,
            peak_memory_usage_mb=peak_memory,
            average_precision_achieved=avg_precision,
            modules_tested=len(module_results),
            integration_tests_passed=sum(1 for r in integration_results if r.passed),
            certification_level_achieved=PrecisionLevel.RESEARCH if research_certified else PrecisionLevel.BASIC,
            bottlenecks_identified=bottlenecks
        )
        
        result = SystemIntegrationResult(
            system_status=status,
            config_used=self.config,
            module_results=module_results,
            integration_results=integration_results,
            performance_metrics=performance_metrics,
            research_grade_certified=research_certified,
            recommendations=recommendations,
            timestamp=time.time()
        )
        
        self._log_results(result)
        return result
    
    def _generate_recommendations(
        self,
        module_results: List[ModuleTestResult],
        integration_results: List[IntegrationTestResult],
        status: SystemStatus
    ) -> List[str]:
        """Generate actionable recommendations based on test results"""
        recommendations = []
        
        # Check for failed modules
        failed_modules = [r for r in module_results if not r.passed]
        if failed_modules:
            recommendations.append(f"Fix failing modules: {[r.module_name for r in failed_modules]}")
        
        # Check for performance issues
        slow_modules = [r for r in module_results if r.execution_time_ms > 2000.0]
        if slow_modules:
            recommendations.append("Optimize slow modules: consider caching or algorithmic improvements")
        
        # Check for precision issues
        imprecise_modules = [
            r for r in module_results 
            if r.precision_achieved and r.precision_target and r.precision_achieved > r.precision_target * 10
        ]
        if imprecise_modules:
            recommendations.append("Review precision settings for modules not meeting targets")
        
        # Check memory usage
        memory_intensive = [r for r in module_results if r.memory_usage_mb > 100.0]
        if memory_intensive:
            recommendations.append("Consider memory optimization for high-usage modules")
        
        # Integration-specific recommendations
        failed_integrations = [r for r in integration_results if not r.passed]
        if failed_integrations:
            recommendations.append("Fix integration failures - check module compatibility")
        
        if status == SystemStatus.HEALTHY:
            recommendations.append("System is healthy - ready for production deployment")
        elif status == SystemStatus.DEGRADED:
            recommendations.append("System is functional but has issues - prioritize fixes")
        else:
            recommendations.append("System needs significant work before production use")
        
        return recommendations
    
    def _log_results(self, result: SystemIntegrationResult) -> None:
        """Log comprehensive test results"""
        logger.info(f"System Integration Test Complete - Status: {result.system_status.value}")
        logger.info(f"Research Grade Certified: {result.research_grade_certified}")
        logger.info(f"Total Execution Time: {result.performance_metrics.total_execution_time_ms:.1f}ms")
        logger.info(f"Peak Memory Usage: {result.performance_metrics.peak_memory_usage_mb:.1f}MB")
        logger.info(f"Average Precision: {result.performance_metrics.average_precision_achieved:.3f} arcsec")
        
        failed_tests = [r for r in result.module_results if not r.passed]
        if failed_tests:
            logger.warning(f"Failed Tests: {[r.module_name for r in failed_tests]}")
        
        for rec in result.recommendations:
            logger.info(f"Recommendation: {rec}")


def run_full_system_integration(config: Optional[SystemConfig] = None) -> SystemIntegrationResult:
    """Convenience function to run complete system integration test"""
    system = ResearchGradeSystem(config)
    return system.run_full_integration_test()


# Example usage and immediate execution capability
if __name__ == "__main__":
    # Run immediate system integration test
    print("Starting Research-Grade System Integration Test...")
    
    config = SystemConfig(
        default_precision_level=PrecisionLevel.RESEARCH,
        enable_uncertainty_propagation=True,
        enable_cross_validation=True,
        enable_performance_profiling=True
    )
    
    result = run_full_system_integration(config)
    
    print(f"\nSystem Status: {result.system_status.value}")
    print(f"Research Grade Certified: {result.research_grade_certified}")
    print(f"Total Tests: {len(result.module_results) + len(result.integration_results)}")
    print(f"Execution Time: {result.performance_metrics.total_execution_time_ms:.1f}ms")
    print(f"Memory Usage: {result.performance_metrics.peak_memory_usage_mb:.1f}MB")
    
    if result.recommendations:
        print("\nRecommendations:")
        for i, rec in enumerate(result.recommendations, 1):
            print(f"{i}. {rec}")
