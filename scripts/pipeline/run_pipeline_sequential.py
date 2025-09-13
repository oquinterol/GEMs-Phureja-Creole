#!/usr/bin/env python3
"""
Sequential model improvement pipeline script.
Runs the complete improvement pipeline from phase2 to final v1.3 model.

Author: Claude Code
Date: 2025-09-13
"""

import os
import sys
import time
import logging
import argparse
import subprocess
from pathlib import Path
from typing import List, Tuple

class SequentialPipeline:
    """Sequential model improvement pipeline manager."""

    def __init__(self, base_model: str, final_output: str, cache_dir: str = "cache"):
        self.base_model = Path(base_model)
        self.final_output = Path(final_output)
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)

        # Pipeline stages and their outputs
        self.stages = [
            {
                'name': 'Balance Mass and Charge',
                'script': 'scripts/balance_model_reactions.py',
                'input': self.base_model,
                'output': Path('models/curated/creole_v1.1_balanced.xml'),
                'description': 'Fix stoichiometric imbalances'
            },
            {
                'name': 'API Enrichment',
                'script': 'scripts/enrich_model_apis.py',
                'input': Path('models/curated/creole_v1.1_balanced.xml'),
                'output': Path('models/curated/creole_v1.2_annotated.xml'),
                'description': 'Add MIRIAM annotations via APIs'
            },
            {
                'name': 'Biomass Optimization',
                'script': 'scripts/optimize_biomass_ngam.py',
                'input': Path('models/curated/creole_v1.2_annotated.xml'),
                'output': Path('models/curated/creole_v1.3_optimized.xml'),
                'description': 'Optimize biomass and NGAM'
            },
            {
                'name': 'Final Validation',
                'script': 'scripts/validate_final_model.py',
                'input': Path('models/curated/creole_v1.3_optimized.xml'),
                'output': self.final_output,
                'description': 'Final MEMOTE validation and promotion'
            }
        ]

        self.setup_logging()

    def setup_logging(self):
        """Setup pipeline logging."""
        log_file = f"pipeline_{int(time.time())}.log"

        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler(sys.stdout)
            ]
        )
        self.logger = logging.getLogger(__name__)
        self.logger.info("Sequential Pipeline Manager initialized")

    def check_dependencies(self) -> bool:
        """Check if all required scripts and dependencies exist."""
        self.logger.info("Checking pipeline dependencies...")

        missing = []

        for stage in self.stages:
            script_path = Path(stage['script'])
            if not script_path.exists():
                missing.append(str(script_path))

        if missing:
            self.logger.error(f"Missing scripts: {missing}")
            return False

        # Check Python dependencies
        try:
            import cobra
            import requests
            self.logger.info("All dependencies satisfied")
            return True
        except ImportError as e:
            self.logger.error(f"Missing Python dependency: {e}")
            return False

    def run_stage(self, stage: dict) -> Tuple[bool, str]:
        """Run a single pipeline stage."""
        stage_name = stage['name']
        script = stage['script']
        input_file = stage['input']
        output_file = stage['output']

        self.logger.info(f"Starting stage: {stage_name}")
        self.logger.info(f"Input: {input_file}")
        self.logger.info(f"Output: {output_file}")

        # Check input exists
        if not input_file.exists():
            error_msg = f"Input file {input_file} does not exist"
            self.logger.error(error_msg)
            return False, error_msg

        # Create output directory
        output_file.parent.mkdir(parents=True, exist_ok=True)

        # Prepare command
        cmd = [
            sys.executable,
            script,
            '--model', str(input_file),
            '--output', str(output_file)
        ]

        # Add cache directory for API scripts
        if 'enrich' in script or 'api' in script:
            cmd.extend(['--cache', str(self.cache_dir / 'api_enrichment')])

        start_time = time.time()

        try:
            # Run the script
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=3600 * 6  # 6 hour timeout
            )

            elapsed = time.time() - start_time

            if result.returncode == 0:
                self.logger.info(f"Stage '{stage_name}' completed successfully in {elapsed/60:.1f} minutes")
                self.logger.info(f"Output saved to: {output_file}")

                # Log script output if verbose
                if result.stdout:
                    self.logger.debug(f"Script output:\n{result.stdout}")

                return True, f"Completed in {elapsed/60:.1f} minutes"
            else:
                error_msg = f"Stage failed with return code {result.returncode}"
                if result.stderr:
                    error_msg += f"\nError output:\n{result.stderr}"
                self.logger.error(error_msg)
                return False, error_msg

        except subprocess.TimeoutExpired:
            error_msg = f"Stage '{stage_name}' timed out after 6 hours"
            self.logger.error(error_msg)
            return False, error_msg
        except Exception as e:
            error_msg = f"Stage '{stage_name}' failed with exception: {e}"
            self.logger.error(error_msg)
            return False, error_msg

    def run_memote_check(self, model_path: Path) -> Tuple[bool, float]:
        """Run MEMOTE validation on a model."""
        self.logger.info(f"Running MEMOTE validation on {model_path}")

        if not model_path.exists():
            self.logger.error(f"Model file {model_path} does not exist")
            return False, 0.0

        # Run MEMOTE
        report_path = model_path.parent / f"{model_path.stem}_memote_report.html"

        cmd = [
            'memote', 'report', 'snapshot',
            str(model_path),
            '--filename', str(report_path)
        ]

        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)

            if result.returncode == 0:
                self.logger.info(f"MEMOTE report generated: {report_path}")

                # Try to extract score (this is simplified - real extraction would parse HTML/JSON)
                # For now, just log success
                return True, 0.0  # Placeholder score
            else:
                self.logger.warning(f"MEMOTE failed: {result.stderr}")
                return False, 0.0

        except subprocess.TimeoutExpired:
            self.logger.warning("MEMOTE timed out")
            return False, 0.0
        except Exception as e:
            self.logger.warning(f"MEMOTE error: {e}")
            return False, 0.0

    def run_fba_check(self, model_path: Path) -> bool:
        """Run FBA smoke test on a model."""
        self.logger.info(f"Running FBA smoke test on {model_path}")

        if not model_path.exists():
            return False

        cmd = [
            sys.executable,
            'scripts/fba_smoke.py',
            '--model', str(model_path),
            '--out', str(model_path.parent / f"{model_path.stem}_fba_test.log")
        ]

        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
            success = result.returncode == 0

            if success:
                self.logger.info("FBA smoke test passed")
            else:
                self.logger.warning(f"FBA smoke test failed: {result.stderr}")

            return success

        except Exception as e:
            self.logger.warning(f"FBA test error: {e}")
            return False

    def run_pipeline(self, skip_validation: bool = False) -> bool:
        """Run the complete sequential pipeline."""
        try:
            self.logger.info("="*60)
            self.logger.info("STARTING SEQUENTIAL IMPROVEMENT PIPELINE")
            self.logger.info(f"Base model: {self.base_model}")
            self.logger.info(f"Final output: {self.final_output}")
            self.logger.info("="*60)

            # Check dependencies
            if not self.check_dependencies():
                return False

            pipeline_start = time.time()
            stage_results = []

            # Run each stage
            for i, stage in enumerate(self.stages, 1):
                self.logger.info(f"\n--- STAGE {i}/{len(self.stages)}: {stage['name']} ---")
                self.logger.info(f"Description: {stage['description']}")

                success, message = self.run_stage(stage)
                stage_results.append((stage['name'], success, message))

                if not success:
                    self.logger.error(f"Pipeline failed at stage: {stage['name']}")
                    self.logger.error(f"Error: {message}")
                    return False

                # Run validation after each stage (if not skipped)
                if not skip_validation and stage['output'].exists():
                    self.logger.info("Running post-stage validation...")

                    # FBA test
                    fba_success = self.run_fba_check(stage['output'])
                    if not fba_success:
                        self.logger.warning("FBA test failed - continuing anyway")

                    # MEMOTE test (optional)
                    memote_success, score = self.run_memote_check(stage['output'])
                    if memote_success:
                        self.logger.info(f"MEMOTE validation completed")

                self.logger.info(f"Stage {i} completed successfully")

            # Pipeline summary
            total_time = time.time() - pipeline_start

            self.logger.info("\n" + "="*60)
            self.logger.info("PIPELINE COMPLETED SUCCESSFULLY")
            self.logger.info(f"Total time: {total_time/3600:.2f} hours")
            self.logger.info(f"Final model: {self.final_output}")

            self.logger.info("\nStage Summary:")
            for stage_name, success, message in stage_results:
                status = "✓" if success else "✗"
                self.logger.info(f"  {status} {stage_name}: {message}")

            self.logger.info("="*60)

            return True

        except KeyboardInterrupt:
            self.logger.info("Pipeline interrupted by user")
            return False
        except Exception as e:
            self.logger.error(f"Pipeline failed with exception: {e}")
            return False

def main():
    parser = argparse.ArgumentParser(description="Sequential model improvement pipeline")
    parser.add_argument("--model", required=True, help="Input model (phase2_emapper_expanded.xml)")
    parser.add_argument("--output", required=True, help="Final output model path")
    parser.add_argument("--cache", default="cache", help="Cache directory for API calls")
    parser.add_argument("--skip-validation", action="store_true", help="Skip MEMOTE validation steps")

    args = parser.parse_args()

    # Validate input
    if not Path(args.model).exists():
        print(f"Error: Input model {args.model} does not exist")
        sys.exit(1)

    # Create output directory
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)

    # Run pipeline
    pipeline = SequentialPipeline(args.model, args.output, args.cache)
    success = pipeline.run_pipeline(skip_validation=args.skip_validation)

    sys.exit(0 if success else 1)

if __name__ == "__main__":
    main()