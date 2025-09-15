#!/usr/bin/env python3
"""
Simple MEMOTE Parser - Extract section scores and calculate total

Extracts the 5 main sections from MEMOTE HTML reports:
- Consistency
- Annotation Metabolites
- Annotation Reactions
- Annotation Genes
- Annotation SBO Terms

Calculates weighted total: (3×Consistency + 1×Met + 1×Rxn + 1×Gene + 2×SBO) ÷ 800

Usage:
  python scripts/validation/simple_memote_parser.py --report reports/memote_v1.3.html
"""
import argparse
import json
import re
import sys
from typing import Dict, Optional


def extract_memote_sections(html_file: str) -> Dict[str, Optional[float]]:
    """Extract the 5 main section scores from MEMOTE HTML report."""
    try:
        with open(html_file, 'r', encoding='utf-8') as f:
            html_content = f.read()
    except Exception as e:
        print(f"Error reading HTML file: {e}")
        return {}

    # Extract window.data JSON
    window_data_pattern = r'window\.data\s*=\s*({.*?});'
    window_match = re.search(window_data_pattern, html_content, re.DOTALL)

    sections = {}

    if window_match:
        try:
            json_str = window_match.group(1)
            window_data = json.loads(json_str)

            # Look for sections in window.data.score.sections
            if 'score' in window_data and 'sections' in window_data['score']:
                section_data = window_data['score']['sections']

                for section in section_data:
                    if isinstance(section, dict) and 'section' in section and 'score' in section:
                        section_name = section['section']
                        score_value = section['score']

                        if isinstance(score_value, (int, float)) and 0 <= score_value <= 1:
                            # Convert to percentage
                            sections[section_name] = score_value * 100

        except (json.JSONDecodeError, KeyError) as e:
            print(f"Error parsing window.data: {e}")

    return sections


def calculate_memote_total(sections: Dict[str, float]) -> Optional[float]:
    """Calculate MEMOTE weighted total score."""
    # MEMOTE section weights
    weights = {
        "consistency": 3,
        "annotation_met": 1,
        "annotation_rxn": 1,
        "annotation_gene": 1,
        "annotation_sbo": 2
    }

    # Check if we have all required sections
    required_sections = list(weights.keys())
    missing_sections = [s for s in required_sections if s not in sections]

    if missing_sections:
        print(f"Warning: Missing sections: {missing_sections}")

    # Calculate weighted sum
    numerator = 0
    denominator = 0

    for section, weight in weights.items():
        if section in sections:
            numerator += weight * sections[section]
            denominator += weight * 100

    if denominator > 0:
        return numerator / denominator * 100
    else:
        return None


def format_results(html_file: str, sections: Dict[str, float], total_score: Optional[float]):
    """Format and print results in clean format."""
    print("=" * 60)
    print("MEMOTE SCORE EXTRACTION RESULTS")
    print("=" * 60)
    print(f"📁 File: {html_file}")
    print()

    if sections:
        print("📊 SECTION SCORES:")
        section_order = ["consistency", "annotation_met", "annotation_rxn", "annotation_gene", "annotation_sbo"]
        section_labels = {
            "consistency": "Consistency",
            "annotation_met": "Annotation - Metabolites",
            "annotation_rxn": "Annotation - Reactions",
            "annotation_gene": "Annotation - Genes",
            "annotation_sbo": "Annotation - SBO Terms"
        }

        for section in section_order:
            if section in sections:
                label = section_labels[section]
                value = sections[section]
                print(f"   {label:<25}: {value:6.2f}%")
            else:
                label = section_labels[section]
                print(f"   {label:<25}: NOT FOUND")
        print()

        if total_score is not None:
            print(f"🏆 TOTAL SCORE (CALCULATED): {total_score:.2f}%")

            # Show formula
            weights = [3, 1, 1, 1, 2]
            values = [sections.get(s, 0) for s in section_order]
            formula_parts = [f"{w}×{v:.2f}" for w, v in zip(weights, values)]
            formula = f"({' + '.join(formula_parts)}) ÷ 800"
            print(f"📐 Formula: {formula}")
        else:
            print("❌ Could not calculate total score (missing sections)")
    else:
        print("❌ No sections found in HTML report")

    print("=" * 60)


def main():
    parser = argparse.ArgumentParser(description="Simple MEMOTE section score extractor")
    parser.add_argument("--report", required=True, help="Path to MEMOTE HTML report")
    parser.add_argument("--json", action="store_true", help="Output as JSON instead of formatted text")

    args = parser.parse_args()

    # Extract sections
    sections = extract_memote_sections(args.report)

    if not sections:
        print("❌ Failed to extract sections from HTML file", file=sys.stderr)
        return 1

    # Calculate total score
    total_score = calculate_memote_total(sections)

    # Output results
    if args.json:
        result = {
            "file": args.report,
            "sections": sections,
            "total_score": total_score
        }
        print(json.dumps(result, indent=2))
    else:
        format_results(args.report, sections, total_score)

    return 0


if __name__ == "__main__":
    exit(main())