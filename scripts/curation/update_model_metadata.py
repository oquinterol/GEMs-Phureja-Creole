#!/usr/bin/env python3
"""
Update SBML model metadata (id, metaid, name, RDF annotations) to a new version.

Uses python-libsbml for direct XML-level manipulation to preserve all existing
annotations and structure while updating version-related fields.

Usage:
    python scripts/curation/update_model_metadata.py \
        --model models/current/creole_v1.9_consistency_fixed.xml \
        --new-id CreolePotato_v2_0 \
        --new-name "Genome-scale metabolic model of Solanum tuberosum Group Phureja (Cultivar Criolla Colombia) v2.0" \
        --new-version 2.0 \
        --out models/current/CreolePotato.xml
"""
from __future__ import annotations

import argparse
import sys

import libsbml


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Update SBML model metadata to new version")
    p.add_argument("--model", required=True, help="Input SBML file")
    p.add_argument("--new-id", required=True, help="New model id (e.g. CreolePotato_v2_0)")
    p.add_argument("--new-name", required=True, help="New model name")
    p.add_argument("--new-version", required=True, help="New version string for RDF annotation")
    p.add_argument("--out", required=True, help="Output SBML file")
    return p.parse_args()


def update_rdf_version(model: libsbml.Model, old_metaid: str, new_metaid: str, new_version: str) -> None:
    """Update RDF annotation: rdf:about attribute and version identifier."""
    annotation = model.getAnnotationString()
    if not annotation:
        return

    # Update rdf:about to point to new metaid
    annotation = annotation.replace(
        f'rdf:about="#{old_metaid}"',
        f'rdf:about="#{new_metaid}"',
    )

    # Update version identifier
    import re
    annotation = re.sub(
        r'rdf:resource="https://identifiers\.org/version/[^"]*"',
        f'rdf:resource="https://identifiers.org/version/{new_version}"',
        annotation,
    )

    # Update created date to today
    from datetime import date
    annotation = re.sub(
        r'rdf:resource="https://identifiers\.org/created/[^"]*"',
        f'rdf:resource="https://identifiers.org/created/{date.today().isoformat()}"',
        annotation,
    )

    model.setAnnotation(annotation)


def main() -> None:
    args = parse_args()

    reader = libsbml.SBMLReader()
    doc = reader.readSBML(args.model)

    if doc.getNumErrors(libsbml.LIBSBML_SEV_ERROR) > 0:
        print(f"SBML read errors:\n{doc.getErrorLog().toString()}", file=sys.stderr)
        sys.exit(1)

    model = doc.getModel()
    if model is None:
        print("Error: no model found in SBML document", file=sys.stderr)
        sys.exit(1)

    old_metaid = model.getMetaId()
    new_metaid = f"meta_{args.new_id}"

    print(f"Updating model metadata:")
    print(f"  id:     {model.getId()} -> {args.new_id}")
    print(f"  metaid: {old_metaid} -> {new_metaid}")
    print(f"  name:   {model.getName()} -> {args.new_name}")

    model.setId(args.new_id)
    model.setMetaId(new_metaid)
    model.setName(args.new_name)

    update_rdf_version(model, old_metaid, new_metaid, args.new_version)

    writer = libsbml.SBMLWriter()
    success = writer.writeSBML(doc, args.out)
    if not success:
        print("Error writing SBML file", file=sys.stderr)
        sys.exit(1)

    print(f"Written: {args.out}")

    # Verify
    verify_doc = reader.readSBML(args.out)
    verify_model = verify_doc.getModel()
    assert verify_model.getId() == args.new_id, "ID verification failed"
    assert verify_model.getMetaId() == new_metaid, "MetaId verification failed"
    print("Verification passed.")


if __name__ == "__main__":
    main()
