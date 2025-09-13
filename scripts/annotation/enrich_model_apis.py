#!/usr/bin/env python3
"""
Long-duration API enrichment script for SBML metabolic models.
Adds MIRIAM annotations using multiple APIs with proper rate limiting.

Author: Claude Code
Date: 2025-09-13
"""

import os
import sys
import time
import json
import logging
import argparse
import requests
import cobra
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Set
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
import xml.etree.ElementTree as ET

# Rate limiting and retry configuration
RATE_LIMITS = {
    'kegg': 0.34,      # 3 req/sec = 0.33s interval + buffer
    'chebi': 0.5,      # Conservative rate for ChEBI REST
    'pubchem': 0.2,    # PubChem is generally more permissive
    'rhea': 0.5,       # Rhea conservative rate
    'uniprot': 0.3,    # UniProt rate limit
}

BATCH_SIZES = {
    'kegg': 10,        # KEGG max batch size
    'chebi': 5,        # Conservative ChEBI batch
    'pubchem': 20,     # PubChem allows larger batches
    'rhea': 5,         # Conservative Rhea batch
    'uniprot': 10,     # UniProt batch size
}

class APIEnricher:
    """Main class for API-based model enrichment with rate limiting."""

    def __init__(self, model_path: str, output_path: str, cache_dir: str = "cache"):
        self.model_path = Path(model_path)
        self.output_path = Path(output_path)
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(exist_ok=True)

        # Load model
        self.model = cobra.io.read_sbml_model(str(self.model_path))

        # Setup logging
        self.setup_logging()

        # Setup HTTP session with retries
        self.session = self.setup_session()

        # Tracking variables
        self.stats = {
            'metabolites_processed': 0,
            'reactions_processed': 0,
            'genes_processed': 0,
            'api_calls_made': 0,
            'annotations_added': 0,
            'errors': 0,
            'start_time': time.time()
        }

        # Load existing cache
        self.load_caches()

    def setup_logging(self):
        """Setup comprehensive logging."""
        log_file = self.output_path.parent / f"enrichment_{int(time.time())}.log"

        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler(sys.stdout)
            ]
        )
        self.logger = logging.getLogger(__name__)
        self.logger.info(f"Starting API enrichment for {self.model_path}")
        self.logger.info(f"Output will be saved to {self.output_path}")

    def setup_session(self) -> requests.Session:
        """Setup HTTP session with retry strategy."""
        session = requests.Session()

        retry_strategy = Retry(
            total=3,
            backoff_factor=1,
            status_forcelist=[429, 500, 502, 503, 504],
        )

        adapter = HTTPAdapter(max_retries=retry_strategy)
        session.mount("http://", adapter)
        session.mount("https://", adapter)

        # Set user agent
        session.headers.update({
            'User-Agent': 'GEMs-Phureja-Creole/1.0 (https://github.com/user/repo) Contact: researcher@email.com'
        })

        return session

    def load_caches(self):
        """Load existing API response caches."""
        self.caches = {}
        for api in RATE_LIMITS.keys():
            cache_file = self.cache_dir / f"{api}_cache.json"
            if cache_file.exists():
                try:
                    with open(cache_file, 'r') as f:
                        self.caches[api] = json.load(f)
                    self.logger.info(f"Loaded {len(self.caches[api])} cached {api} responses")
                except Exception as e:
                    self.logger.warning(f"Failed to load {api} cache: {e}")
                    self.caches[api] = {}
            else:
                self.caches[api] = {}

    def save_cache(self, api: str):
        """Save API cache to disk."""
        cache_file = self.cache_dir / f"{api}_cache.json"
        try:
            with open(cache_file, 'w') as f:
                json.dump(self.caches[api], f, indent=2)
        except Exception as e:
            self.logger.error(f"Failed to save {api} cache: {e}")

    def rate_limit(self, api: str):
        """Apply rate limiting for specific API."""
        time.sleep(RATE_LIMITS[api])

    def api_request(self, api: str, url: str, params: dict = None) -> Optional[dict]:
        """Make rate-limited API request with caching."""
        # Create cache key
        cache_key = f"{url}_{str(sorted(params.items()) if params else '')}"

        # Check cache first
        if cache_key in self.caches[api]:
            return self.caches[api][cache_key]

        # Apply rate limiting
        self.rate_limit(api)

        try:
            response = self.session.get(url, params=params, timeout=30)
            self.stats['api_calls_made'] += 1

            if response.status_code == 200:
                data = response.json() if 'json' in response.headers.get('content-type', '') else response.text

                # Cache the response
                self.caches[api][cache_key] = data

                # Periodic cache saving
                if self.stats['api_calls_made'] % 50 == 0:
                    self.save_cache(api)

                return data
            elif response.status_code == 429:
                # Rate limit hit, wait longer
                self.logger.warning(f"Rate limit hit for {api}, waiting 60s")
                time.sleep(60)
                return self.api_request(api, url, params)  # Retry
            else:
                self.logger.warning(f"API request failed: {response.status_code} for {url}")
                return None

        except Exception as e:
            self.logger.error(f"API request error for {api}: {e}")
            self.stats['errors'] += 1
            return None

    def extract_compound_id(self, metabolite_id: str) -> str:
        """Extract base compound ID from SBML metabolite ID."""
        # Remove compartment suffix (e.g., cpd00001_c0 -> cpd00001)
        return metabolite_id.split('_')[0]

    def enrich_metabolite_kegg(self, metabolite) -> Dict[str, str]:
        """Enrich metabolite with KEGG compound data."""
        annotations = {}
        compound_id = self.extract_compound_id(metabolite.id)

        # Try multiple KEGG formats
        kegg_ids = [compound_id, compound_id.replace('cpd', 'C')]

        for kegg_id in kegg_ids:
            # Get compound info
            data = self.api_request('kegg', f'https://rest.kegg.jp/get/{kegg_id}')
            if data and isinstance(data, str) and 'ENTRY' in data:
                annotations['kegg.compound'] = kegg_id

                # Parse additional identifiers from KEGG entry
                lines = data.split('\n')
                for line in lines:
                    if line.startswith('FORMULA'):
                        # Could extract formula
                        pass
                break

        return annotations

    def enrich_metabolite_chebi(self, metabolite) -> Dict[str, str]:
        """Enrich metabolite with ChEBI data."""
        annotations = {}

        # Search ChEBI by name/synonyms
        search_terms = [metabolite.name] if metabolite.name else []

        for term in search_terms:
            if not term or len(term) < 3:
                continue

            data = self.api_request(
                'chebi',
                'https://www.ebi.ac.uk/chebi/backend/api/public/search',
                {'term': term, 'category': 'name', 'maximum_results': 5}
            )

            if data and isinstance(data, dict) and 'compounds' in data:
                compounds = data['compounds']
                if compounds:
                    # Take the first match (best score)
                    chebi_id = compounds[0].get('chebiId')
                    if chebi_id:
                        annotations['chebi'] = f"CHEBI:{chebi_id}"
                        break

        return annotations

    def enrich_metabolite_pubchem(self, metabolite) -> Dict[str, str]:
        """Enrich metabolite with PubChem data."""
        annotations = {}

        # Search PubChem by name
        if metabolite.name and len(metabolite.name) > 2:
            data = self.api_request(
                'pubchem',
                'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/json',
                {'name': metabolite.name}
            )

            if data and isinstance(data, dict) and 'PC_Compounds' in data:
                compounds = data['PC_Compounds']
                if compounds:
                    cid = compounds[0].get('id', {}).get('id', {}).get('cid')
                    if cid:
                        annotations['pubchem.compound'] = str(cid)

        return annotations

    def enrich_reaction_rhea(self, reaction) -> Dict[str, str]:
        """Enrich reaction with Rhea data."""
        annotations = {}

        # Search Rhea by reaction equation
        if reaction.reaction:
            # Simplify equation for search
            equation = str(reaction.reaction).replace(' <=> ', ' = ')

            data = self.api_request(
                'rhea',
                'https://www.rhea-db.org/rhea/rest/1.0/ws/reaction/search',
                {'query': equation, 'limit': 5}
            )

            if data and isinstance(data, dict) and 'results' in data:
                results = data['results']
                if results:
                    rhea_id = results[0].get('id')
                    if rhea_id:
                        annotations['rhea'] = rhea_id

        return annotations

    def enrich_gene_uniprot(self, gene) -> Dict[str, str]:
        """Enrich gene with UniProt data."""
        annotations = {}

        # Search UniProt by gene name for Solanum tuberosum
        gene_name = gene.name or gene.id

        data = self.api_request(
            'uniprot',
            'https://rest.uniprot.org/uniprotkb/search',
            {
                'query': f'gene:{gene_name} AND organism_id:4113',  # 4113 = Solanum tuberosum
                'format': 'json',
                'limit': 5
            }
        )

        if data and isinstance(data, dict) and 'results' in data:
            results = data['results']
            if results:
                uniprot_id = results[0].get('primaryAccession')
                if uniprot_id:
                    annotations['uniprot'] = uniprot_id

        return annotations

    def enrich_metabolites(self):
        """Enrich all metabolites with API data."""
        self.logger.info(f"Starting metabolite enrichment for {len(self.model.metabolites)} metabolites")

        for i, metabolite in enumerate(self.model.metabolites):
            if i % 100 == 0:
                self.logger.info(f"Processing metabolite {i+1}/{len(self.model.metabolites)}")

            # Skip if already well-annotated
            existing_annotations = getattr(metabolite, 'annotation', {})
            if len(existing_annotations) > 5:  # Already has many annotations
                continue

            new_annotations = {}

            # KEGG enrichment
            kegg_data = self.enrich_metabolite_kegg(metabolite)
            new_annotations.update(kegg_data)

            # ChEBI enrichment
            chebi_data = self.enrich_metabolite_chebi(metabolite)
            new_annotations.update(chebi_data)

            # PubChem enrichment
            pubchem_data = self.enrich_metabolite_pubchem(metabolite)
            new_annotations.update(pubchem_data)

            # Add new annotations to model
            if new_annotations:
                if not hasattr(metabolite, 'annotation'):
                    metabolite.annotation = {}
                metabolite.annotation.update(new_annotations)
                self.stats['annotations_added'] += len(new_annotations)

            self.stats['metabolites_processed'] += 1

            # Save progress periodically
            if i % 200 == 0 and i > 0:
                self.save_progress()

    def enrich_reactions(self):
        """Enrich all reactions with API data."""
        self.logger.info(f"Starting reaction enrichment for {len(self.model.reactions)} reactions")

        for i, reaction in enumerate(self.model.reactions):
            if i % 50 == 0:
                self.logger.info(f"Processing reaction {i+1}/{len(self.model.reactions)}")

            # Skip if already well-annotated
            existing_annotations = getattr(reaction, 'annotation', {})
            if len(existing_annotations) > 3:
                continue

            new_annotations = {}

            # Rhea enrichment
            rhea_data = self.enrich_reaction_rhea(reaction)
            new_annotations.update(rhea_data)

            # Add new annotations to model
            if new_annotations:
                if not hasattr(reaction, 'annotation'):
                    reaction.annotation = {}
                reaction.annotation.update(new_annotations)
                self.stats['annotations_added'] += len(new_annotations)

            self.stats['reactions_processed'] += 1

    def enrich_genes(self):
        """Enrich all genes with API data."""
        self.logger.info(f"Starting gene enrichment for {len(self.model.genes)} genes")

        for i, gene in enumerate(self.model.genes):
            if i % 100 == 0:
                self.logger.info(f"Processing gene {i+1}/{len(self.model.genes)}")

            # Skip if already annotated
            existing_annotations = getattr(gene, 'annotation', {})
            if existing_annotations:
                continue

            new_annotations = {}

            # UniProt enrichment
            uniprot_data = self.enrich_gene_uniprot(gene)
            new_annotations.update(uniprot_data)

            # Add new annotations to model
            if new_annotations:
                if not hasattr(gene, 'annotation'):
                    gene.annotation = {}
                gene.annotation.update(new_annotations)
                self.stats['annotations_added'] += len(new_annotations)

            self.stats['genes_processed'] += 1

    def save_progress(self):
        """Save current progress to output file."""
        try:
            cobra.io.write_sbml_model(self.model, str(self.output_path))

            # Save all caches
            for api in self.caches:
                self.save_cache(api)

            # Log progress
            elapsed = time.time() - self.stats['start_time']
            self.logger.info(f"Progress saved. Elapsed: {elapsed/3600:.1f}h, "
                           f"API calls: {self.stats['api_calls_made']}, "
                           f"Annotations added: {self.stats['annotations_added']}")
        except Exception as e:
            self.logger.error(f"Failed to save progress: {e}")

    def run_enrichment(self):
        """Run complete enrichment pipeline."""
        try:
            self.logger.info("Starting complete API enrichment pipeline")

            # Phase 1: Metabolites (most important for MEMOTE)
            self.enrich_metabolites()

            # Phase 2: Reactions
            self.enrich_reactions()

            # Phase 3: Genes
            self.enrich_genes()

            # Final save
            self.save_progress()

            # Final statistics
            elapsed = time.time() - self.stats['start_time']
            self.logger.info("="*60)
            self.logger.info("ENRICHMENT COMPLETED")
            self.logger.info(f"Total time: {elapsed/3600:.2f} hours")
            self.logger.info(f"Metabolites processed: {self.stats['metabolites_processed']}")
            self.logger.info(f"Reactions processed: {self.stats['reactions_processed']}")
            self.logger.info(f"Genes processed: {self.stats['genes_processed']}")
            self.logger.info(f"Total API calls made: {self.stats['api_calls_made']}")
            self.logger.info(f"Total annotations added: {self.stats['annotations_added']}")
            self.logger.info(f"Errors encountered: {self.stats['errors']}")
            self.logger.info(f"Output saved to: {self.output_path}")
            self.logger.info("="*60)

        except KeyboardInterrupt:
            self.logger.info("Enrichment interrupted by user. Saving progress...")
            self.save_progress()
            sys.exit(0)
        except Exception as e:
            self.logger.error(f"Enrichment failed: {e}")
            self.save_progress()
            raise


def main():
    parser = argparse.ArgumentParser(description="API-based SBML model enrichment")
    parser.add_argument("--model", required=True, help="Input SBML model path")
    parser.add_argument("--output", required=True, help="Output SBML model path")
    parser.add_argument("--cache", default="cache/api_enrichment", help="Cache directory")

    args = parser.parse_args()

    # Validate inputs
    if not Path(args.model).exists():
        print(f"Error: Input model {args.model} does not exist")
        sys.exit(1)

    # Create output directory
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)

    # Run enrichment
    enricher = APIEnricher(args.model, args.output, args.cache)
    enricher.run_enrichment()


if __name__ == "__main__":
    main()