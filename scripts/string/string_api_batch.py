#!/usr/bin/env python3
"""
Query STRING database API for protein-protein interaction networks in batches.

This script uses the STRING API to:
1. Map protein identifiers to STRING database IDs
2. Retrieve protein-protein interaction networks
3. Get functional annotations and enrichment data
4. Export results in TSV format

STRING API documentation: https://string-db.org/help/api/

Author: GEM_Creole Team
Date: 2025-10-19
"""

import argparse
import json
import time
from pathlib import Path
from typing import Dict, List, Set, Tuple
from urllib.parse import quote
import sys

try:
    import requests
except ImportError:
    print("Error: 'requests' library not found. Install with: pip install requests")
    sys.exit(1)


class StringAPI:
    """
    Wrapper for STRING database REST API.

    Rate limit: Be conservative with requests to avoid blocking.
    """

    BASE_URL = "https://string-db.org/api"
    API_VERSION = "json"

    def __init__(self, species: int = 4113, min_score: int = 400):
        """
        Initialize STRING API client.

        Args:
            species: NCBI taxonomy ID (4113 for Solanum tuberosum)
            min_score: Minimum interaction score (0-1000, default 400 = medium confidence)
        """
        self.species = species
        self.min_score = min_score
        self.session = requests.Session()
        self.request_count = 0
        self.last_request_time = 0

    def _rate_limit(self, min_interval: float = 0.5):
        """
        Enforce rate limiting between API requests.

        Args:
            min_interval: Minimum seconds between requests
        """
        current_time = time.time()
        elapsed = current_time - self.last_request_time

        if elapsed < min_interval:
            sleep_time = min_interval - elapsed
            time.sleep(sleep_time)

        self.last_request_time = time.time()

    def _make_request(self, endpoint: str, params: Dict) -> Dict:
        """
        Make API request with error handling and rate limiting.

        Args:
            endpoint: API endpoint (e.g., 'network', 'enrichment')
            params: Query parameters

        Returns:
            Parsed JSON response
        """
        self._rate_limit()

        url = f"{self.BASE_URL}/{self.API_VERSION}/{endpoint}"

        try:
            response = self.session.get(url, params=params, timeout=30)
            response.raise_for_status()
            self.request_count += 1

            return response.json()

        except requests.exceptions.HTTPError as e:
            if e.response.status_code == 429:
                print(f"Rate limit exceeded. Waiting 60 seconds...")
                time.sleep(60)
                return self._make_request(endpoint, params)
            else:
                print(f"HTTP Error: {e}")
                return None

        except requests.exceptions.RequestException as e:
            print(f"Request failed: {e}")
            return None

    def get_string_ids(self, identifiers: List[str], limit: int = 1) -> Dict[str, str]:
        """
        Map protein identifiers to STRING database IDs.

        Args:
            identifiers: List of protein identifiers (gene names, UniProt IDs, etc.)
            limit: Maximum number of matches per identifier

        Returns:
            Dict mapping input_id -> STRING ID
        """
        # STRING API can handle multiple identifiers
        # But we batch in groups of 100 to be safe
        batch_size = 100
        all_mappings = {}

        for i in range(0, len(identifiers), batch_size):
            batch = identifiers[i:i+batch_size]
            identifiers_str = "\n".join(batch)

            params = {
                'identifiers': identifiers_str,
                'species': self.species,
                'limit': limit,
                'echo_query': 1
            }

            print(f"  Mapping batch {i//batch_size + 1} ({len(batch)} identifiers)...")
            response = self._make_request('get_string_ids', params)

            if response:
                for item in response:
                    query_id = item.get('queryItem', '')
                    string_id = item.get('stringId', '')

                    if string_id:
                        all_mappings[query_id] = string_id

        return all_mappings

    def get_network(self, string_ids: List[str], network_type: str = "functional") -> List[Dict]:
        """
        Get protein-protein interaction network.

        Args:
            string_ids: List of STRING protein IDs
            network_type: Type of network ('functional' or 'physical')

        Returns:
            List of interactions with scores
        """
        # STRING network endpoint can handle multiple proteins
        # Limit to reasonable batch size
        max_proteins = 200

        if len(string_ids) > max_proteins:
            print(f"  Warning: {len(string_ids)} proteins requested, limiting to {max_proteins}")
            string_ids = string_ids[:max_proteins]

        identifiers_str = "%0d".join(string_ids)

        params = {
            'identifiers': identifiers_str,
            'species': self.species,
            'required_score': self.min_score,
            'network_type': network_type
        }

        print(f"  Retrieving network for {len(string_ids)} proteins...")
        response = self._make_request('network', params)

        return response if response else []

    def get_enrichment(self, string_ids: List[str]) -> Dict:
        """
        Get functional enrichment analysis (GO, KEGG, etc.).

        Args:
            string_ids: List of STRING protein IDs

        Returns:
            Enrichment results
        """
        identifiers_str = "%0d".join(string_ids)

        params = {
            'identifiers': identifiers_str,
            'species': self.species
        }

        print(f"  Retrieving enrichment for {len(string_ids)} proteins...")
        response = self._make_request('enrichment', params)

        return response if response else {}

    def get_interaction_partners(self, string_id: str, limit: int = 10) -> List[Dict]:
        """
        Get interaction partners for a single protein.

        Args:
            string_id: STRING protein ID
            limit: Maximum number of partners to return

        Returns:
            List of interaction partners with scores
        """
        params = {
            'identifiers': string_id,
            'species': self.species,
            'required_score': self.min_score,
            'limit': limit
        }

        response = self._make_request('interaction_partners', params)

        return response if response else []


def parse_fasta_headers(fasta_file: Path) -> List[str]:
    """
    Extract protein identifiers from FASTA file.

    Returns:
        List of protein IDs
    """
    identifiers = []

    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                # Extract first token (protein ID)
                header = line[1:].strip()
                protein_id = header.split()[0]
                identifiers.append(protein_id)

    return identifiers


def write_network_tsv(interactions: List[Dict], output_file: Path, id_mapping: Dict[str, str] = None):
    """
    Write protein-protein interactions to TSV file.

    Args:
        interactions: List of interaction dictionaries from STRING API
        output_file: Output TSV file path
        id_mapping: Optional mapping of STRING IDs back to original IDs
    """
    # Create reverse mapping if provided
    reverse_map = {}
    if id_mapping:
        reverse_map = {v: k for k, v in id_mapping.items()}

    with open(output_file, 'w') as f:
        # Write header
        f.write("protein1\tprotein2\tstring_id1\tstring_id2\tscore\n")

        for interaction in interactions:
            string_id1 = interaction.get('stringId_A', '')
            string_id2 = interaction.get('stringId_B', '')
            score = interaction.get('score', 0)

            # Map back to original IDs if possible
            orig_id1 = reverse_map.get(string_id1, string_id1)
            orig_id2 = reverse_map.get(string_id2, string_id2)

            f.write(f"{orig_id1}\t{orig_id2}\t{string_id1}\t{string_id2}\t{score}\n")

    print(f"Wrote {len(interactions)} interactions to {output_file}")


def write_enrichment_tsv(enrichment_data: List[Dict], output_file: Path):
    """
    Write functional enrichment results to TSV file.

    Args:
        enrichment_data: Enrichment results from STRING API
        output_file: Output TSV file path
    """
    with open(output_file, 'w') as f:
        # Write header
        f.write("category\tterm\tdescription\tnumber_of_genes\tncbiTaxonId\tinputGenes\tpreferredNames\tp_value\tfdr\n")

        for item in enrichment_data:
            category = item.get('category', '')
            term_id = item.get('term', '')
            description = item.get('description', '')
            num_genes = item.get('number_of_genes', 0)
            taxon = item.get('ncbiTaxonId', '')
            input_genes = item.get('inputGenes', 0)
            pref_names = ','.join(item.get('preferredNames', []))
            p_value = item.get('p_value', 1.0)
            fdr = item.get('fdr', 1.0)

            f.write(f"{category}\t{term_id}\t{description}\t{num_genes}\t{taxon}\t{input_genes}\t{pref_names}\t{p_value}\t{fdr}\n")

    print(f"Wrote {len(enrichment_data)} enrichment terms to {output_file}")


def process_batch_file(
    fasta_file: Path,
    output_dir: Path,
    api_client: StringAPI,
    get_network: bool = True,
    get_enrichment: bool = False
):
    """
    Process a single batch FASTA file.

    Args:
        fasta_file: Input FASTA file
        output_dir: Output directory for results
        api_client: Initialized StringAPI client
        get_network: Whether to retrieve PPI network
        get_enrichment: Whether to retrieve functional enrichment
    """
    print(f"\nProcessing {fasta_file.name}...")

    # Parse identifiers from FASTA
    identifiers = parse_fasta_headers(fasta_file)
    print(f"  Found {len(identifiers)} protein sequences")

    # Map to STRING IDs
    print("  Mapping identifiers to STRING database...")
    id_mapping = api_client.get_string_ids(identifiers)
    print(f"  Mapped {len(id_mapping)}/{len(identifiers)} proteins to STRING IDs")

    if len(id_mapping) == 0:
        print("  Warning: No proteins mapped. Check species ID and identifier format.")
        return

    # Get STRING IDs
    string_ids = list(id_mapping.values())

    # Get network if requested
    if get_network and len(string_ids) > 0:
        print("\nRetrieving protein-protein interaction network...")
        interactions = api_client.get_network(string_ids)

        if interactions:
            network_file = output_dir / f"{fasta_file.stem}_network.tsv"
            write_network_tsv(interactions, network_file, id_mapping)
        else:
            print("  No interactions retrieved")

    # Get enrichment if requested
    if get_enrichment and len(string_ids) > 0:
        print("\nRetrieving functional enrichment...")
        enrichment = api_client.get_enrichment(string_ids)

        if enrichment:
            enrichment_file = output_dir / f"{fasta_file.stem}_enrichment.tsv"
            write_enrichment_tsv(enrichment, enrichment_file)
        else:
            print("  No enrichment data retrieved")

    # Save ID mapping
    mapping_file = output_dir / f"{fasta_file.stem}_string_ids.tsv"
    with open(mapping_file, 'w') as f:
        f.write("query_id\tstring_id\n")
        for query_id, string_id in id_mapping.items():
            f.write(f"{query_id}\t{string_id}\n")

    print(f"  Saved ID mapping to {mapping_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Query STRING database API for protein-protein interactions",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process a single batch file
  python string_api_batch.py \\
    --input data/string/batches/creole_metabolic_batch01.faa \\
    --out data/string/results \\
    --species 4113

  # Process all batch files in a directory
  python string_api_batch.py \\
    --input data/string/batches/*.faa \\
    --out data/string/results \\
    --species 4113 \\
    --get-enrichment

  # Use strict confidence threshold
  python string_api_batch.py \\
    --input data/string/batches/creole_metabolic_batch01.faa \\
    --out data/string/results \\
    --min-score 700
        """
    )

    parser.add_argument(
        '--input',
        type=Path,
        required=True,
        nargs='+',
        help='Input FASTA file(s) to process'
    )
    parser.add_argument(
        '--out',
        type=Path,
        required=True,
        help='Output directory for results'
    )
    parser.add_argument(
        '--species',
        type=int,
        default=4113,
        help='NCBI taxonomy ID (default: 4113 for Solanum tuberosum)'
    )
    parser.add_argument(
        '--min-score',
        type=int,
        default=400,
        choices=range(0, 1001),
        metavar='[0-1000]',
        help='Minimum interaction confidence score (default: 400 = medium confidence)'
    )
    parser.add_argument(
        '--get-enrichment',
        action='store_true',
        help='Retrieve functional enrichment analysis (GO, KEGG, etc.)'
    )
    parser.add_argument(
        '--network-type',
        type=str,
        default='functional',
        choices=['functional', 'physical'],
        help='Type of network to retrieve (default: functional)'
    )

    args = parser.parse_args()

    # Create output directory
    args.out.mkdir(parents=True, exist_ok=True)

    # Initialize API client
    print("="*70)
    print("STRING Database API Query Tool")
    print("="*70)
    print(f"\nSpecies: {args.species} (Solanum tuberosum)")
    print(f"Minimum score: {args.min_score}")
    print(f"Network type: {args.network_type}")
    print(f"Output directory: {args.out}")

    api_client = StringAPI(species=args.species, min_score=args.min_score)

    # Process each input file
    for fasta_file in args.input:
        if not fasta_file.exists():
            print(f"Warning: File not found: {fasta_file}")
            continue

        try:
            process_batch_file(
                fasta_file,
                args.out,
                api_client,
                get_network=True,
                get_enrichment=args.get_enrichment
            )
        except Exception as e:
            print(f"Error processing {fasta_file}: {e}")
            continue

    print("\n" + "="*70)
    print(f"Processing complete! Made {api_client.request_count} API requests")
    print(f"Results saved to: {args.out}")
    print("="*70)


if __name__ == '__main__':
    main()
