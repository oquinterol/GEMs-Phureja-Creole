#!/usr/bin/env python3
"""
KEGG API Client with proper rate limiting and best practices.

Based on KEGG API usage guidelines:
- Maximum 3 requests per second
- Maximum 10 IDs per batch in get/list operations
- Exponential backoff for 429/503 errors
- Local caching to avoid repeated calls
- Solanum tuberosum organism code: 'sot'

Usage:
  from scripts.kegg_api_client import KEGGClient
  
  client = KEGGClient()
  pathways = client.list_pathways('sot')
  reactions = client.get_reactions(['R00001', 'R00002'])
"""
from __future__ import annotations

import json
import time
import itertools
from pathlib import Path
from typing import Dict, List, Optional, Iterable
import requests
import requests_cache


class KEGGClient:
    """KEGG REST API client with rate limiting and best practices."""
    
    def __init__(self, max_rps: float = 3.0, cache_dir: Optional[str] = None, 
                 max_retries: int = 5, backoff_factor: float = 1.5):
        """
        Initialize KEGG client.
        
        Args:
            max_rps: Maximum requests per second (KEGG limit: 3)
            cache_dir: Directory for local caching (None = no cache)
            max_retries: Maximum retry attempts for failed requests
            backoff_factor: Exponential backoff factor
        """
        self.base_url = "https://rest.kegg.jp"
        self.max_rps = max_rps
        self.max_retries = max_retries
        self.backoff_factor = backoff_factor
        self.last_request_time = 0.0
        self.request_count = 0
        
        # Setup caching if requested
        if cache_dir:
            cache_path = Path(cache_dir)
            cache_path.mkdir(parents=True, exist_ok=True)
            self.session = requests_cache.CachedSession(
                cache_name=str(cache_path / "kegg_cache"),
                expire_after=86400  # 24 hours
            )
        else:
            self.session = requests.Session()
    
    def _rate_limit(self) -> None:
        """Enforce rate limiting to respect KEGG's 3 req/s limit."""
        min_interval = 1.0 / self.max_rps
        current_time = time.time()
        time_since_last = current_time - self.last_request_time
        
        if time_since_last < min_interval:
            sleep_time = min_interval - time_since_last
            time.sleep(sleep_time)
        
        self.last_request_time = time.time()
        self.request_count += 1
    
    def _request(self, endpoint: str, params: Optional[Dict] = None) -> str:
        """Make rate-limited request with exponential backoff."""
        url = f"{self.base_url}{endpoint}"
        
        for attempt in range(self.max_retries):
            self._rate_limit()
            
            try:
                response = self.session.get(url, params=params, timeout=30)
                
                if response.status_code == 200:
                    return response.text
                elif response.status_code in (429, 503):
                    # Rate limit or service busy - apply exponential backoff
                    wait_time = self.backoff_factor ** attempt
                    print(f"Rate limited (HTTP {response.status_code}), waiting {wait_time:.1f}s...")
                    time.sleep(wait_time)
                    continue
                else:
                    response.raise_for_status()
                    
            except requests.RequestException as e:
                if attempt == self.max_retries - 1:
                    raise
                wait_time = self.backoff_factor ** attempt
                print(f"Request failed (attempt {attempt + 1}): {e}, retrying in {wait_time:.1f}s...")
                time.sleep(wait_time)
        
        raise RuntimeError(f"Failed to fetch {endpoint} after {self.max_retries} attempts")
    
    def _chunk_ids(self, ids: Iterable[str], chunk_size: int = 10) -> Iterable[List[str]]:
        """Split IDs into chunks respecting KEGG's 10 ID limit."""
        iterator = iter(ids)
        while True:
            chunk = list(itertools.islice(iterator, chunk_size))
            if not chunk:
                break
            yield chunk
    
    def _parse_tsv(self, text: str) -> List[List[str]]:
        """Parse tab-separated values from KEGG response."""
        rows = []
        for line in text.strip().splitlines():
            if line:
                parts = line.split('\t')
                if len(parts) >= 2:
                    rows.append(parts)
        return rows
    
    def _parse_kegg_flat_file(self, text: str) -> Dict[str, str]:
        """Parse KEGG flat file format."""
        data = {}
        current_key = None
        
        for line in text.splitlines():
            if line and line[0] != ' ':
                # New field
                if ' ' in line:
                    key = line.split(' ', 1)[0]
                    value = line.split(' ', 1)[1] if len(line.split(' ', 1)) > 1 else ''
                    data[key] = value.strip()
                    current_key = key
            elif current_key and line.startswith('            '):
                # Continuation line
                data[current_key] += ' ' + line.strip()
        
        return data
    
    # Core KEGG operations
    
    def list_pathways(self, organism: str = 'sot') -> List[Dict[str, str]]:
        """List pathways for organism (default: Solanum tuberosum 'sot')."""
        response = self._request(f"/list/pathway/{organism}")
        rows = self._parse_tsv(response)
        
        pathways = []
        for row in rows:
            pathways.append({
                'id': row[0].split(':')[-1] if ':' in row[0] else row[0],
                'name': row[1] if len(row) > 1 else '',
                'full_id': row[0]
            })
        
        return pathways
    
    def list_compounds(self) -> List[Dict[str, str]]:
        """List all KEGG compounds."""
        response = self._request("/list/cpd")
        rows = self._parse_tsv(response)
        
        compounds = []
        for row in rows:
            compounds.append({
                'id': row[0].split(':')[-1] if ':' in row[0] else row[0],
                'name': row[1] if len(row) > 1 else '',
                'full_id': row[0]
            })
        
        return compounds
    
    def list_reactions(self) -> List[Dict[str, str]]:
        """List all KEGG reactions."""
        response = self._request("/list/rn")
        rows = self._parse_tsv(response)
        
        reactions = []
        for row in rows:
            reactions.append({
                'id': row[0].split(':')[-1] if ':' in row[0] else row[0],
                'name': row[1] if len(row) > 1 else '',
                'full_id': row[0]
            })
        
        return reactions
    
    def get_reactions(self, reaction_ids: List[str]) -> List[Dict[str, str]]:
        """Get detailed information for reactions (batched in groups of 10)."""
        all_reactions = []
        
        for chunk in self._chunk_ids(reaction_ids, 10):
            # Format IDs for batch request
            batch_ids = '+'.join(chunk)
            response = self._request(f"/get/{batch_ids}")
            
            # Split response by /// delimiter
            reaction_blocks = response.split('///')
            for block in reaction_blocks:
                block = block.strip()
                if not block:
                    continue
                    
                parsed = self._parse_kegg_flat_file(block)
                if 'ENTRY' in parsed:
                    all_reactions.append(parsed)
        
        return all_reactions
    
    def get_compounds(self, compound_ids: List[str]) -> List[Dict[str, str]]:
        """Get detailed information for compounds (batched in groups of 10)."""
        all_compounds = []
        
        for chunk in self._chunk_ids(compound_ids, 10):
            batch_ids = '+'.join(chunk)
            response = self._request(f"/get/{batch_ids}")
            
            compound_blocks = response.split('///')
            for block in compound_blocks:
                block = block.strip()
                if not block:
                    continue
                    
                parsed = self._parse_kegg_flat_file(block)
                if 'ENTRY' in parsed:
                    all_compounds.append(parsed)
        
        return all_compounds
    
    def link_pathway_reactions(self, pathway_id: str) -> List[str]:
        """Get reaction IDs linked to a pathway."""
        response = self._request(f"/link/rn/{pathway_id}")
        rows = self._parse_tsv(response)
        
        reaction_ids = []
        for row in rows:
            # Response format: "pathway_id \t rn:reaction_id"
            if len(row) >= 2 and row[1].startswith('rn:'):
                reaction_ids.append(row[1].split(':')[1])
        
        return reaction_ids
    
    def link_pathway_compounds(self, pathway_id: str) -> List[str]:
        """Get compound IDs linked to a pathway."""
        response = self._request(f"/link/cpd/{pathway_id}")
        rows = self._parse_tsv(response)
        
        compound_ids = []
        for row in rows:
            if len(row) >= 2 and row[1].startswith('cpd:'):
                compound_ids.append(row[1].split(':')[1])
        
        return compound_ids
    
    def convert_ids(self, target_db: str, source_db: str, ids: List[str]) -> Dict[str, List[str]]:
        """Convert IDs between databases."""
        conversion_map = {}
        
        for chunk in self._chunk_ids(ids, 10):
            source_list = '+'.join(f"{source_db}:{id_}" for id_ in chunk)
            response = self._request(f"/conv/{target_db}/{source_list}")
            rows = self._parse_tsv(response)
            
            for row in rows:
                if len(row) >= 2:
                    source_id = row[0].split(':')[-1]
                    target_id = row[1].split(':')[-1]
                    
                    if source_id not in conversion_map:
                        conversion_map[source_id] = []
                    conversion_map[source_id].append(target_id)
        
        return conversion_map
    
    # Convenience methods for Solanum tuberosum
    
    def get_potato_pathways(self) -> List[Dict[str, str]]:
        """Get all pathways for Solanum tuberosum."""
        return self.list_pathways('sot')
    
    def get_potato_pathway_details(self, max_pathways: Optional[int] = None) -> List[Dict]:
        """Get detailed pathway information for potato."""
        pathways = self.get_potato_pathways()
        
        if max_pathways:
            pathways = pathways[:max_pathways]
        
        detailed_pathways = []
        for pathway in pathways:
            pathway_id = pathway['full_id']
            
            # Get reactions and compounds for this pathway
            reactions = self.link_pathway_reactions(pathway_id)
            compounds = self.link_pathway_compounds(pathway_id)
            
            detailed_pathways.append({
                'pathway': pathway,
                'reaction_ids': reactions,
                'compound_ids': compounds,
                'reaction_count': len(reactions),
                'compound_count': len(compounds)
            })
        
        return detailed_pathways
    
    def get_statistics(self) -> Dict[str, int]:
        """Get client usage statistics."""
        return {
            'total_requests': self.request_count,
            'cache_hits': getattr(self.session.cache, 'hits', 0) if hasattr(self.session, 'cache') else 0,
            'cache_misses': getattr(self.session.cache, 'misses', 0) if hasattr(self.session, 'cache') else 0
        }


def main():
    """Example usage of KEGGClient."""
    print("🔬 KEGG API Client Example")
    print("=" * 40)
    
    # Initialize client with caching
    client = KEGGClient(cache_dir="data/kegg_cache")
    
    # Get potato pathways
    print("Getting Solanum tuberosum pathways...")
    pathways = client.get_potato_pathways()
    print(f"Found {len(pathways)} pathways")
    
    # Show first 5 pathways
    for i, pathway in enumerate(pathways[:5], 1):
        print(f"  {i}. {pathway['id']} - {pathway['name']}")
    
    # Get detailed info for first pathway
    if pathways:
        first_pathway = pathways[0]
        print(f"\nGetting details for: {first_pathway['name']}")
        
        reactions = client.link_pathway_reactions(first_pathway['full_id'])
        compounds = client.link_pathway_compounds(first_pathway['full_id'])
        
        print(f"  Reactions: {len(reactions)}")
        print(f"  Compounds: {len(compounds)}")
        
        # Get reaction details (limit to first 3)
        if reactions:
            print(f"\nGetting details for first 3 reactions...")
            reaction_details = client.get_reactions(reactions[:3])
            
            for rxn in reaction_details:
                print(f"  {rxn.get('ENTRY', 'Unknown')}: {rxn.get('NAME', 'No name')}")
    
    # Show statistics
    stats = client.get_statistics()
    print(f"\nClient Statistics:")
    print(f"  Total requests: {stats['total_requests']}")
    print(f"  Cache hits: {stats['cache_hits']}")
    print(f"  Cache misses: {stats['cache_misses']}")


if __name__ == "__main__":
    main()