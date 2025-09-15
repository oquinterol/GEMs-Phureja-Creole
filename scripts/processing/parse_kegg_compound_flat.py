#!/usr/bin/env python3
import re, sys, argparse, pathlib

ENTRY = re.compile(r'^ENTRY\s+CPD:\s*(C\d{5})', re.M)
FIELD = re.compile(r'^([A-Z][A-Z0-9_]+)\s{2,}(.*)$')

def parse_blocks(txt):
    return [b.strip() for b in txt.split('///') if b.strip()]

def parse_fields(block):
    fields = {}
    curr = None
    for line in block.splitlines():
        m = FIELD.match(line)
        if m:
            curr = m.group(1)
            fields.setdefault(curr, [])
            fields[curr].append(m.group(2).rstrip())
        else:
            if curr:
                fields[curr].append(line.strip())
    return {k: ' '.join(v).strip() for k,v in fields.items()}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('compound_flat')
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    txt = pathlib.Path(args.compound_flat).read_text(encoding='utf-8', errors='ignore')
    blocks = parse_blocks(txt)
    rows = [('compound_id','name','formula','exact_mass')]

    for b in blocks:
        m = ENTRY.search(b)
        if not m: continue
        cid = m.group(1)
        fields = parse_fields(b)
        name = fields.get('NAME','').split(';')[0].strip()
        formula = fields.get('FORMULA','')
        mass = fields.get('EXACT_MASS','')
        rows.append((cid, name, formula, mass))

    with open(args.out, 'w', encoding='utf-8') as f:
        for r in rows:
            f.write(','.join(f'"{c}"' if (',' in c or ' ' in c) else c for c in r) + '\n')

if __name__ == '__main__':
    main()
