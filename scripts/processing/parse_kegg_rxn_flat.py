#!/usr/bin/env python3
import re, sys, argparse, pathlib

FIELD_RE = re.compile(r'^([A-Z][A-Z0-9_]+)\s{2,}(.*)$')
ENTRY_RE = re.compile(r'^ENTRY\s+RN:\s*(R\d{5})', re.M)

def parse_blocks(txt:str):
    return [b.strip() for b in txt.split('///') if b.strip()]

def parse_fields(block:str):
    fields = {}
    curr = None
    for line in block.splitlines():
        m = FIELD_RE.match(line)
        if m:
            curr = m.group(1)
            fields.setdefault(curr, [])
            fields[curr].append(m.group(2).rstrip())
        else:
            if curr:
                fields[curr].append(line.strip())
    return {k: ' '.join(v).strip() for k,v in fields.items()}

def parse_equation(eq:str):
    # Detecta reversibilidad
    reversible = '<=>' in eq
    arrow = '<=>' if reversible else '='
    left, right = [s.strip() for s in eq.split(arrow, 1)]
    def parse_side(s):
        items = []
        # tokens separados por '+'
        for t in s.split('+'):
            t = t.strip()
            m = re.match(r'^(\d+(?:\.\d+)?)\s+(C\d{5})$', t)
            if m:
                coeff, cid = float(m.group(1)), m.group(2)
            else:
                m2 = re.match(r'^(C\d{5})$', t)
                if not m2:
                    # ignorar nombres sin CIDs
                    continue
                coeff, cid = 1.0, m2.group(1)
            items.append((cid, coeff))
        return items
    return reversible, parse_side(left), parse_side(right)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('rxn_flat')
    ap.add_argument('--rxn-equation', required=True)
    ap.add_argument('--rxn-ec', required=True)
    ap.add_argument('--rxn-compound-edges', required=True)
    ap.add_argument('--rclass-map', required=True)
    args = ap.parse_args()

    txt = pathlib.Path(args.rxn_flat).read_text(encoding='utf-8', errors='ignore')
    blocks = parse_blocks(txt)

    eq_rows = [('reaction_id','equation','reversible')]
    ec_rows = [('reaction_id','ec_list')]
    edge_rows = [('reaction_id','role','compound_id','stoich')]
    rc_rows = [('reaction_id','rclass_list')]

    for b in blocks:
        m = ENTRY_RE.search(b)
        if not m: continue
        rid = m.group(1)
        fields = parse_fields(b)
        # EQUATION
        eq = fields.get('EQUATION','')
        if eq:
            try:
                reversible, left, right = parse_equation(eq)
                eq_rows.append((rid, eq, '1' if reversible else '0'))
                for cid, coeff in left:
                    edge_rows.append((rid, 'substrate', cid, str(coeff)))
                for cid, coeff in right:
                    edge_rows.append((rid, 'product', cid, str(coeff)))
            except Exception:
                # si falla el parseo, al menos registra la ecuación cruda
                eq_rows.append((rid, eq, ''))
        # ENZYME -> ECs separados por espacios o comas
        ecs = fields.get('ENZYME','').replace(',', ' ').split()
        if ecs:
            ec_rows.append((rid, ';'.join(ecs)))
        # RCLASS (puede aparecer en varias líneas)
        rclass = fields.get('RCLASS','').replace(' ', '')
        if rclass:
            # típicamente: RCxxxx RCyyyy ...
            rc = re.findall(r'RC\d{5}', rclass)
            if rc:
                rc_rows.append((rid, ';'.join(sorted(set(rc)))))

    # escribe CSVs
    for path, rows in [
        (args.rxn_equation, eq_rows),
        (args.rxn_ec, ec_rows),
        (args.rxn_compound_edges, edge_rows),
        (args.rclass_map, rc_rows),
    ]:
        with open(path, 'w', encoding='utf-8') as f:
            for r in rows:
                f.write(','.join(f'"{c}"' if (',' in c or ' ' in c) else c for c in r) + '\n')

if __name__ == '__main__':
    main()
