#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MODPlus 결과 비교 (헤더 유연 처리 버전)
- 블록 시작: '>>' 로 시작하는 모든 줄
  * 헤더에 MGF 경로가 있거나(예: '>>E:\\...\\file.mgf 1 745.35 ...')
    없거나(예: '>> 1 745.35 ...') 모두 허용
  * 헤더 내용은 비교에서 제외, 'scan=NNNN' 또는 파일명의 '.NNNN.NNNN.k' 패턴,
    혹은 단순 '>> 17694' 숫자를 우선적으로 블록 ID로 사용
  * 위 패턴이 없으면 파일 내 등장 순번을 블록 ID로 사용 (예: '#1', '#2', ...)
- 블록 내부 라인(펩타이드 결과)만 비교 (헤더 무시)
- 블록 내부 순서는 무시(멀티셋 비교)
- 수치 컬럼 정규화: 0,1,3번째 컬럼(질량, 델타, 확률)을 반올림 (기본 4자리)
- 결과 전체를 지정 파일(기본: ./compare/cmp_result.txt)에 저장
"""

import re
from collections import Counter, defaultdict
from pathlib import Path
import argparse
import sys
from typing import Dict, Tuple, List

# 패턴들
RE_SIMPLE_NUM   = re.compile(r'^\s*>>\s*(\d+)\b')                    # >> 17694
RE_SCAN_EQ      = re.compile(r'scan=(\d+)\b', re.IGNORECASE)         # scan=1416_
RE_TRIPLE_DOT   = re.compile(r'\.(\d+)\.\1\.\d+\b')                  # ... .1416.1416.2
RE_HEADER_START = re.compile(r'^\s*>>')                              # 헤더 시작

def natural_key(s: str):
    """정렬을 위한 키: 숫자는 숫자로, 그 외는 문자열로"""
    try:
        return (0, int(s))
    except ValueError:
        return (1, s)

def parse_header_id(header_line: str, fallback_idx: int) -> str:
    """
    헤더 문자열에서 가능한 한 안정적인 블록 ID를 추출.
    우선순위: scan=NNNN > '.NNNN.NNNN.k' > '>> NNNN' > 순번 fallback
    """
    m = RE_SCAN_EQ.search(header_line)
    if m:
        return m.group(1)

    m = RE_TRIPLE_DOT.search(header_line)
    if m:
        return m.group(1)

    m = RE_SIMPLE_NUM.match(header_line)
    if m:
        return m.group(1)

    # 어떤 것도 못 찾으면 순번으로
    return f"#{fallback_idx}"

def normalize_line(line: str, decimals: int = 4) -> str:
    """
    펩타이드 결과 라인 정규화:
    - 공백/탭 구분 모두 허용
    - [0]=mass, [1]=delta, [2]=rank/int 그대로, [3]=prob 을 반올림
    """
    s = line.strip()
    if not s:
        return ""
    # 탭 우선, 없으면 공백 스플릿
    cols = re.split(r'\t+', s) if '\t' in s else s.split()

    def fmt_float(tok: str) -> str:
        try:
            return f"{float(tok):.{decimals}f}"
        except Exception:
            return tok

    if len(cols) >= 4:
        cols[0] = fmt_float(cols[0])  # mass
        cols[1] = fmt_float(cols[1])  # delta
        # cols[2] 그대로(보통 rank/int)
        cols[3] = fmt_float(cols[3])  # probability

    return "\t".join(cols)

def parse_blocks(path: Path, decimals: int = 4):
    """
    파일을 블록 단위로 파싱.
    - 헤더(>>...)를 만나면 새 블록 시작
    - 헤더 내용은 비교에서 제외
    - 이후 빈 줄/주석은 건너뛰고, 펩타이드 결과 라인을 정규화하여 카운팅
    """
    blocks: Dict[str, Counter] = {}
    order: List[str] = []
    current_id = None

    # 파일 내 ID 중복 시 충돌 방지를 위해 접미사 부여
    id_counts = defaultdict(int)
    fallback_seq = 0

    with path.open('r', encoding='utf-8', errors='replace') as f:
        for raw in f:
            if RE_HEADER_START.match(raw):
                fallback_seq += 1
                sid = parse_header_id(raw.rstrip("\n"), fallback_seq)
                id_counts[sid] += 1
                # 같은 ID가 이미 있다면 '#2' 같은 접미사로 구분
                blk_id = sid if id_counts[sid] == 1 else f"{sid}#{id_counts[sid]}"
                current_id = blk_id
                if current_id not in blocks:
                    blocks[current_id] = Counter()
                    order.append(current_id)
                continue

            if current_id is None:
                # 헤더 전에 등장하는 라인은 무시
                continue

            s = raw.strip()
            if not s or s.startswith("#"):
                continue

            norm = normalize_line(s, decimals=decimals)
            if norm:
                blocks[current_id][norm] += 1

    return blocks, order

def compare_files(ref_path: Path, cand_path: Path, decimals: int = 4):
    ref_blocks, _ = parse_blocks(ref_path, decimals=decimals)
    cand_blocks, _ = parse_blocks(cand_path, decimals=decimals)

    ref_ids = set(ref_blocks.keys())
    cand_ids = set(cand_blocks.keys())

    # 동일 ID 기준 비교 (ID가 순번/scan으로 매칭됨)
    missing_scans = sorted(ref_ids - cand_ids, key=natural_key)
    extra_scans   = sorted(cand_ids - ref_ids, key=natural_key)
    common_scans  = sorted(ref_ids & cand_ids, key=natural_key)

    per_scan_diff: Dict[str, Tuple[Counter, Counter]] = {}
    total_missing_lines = 0
    total_extra_lines = 0

    for sid in common_scans:
        r = ref_blocks[sid]
        c = cand_blocks[sid]
        miss = r - c
        extra = c - r
        if miss or extra:
            per_scan_diff[sid] = (miss, extra)
            total_missing_lines += sum(miss.values())
            total_extra_lines   += sum(extra.values())

    return {
        "missing_scans": missing_scans,
        "extra_scans": extra_scans,
        "per_scan_diff": per_scan_diff,
        "totals": {
            "missing_scans": len(missing_scans),
            "extra_scans": len(extra_scans),
            "missing_lines": total_missing_lines,
            "extra_lines": total_extra_lines
        }
    }

def write_full_report(report, ref_path: Path, cand_path: Path, out_path: Path):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open('w', encoding='utf-8', newline='\n') as w:
        t = report["totals"]
        w.write(f"[Reference] {ref_path}\n")
        w.write(f"[Candidate] {cand_path}\n\n")
        w.write("=== SUMMARY ===\n")
        w.write(f"- Missing scan blocks  : {t['missing_scans']}\n")
        w.write(f"- Extra scan blocks    : {t['extra_scans']}\n")
        w.write(f"- Missing lines (total): {t['missing_lines']}\n")
        w.write(f"- Extra lines (total)  : {t['extra_lines']}\n\n")

        if report["missing_scans"]:
            w.write("Missing scan IDs (in ref, not in cand):\n")
            w.write(", ".join(report["missing_scans"]) + "\n\n")

        if report["extra_scans"]:
            w.write("Extra scan IDs (in cand, not in ref):\n")
            w.write(", ".join(report["extra_scans"]) + "\n\n")

        if report["per_scan_diff"]:
            w.write("=== PER-SCAN DIFFS (full) ===\n")
            for sid in sorted(report["per_scan_diff"].keys(), key=natural_key):
                miss, extra = report["per_scan_diff"][sid]
                w.write(f">> {sid}\n")
                if miss:
                    w.write("  - Missing lines:\n")
                    for ln in sorted(miss):
                        cnt = miss[ln]
                        w.write(f"    x{cnt}\t{ln}\n")
                if extra:
                    w.write("  - Extra lines:\n")
                    for ln in sorted(extra):
                        cnt = extra[ln]
                        w.write(f"    x{cnt}\t{ln}\n")
                w.write("\n")
        else:
            w.write("No line-level diffs inside common scans.\n")

def main():
    ap = argparse.ArgumentParser(
        description="Compare MODPlus-like outputs ignoring within-block order (robust header)."
    )
    ap.add_argument("--ref", required=True, help="Reference file (e.g., Single.txt).")
    ap.add_argument("--cand", required=True, help="Candidate file (parallel result).")
    ap.add_argument("--decimals", type=int, default=4,
                    help="Rounding decimals for mass/delta/prob columns (default: 4).")
    ap.add_argument("--out", default=str(Path("compare") / "cmp_result.txt"),
                    help="Output report path (default: ./compare/cmp_result.txt)")
    args = ap.parse_args()

    ref_path = Path(args.ref)
    cand_path = Path(args.cand)
    out_path = Path(args.out)

    if not ref_path.exists():
        print(f"Error: reference file not found: {ref_path}", file=sys.stderr)
        sys.exit(2)
    if not cand_path.exists():
        print(f"Error: candidate file not found: {cand_path}", file=sys.stderr)
        sys.exit(2)

    report = compare_files(ref_path, cand_path, decimals=args.decimals)
    write_full_report(report, ref_path, cand_path, out_path)

    print(f"Full report saved to: {out_path}")

if __name__ == "__main__":
    main()
