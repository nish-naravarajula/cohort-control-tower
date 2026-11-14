import json, sys, time
from detect import compare_releases
from impact import run_analyses, events, add_impact
from report import build

if len(sys.argv) < 3:
    sys.exit("Usage: python control_tower/run.py <release_A_dir> <release_B_dir> [out.json]")
dir_A, dir_B = sys.argv[1], sys.argv[2]
out = sys.argv[3] if len(sys.argv) > 3 else "findings.json"

t = time.time()
res_A, res_B = run_analyses(dir_A), run_analyses(dir_B)
evs = add_impact(events(compare_releases(dir_A, dir_B)), res_A, res_B)
doc, notes = build(dir_A, dir_B, evs, res_A, res_B)
with open(out, "w") as f:
    json.dump(doc, f, indent=2)

print(f"Verdict: {doc['verdict']}   ({len(doc['findings'])} findings, {time.time() - t:.1f}s)\n")
for fd in doc["findings"]:
    print(f"{fd['id']} [{fd['severity']}] {fd['summary']}")
    print(f"    affects: {fd['affects_analyses'] or 'nothing'}   not affected: {fd['not_affected_analyses']}")
    print(f"    limitation: {fd['limitation']}")
    if fd["needs_human_review"]:
        print("    NEEDS HUMAN REVIEW before release B is used")
for a in doc["analyses"]:
    print(f"{a['id']}: {a['status']}")
for n in notes:
    print(f"note: {n}")
print(f"\nWrote {out}  (research use only - not a clinical result)")
