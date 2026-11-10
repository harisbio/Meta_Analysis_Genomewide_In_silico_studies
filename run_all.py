"""Regenerate every figure and the reported values.

python run_all.py              main article figures (2-11) and values
python run_all.py --global     also the Supplementary File figures (needs the global abstract dataset)

Each script runs in its own process, so one failure does not stop the others. A summary is printed at the end.
"""
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
MAIN = ["fig2_country", "fig3_trends_network", "fig4_team_size", "fig5_radial", "fig6_organisms", "fig7_families",
        "fig8_genes", "fig9_pipeline", "fig10_phylo_dup", "fig11_validation"]


def run(script):
    print(f"--- {script}", flush=True)
    r = subprocess.run([sys.executable, "-X", "utf8", str(ROOT / script)], cwd=ROOT)
    return script, r.returncode


def main():
    jobs = [f"figures/main/{m}.py" for m in MAIN] + ["analysis/manuscript_values.py"]
    if "--global" in sys.argv:
        jobs += ["figures/supplementary/global_analysis.py", "figures/supplementary/extra_figures.py"]
    results = [run(j) for j in jobs]
    failed = [s for s, rc in results if rc]
    print(f"\n{len(results) - len(failed)} of {len(results)} scripts finished" + (f"; failed: {', '.join(failed)}" if failed else ""))
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
