# AgriPAM-AI

Open software that takes a tested microbial bank to a designed community, an editing target and a screened construct.

- **Bank → community.** Score bacteria and fungi on agronomic traits and pairwise compatibility, behind a biosafety gate; find the functions the community lacks.
- **Community → chassis.** Rank each member as an editing chassis by fit, dispensability and editing precedent.
- **Chassis → target.** Infer a native PAM from the isolate's own spacers (all 64 three-base motifs are tested; no default PAM) and scan introduced editors (SpCas9 NGG, Cas12a TTTV, dCas9).
- **Target → construct.** Assemble BioBrick parts into a screened reporter or CRISPRi construct and follow the experiment in a tracker.

All outputs are computational hypotheses, not measured editing efficiency.

## For the judge: try it online

Live application: **https://agripam-ai.streamlit.app** (no installation or sign-in needed; if the page says the app is asleep, press the wake-up button and wait about a minute).

## For the judge: run it locally in three minutes

The application has seven windows. Start with **SynCom candidate bank → Analyze your own bank**, where you can download the filled example workbook, upload it and press *Run SynCom analysis*.

**Option 1: Python 3.9 or later**

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
streamlit run app/streamlit_app.py
```

Then open http://localhost:8501.

**Option 2: Docker**

```bash
docker build -t agripam-ai .
docker run -p 8501:8501 agripam-ai
```

**Run the tests** (about one second):

```bash
python -m unittest discover -s tests
```

## What works out of the box, and what needs extra tools

- Works with the packages in `requirements.txt`: the bank analysis, community design and chassis ranking, the parts and constructs tools, the knowledgebase, the reference panels and the external-validation tools.
- The full genome workflow calls specialist tools (CRISPRCasTyper, DefenseFinder, FastANI, Parsnp, IQ-TREE). If one is not installed, its module reports **"not analysed"**; this is intended and is never shown as "not present". `environment-full.yml` installs them.
- The iGEM Registry fetch needs internet access.

## Documentation

- `docs/WINDOW_WORKFLOW.md`: the software window by window.
- `docs/HOW_IT_WORKS.md`: step by step, and `docs/SIDEBAR_SUMMARY.md`: the in-app summary.
- `docs/EXTERNAL_VALIDATION.md`: the validation benchmarks and what they do not show.
- `docs/README_engine_details.md`: reproducibility rules and engine details.

## Reproducibility

Random procedures use the fixed seed 42. Outputs carry versions and SHA-256 checksums. Isolates are identified by species and strain code.

## Licence and contact

MIT licence (see `LICENSE`). Contact: rhizoforgeai@gmail.com
