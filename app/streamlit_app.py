"""Judge-facing AgriPAM-AI demonstration application."""

from __future__ import annotations

import csv
import io
import json
import platform
import re
import tempfile
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import streamlit as st

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from agripam.core import analyze_ttc_genome, parse_fasta, read_pam_scores, scan_ttc_targets, score_pam  # noqa: E402
from agripam.editor_targeting import EDITOR_PRESETS, scan_editor_targets  # noqa: E402
from agripam.knowledgebase import delete_records, merge_records, with_record_ids  # noqa: E402
from agripam.ncbi import download_genome_fasta, get_assembly_identity, validate_assembly_accession  # noqa: E402
from agripam.mobile_references import build_tiered_mobile_reference  # noqa: E402
from agripam.workflow import assemble_reads, read_gff, run_genome_workflow, tool_inventory, zip_results  # noqa: E402
from agripam.validation import ValidationConfig, evaluate_ablations, evaluate_predictions, permuted_outcome_control  # noqa: E402
from agripam.visualization import build_editor_target_explorer, build_introduced_editor_funnel, build_strategy_map, build_target_explorer  # noqa: E402


st.set_page_config(page_title="AgriPAM-AI", page_icon="🧬", layout="wide")
st.title("AgriPAM-AI")
st.caption("From a tested microbial bank to a designed community and its first edit")
st.markdown(
    "<div style='display:flex; flex-wrap:wrap; gap:0.5rem; margin:0.2rem 0 0.6rem 0;'>"
    + "".join(
        f"<span style='background:#E6F3EC; color:#0E3B36; border:1px solid #BFE3D0; border-radius:999px; "
        f"padding:0.3rem 0.9rem; font-size:1.05rem; font-weight:650;'>{step}</span>"
        for step in ["1 Bank", "2 Community", "3 Chassis", "4 Target", "5 Construct", "6 Tracker"]
    )
    + "</div>",
    unsafe_allow_html=True,
)
st.markdown(
    "<p style='font-size:1.35rem; font-weight:650; color:#135EA8; margin:0.35rem 0 1rem 0;'>"
    "AI-guided, chassis-aware CRISPR genome-editing design for agricultural microbial synthetic communities"
    "</p>",
    unsafe_allow_html=True,
)
st.markdown(
    """
    <style>
      h1 { font-size: 2.45rem !important; line-height: 1.15 !important; }
      h2 { font-size: 1.85rem !important; line-height: 1.2 !important; }
      h3 { font-size: 1.35rem !important; line-height: 1.25 !important; }
      [data-testid="stCaptionContainer"] { font-size: 1.02rem !important; }
      [data-testid="stDataFrame"] { font-size: 1.02rem !important; }
      /* use most of the screen: no width cap, slim side margins */
      .block-container, [data-testid="stMainBlockContainer"] {
        padding-top: 2.2rem !important; padding-left: 1.6rem !important; padding-right: 1.6rem !important;
        max-width: 100% !important; }
      /* larger, easier-to-scan tabs */
      button[data-baseweb="tab"] { padding: 0.7rem 1.1rem !important; }
      button[data-baseweb="tab"] p { font-size: 1.1rem !important; font-weight: 650 !important; }
      /* collapsible panels: clear headers, soft cards */
      [data-testid="stExpander"] { border: 1px solid #cfe3da !important; border-radius: 14px !important;
        background: #ffffff; margin-bottom: 0.9rem; }
      [data-testid="stExpander"] summary { padding: 0.85rem 1.1rem !important; }
      [data-testid="stExpander"] summary p { font-size: 1.18rem !important; font-weight: 650 !important; color: #16324A; }
      [data-testid="stExpander"] summary:hover { background: #F3F8F6; border-radius: 14px; }
      [data-testid="stExpanderDetails"] { padding: 0.4rem 1.1rem 1.1rem 1.1rem !important; }
      /* readable metrics, buttons and sidebar */
      [data-testid="stMetricValue"] { font-size: 2.1rem !important; }
      [data-testid="stMetricLabel"] p { font-size: 1.02rem !important; }
      .stButton button, .stDownloadButton button { font-size: 1.05rem !important; padding: 0.55rem 1.1rem !important; }
      [data-testid="stSidebar"] { font-size: 1.02rem !important; }
      [data-testid="stSidebar"][aria-expanded="true"] { min-width: 320px !important; }
      [data-testid="stSidebar"] [data-testid="stExpander"] summary { padding: 0.55rem 0.75rem !important; }
      [data-testid="stSidebar"] [data-testid="stExpander"] summary p { font-size: 1.02rem !important; }
      [data-testid="stSidebar"] [data-testid="stExpanderDetails"] { padding: 0.2rem 0.75rem 0.8rem 0.75rem !important; }
      [data-testid="stAlert"] p { font-size: 1.04rem !important; }
    </style>
    """,
    unsafe_allow_html=True,
)
st.info("AgriPAM-AI links the choice of community, the member to edit, the editing target and the construct in one auditable workflow.")

with st.sidebar:
    text_size = st.radio("Text size", ["Normal", "Large", "Extra large"], index=1, horizontal=True, key="text_size",
                         help="Makes all text larger or smaller. Tables follow your browser zoom (Ctrl/Cmd and +).")
    _scale = {"Normal": 1.0, "Large": 1.11, "Extra large": 1.25}[text_size]
    st.markdown(
        f"""
        <style>
          html {{ font-size: {18 * _scale:.1f}px !important; }}
          [data-testid="stCaptionContainer"], [data-testid="stCaptionContainer"] p {{ color: #3F5566 !important; }}
          [data-testid="stSidebar"] [data-testid="stCaptionContainer"] p {{ color: #3F5566 !important; }}
        </style>
        """,
        unsafe_allow_html=True,
    )
    with st.expander("How this software works", expanded=False):
        _summary_file = ROOT / "docs" / "SIDEBAR_SUMMARY.md"
        st.markdown(_summary_file.read_text() if _summary_file.exists() else "The summary file is missing (docs/SIDEBAR_SUMMARY.md).")
        _guide_file = ROOT / "docs" / "WINDOW_WORKFLOW.md"
        if _guide_file.exists():
            st.download_button("Download the window-by-window guide", _guide_file.read_bytes(), _guide_file.name,
                               mime="text/markdown", key="dl_window_guide", width="stretch")
    with st.expander("Analysis contract (genome workflow)", expanded=False):
        st.markdown("**Inputs**\n\nFASTA or versioned NCBI accession, target length and mobile-context assumption")
        st.markdown("**Outputs**\n\nGenome annotation, PAM-to-gene consequences, defence and mobility screens, candidate targets and a machine-readable manifest")
    st.caption("Predictions remain hypotheses until wet-lab validation.")

score_path = ROOT / "machine_learning" / "type_ic_all_64_pam_scores.tsv"
context_path = ROOT / "results" / "type_ic_context_level_evidence.tsv"
summary_path = ROOT / "figures" / "prism_type_ic_summary_data.tsv"
species_path = ROOT / "figures" / "prism_three_species_crispr_comparison.tsv"
model_path = ROOT / "machine_learning" / "model_cross_validation_summary.tsv"
defence_summary_path = ROOT / "results" / "defense_atlas_species_summary.tsv"
defence_prev_path = ROOT / "results" / "defense_atlas_system_prevalence.tsv"
chassis_path = ROOT / "results" / "defense_atlas_chassis_prioritization.tsv"
agricultural_kb_path = ROOT / "data" / "agricultural_editing_knowledgebase.tsv"
native_cohort_path = ROOT / "data" / "agricultural_isolate_native_system_cohort.tsv"
native_validation_path = ROOT / "data" / "native_isolate_validation_status.tsv"
external_validation_path = ROOT / "data" / "external_validation_datasets.tsv"
gse196911_outcomes_path = ROOT / "data" / "external_validation" / "processed" / "YU2024_GSE196911_PURINE" / "heldout_outcomes.tsv"
gse196911_manifest_path = ROOT / "data" / "external_validation" / "processed" / "YU2024_GSE196911_PURINE" / "provenance_manifest.json"
gse196911_mapping_manifest_path = ROOT / "data" / "external_validation" / "processed" / "YU2024_GSE196911_PURINE" / "guide_sequences_and_outcomes.manifest.json"
published_validation_summary_path = ROOT / "data" / "external_validation" / "published_validation_summary.tsv"
gse196911_validation_metrics_path = ROOT / "data" / "external_validation" / "processed" / "YU2024_GSE196911_PURINE" / "validation_bundle" / "validation_metrics.json"
gse196911_validation_table_path = ROOT / "data" / "external_validation" / "processed" / "YU2024_GSE196911_PURINE" / "validation_bundle" / "frozen_predictions_and_outcomes.tsv"
hawkins_metrics_path = ROOT / "data" / "external_validation" / "processed" / "HAWKINS2020_GFP_TRANSFER" / "validation_metrics.json"
hawkins_table_path = ROOT / "data" / "external_validation" / "processed" / "HAWKINS2020_GFP_TRANSFER" / "frozen_scores_and_heldout_outcomes.tsv"
crisprhal_metrics_path = ROOT / "data" / "external_validation" / "processed" / "CRISPRHAL_TEV_RELEASED_HOLDOUT" / "validation_metrics.json"
crisprhal_table_path = ROOT / "data" / "external_validation" / "processed" / "CRISPRHAL_TEV_RELEASED_HOLDOUT" / "frozen_predictions_and_outcomes.tsv"
challenge_manifest_path = ROOT / "results" / "unseen_genome_challenge" / "challenge_manifest.json"
challenge_candidates_path = ROOT / "results" / "unseen_genome_challenge" / "frozen_reference_ttc_candidates.tsv"
public_native_survey_path = ROOT / "results" / "public_native_system_survey" / "public_native_system_survey.tsv"
all_ncbi_inventory_path = ROOT / "results" / "public_native_system_survey" / "all_ncbi_accession_inventory.tsv"
all_ncbi_inventory_manifest_path = ROOT / "results" / "public_native_system_survey" / "all_ncbi_accession_inventory_manifest.json"
all_accession_results_path = ROOT / "results" / "public_native_system_survey" / "all_accession_crispr_cas_results.tsv"
fungal_reference_screen_path = ROOT / "results" / "public_native_system_survey" / "fungal_reference_screen.tsv"

score_rows = read_pam_scores(score_path)
contexts = pd.read_csv(context_path, sep="\t")
summary = pd.read_csv(summary_path, sep="\t")
species = pd.read_csv(species_path, sep="\t")
models = pd.read_csv(model_path, sep="\t")
defence_summary = pd.read_csv(defence_summary_path, sep="\t") if defence_summary_path.exists() else pd.DataFrame()
defence_prev = pd.read_csv(defence_prev_path, sep="\t") if defence_prev_path.exists() else pd.DataFrame()
chassis = pd.read_csv(chassis_path, sep="\t") if chassis_path.exists() else pd.DataFrame()
agricultural_kb = pd.read_csv(agricultural_kb_path, sep="\t", keep_default_na=False) if agricultural_kb_path.exists() else pd.DataFrame()
native_cohort = pd.read_csv(native_cohort_path, sep="\t", keep_default_na=False) if native_cohort_path.exists() else pd.DataFrame()
native_validation = pd.read_csv(native_validation_path, sep="\t", keep_default_na=False) if native_validation_path.exists() else pd.DataFrame()
external_validation = pd.read_csv(external_validation_path, sep="\t", keep_default_na=False) if external_validation_path.exists() else pd.DataFrame()
gse196911_outcomes = pd.read_csv(gse196911_outcomes_path, sep="\t") if gse196911_outcomes_path.exists() else pd.DataFrame()
gse196911_manifest = json.loads(gse196911_manifest_path.read_text()) if gse196911_manifest_path.exists() else {}
gse196911_mapping_manifest = json.loads(gse196911_mapping_manifest_path.read_text()) if gse196911_mapping_manifest_path.exists() else {}
published_validation_summary = pd.read_csv(published_validation_summary_path, sep="\t") if published_validation_summary_path.exists() else pd.DataFrame()
hawkins_metrics = json.loads(hawkins_metrics_path.read_text()) if hawkins_metrics_path.exists() else {}
crisprhal_metrics = json.loads(crisprhal_metrics_path.read_text()) if crisprhal_metrics_path.exists() else {}
crisprhal_table = pd.read_csv(crisprhal_table_path, sep="\t") if crisprhal_table_path.exists() else pd.DataFrame()
challenge_manifest = json.loads(challenge_manifest_path.read_text()) if challenge_manifest_path.exists() else {}
challenge_candidates = pd.read_csv(challenge_candidates_path, sep="\t") if challenge_candidates_path.exists() else pd.DataFrame()
hawkins_table = pd.read_csv(hawkins_table_path, sep="\t") if hawkins_table_path.exists() else pd.DataFrame()
gse196911_validation_metrics = json.loads(gse196911_validation_metrics_path.read_text()) if gse196911_validation_metrics_path.exists() else {}
gse196911_validation_table = pd.read_csv(gse196911_validation_table_path, sep="\t") if gse196911_validation_table_path.exists() else pd.DataFrame()
public_native_survey = pd.read_csv(public_native_survey_path, sep="\t", keep_default_na=False) if public_native_survey_path.exists() else pd.DataFrame()
all_ncbi_inventory = pd.read_csv(all_ncbi_inventory_path, sep="\t", keep_default_na=False) if all_ncbi_inventory_path.exists() else pd.DataFrame()
all_ncbi_inventory_manifest = json.loads(all_ncbi_inventory_manifest_path.read_text()) if all_ncbi_inventory_manifest_path.exists() else {}
expanded_result_files = ([all_accession_results_path] if all_accession_results_path.exists() else []) + sorted(
    (ROOT / "results" / "public_native_system_survey").glob("expanded_*.tsv"))
all_accession_results = (pd.concat([pd.read_csv(path, sep="\t", keep_default_na=False) for path in expanded_result_files], ignore_index=True)
                         .drop_duplicates("canonical_assembly", keep="last") if expanded_result_files else pd.DataFrame())
if not all_accession_results.empty:
    all_accession_results["pipeline_validation"] = all_accession_results.get(
        "parser_revision", pd.Series("", index=all_accession_results.index)).map(
            lambda value: "Validated parser v2" if value == "confirmed_cas_operons_tab_v2"
            else "WITHDRAWN — legacy parser read the wrong output table")
fungal_reference_screen = pd.read_csv(fungal_reference_screen_path, sep="\t", keep_default_na=False) if fungal_reference_screen_path.exists() else pd.DataFrame()
priority_result_files = sorted((ROOT / "results" / "public_native_system_survey").glob("priority_results_*.tsv"))
priority_results = (pd.concat([pd.read_csv(path, sep="\t", keep_default_na=False) for path in priority_result_files], ignore_index=True)
                    .drop_duplicates("canonical_assembly", keep="last") if priority_result_files else pd.DataFrame())
wave_result_files = sorted((ROOT / "results" / "public_native_system_survey").glob("wave*.tsv"))
wave_results = (pd.concat([pd.read_csv(path, sep="\t", keep_default_na=False) for path in wave_result_files], ignore_index=True)
                .drop_duplicates("canonical_assembly", keep="last") if wave_result_files else pd.DataFrame())
if not wave_results.empty and "related_submitted_isolates" in wave_results.columns:
    wave_results = wave_results.rename(columns={"related_submitted_isolates": "isolate"})
corrected_bank_results = pd.concat([frame for frame in (priority_results, wave_results) if not frame.empty],
                                   ignore_index=True) if not priority_results.empty or not wave_results.empty else pd.DataFrame()
if not corrected_bank_results.empty:
    corrected_bank_results = corrected_bank_results.drop_duplicates("canonical_assembly", keep="last")


@st.cache_data(show_spinner=False)
def load_mobile_reference_collection(paths: tuple[str, ...]) -> str:
    """Load a versioned local mobile-element collection once per app session."""
    return "\n".join(Path(path).read_text() for path in paths)


@st.cache_data(show_spinner="Building a host-matched mobile-element collection…")
def load_host_aware_reference(organism: str) -> tuple[str | None, dict]:
    return build_tiered_mobile_reference(ROOT, organism)

# Code below refers to tabs by their original position; DISPLAY_ORDER only changes where each one appears.
TAB_NAMES = [
    "Start here",
    "Genome evaluation",
    "SynCom candidate bank",
    "Agricultural editing knowledgebase",
    "Software & reproducibility",
    "External validation",
    "Parts & constructs",
]
DISPLAY_ORDER = [
    "Start here",
    "Agricultural editing knowledgebase",
    "SynCom candidate bank",
    "Genome evaluation",
    "Parts & constructs",
    "External validation",
    "Software & reproducibility",
]
_displayed = st.tabs(DISPLAY_ORDER)
tabs = [_displayed[DISPLAY_ORDER.index(name)] for name in TAB_NAMES]

with tabs[0]:
    st.header("From target organism to validated edit")
    st.info(
        "Begin here. This roadmap explains everything required for a defensible microbial editing project and "
        "shows which decisions belong to AgriPAM-AI and which require laboratory confirmation.",
        icon=":material/info:",
    )
    st.subheader("Interactive CRISPR learning model")
    st.caption(
        "Change the assumptions to see why a PAM alone is insufficient. This is an educational simulation; "
        "it does not design an experimental guide or predict editing efficiency."
    )
    nuclease_models = {
        "SpCas9 (introduced tool)": {"target_rule": "NGG", "action": "DNA double-strand cleavage"},
        "Cas12a/Cpf1 (introduced tool)": {"target_rule": "TTTV", "action": "Staggered DNA cleavage"},
        "CRISPRi with dCas9": {"target_rule": "NGG", "action": "DNA binding and transcriptional repression"},
        "Native strain CRISPR system": {"target_rule": "Must be inferred for the strain", "action": "Depends on detected subtype and activity"},
    }
    learning_a, learning_b, learning_c = st.columns(3)
    with learning_a:
        learning_nuclease = st.selectbox("Editing system", list(nuclease_models), key="learning_nuclease")
        learning_goal = st.selectbox(
            "Educational objective",
            ["Precise DNA change", "Gene disruption", "Gene repression without changing DNA"],
            key="learning_goal",
        )
    with learning_b:
        learning_pam = st.selectbox(
            "Target-site evidence",
            ["Compatible PAM/target rule present", "PAM/target rule unknown", "No compatible PAM/target rule"],
            key="learning_pam",
        )
        learning_guide = st.selectbox(
            "Guide specificity",
            ["Unique candidate", "Repeated or off-target-prone candidate", "Not evaluated"],
            key="learning_guide",
        )
    with learning_c:
        learning_delivery = st.toggle("Compatible delivery demonstrated", key="learning_delivery")
        learning_repair = st.toggle(
            "Compatible repair strategy available",
            key="learning_repair",
            help="A repair template or repair route is generally needed for a defined sequence change; CRISPRi acts without changing DNA.",
        )

    selected_nuclease = nuclease_models[learning_nuclease]
    pam_supported = learning_pam == "Compatible PAM/target rule present"
    pam_unknown = learning_pam == "PAM/target rule unknown"
    guide_supported = learning_guide == "Unique candidate"
    repression_route = learning_nuclease == "CRISPRi with dCas9" or learning_goal == "Gene repression without changing DNA"
    repair_supported = learning_repair or repression_route
    stage_states = [
        ("1", "Target rule", "Ready" if pam_supported else "Unknown" if pam_unknown else "Blocked"),
        ("2", "Guide", "Ready" if guide_supported else "Unknown" if learning_guide == "Not evaluated" else "Blocked"),
        ("3", "Delivery", "Ready" if learning_delivery else "Unverified"),
        ("4", "Repair/action", "Ready" if repair_supported else "Missing"),
        ("5", "Validation", "Required"),
    ]
    stage_columns = st.columns(5)
    state_icons = {"Ready": "✅", "Required": "🧪", "Unknown": "❓", "Unverified": "⚠️", "Blocked": "⛔", "Missing": "⛔"}
    for stage_column, (number, label, state) in zip(stage_columns, stage_states):
        stage_column.markdown(f"**{number}. {label}**")
        stage_column.markdown(f"{state_icons[state]} {state}")

    pam_color = "#1B8A5A" if pam_supported else "#D08A00" if pam_unknown else "#C83E4D"
    guide_color = "#1B8A5A" if guide_supported else "#D08A00" if learning_guide == "Not evaluated" else "#C83E4D"
    delivery_color = "#1B8A5A" if learning_delivery else "#D08A00"
    repair_color = "#1B8A5A" if repair_supported else "#C83E4D"
    nuclease_short = "dCas9" if "dCas9" in learning_nuclease else "Cas12a" if "Cas12a" in learning_nuclease else "Native Cas" if "Native" in learning_nuclease else "Cas9"
    display_target_rule = "PAM?" if learning_nuclease == "Native strain CRISPR system" else selected_nuclease["target_rule"]
    molecular_result = "REPRESS" if repression_route else "EDIT" if repair_supported else "UNRESOLVED"
    cut_symbol = "⊣" if repression_route else "✂"
    action_caption = "DNA binding blocks transcription" if repression_route else "nuclease cleavage at the target"
    outcome_caption = "reduced transcription" if repression_route else "designed sequence change" if repair_supported else "repair outcome unknown"
    repair_template_title = "NO DNA DONOR · CRISPRi REGULATION" if repression_route else "DONOR TEMPLATE / REPAIR"
    repair_template_detail = "expression changes without sequence replacement" if repression_route else "left homology arm  →  intended change  ←  right homology arm"
    st.components.v1.html(
        f"""
        <style>
          html, body {{ width:100%; height:100%; margin:0; padding:0; overflow:hidden; background:#F7F9FC; font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif; }}
          .figure-frame {{ width:100%; height:100vh; padding:12px; box-sizing:border-box; display:flex; align-items:center; justify-content:center; overflow:hidden; }}
          svg {{ display:block; width:100%; height:100%; max-width:1200px; max-height:470px; }}
          .label {{ fill:#26334A; font-size:16px; font-weight:600; }}
          .small {{ fill:#53627A; font-size:13px; }}
          .step {{ fill:#135EA8; font-size:14px; font-weight:700; letter-spacing:.5px; }}
        </style>
        <div class="figure-frame">
          <svg viewBox="0 0 1200 470" preserveAspectRatio="xMidYMid meet" role="img" aria-labelledby="crispr-title crispr-desc">
            <title id="crispr-title">Dynamic CRISPR editing mechanism</title>
            <desc id="crispr-desc">Classical CRISPR diagram showing delivery of a Cas guide complex, guide pairing to double-stranded DNA beside a PAM, nuclease cleavage or repression, and validation of the outcome.</desc>
            <defs>
              <marker id="rf-arrow" markerWidth="10" markerHeight="10" refX="9" refY="5" orient="auto">
                <path d="M0,0 L10,5 L0,10 Z" fill="#65738A"></path>
              </marker>
              <filter id="soft-shadow" x="-20%" y="-20%" width="140%" height="140%">
                <feDropShadow dx="0" dy="3" stdDeviation="4" flood-color="#26334A" flood-opacity="0.16"/>
              </filter>
            </defs>
            <text x="600" y="30" text-anchor="middle" font-size="22" font-weight="700" fill="#26334A">How CRISPR reaches and acts on a bacterial DNA target</text>

            <text x="105" y="70" text-anchor="middle" class="step">1 · DELIVERY</text>
            <rect x="30" y="88" width="150" height="74" rx="37" fill="#EAF3FF" stroke="{delivery_color}" stroke-width="4" filter="url(#soft-shadow)"></rect>
            <circle cx="77" cy="125" r="25" fill="#DCD3F2" stroke="#7257A8" stroke-width="3"></circle>
            <path d="M111 139 C120 107,136 107,145 139 C151 158,162 150,166 132" fill="none" stroke="{guide_color}" stroke-width="5" stroke-linecap="round"></path>
            <text x="77" y="130" text-anchor="middle" font-size="12" font-weight="700" fill="#47366B">{nuclease_short}</text>
            <text x="105" y="184" text-anchor="middle" font-size="13" font-weight="600" fill="{delivery_color}">{'delivery supported' if learning_delivery else 'delivery unverified'}</text>
            <line x1="190" y1="125" x2="255" y2="125" stroke="#65738A" stroke-width="3" marker-end="url(#rf-arrow)"></line>

            <text x="555" y="70" text-anchor="middle" class="step">2 · TARGET RECOGNITION AND ACTION</text>
            <text x="265" y="221" class="small">5′</text>
            <text x="265" y="286" class="small">3′</text>
            <text x="840" y="221" class="small">3′</text>
            <text x="840" y="286" class="small">5′</text>
            <line x1="290" y1="215" x2="825" y2="215" stroke="#315C9B" stroke-width="8" stroke-linecap="round"></line>
            <line x1="290" y1="280" x2="825" y2="280" stroke="#78A6DB" stroke-width="8" stroke-linecap="round"></line>
            <g stroke="#AFC4DF" stroke-width="2">
              <line x1="310" y1="222" x2="310" y2="273"></line><line x1="335" y1="222" x2="335" y2="273"></line>
              <line x1="360" y1="222" x2="360" y2="273"></line><line x1="385" y1="222" x2="385" y2="273"></line>
              <line x1="410" y1="222" x2="410" y2="273"></line><line x1="435" y1="222" x2="435" y2="273"></line>
              <line x1="460" y1="222" x2="460" y2="273"></line><line x1="485" y1="222" x2="485" y2="273"></line>
              <line x1="510" y1="222" x2="510" y2="273"></line><line x1="535" y1="222" x2="535" y2="273"></line>
              <line x1="560" y1="222" x2="560" y2="273"></line><line x1="585" y1="222" x2="585" y2="273"></line>
              <line x1="610" y1="222" x2="610" y2="273"></line><line x1="635" y1="222" x2="635" y2="273"></line>
            </g>
            <rect x="350" y="207" width="292" height="81" rx="8" fill="#DDF3E9" fill-opacity="0.72" stroke="{guide_color}" stroke-width="3"></rect>
            <text x="496" y="310" text-anchor="middle" font-size="14" font-weight="600" fill="{guide_color}">protospacer · {learning_guide.lower()}</text>
            <rect x="650" y="207" width="96" height="81" rx="8" fill="#FFF1C7" stroke="{pam_color}" stroke-width="4"></rect>
            <text x="698" y="239" text-anchor="middle" class="label">PAM</text>
            <text x="698" y="266" text-anchor="middle" font-size="17" font-weight="700" fill="{pam_color}">{display_target_rule}</text>

            <path d="M353 198 L353 169 C353 142,327 142,327 169 C327 190,304 190,304 166 C304 132,278 134,278 161" fill="none" stroke="{guide_color}" stroke-width="6" stroke-linecap="round"></path>
            <line x1="365" y1="197" x2="635" y2="197" stroke="{guide_color}" stroke-width="6" stroke-linecap="round"></line>
            <text x="458" y="143" text-anchor="middle" font-size="14" font-weight="700" fill="{guide_color}">guide RNA</text>

            <path d="M390 105 C450 72,574 78,626 120 C670 157,668 210,630 235 C585 264,474 256,411 224 C357 196,342 132,390 105 Z" fill="#DCD3F2" fill-opacity="0.82" stroke="#7257A8" stroke-width="5" filter="url(#soft-shadow)"></path>
            <text x="510" y="132" text-anchor="middle" font-size="23" font-weight="700" fill="#47366B">{nuclease_short}</text>
            <text x="510" y="155" text-anchor="middle" font-size="13" fill="#47366B">guide-directed recognition</text>
            <line x1="635" y1="190" x2="635" y2="296" stroke="#C83E4D" stroke-width="3" stroke-dasharray="7 5"></line>
            <text x="632" y="185" text-anchor="middle" font-size="31" fill="#C83E4D">{cut_symbol}</text>
            <text x="755" y="330" text-anchor="middle" class="small">{action_caption}</text>

            <line x1="635" y1="299" x2="635" y2="334" stroke="#65738A" stroke-width="2" stroke-dasharray="5 4" marker-end="url(#rf-arrow)"></line>
            <rect x="300" y="338" width="600" height="38" rx="9" fill="#FFF8E6" stroke="{repair_color}" stroke-width="2"></rect>
            <text x="600" y="354" text-anchor="middle" font-size="12" font-weight="700" fill="#7A5200">{repair_template_title}</text>
            <text x="600" y="370" text-anchor="middle" font-size="11" fill="#53627A">{repair_template_detail}</text>

            <line x1="855" y1="248" x2="920" y2="248" stroke="#65738A" stroke-width="3" marker-end="url(#rf-arrow)"></line>
            <text x="1035" y="70" text-anchor="middle" class="step">3 · OUTCOME</text>
            <rect x="930" y="90" width="240" height="250" rx="18" fill="#FFFFFF" stroke="{repair_color}" stroke-width="4" filter="url(#soft-shadow)"></rect>
            <text x="1050" y="126" text-anchor="middle" font-size="20" font-weight="700" fill="{repair_color}">{molecular_result}</text>
            <text x="1050" y="152" text-anchor="middle" class="small">{outcome_caption}</text>
            <line x1="970" y1="205" x2="1130" y2="205" stroke="#315C9B" stroke-width="7" stroke-linecap="round"></line>
            <line x1="970" y1="245" x2="1130" y2="245" stroke="#78A6DB" stroke-width="7" stroke-linecap="round"></line>
            <rect x="1032" y="194" width="40" height="62" rx="5" fill="{repair_color}" fill-opacity="0.25" stroke="{repair_color}" stroke-width="2"></rect>
            <text x="1050" y="290" text-anchor="middle" class="label">must be measured</text>
            <text x="1050" y="315" text-anchor="middle" class="small">not predicted as success</text>

            <line x1="600" y1="378" x2="600" y2="380" stroke="#65738A" stroke-width="3" marker-end="url(#rf-arrow)"></line>
            <rect x="345" y="386" width="510" height="58" rx="12" fill="#EAF3FF" stroke="#315C9B" stroke-width="2"></rect>
            <text x="600" y="412" text-anchor="middle" class="label">CONFIRMATION AND VALIDATION</text>
            <text x="600" y="434" text-anchor="middle" class="small">sequence → phenotype → preserved function → non-target safety</text>
          </svg>
        </div>
        """,
        width="stretch",
        # Keep the title and validation strip visible on narrow layouts.
        height=620,
    )

    if pam_supported and guide_supported and learning_delivery and repair_supported:
        st.success(
            "The simulated design chain is complete enough to proceed to controlled experimental validation. "
            "This does not predict that editing will succeed."
        )
    elif learning_pam == "No compatible PAM/target rule":
        st.error(
            "This nuclease cannot use the simulated target. Choose another target, another nuclease or a non-PAM-dependent route."
        )
    else:
        missing_learning = [label for _, label, state in stage_states[:-1] if state != "Ready"]
        st.warning("Evidence still needed for: " + ", ".join(missing_learning) + ".")

    with st.expander("Explain the current simulation", expanded=False):
        st.markdown(f"**Selected system:** {learning_nuclease}")
        st.markdown(f"**Target requirement:** {selected_nuclease['target_rule']}")
        st.markdown(f"**Molecular action:** {selected_nuclease['action']}")
        st.markdown(
            "**How to interpret it:** the guide supplies sequence specificity; the PAM or target rule permits recognition; "
            "the nuclease binds, cuts or represses; delivery determines whether the machinery reaches the cell; repair "
            "determines the resulting DNA change; and physical measurements establish whether the intended result occurred."
        )
        st.info(
            "For an introduced SpCas9 or Cas12a tool, the displayed PAM is a property of that tool. For a native "
            "CRISPR system, its strain-specific targeting rule must be inferred and experimentally validated."
        )

    with st.expander("How to read PAM notation (IUPAC legend)", expanded=False):
        st.markdown(
            """
            **Common IUPAC nucleotide symbols**

            - **N** = A, C, G or T
            - **R** = A or G
            - **V** = A, C or G

            For a model that evaluates every exact three-base sequence, there are **4³ = 64 possible triplets**.

            These 64 triplets are not universally valid PAMs. Each nuclease recognizes only a subset. For example,
            **NGG** represents AGG, CGG, GGG and TGG, whereas **TTTV** represents TTTA, TTTC and TTTG.
            """
        )
        st.info(
            f"The figure currently shows **{display_target_rule}** because **{learning_nuclease}** is selected. "
            "For an introduced editing tool, the PAM belongs to that nuclease. It is not automatically the PAM of "
            "the bacterium's native CRISPR system."
        )

    st.subheader("Editing route selector")
    st.write(
        "AgriPAM-AI evaluates several routes separately. The preferred route must be able to create the intended "
        "effect while preserving essential agricultural, developmental and community-complementarity functions. "
        "A PAM elsewhere in the genome is irrelevant unless editing at that site can produce the intended effect."
    )
    start_route_table = pd.DataFrame([
        {"Route": "Native Type I-C candidate", "PAM or target rule": "Usually a 5′ PAM; TTC is the current P. polymyxa hypothesis", "Action": "Cascade recognition followed by Cas3 interference", "When it can be preferred": "The isolate has an intact, active and deliverable strain-specific locus"},
        {"Route": "SpCas9", "PAM or target rule": "NGG adjacent to the intended DNA target", "Action": "DNA cleavage", "When it can be preferred": "A unique target and compatible repair route are available"},
        {"Route": "dCas9 / CRISPRi", "PAM or target rule": "NGG in a useful promoter or transcriptional region", "Action": "Repression without DNA cleavage", "When it can be preferred": "Reversible regulation is safer or more appropriate than a permanent edit"},
        {"Route": "Cas12a / Cpf1", "PAM or target rule": "Typically TTTV; variant-specific", "Action": "Staggered DNA cleavage", "When it can be preferred": "A suitable target and supported delivery route are available"},
        {"Route": "Alternative or insertion route", "PAM or target rule": "Variant-specific PAM, regulatory site or validated neutral locus", "Action": "Alternative targeting, regulation or function insertion", "When it can be preferred": "The intended locus lacks a useful PAM or the objective is to add a missing benefit"},
    ])
    st.dataframe(start_route_table, hide_index=True, width="stretch")
    st.caption(
        "Ranking order: preserve required functions; match the biological objective; require a correctly positioned "
        "target; evaluate guide uniqueness and off-target risk; then assess delivery, defence, repair and controls. "
        "PAM compatibility indicates computational targetability, not successful editing."
    )

    st.divider()
    with st.expander("Define the intended edit and success measurement", expanded=False):
        st.caption(
            "These choices are carried into Genome evaluation, route prioritization, construct design and validation."
        )
        objective_col, readout_col = st.columns(2)
        with objective_col:
            editing_objective = st.selectbox(
                "Intended modification",
                ["Not yet defined", "Gene deletion", "Gene insertion", "Sequence replacement", "Point mutation", "CRISPR interference (repression)", "Gene activation"],
                key="editing_objective",
                help="Defines which editing routes and construct types are prioritized downstream.",
            )
        with readout_col:
            phenotype_readout_choice = st.selectbox(
                "Primary measurable readout",
                [
                    "Not yet selected", "Fluorescence reporter", "Antifungal inhibition",
                    "Antibacterial inhibition", "Phosphate solubilization", "Nitrogen transformation",
                    "Siderophore production", "Lipopeptide or metabolite production",
                    "Hydrolytic-enzyme activity", "Root colonization", "Growth or fitness",
                    "Transformation or editing efficiency", "Other — specify",
                ],
                key="phenotype_readout_choice",
                help="Defines the primary validation outcome carried into later design steps.",
            )
            if phenotype_readout_choice == "Other — specify":
                phenotype_readout = st.text_input(
                    "Describe the other measurable readout",
                    key="phenotype_readout_other",
                    placeholder="Enter a measurable phenotype and, if possible, its unit",
                ).strip()
            elif phenotype_readout_choice == "Not yet selected":
                phenotype_readout = ""
            else:
                phenotype_readout = phenotype_readout_choice
        if editing_objective != "Not yet defined" or phenotype_readout:
            st.info(
                f"Active project definition: **{editing_objective}**; success measured primarily by "
                f"**{phenotype_readout or 'a readout still to be selected'}**."
            )
    st.session_state["project_brief"] = {
        "modification": editing_objective,
        "readout": phenotype_readout,
    }
    genome_ready = bool(st.session_state.get("full_workflow"))
    objective_ready = editing_objective != "Not yet defined" and bool(phenotype_readout.strip())
    roadmap_steps = [
        ("1. Define the objective", "Specify deletion, insertion, replacement, repression or activation and its quantitative phenotype.", "Records the objective and connects it to candidate ranking.", "Specified" if objective_ready else "User input required"),
        ("2. Verify the strain and genome", "Confirm strain identity, assembly quality, annotation and the exact target sequence.", "Processes FASTA, reads or an NCBI assembly into an auditable genome package.", "Genome analysed" if genome_ready else "Run Genome evaluation"),
        ("3. Choose a compatible editing system", "Assess a native or introduced nuclease, recombinase, integrase or alternative editing route.", "Screens system components and compares literature-supported alternatives.", "Evaluate after genome analysis"),
        ("4. Identify a compatible target", "Establish the PAM or targeting rule, guide uniqueness, position and gene consequence.", "Ranks supported sites and maps them to annotated genes.", "Evaluate after PAM evidence"),
        ("5. Define the repair strategy", "Choose an appropriate repair template or a regulatory strategy such as CRISPRi.", "Reports repair-machinery indicators, compatibility and limitations.", "Experimental design required"),
        ("6. Establish delivery", "Select a strain-compatible delivery, expression and selection architecture.", "Flags restriction, defence and plasmid-compatibility barriers.", "Strain-specific SOP required"),
        ("7. Evaluate defence barriers", "Assess restriction–modification and other systems acting against incoming DNA.", "Runs available defence and restriction-system analyses.", "Evaluate after genome analysis"),
        ("8. Define experimental controls", "Include parental, no-guide/delivery-only, positive and lower-ranked comparator controls.", "Generates a control-oriented results template.", "Approve before experiment"),
        ("9. Confirm the genotype", "Verify the target locus, edit junctions, sequence and construct loss or persistence.", "Links confirmation records to the selected candidate.", "Wet-lab confirmation required"),
        ("10. Validate function and safety", "Measure fitness, intended phenotype, preservation of beneficial functions and non-target effects.", "Links measurements to predictions for evidence updating.", "Wet-lab validation required"),
    ]
    completed_foundations = int(objective_ready) + int(genome_ready)
    st.progress(completed_foundations / 2, text=f"Computational preparation foundations completed: {completed_foundations}/2")
    with st.expander("View the complete project roadmap", expanded=False):
        st.dataframe(
            pd.DataFrame(roadmap_steps, columns=["Stage", "What is required", "AgriPAM-AI contribution", "Status"]),
            hide_index=True,
            width="stretch",
        )
    st.caption(
        "AgriPAM-AI supports evidence collection and design prioritization. Delivery, successful editing, "
        "phenotype and safety require an approved strain-specific procedure and physical validation."
    )

with tabs[1]:
    st.header("Genome-to-design workspace")
    active_design = st.session_state.get("design_handoff")
    project_brief = st.session_state.get("project_brief", {})
    if project_brief.get("modification") != "Not yet defined" or project_brief.get("readout"):
        st.info(
            f"**Project definition carried forward:** {project_brief.get('modification', 'not defined')}; "
            f"primary validation readout: **{project_brief.get('readout') or 'not selected'}**."
        )
    if active_design:
        st.success(
            f"**Active SynCom design:** community {', '.join(active_design.get('community', []))}; "
            f"selected chassis **{active_design.get('chassis', 'not selected')}**; "
            f"objective **{active_design.get('objective', 'not defined')}**."
        )
        st.caption(
            "Provide the exact selected-isolate genome below. AgriPAM-AI will compare native evidence, "
            "SpCas9/dCas9 NGG sites and Cas12a TTTV sites, map candidates to genes and neutral regions, "
            "and preserve the functions marked as protected. A related-species reference cannot confirm isolate targets."
        )
        if st.button("Clear SynCom design handoff", key="clear_design_handoff"):
            del st.session_state["design_handoff"]
            st.rerun()
    st.info(
        "This workspace analyzes the genome supplied below. Its outputs are separate "
        "from the fixed comparative P. polymyxa reference study shown in the other tabs.",
        icon=":material/info:",
    )
    st.write(
        "Start from an assembled genome, raw sequencing reads, or a versioned NCBI "
        "assembly accession. The workflow creates an auditable annotation package and "
        "separates database-supported calls from transparent preliminary screens."
    )
    inventory = pd.DataFrame(tool_inventory(ROOT))
    with st.expander("Pipeline capability check", icon=":material/fact_check:"):
        st.dataframe(inventory, hide_index=True, width="stretch")
        st.caption("Unavailable modules are reported as unavailable—not as absent biological systems.")

    source_mode = st.segmented_control(
        "Input source",
        ["Assembled FASTA", "Raw sequencing reads", "NCBI accession"],
        default="Assembled FASTA",
        key="genome_input_source",
        help=(
            "FASTA: a completed genome assembly. Raw reads: FASTQ files that must first be assembled. "
            "NCBI: a public, versioned GCF/GCA genome-assembly accession. Strain-specific PAM inference "
            "also requires the optional mobile-sequence reference FASTA."
        ),
    )
    source_guidance = {
        "Assembled FASTA": ("Upload a completed genome assembly in FA, FASTA or FNA format. "
                            "Use this for your native strain after genome assembly."),
        "Raw sequencing reads": ("Upload paired Illumina reads, one long-read file, or paired short reads "
                                 "plus long reads for hybrid assembly."),
        "NCBI accession": ("Enter a public, versioned assembly accession—for example GCF_000597985.1 or "
                           "GCA_000597985.1. BioSample, BioProject, gene and protein accessions are not accepted."),
    }
    # During a hot reload Streamlit may briefly return None before the
    # segmented control has initialized. Render a safe default instead of
    # raising KeyError: None and breaking the whole page.
    if source_mode is None:
        source_mode = "Assembled FASTA"
    st.info(source_guidance.get(source_mode, source_guidance["Assembled FASTA"]), icon=":material/info:")
    if source_mode == "NCBI accession":
        st.warning(
            "Important: an NCBI assembly accession downloads only the uploaded organism’s "
            "assembled chromosomes/contigs. NCBI does not automatically supply the separate "
            "phage/plasmid/mobile-element sequences used for strain-specific PAM discovery. "
            "AgriPAM-AI reads the organism and strain from NCBI metadata, then uses an "
            "installed mobile collection only when its host scope matches. Otherwise it requests "
            "a matching collection and reports PAM discovery as not run.",
            icon=":material/database:",
        )

    host_aware_registry_paths = [
        ROOT / "metadata" / "caudoviricetes_refseq_sequence_manifest.tsv",
        ROOT / "data" / "phages" / "processed" / "caudoviricetes_refseq_exact_deduplicated.fna",
        ROOT / "data" / "plasmids" / "ncbi_raw_chunks",
    ]
    mobile_reference_mode = st.selectbox(
        "PAM-inference reference collection",
        ["Automatic host-aware selection", "Upload a reference FASTA", "No external collection"],
        index=0,
        key="mobile_reference_mode",
        help=("Automatic mode first identifies the submitted target organism, then selects locally installed "
              "phage and plasmid records associated with that species, its defined species complex, or its genus. "
              "It does not substitute Paenibacillus records for another target organism."),
    )
    st.caption(
        "**What this means:** The reference collection is separate from the uploaded bacterial genome. "
        "It contains comparative phage, plasmid and mobile-element sequences used to infer PAM preferences. "
        "**Automatic host-aware selection** chooses records associated with the submitted organism, species "
        "complex or genus; **Upload a reference FASTA** uses your own collection; **No external collection** "
        "analyzes only the bacterial genome and does not run strain-specific PAM inference. All inferred PAMs "
        "are hypotheses that require experimental validation."
    )
    if mobile_reference_mode == "Automatic host-aware selection":
        if all(path.exists() for path in host_aware_registry_paths):
            st.success("Host-aware reference registry ready.", icon=":material/library_books:")
        else:
            st.warning("The host-aware reference registry is incomplete; upload an appropriate reference FASTA instead.")

    with st.form("complete_genome_workflow"):
        organism_name = st.text_input(
            "Target organism name",
            placeholder="Paenibacillus polymyxa",
            disabled=source_mode == "NCBI accession",
            help="Required for uploaded FASTA or raw reads. For an NCBI assembly, the organism is retrieved automatically from NCBI metadata.",
        )
        strain_name = st.text_input(
            "Strain or isolate name",
            placeholder="B34 (optional)",
            disabled=source_mode == "NCBI accession",
            help="Use your laboratory isolate identifier. If omitted, the uploaded filename is used as a provisional label.",
        )
        genome_upload = st.file_uploader(
            "Genome assembly (FASTA)", type=["fa", "fasta", "fna"],
            disabled=source_mode != "Assembled FASTA",
        )
        read_technology = st.selectbox(
            "Sequencing technology",
            ["Illumina paired-end", "Long reads", "Hybrid"],
            disabled=source_mode != "Raw sequencing reads",
        )
        read_uploads = st.file_uploader(
            "Raw reads (FASTQ or FASTQ.GZ)",
            type=["fastq", "fq", "gz"], accept_multiple_files=True,
            disabled=source_mode != "Raw sequencing reads",
            help="Paired-end: R1 and R2. Long reads: one file. Hybrid: R1, R2, then long reads.",
        )
        ncbi_accession = st.text_input(
            "Versioned NCBI assembly accession",
            placeholder="GCF_000597985.1",
            disabled=source_mode != "NCBI accession",
        )
        specialist = st.checkbox(
            "Run installed specialist analyses (slower)", value=True,
            help="Runs CRISPRCasTyper, Bakta/Prokka, DefenseFinder and MOB-suite when available.",
        )
        mobile_reference_upload = st.file_uploader(
            "Separate mobile-sequence reference collection (FASTA, optional)",
            type=["fa", "fasta", "fna"],
            help=("Versioned phage, plasmid or mobile-element sequences used to find "
                  "protospacer matches. This file is optional for general genome annotation, but required "
                  "for strain-specific PAM inference. Without it, no strain-specific PAM is claimed."),
            disabled=mobile_reference_mode != "Upload a reference FASTA",
        )
        submitted = st.form_submit_button("Run genome workflow", type="primary", icon=":material/play_arrow:")

    if submitted:
        try:
            with st.status("Running genome-to-design workflow", expanded=True) as status:
                if source_mode == "Assembled FASTA":
                    if genome_upload is None:
                        raise ValueError("Upload an assembled FASTA file.")
                    fasta_text = genome_upload.getvalue().decode("utf-8")
                    supplied_organism = organism_name.strip()
                    filename_accession = re.search(r"GC[AF]_\d{9}\.\d+", genome_upload.name, re.I)
                    if supplied_organism:
                        identity = {"organism": supplied_organism,
                                    "strain": strain_name.strip() or Path(genome_upload.name).stem}
                    elif filename_accession:
                        inferred_accession = filename_accession.group(0).upper()
                        st.write(f"Recognized {inferred_accession} in the FASTA filename; retrieving its organism identity")
                        identity = get_assembly_identity(inferred_accession, ROOT)
                        if strain_name.strip():
                            identity["strain"] = strain_name.strip()
                    else:
                        identity = {"organism": "Unspecified uploaded organism",
                                    "strain": strain_name.strip() or Path(genome_upload.name).stem}
                        st.warning(
                            "No organism name or versioned NCBI accession was available. Genome annotation will "
                            "continue under a provisional identity, but automatic host-aware mobile-reference "
                            "selection and taxonomic PAM inference cannot be performed reliably."
                        )
                    source_name = f"{identity['organism']} {identity['strain']}"
                elif source_mode == "NCBI accession":
                    clean = validate_assembly_accession(ncbi_accession)
                    st.write(f"Retrieving {clean} from NCBI Datasets")
                    identity = get_assembly_identity(clean, ROOT)
                    fasta_text = download_genome_fasta(clean, ROOT)
                    source_name = f"{identity['organism']} {identity['strain']} ({clean})"
                    st.write(f"NCBI identity: {identity['organism']} • strain/isolate: {identity['strain']}")
                else:
                    if not read_uploads:
                        raise ValueError("Upload the required FASTQ files.")
                    if not organism_name.strip():
                        raise ValueError("Enter the target organism name for the raw reads.")
                    with tempfile.TemporaryDirectory(prefix="rhizoforge_reads_") as tmp:
                        read_dir = Path(tmp)
                        paths = []
                        for upload in read_uploads:
                            path = read_dir / Path(upload.name).name
                            path.write_bytes(upload.getvalue())
                            paths.append(path)
                        assembly, _ = assemble_reads(paths, read_technology, ROOT, read_dir)
                        fasta_text = assembly.read_text()
                    identity = {"organism": organism_name.strip(), "strain": strain_name.strip() or "provisional isolate"}
                    source_name = f"{identity['organism']} {identity['strain']} — raw reads ({read_technology})"
                st.write("Running available annotation and defense modules")
                reference_label = "none"
                reference_metadata = {}
                if mobile_reference_mode == "Automatic host-aware selection":
                    mobile_reference_text, reference_metadata = load_host_aware_reference(identity["organism"])
                    reference_label = (f"automatic tiered collection: {reference_metadata.get('records', 0)} records "
                                       f"(T2={reference_metadata.get('tier_2_same_species', 0)}, "
                                       f"T3={reference_metadata.get('tier_3_species_complex', 0)}, "
                                       f"T4={reference_metadata.get('tier_4_same_genus', 0)})")
                elif mobile_reference_mode == "Upload a reference FASTA" and mobile_reference_upload:
                    mobile_reference_text = mobile_reference_upload.getvalue().decode("utf-8")
                    reference_label = mobile_reference_upload.name
                    reference_metadata = {"scope": "user supplied", "records": "reported in PAM summary"}
                else:
                    mobile_reference_text = None
                if mobile_reference_mode == "Automatic host-aware selection" and mobile_reference_text is None:
                    st.warning(
                        f"No host-matched mobile records were found at the same-species, related-species-complex, "
                        f"or same-genus level for {identity['organism']}. General annotation will continue, but "
                        "strain-specific PAM discovery will remain not run. If this organism has a former or "
                        "alternative genus name, upload a curated reference FASTA or add the synonym to the registry."
                    )
                result_dir, workflow_manifest = run_genome_workflow(
                    fasta_text, source_name, ROOT, include_heavy=specialist,
                    mobile_reference_fasta=mobile_reference_text,
                )
                tables = {}
                for table_name in ["pipeline_status.tsv", "genome_qc.tsv", "uncertainty_summary.tsv", "agricultural_function_genes.tsv",
                                   "restriction_modification_systems.tsv", "mobile_elements.tsv",
                                   "plasmid_predictions.tsv", "candidate_neutral_regions.tsv",
                                   "pam_gene_context.tsv", "crispr_arrays.tsv", "crispr_spacers.tsv",
                                   "introduced_editor_targets.tsv",
                                   "effector_chassis_compatibility.tsv",
                                   "pam_protospacer_observations.tsv", "pam_64_triplet_ranking.tsv",
                                   "pam_logo_frequencies.tsv", "editing_systems_screen.tsv",
                                   "editing_system_assessments.tsv",
                                   "multiobjective_design_portfolio.tsv", "experimental_results_template.tsv"]:
                    path = result_dir / table_name
                    if path.exists():
                        # Keep the interface responsive; the ZIP always contains the complete table.
                        tables[table_name] = pd.read_csv(path, sep="\t", nrows=25_000).to_dict(orient="records")
                st.session_state["full_workflow"] = {
                    "manifest": workflow_manifest, "tables": tables,
                    "result_dir": str(result_dir),
                    "bundle": zip_results(result_dir), "source": source_name,
                    "identity": identity, "mobile_reference": reference_label,
                    "mobile_reference_metadata": reference_metadata,
                }
                status.update(label="Genome workflow complete", state="complete", expanded=False)
        except (UnicodeDecodeError, ValueError, RuntimeError, TimeoutError) as error:
            st.error(str(error), icon=":material/error:")

    full = st.session_state.get("full_workflow")
    if full:
        summary_full = full["manifest"]["summary"]
        st.success(f"Analysis package ready for {full['source']}", icon=":material/check_circle:")
        st.caption(f"Current uploaded-strain result • source: {full['source']}")
        if full.get("identity"):
            st.caption(f"Target identity: {full['identity']['organism']} • strain/isolate: {full['identity']['strain']} • mobile reference: {full.get('mobile_reference', 'none')}")
        # Competition-facing design summary
        pam_context_df = pd.DataFrame(
            full["tables"].get("pam_gene_context.tsv", [])
        )
        portfolio_df = pd.DataFrame(
            full["tables"].get("multiobjective_design_portfolio.tsv", [])
        )

        annotated_targets = int(
            summary_full.get("pam_targets_annotated", 0) or 0
        )

        if not portfolio_df.empty and "decision" in portfolio_df.columns:
            retained_designs = int(
                (portfolio_df["decision"].astype(str) == "retain").sum()
            )
        else:
            retained_designs = 0

        if not portfolio_df.empty and "pareto_optimal" in portfolio_df.columns:
            pareto_designs = int(
                (portfolio_df["pareto_optimal"].astype(str) == "yes").sum()
            )
        else:
            pareto_designs = 0

        st.markdown("## Design Summary")
        st.caption(
            "Target-space summary for strain-specific PAM-supported designs. "
            "Editing-system evidence is assessed separately below. Counts are "
            "computational decision support, not measured editing efficiency."
        )

        # Judge-facing results status strip
        native_pam_status = (
            "Supported strain-specific PAM evidence"
            if annotated_targets > 0
            else "No supported strain-specific PAM evidence"
        )

        st.markdown(
            f"""
            <div style="
                display:grid;
                grid-template-columns:repeat(4,1fr);
                gap:12px;
                margin:0.5rem 0 1.2rem 0;
            ">
                <div style="
                    border:1px solid #cfe0ef;
                    border-radius:12px;
                    padding:14px 16px;
                    background:#f7fbff;
                ">
                    <div style="font-size:12px;color:#536273;font-weight:600;">
                        ANALYSIS STATUS
                    </div>
                    <div style="font-size:17px;font-weight:700;color:#135EA8;">
                        ✓ Genome analyzed
                    </div>
                </div>

                <div style="
                    border:1px solid #cfe0ef;
                    border-radius:12px;
                    padding:14px 16px;
                    background:#f7fbff;
                ">
                    <div style="font-size:12px;color:#536273;font-weight:600;">
                        EDITOR TARGETABILITY
                    </div>
                    <div style="font-size:17px;font-weight:700;color:#135EA8;">
                        ✓ Introduced editors evaluated
                    </div>
                </div>

                <div style="
                    border:1px solid #e3bd70;
                    border-radius:12px;
                    padding:14px 16px;
                    background:#fffaf0;
                ">
                    <div style="font-size:12px;color:#76520b;font-weight:600;">
                        NATIVE PAM EVIDENCE
                    </div>
                    <div style="font-size:17px;font-weight:700;color:#76520b;">
                        ⚠ {native_pam_status}
                    </div>
                </div>

                <div style="
                    border:1px solid #d7dee7;
                    border-radius:12px;
                    padding:14px 16px;
                    background:#fafbfc;
                ">
                    <div style="font-size:12px;color:#536273;font-weight:600;">
                        EVIDENCE SCOPE
                    </div>
                    <div style="font-size:17px;font-weight:700;color:#34495e;">
                        Computational decision support
                    </div>
                </div>
            </div>
            """,
            unsafe_allow_html=True,
        )

        d1, d2, d3, d4 = st.columns(4)

        d1.metric(
            "Genome size",
            f"{int(summary_full['total_length_bp']):,} bp",
        )
        d2.metric(
            "Strain-specific PAM targets",
            f"{annotated_targets:,}",
        )
        d3.metric(
            "Retained designs",
            f"{retained_designs:,}",
        )
        d4.metric(
            "Pareto designs",
            f"{pareto_designs:,}",
        )

        target_funnel_df = pd.DataFrame(
            [
                {
                    "Stage": "Strain-specific PAM-supported targets",
                    "Candidates": annotated_targets,
                    "Interpretation": (
                        "Supported PAM motifs mapped to annotated genomic coordinates"
                    ),
                },
                {
                    "Stage": "Retained multi-objective designs",
                    "Candidates": retained_designs,
                    "Interpretation": (
                        "Targets passing the current transparent retain/reject criteria"
                    ),
                },
                {
                    "Stage": "Pareto-optimal designs",
                    "Candidates": pareto_designs,
                    "Interpretation": (
                        "Non-dominated retained designs across the decision objectives"
                    ),
                },
            ]
        )

        st.dataframe(
            target_funnel_df,
            hide_index=True,
            width="stretch",
            column_config={
                "Candidates": st.column_config.NumberColumn(
                    "Candidates",
                    format="%d",
                ),
            },
        )

        if annotated_targets == 0:
            st.warning(
                "No strain-specific PAM motif currently passes the evidence "
                "threshold, so no native-PAM target funnel is generated. "
                "This does NOT mean that the genome is uneditable. Introduced "
                "editors and reference-supported hypotheses are evaluated "
                "separately below."
            )

        # Editing Strategy Overview
        st.markdown("## Editing Strategy Overview")
        st.caption(
            "AgriPAM-AI evaluates three distinct evidence layers: "
            "genome-derived editing machinery, strain-specific targeting evidence, "
            "and introduced-editor targetability. None of these alone demonstrates "
            "experimental editing efficiency."
        )

        native_system_df = pd.DataFrame(
            full["tables"].get("editing_system_assessments.tsv", [])
        )
        compatibility_df = pd.DataFrame(
            full["tables"].get("effector_chassis_compatibility.tsv", [])
        )
        pam_ranking_df = pd.DataFrame(
            full["tables"].get("pam_ranking.tsv", [])
        )

        # 1. Genome-derived machinery
        st.markdown("### 1. Genome-derived machinery")

        native_supported = 0
        native_partial = 0

        if not native_system_df.empty:
            native_view = native_system_df.copy()

            def _native_evidence(value):
                value = str(value)
                if value == "minimum component model satisfied":
                    return "Genome-supported candidate"
                if value == "partial candidate":
                    return "Partial genomic evidence"
                if value == "not detected":
                    return "Not detected"
                return "Unresolved"

            native_view["Evidence"] = native_view["completeness"].map(
                _native_evidence
            )

            native_supported = int(
                (native_view["Evidence"] == "Genome-supported candidate").sum()
            )
            native_partial = int(
                (native_view["Evidence"] == "Partial genomic evidence").sum()
            )

            n1, n2, n3 = st.columns(3)
            n1.metric("Genome-supported routes", native_supported)
            n2.metric("Partial routes", native_partial)
            n3.metric("Models assessed", len(native_view))

            native_display = native_view[
                [
                    "system_model",
                    "Evidence",
                    "detected_component_groups",
                    "locus_organization",
                    "pam_or_target_requirement",
                    "confidence",
                ]
            ].rename(
                columns={
                    "system_model": "Strategy",
                    "detected_component_groups": "Detected components",
                    "locus_organization": "Locus organization",
                    "pam_or_target_requirement": "Target/PAM requirement",
                    "confidence": "Confidence",
                }
            )

            st.dataframe(
                native_display,
                hide_index=True,
                width="stretch",
            )
        else:
            st.info(
                "Genome-derived machinery could not be assessed in this run."
            )

        # 2. Strain-specific targeting evidence
        st.markdown("### 2. Strain-specific targeting evidence")

        supported_pam_rows = pd.DataFrame()

        if not pam_ranking_df.empty:
            observations = pd.to_numeric(
                pam_ranking_df.get("observations", 0),
                errors="coerce",
            ).fillna(0)

            q_values = pd.to_numeric(
                pam_ranking_df.get("bh_q_value", 1),
                errors="coerce",
            ).fillna(1)

            supported_mask = (observations >= 3) & (q_values < 0.10)
            supported_pam_rows = pam_ranking_df.loc[supported_mask].copy()

        p1, p2, p3 = st.columns(3)
        p1.metric(
            "Supported strain-specific PAMs",
            len(supported_pam_rows),
        )
        p2.metric(
            "Annotated strain-specific targets",
            annotated_targets,
        )
        p3.metric(
            "Reference collection",
            (
                f"{len(full.get('mobile_reference_metadata', {}).get('records', [])):,}"
                if isinstance(
                    full.get("mobile_reference_metadata", {}).get("records", []),
                    list,
                )
                else "available"
            ),
        )

        if not supported_pam_rows.empty:
            pam_columns = [
                column
                for column in (
                    "pam3",
                    "observations",
                    "bh_q_value",
                    "enrichment",
                )
                if column in supported_pam_rows.columns
            ]

            st.dataframe(
                supported_pam_rows[pam_columns],
                hide_index=True,
                width="stretch",
            )

            st.success(
                "At least one strain-specific PAM candidate passes the current "
                "computational evidence threshold. These motifs remain hypotheses "
                "until experimentally validated."
            )
        else:
            st.info(
                "No strain-specific PAM currently passes the evidence threshold. "
                "AgriPAM-AI therefore does not substitute TTC or another "
                "reference PAM as strain-specific evidence."
            )

        # 3. Introduced editors
        st.markdown("### 3. Introduced editor options")

        if not compatibility_df.empty:
            st.components.v1.html(
                build_introduced_editor_funnel(
                    compatibility_df.to_dict("records")
                ),
                height=560,
                scrolling=False,
            )

        if not compatibility_df.empty:
            introduced_view = compatibility_df.copy()

            for column in (
                "pam_compatible_sites",
                "exact_unique_sites",
                "approximately_screened_top_sites",
            ):
                if column in introduced_view.columns:
                    introduced_view[column] = pd.to_numeric(
                        introduced_view[column],
                        errors="coerce",
                    ).fillna(0).astype(int)

            targetable_editors = int(
                (
                    introduced_view["pam_compatible_sites"] > 0
                ).sum()
            )

            total_editor_sites = int(
                introduced_view["pam_compatible_sites"].sum()
            )

            unique_editor_sites = int(
                introduced_view["exact_unique_sites"].sum()
            )

            e1, e2, e3 = st.columns(3)
            e1.metric("Targetable introduced editors", targetable_editors)
            e2.metric("PAM-compatible sites", f"{total_editor_sites:,}")
            e3.metric("Exact-unique sites", f"{unique_editor_sites:,}")
            def screening_display(row):
                status = str(
                    row.get(
                        "off_target_screening_status",
                        "",
                    )
                )

                screened_value = row.get(
                    "approximately_screened_top_sites",
                    None,
                )
                clean_value = row.get(
                    "screened_sites_without_1_or_2_mismatch_hit",
                    None,
                )

                screened = (
                    int(screened_value)
                    if pd.notna(screened_value)
                    else 0
                )

                if status == "screened" and screened > 0:
                    clean = (
                        int(clean_value)
                        if pd.notna(clean_value)
                        else None
                    )
                    if clean is not None:
                        return f"{clean:,} / {screened:,} screened"
                    return f"{screened:,} screened; clean count unresolved"

                if status == "not_screened_reference_limit":
                    pam_sites = int(
                        row.get("pam_compatible_sites", 0) or 0
                    )
                    return (
                        f"Not screened — {pam_sites:,} PAM-compatible sites "
                        "exceed the 500,000 reference-space limit"
                    )

                if status == "not_screened":
                    return "Not screened — off-target status unresolved"

                pam_sites = int(
                    row.get("pam_compatible_sites", 0) or 0
                )

                if pam_sites > 500_000:
                    return (
                        f"Not screened — {pam_sites:,} PAM-compatible sites "
                        "exceed the 500,000 reference-space limit"
                    )

                return "Not screened — off-target status unresolved"

            introduced_view["screening_summary"] = introduced_view.apply(
                screening_display,
                axis=1,
            )
            editor_columns = [
                column
                for column in (
                    "editor",
                    "pam_rule",
                    "pam_compatible_sites",
                    "exact_unique_sites",
                    "screening_summary",
                    "compatibility_tier",
                    "repair_assessment",
                    "delivery_assessment",
                )
                if column in introduced_view.columns
            ]

            editor_display = introduced_view[editor_columns].rename(
                columns={
                    "editor": "Editor",
                    "pam_rule": "PAM rule",
                    "pam_compatible_sites": "PAM-compatible sites",
                    "exact_unique_sites": "Exact-unique sites",
                    "screening_summary": "Off-target screening",
                    "compatibility_tier": "Compatibility assessment",
                    "repair_assessment": "Repair evidence",
                    "delivery_assessment": "Delivery / defense",
                }
            )

            st.dataframe(
                editor_display,
                hide_index=True,
                width="stretch",
            )

            st.caption(
                "Introduced-editor targetability is independent of strain-specific "
                "native PAM inference. PAM-compatible sites establish sequence "
                "targetability only; chassis function, delivery, toxicity and editing "
                "efficiency still require experimental validation."
            )
        else:
            st.info(
                "Introduced-editor targetability was not available in this run."
            )

        # Interactive introduced-editor Candidate Explorer
        introduced_targets_df = pd.DataFrame(
            full["tables"].get("introduced_editor_targets.tsv", [])
        )

        result_dir = Path(full["result_dir"])
        annotation_path = result_dir / "annotation.gff"
        annotation_features = (
            read_gff(annotation_path)
            if annotation_path.exists()
            else []
        )

        if not introduced_targets_df.empty:
            project_modification = str(project_brief.get("modification", "")).strip()
            route_context_ready = bool(active_design) or project_modification not in ("", "Not yet defined")
            if route_context_ready and not compatibility_df.empty:
                st.markdown("### Objective-aware editing-route recommendation")
                st.write(
                    "This step connects the intended modification, selected community and chassis to the exact genome analysis. "
                    "It compares introduced editors, removes candidates whose annotations overlap a function "
                    "you asked to preserve, and prioritizes a route appropriate for the stated objective."
                )

                syncom_objective = str((active_design or {}).get("objective", "")).strip()
                objective_text = "; ".join(
                    value for value in (project_modification, syncom_objective) if value and value != "Not yet defined"
                )
                protected_functions = [
                    str(value).strip()
                    for value in (active_design or {}).get("protected_functions", [])
                    if str(value).strip()
                ]
                ignored_terms = {
                    "activity", "agricultural", "benefit", "function", "growth",
                    "production", "plant", "preserve", "existing", "other",
                }

                def _design_terms(values):
                    terms = set()
                    for value in values:
                        terms.update(
                            token
                            for token in re.findall(r"[a-z0-9]+", value.lower())
                            if len(token) >= 4 and token not in ignored_terms
                        )
                    return terms

                protected_terms = _design_terms(protected_functions)
                objective_terms = _design_terms([objective_text])
                route_candidates = introduced_targets_df.copy()
                product_text = route_candidates.get(
                    "product", pd.Series("", index=route_candidates.index)
                ).fillna("").astype(str).str.lower()
                route_candidates["protected_annotation_overlap"] = product_text.apply(
                    lambda value: any(term in value for term in protected_terms)
                )
                route_candidates["objective_annotation_match"] = product_text.apply(
                    lambda value: any(term in value for term in objective_terms)
                )
                route_candidates["is_intergenic_candidate"] = route_candidates.get(
                    "gene_relation", pd.Series("", index=route_candidates.index)
                ).fillna("").astype(str).str.lower().str.contains("intergenic|uncalled", regex=True)
                route_candidates["is_exact_unique"] = route_candidates.get(
                    "exact_uniqueness", pd.Series("", index=route_candidates.index)
                ).fillna("").astype(str).str.lower().eq("unique")
                route_score_column = (
                    "offtarget_adjusted_priority_score"
                    if "offtarget_adjusted_priority_score" in route_candidates.columns
                    else "guide_priority_score"
                )
                route_candidates["_route_score"] = pd.to_numeric(
                    route_candidates.get(route_score_column, 0), errors="coerce"
                ).fillna(0)

                regulation_objective = any(
                    term in objective_text.lower()
                    for term in ("regulat", "repress", "activate", "expression", "preserve")
                )
                addition_objective = any(
                    term in objective_text.lower()
                    for term in ("add", "insert", "introduc", "missing", "benefit")
                )
                compatibility_lookup = {
                    str(row.get("editor", "")): row
                    for row in compatibility_df.to_dict("records")
                }
                route_rows = []
                safe_candidate_sets = {}
                for editor, group in route_candidates.groupby("editor", sort=False):
                    editor = str(editor)
                    unprotected = group.loc[~group["protected_annotation_overlap"]].copy()
                    exact_unique = unprotected.loc[unprotected["is_exact_unique"]].copy()
                    preferred_pool = exact_unique if not exact_unique.empty else unprotected
                    if addition_objective:
                        neutral_pool = preferred_pool.loc[preferred_pool["is_intergenic_candidate"]]
                        if not neutral_pool.empty:
                            preferred_pool = neutral_pool
                    elif objective_terms:
                        objective_pool = preferred_pool.loc[preferred_pool["objective_annotation_match"]]
                        if not objective_pool.empty:
                            preferred_pool = objective_pool
                    preferred_pool = preferred_pool.sort_values("_route_score", ascending=False)
                    safe_candidate_sets[editor] = preferred_pool
                    compatibility = compatibility_lookup.get(editor, {})
                    route_priority = float(preferred_pool["_route_score"].max()) if not preferred_pool.empty else -1
                    if regulation_objective and editor == "dCas9/CRISPRi":
                        route_priority += 25
                    if not regulation_objective and editor == "dCas9/CRISPRi":
                        route_priority -= 10
                    if addition_objective:
                        route_priority += min(20, 2 * int(preferred_pool["is_intergenic_candidate"].sum()))
                    if str(compatibility.get("compatibility_tier", "")).startswith("computationally compatible"):
                        route_priority += 10
                    route_rows.append({
                        "Editing route": editor,
                        "PAM": compatibility.get("pam_rule", "not reported"),
                        "Genome-wide sites": int(compatibility.get("pam_compatible_sites", len(group)) or 0),
                        "Protected-overlap warnings": int(group["protected_annotation_overlap"].sum()),
                        "Remaining exact-unique": int(len(exact_unique)),
                        "Candidate neutral regions": int(unprotected["is_intergenic_candidate"].sum()),
                        "Objective-annotation matches": int(unprotected["objective_annotation_match"].sum()),
                        "Chassis assessment": compatibility.get("compatibility_tier", "unresolved"),
                        "_priority": route_priority,
                    })

                route_summary = pd.DataFrame(route_rows).sort_values(
                    ["_priority", "Remaining exact-unique"], ascending=[False, False]
                ).reset_index(drop=True)
                route_summary.insert(0, "Priority", range(1, len(route_summary) + 1))
                best_editor = str(route_summary.iloc[0]["Editing route"])
                best_pool = safe_candidate_sets.get(best_editor, pd.DataFrame())
                st.success(
                    f"Highest-priority computational route for **{(active_design or {}).get('chassis', 'the analysed chassis')}**: "
                    f"**{best_editor}** for the objective **{objective_text or 'not defined'}**. "
                    "This is a design recommendation, not a prediction of successful editing."
                )
                st.dataframe(
                    route_summary.drop(columns=["_priority"]), hide_index=True, width="stretch"
                )
                if not best_pool.empty:
                    top = best_pool.iloc[0]
                    target_context = str(top.get("gene_relation", "unresolved"))
                    target_product = str(top.get("product", "")).strip() or "no annotated product"
                    st.markdown(
                        f"**Top loaded candidate for this route:** `{top.get('contig', '')}:"
                        f"{int(top.get('start', 0))}-{int(top.get('end', 0))}`; PAM "
                        f"**{top.get('pam', 'NA')}**; context **{target_context}**; annotation "
                        f"**{target_product}**. Inspect it in the candidate designer below before export."
                    )
                if protected_functions:
                    st.caption(
                        "Protected functions requested: " + ", ".join(protected_functions) + ". "
                        "Exclusion here is based on annotation-keyword overlap. It cannot replace essential-gene, "
                        "operon, promoter, polar-effect or phenotype analysis."
                    )
                st.warning(
                    "A candidate intergenic site is not automatically a validated neutral insertion site. "
                    "Confirm regulatory context, neighboring genes, operon structure, essentiality, off-targets, "
                    "delivery and repair in the exact isolate before construct design."
                )

            st.markdown("### Editing Strategy Candidate Designer")
            st.caption(
                "Genome-wide target discovery → PAM compatibility → exact uniqueness → "
                "bounded specificity assessment → candidate prioritization. Genomic context "
                "is drawn only from the workflow annotation and real candidate coordinates."
            )

            available_editors = (
                introduced_targets_df["editor"]
                .dropna()
                .astype(str)
                .drop_duplicates()
                .tolist()
            )

            selected_editor = st.selectbox(
                "Editing strategy",
                available_editors,
                key="introduced_candidate_editor",
            )

            editor_candidates = introduced_targets_df[
                introduced_targets_df["editor"].astype(str) == selected_editor
            ].copy()

            # Editor-specific genome-wide funnel. Counts come from the complete workflow summary,
            # not from the bounded candidate table loaded into the interface.
            editor_summary = introduced_view[
                introduced_view["editor"].astype(str) == selected_editor
            ]

            if not editor_summary.empty:
                summary_row = editor_summary.iloc[0]
                pam_site_count = int(summary_row.get("pam_compatible_sites", 0) or 0)
                unique_site_count = int(summary_row.get("exact_unique_sites", 0) or 0)
                screening_summary = str(summary_row.get("screening_summary", "Not screened"))

                f1, f2, f3, f4 = st.columns(4)
                f1.metric("PAM-compatible", f"{pam_site_count:,}")
                f2.metric("Exact-unique", f"{unique_site_count:,}")
                f3.metric("Off-target screening", screening_summary)
                f4.metric("Loaded candidates", f"{len(editor_candidates):,}")

                st.caption(
                    "PAM-compatible and exact-unique values are genome-wide workflow counts. "
                    "Loaded candidates are the bounded records retained in the interactive interface; "
                    "they are not the total genome-wide candidate count."
                )

            # Optional genomic-context filter based only on workflow annotation.
            context_filter = st.radio(
                "Genomic context",
                ["All", "Within CDS", "Intergenic / uncalled"],
                horizontal=True,
                key="introduced_candidate_context",
            )

            if "gene_relation" in editor_candidates.columns:
                relation = editor_candidates["gene_relation"].fillna("").astype(str).str.lower()
                if context_filter == "Within CDS":
                    editor_candidates = editor_candidates[relation.str.contains("within")].copy()
                elif context_filter == "Intergenic / uncalled":
                    editor_candidates = editor_candidates[
                        relation.str.contains("intergenic|uncalled", regex=True)
                    ].copy()

            if editor_candidates.empty:
                st.info(
                    "No loaded candidates match this genomic-context filter. "
                    "This reflects the current workflow annotation and loaded candidate subset."
                )

            # Prioritize exact-unique candidates and the workflow priority score.
            if "exact_uniqueness" in editor_candidates.columns:
                editor_candidates["_unique_rank"] = (
                    editor_candidates["exact_uniqueness"]
                    .astype(str)
                    .eq("unique")
                    .astype(int)
                )
            else:
                editor_candidates["_unique_rank"] = 0

            score_column = (
                "offtarget_adjusted_priority_score"
                if "offtarget_adjusted_priority_score" in editor_candidates.columns
                else "guide_priority_score"
            )

            if score_column in editor_candidates.columns:
                editor_candidates["_display_score"] = pd.to_numeric(
                    editor_candidates[score_column],
                    errors="coerce",
                ).fillna(-1)
            else:
                editor_candidates["_display_score"] = -1

            editor_candidates = editor_candidates.sort_values(
                ["_unique_rank", "_display_score"],
                ascending=[False, False],
            ).reset_index(drop=True)

            # Keep the interactive selector compact.
            explorer_candidates = editor_candidates.head(250).copy()

            candidate_labels = []
            for idx, row in explorer_candidates.iterrows():
                score_value = row.get(score_column, "")
                score_label = (
                    f"{float(score_value):.1f}"
                    if pd.notna(score_value)
                    else "NA"
                )
                candidate_labels.append(
                    f"{idx + 1}. {row.get('contig', '')} "
                    f"{int(row.get('start', 0)):,}-{int(row.get('end', 0)):,} "
                    f"({row.get('strand', '')}) | "
                    f"{row.get('pam', '')} | score {score_label}"
                )

            selected_label = st.selectbox(
                "Candidate target",
                candidate_labels,
                key="introduced_candidate_target",
            )

            selected_index = candidate_labels.index(selected_label)
            selected_candidate = explorer_candidates.iloc[
                selected_index
            ].to_dict()

            c1, c2, c3, c4 = st.columns(4)

            c1.metric(
                "Observed PAM",
                str(selected_candidate.get("pam", "NA")),
            )

            c2.metric(
                "GC",
                (
                    f"{float(selected_candidate.get('gc_percent')):.1f}%"
                    if pd.notna(selected_candidate.get("gc_percent"))
                    else "NA"
                ),
            )

            c3.metric(
                "Exact uniqueness",
                str(
                    selected_candidate.get(
                        "exact_uniqueness",
                        selected_candidate.get(
                            "guide_uniqueness",
                            "NA",
                        ),
                    )
                ),
            )

            selected_score = selected_candidate.get(score_column)
            c4.metric(
                "Priority score",
                (
                    f"{float(selected_score):.1f}"
                    if pd.notna(selected_score)
                    else "NA"
                ),
            )

            explorer_html = build_editor_target_explorer(
                candidate=selected_candidate,
                features=annotation_features,
                flank_bp=2500,
            )

            st.components.v1.html(
                explorer_html,
                height=680,
                scrolling=False,
            )

            screen_status = str(
                selected_candidate.get(
                    "off_target_screen_status",
                    "not reported",
                )
            )

            if screen_status == "skipped_reference_limit":
                st.warning(
                    "Approximate off-target screening was not performed for "
                    "this candidate because the PAM-compatible reference space "
                    "exceeded the configured 500,000-sequence limit. "
                    "Off-target risk therefore remains unresolved; a zero "
                    "off-target count must not be inferred."
                )
            elif screen_status == "complete_substitutions_only":
                st.info(
                    "Approximate off-target screening was completed against "
                    "PAM-compatible sites for substitutions up to two "
                    "mismatches. Bulges, indels and noncanonical PAMs remain "
                    "outside this bounded screen."
                )
            elif screen_status == "not_prioritized":
                st.info(
                    "This candidate was not among the prioritized guides used "
                    "for the bounded approximate off-target screen."
                )

            with st.expander("Selected candidate evidence"):
                evidence_fields = [
                    "editor",
                    "application",
                    "pam_pattern",
                    "pam",
                    "protospacer",
                    "strand",
                    "contig",
                    "start",
                    "end",
                    "pam_start",
                    "pam_end",
                    "gene_relation",
                    "feature_id",
                    "product",
                    "gc_percent",
                    "exact_genome_copies",
                    "exact_uniqueness",
                    "guide_priority_score",
                    "off_target_screen_status",
                    "approx_offtargets_le2",
                    "nearest_approx_offtarget_mismatches",
                    "offtarget_adjusted_priority_score",
                    "evidence_level",
                    "validation_required",
                ]

                evidence = {
                    field: selected_candidate.get(field, "")
                    for field in evidence_fields
                    if field in selected_candidate
                }

                st.dataframe(
                    pd.DataFrame(
                        {
                            "Evidence field": list(evidence.keys()),
                            "Result": list(evidence.values()),
                        }
                    ),
                    hide_index=True,
                    width="stretch",
                )

        else:
            st.info(
                "No introduced-editor candidate table is available for this run."
            )


        # Visual genome-informed strategy map
        strategy_map_html = build_strategy_map(
            native_systems=native_system_df.to_dict("records"),
            introduced_editors=compatibility_df.to_dict("records"),
            retained_designs=retained_designs,
            pareto_designs=pareto_designs,
        )
        st.components.v1.html(
            strategy_map_html,
            height=650,
            scrolling=False,
        )

        reference_meta = full.get("mobile_reference_metadata", {})
        if reference_meta.get("records"):
            st.subheader("Automatic mobile-reference evidence")
            tier_table = pd.DataFrame([
                {"Evidence tier": "Tier 2 — same species", "Records": reference_meta.get("tier_2_same_species", 0),
                 "Permitted interpretation": "Strongest taxonomic relevance; still computational"},
                {"Evidence tier": "Tier 3 — related species complex", "Records": reference_meta.get("tier_3_species_complex", 0),
                 "Permitted interpretation": "Transferred supporting evidence; moderate at best"},
                {"Evidence tier": "Tier 4 — same genus", "Records": reference_meta.get("tier_4_same_genus", 0),
                 "Permitted interpretation": "Exploratory evidence only; never treated as strain-specific proof"},
            ])
            st.dataframe(tier_table, hide_index=True, width="stretch")
            inventory = pd.DataFrame(reference_meta.get("reference_inventory", []))
            if not inventory.empty:
                with st.expander("Inspect the phage and plasmid records used in each evidence tier"):
                    st.caption(
                        "These are reference candidates supplied to spacer–protospacer matching. Inclusion in this "
                        "inventory is not itself evidence of infection or PAM recognition by the analyzed isolate."
                    )
                    tier_names = {
                        2: "Tier 2 — same species", 3: "Tier 3 — related species complex",
                        4: "Tier 4 — same genus",
                    }
                    for tier_number, tier_name in tier_names.items():
                        tier_records = inventory[inventory["evidence_tier"] == tier_number]
                        with st.expander(f"{tier_name} ({len(tier_records)} records)"):
                            st.dataframe(
                                tier_records,
                                hide_index=True,
                                width="stretch",
                                column_config={"ncbi_record": st.column_config.LinkColumn("Open NCBI record")},
                            )
                    st.download_button(
                        "Download complete mobile-reference inventory",
                        inventory.to_csv(index=False, sep="\t"),
                        file_name="host_aware_mobile_reference_inventory.tsv",
                        mime="text/tab-separated-values",
                    )
            st.caption(
                "The collection is constructed automatically from locally installed RefSeq virus and NCBI plasmid "
                "records. Taxonomic proximity is recorded for every match. It does not guarantee that a mobile "
                "element infected this isolate."
            )
        stage_rows = full["manifest"].get("stages", [])
        annotation_stage = next((row for row in stage_rows if row.get("name") == "coding_sequence_annotation"), {})
        annotation_ok = annotation_stage.get("status") == "complete"
        product_annotations = annotation_stage.get("method") in {"Bakta", "Prokka"}
        m1, m2, m3, m4 = st.columns(4)
        m1.metric("Genome size", f"{int(summary_full['total_length_bp']):,} bp")
        m2.metric("Coding sequences", f"{int(summary_full['coding_sequences']):,}" if annotation_ok else "Unavailable")
        m3.metric(
            "Functional hits",
            f"{int(summary_full['functional_gene_hits']):,}" if product_annotations else "Unavailable",
            help=("Genes whose names or predicted products matched the agricultural-function categories used by "
                  "AgriPAM-AI, including biocontrol, nutrient transformation, secretion, siderophore and "
                  "colonization functions. These are annotation-supported candidates, not experimentally confirmed activities."),
        )
        m4.metric(
            "Neutral candidates",
            f"{int(summary_full['neutral_region_candidates']):,}" if annotation_ok else "Unavailable",
            help=("Preliminary genomic regions located away from predicted coding and excluded features. "
                  "They may be useful starting points for reporter or construct integration, but they are not "
                  "confirmed safe-harbour sites. Verify nearby regulatory elements, conservation, essentiality, "
                  "growth effects and preservation of agricultural functions experimentally."),
        )
        qc_table = pd.DataFrame(full["tables"].get("genome_qc.tsv", []))
        if not qc_table.empty:
            with st.expander("Genome quality control", expanded=True):
                st.dataframe(qc_table, hide_index=True, width="stretch")
                if (qc_table["status"] == "review").any():
                    st.warning("One or more generic assembly-QC checks require review before treating target absence or uniqueness as reliable.")
                else:
                    st.success("No generic FASTA quality threshold was triggered. This does not replace completeness or contamination assessment.")

        if not product_annotations:
            st.warning(
                "Named functional-gene screening is unavailable because a database-backed Bakta/Prokka annotation "
                "was not completed. A Prodigal fallback can still recover CDS coordinates, proteins and preliminary "
                "neutral-region candidates; it must not be interpreted as functional annotation."
            )
        status_df = pd.DataFrame(full["tables"].get("pipeline_status.tsv", []))
        st.subheader("Module status")
        st.dataframe(status_df, hide_index=True, width="stretch")
        with st.expander("Confidence and uncertainty by module", expanded=True):
            st.dataframe(pd.DataFrame(full["tables"].get("uncertainty_summary.tsv", [])), hide_index=True, width="stretch")
        with st.expander("Editing-readiness assessment", expanded=True):
            status_lookup = dict(zip(status_df.get("name", []), status_df.get("status", []))) if not status_df.empty else {}
            def _state(*names):
                values = [str(status_lookup.get(name, "not_run")) for name in names]
                return "available" if any(v in {"complete", "available", "success"} for v in values) else "not detected / not run"
            readiness = pd.DataFrame([
                {"Evidence layer": "Genome annotation and coding sequences", "Result": "available", "Interpretation": "Required for mapping candidate edits to genes"},
                {"Evidence layer": "CRISPR arrays and Cas operons", "Result": _state("crispr_cas", "crispr_detection", "cctyper"), "Interpretation": "A PAM can be assigned only when a compatible system is supported"},
                {"Evidence layer": "Defense barriers (restriction–modification / other systems)", "Result": _state("defensefinder", "restriction_modification", "defense"), "Interpretation": "Informs DNA-delivery and construct-stability risk"},
                {"Evidence layer": "Mobile elements and spacer matches", "Result": _state("strain_specific_pam_discovery", "mobile_elements"), "Interpretation": "Supports strain-specific PAM inference when reference sequences are supplied"},
                {"Evidence layer": "Alternative editing machinery", "Result": "screened as exploratory", "Interpretation": "Recombinase, retron and transposase signals remain hypotheses requiring validation"},
                {"Evidence layer": "Candidate targets and neutral regions", "Result": "available after annotation", "Interpretation": "Prioritized computational candidates; not proven safe harbors"},
            ])
            st.dataframe(readiness, hide_index=True, width="stretch")

        st.subheader("Introduced effector–chassis compatibility")
        st.dataframe(pd.DataFrame(full["tables"].get("effector_chassis_compatibility.tsv", [])), hide_index=True, width="stretch")
        st.caption("Compatibility tiers combine genome-derived targeting, defense, delivery and repair evidence. They are not editing-efficiency estimates.")
        st.info("This assessment reports computational readiness. It does not establish editing efficiency, viability, or a safe integration site; those require institutionally approved experiments.", icon=":material/science:")
        with st.expander("Introduced editor targetability (NGG / TTTV)", expanded=True):
            st.caption("These are sites compatible with introduced editors, separate from native PAM inference. A match is computationally targetable—not proof of delivery, activity, repair or successful editing.")
            introduced = pd.DataFrame(full["tables"].get("introduced_editor_targets.tsv", []))
            if introduced.empty:
                st.info("No introduced-editor sites were found in the analyzed sequences.")
            else:
                editor_filter = st.multiselect("Editor", sorted(introduced["editor"].unique()), default=sorted(introduced["editor"].unique()))
                visible = introduced[introduced["editor"].isin(editor_filter)]
                if "off_target_screen_status" in visible:
                    screen_filter = st.multiselect("Off-target screen status", sorted(visible["off_target_screen_status"].unique()),
                                                   default=sorted(visible["off_target_screen_status"].unique()))
                    visible = visible[visible["off_target_screen_status"].isin(screen_filter)]
                st.metric("PAM-compatible sites", len(visible)); st.dataframe(visible, hide_index=True, width="stretch")
                st.download_button("Download introduced-editor targets", visible.to_csv(sep="\t", index=False).encode(), "introduced_editor_targets.tsv", mime="text/tab-separated-values")
        result_tabs = st.tabs(["Strain-specific PAM discovery", "PAM → gene consequences", "Agricultural functions", "Defense and mobility", "Neutral regions"])
        with result_tabs[0]:
            pam_ranking = pd.DataFrame(full["tables"].get("pam_64_triplet_ranking.tsv", []))
            pam_observations = pd.DataFrame(full["tables"].get("pam_protospacer_observations.tsv", []))
            pam_logo = pd.DataFrame(full["tables"].get("pam_logo_frequencies.tsv", []))
            if pam_ranking.empty:
                st.info(
                    "No strain-specific PAM ranking is available. A native PAM can be inferred only when both "
                    "usable CRISPR spacers and matching external mobile-element sequences are present. This result "
                    "does not mean that the bacterium cannot be edited with an introduced CRISPR system."
                )
            else:
                supported = pam_ranking[pam_ranking["observations"] > 0]
                st.caption("All 64 triplets are ranked. Zero-count motifs are retained so absence of evidence is visible.")
                st.info(
                    "Evidence tiers: Tier 2 = same species; Tier 3 = related species complex; Tier 4 = same genus. "
                    "Tier-4 matches are exploratory and cannot establish a strain-specific PAM without direct validation."
                )
                st.dataframe(supported.head(20), hide_index=True, width="stretch")
                if not pam_logo.empty:
                    st.bar_chart(pam_logo, x="position", y="frequency", color="nucleotide", stack=True)
                st.write("Auditable spacer–protospacer observations")
                st.dataframe(pam_observations, hide_index=True, width="stretch")
                st.warning(
                    "This is computational inference from exact database matches—not experimental PAM validation. "
                    "Confirm candidate PAMs with an independent search, locus verification, and a PAM-library or "
                    "interference-reporter assay including PAM-mutant and no-guide controls."
                )
        with result_tabs[1]:
            pam_df = pd.DataFrame(full["tables"].get("pam_gene_context.tsv", []))
            strain_specific_status = ""
            if not status_df.empty:
                matches = status_df[status_df["name"] == "strain_specific_pam_discovery"]
                if not matches.empty:
                    strain_specific_status = str(matches.iloc[0]["status"])
            target_basis = ("strain-specific supported PAMs" if strain_specific_status == "complete"
                            else "reference TTC assumption" if strain_specific_status == "not_run"
                            else "strain-specific evidence")
            st.caption(
                f"{int(summary_full.get('pam_targets_annotated', 0)):,} targets were annotated using {target_basis}. "
                "The on-screen preview is limited to 25,000 rows; the downloaded ZIP contains every target."
            )
            if not pam_df.empty:
                context_filter = st.multiselect(
                    "Gene context", sorted(pam_df["gene_relation"].dropna().unique()),
                    default=sorted(pam_df["gene_relation"].dropna().unique()),
                )
                priority_only = st.toggle("Show only growth- or agronomy-flagged targets")
                visible_pam = pam_df[pam_df["gene_relation"].isin(context_filter)]
                if priority_only:
                    visible_pam = visible_pam[
                        (visible_pam["growth_relevance_screen"] != "not flagged") |
                        (visible_pam["agronomic_category_screen"] != "not flagged")
                    ]
                st.dataframe(
                    visible_pam,
                    hide_index=True,
                    width="stretch",
                    column_config={
                        "protospacer": st.column_config.TextColumn("Target sequence", width="large"),
                        "product": st.column_config.TextColumn("Gene product", width="large"),
                        "potential_consequence": st.column_config.TextColumn("Potential consequence", width="large"),
                    },
                )
            st.warning(
                "A growth flag is not proof of essentiality, and an agronomic flag is not proof of phenotype. "
                "Confirm the selected locus with essential-gene evidence, synteny and wet-lab validation."
            )
        with result_tabs[2]:
            st.dataframe(pd.DataFrame(full["tables"].get("agricultural_function_genes.tsv", [])), hide_index=True, width="stretch")
        with result_tabs[3]:
            for label, filename in [("Restriction–modification systems", "restriction_modification_systems.tsv"),
                                    ("Mobile elements", "mobile_elements.tsv"),
                                    ("Plasmid predictions", "plasmid_predictions.tsv")]:
                st.markdown(f"**{label}**")
                st.dataframe(pd.DataFrame(full["tables"].get(filename, [])), hide_index=True, width="stretch")
        with result_tabs[4]:
            st.warning("Candidate neutral regions are computational exclusions, not proven safe-harbor loci.")
            st.dataframe(pd.DataFrame(full["tables"].get("candidate_neutral_regions.tsv", [])), hide_index=True, width="stretch")
        st.download_button(
            "Download complete results package", full["bundle"],
            file_name="rhizoforge_select_results.zip", mime="application/zip",
            type="primary", icon=":material/download:",
        )

with tabs[3]:
    st.header("Agricultural Microorganism Editing Knowledgebase")
    st.info(
        "A strain-resolved evidence catalogue connecting agricultural microorganisms to published genome-editing "
        "routes, PAM requirements, genome records and validation needs. It is separate from the fixed three-species "
        "comparative study.",
        icon=":material/menu_book:",
    )
    st.warning(
        "‘PAM unknown’ or ‘no native CRISPR system reported’ does not mean that the microorganism cannot be edited. "
        "It means that the available evidence does not support that conclusion; an introduced nuclease, "
        "recombineering or another editing route may still be feasible."
    )
    st.markdown("### 1. Add your microbial bank and rank the community")
    st.write(
        "Begin by adding the agricultural characteristics measured for your isolates and their pairwise "
        "compatibility or antagonism results. AgriPAM-AI combines those observations with the current "
        "knowledgebase to rank candidate SynCom members by agricultural contribution, functional coverage, "
        "community compatibility, interaction evidence and project biosafety flags."
    )
    knowledgebase_bank_entry = st.container()

    cohort_size = (
        native_cohort["strain"].astype(str).str.strip().replace("", pd.NA).nunique()
        if "strain" in native_cohort.columns else 0
    )
    with st.expander(
        f"Microbial-bank native-system discovery registry — {cohort_size} isolates",
        expanded=True,
    ):
        st.write(
            "The count reflects the unique isolate codes currently in this discovery registry and updates "
            "when registry records are added. Uploaded bank measurements support SynCom analysis separately; "
            "new isolates require discovery records and analysis before results appear here. Priority helps "
            "select candidates for sequencing and native-system validation; it does not establish that a "
            "CRISPR-Cas system is present or active."
        )
        if native_cohort.empty:
            st.info("The isolate cohort registry is not available.")
        else:
            cohort_priorities = list(native_cohort["native_discovery_priority"].drop_duplicates())
            selected_cohort_priorities = st.multiselect(
                "Native-system discovery priority",
                cohort_priorities,
                default=cohort_priorities,
                key="native_cohort_priority_filter",
            )
            visible_cohort = native_cohort[native_cohort["native_discovery_priority"].isin(selected_cohort_priorities)]
            lead_count = int(native_cohort["native_discovery_priority"].isin(["Very high", "High"]).sum())
            c1, c2, c3 = st.columns(3)
            c1.metric("Isolates retained", len(native_cohort))
            c2.metric("Lead discovery candidates", lead_count)
            c3.metric("Biosafety holds", int((native_cohort["native_discovery_priority"] == "Hold").sum()))
            st.dataframe(
                visible_cohort,
                hide_index=True,
                width="stretch",
                column_config={
                    "strain": st.column_config.TextColumn("Isolate", pinned=True),
                    "identification": st.column_config.TextColumn("Provisional identification", pinned=True),
                    "native_discovery_priority": st.column_config.TextColumn("Discovery priority"),
                    "rationale": st.column_config.TextColumn("Why investigate", width="large"),
                    "required_next_step": st.column_config.TextColumn("Next evidence required", width="large"),
                },
            )
            st.download_button(
                "Download native-system discovery cohort",
                native_cohort.to_csv(sep="\t", index=False).encode("utf-8"),
                "agricultural_isolate_native_system_cohort.tsv",
                mime="text/tab-separated-values",
            )
            st.caption(
                "Recommended first sequencing set: B26 and B35; paired discovery sets: B59/B58 and CC29/B4; "
                "B34 is retained as a comparative reference member. Blank original measurements remain ‘Not measured’."
            )

    st.subheader("Curated evidence database")
    st.write(
        "Upload reviewed CSV or TSV records from papers, patents or genome databases. If the database is empty, "
        "the first valid upload initializes it. Later uploads add new records and update matching records; they "
        "do not replace the rest of the database. A matching record has the same organism, strain, editing route "
        "and source ID. Changes remain in this browser session."
    )
    if "agricultural_kb_records" not in st.session_state:
        st.session_state["agricultural_kb_records"] = agricultural_kb.copy()
    kb = st.session_state["agricultural_kb_records"].copy()

    imported_kb = st.file_uploader(
        "Add reviewed records (TSV or CSV)",
        type=["tsv", "csv"],
        key="agricultural_kb_import",
        help="Use the downloadable table as the schema. Existing matching records are updated; all others remain.",
    )
    upload_token = getattr(imported_kb, "file_id", None) if imported_kb is not None else None
    if imported_kb is not None and upload_token != st.session_state.get("agricultural_kb_last_upload"):
        try:
            incoming_kb = pd.read_csv(imported_kb, sep="\t" if imported_kb.name.lower().endswith(".tsv") else ",", keep_default_na=False)
            kb, merge_summary = merge_records(kb, incoming_kb, agricultural_kb.columns)
            st.session_state["agricultural_kb_records"] = kb
            st.session_state["agricultural_kb_last_upload"] = upload_token
            st.session_state["agricultural_kb_notice"] = (
                f"Upload applied: {merge_summary['added']} added, {merge_summary['updated']} updated and "
                f"{merge_summary['unchanged']} unchanged. The database now contains {len(kb)} records."
            )
            st.rerun()
        except Exception as error:
            st.error(f"The knowledgebase table could not be read: {error}")
    if notice := st.session_state.pop("agricultural_kb_notice", None):
        st.success(notice)

    with st.expander("Delete database records", expanded=False):
        if kb.empty:
            st.info("The database is empty. Upload a valid CSV or TSV file to initialize it.")
        else:
            tagged_kb = with_record_ids(kb)
            delete_labels = {
                f"{row['organism']} — {row['strain']} — {row['editing_route']} — {row['source_id']}": row["record_id"]
                for _, row in tagged_kb.iterrows()
            }
            selected_delete_labels = st.multiselect(
                "Select records to delete",
                list(delete_labels),
                key="agricultural_kb_delete_selection",
            )
            d1, d2 = st.columns(2)
            if d1.button(
                "Delete selected records",
                disabled=not selected_delete_labels,
                key="agricultural_kb_delete_selected",
            ):
                ids_to_delete = [delete_labels[label] for label in selected_delete_labels]
                st.session_state["agricultural_kb_records"] = delete_records(kb, ids_to_delete)
                st.session_state["agricultural_kb_notice"] = f"Deleted {len(ids_to_delete)} selected record(s)."
                st.rerun()
            confirm_clear = d2.checkbox("Confirm full database clear", key="agricultural_kb_confirm_clear")
            if d2.button(
                "Clear full database",
                disabled=not confirm_clear,
                key="agricultural_kb_clear_all",
            ):
                st.session_state["agricultural_kb_records"] = kb.iloc[0:0].copy()
                st.session_state["agricultural_kb_last_upload"] = None
                st.session_state["agricultural_kb_notice"] = "The full database was cleared. The next valid upload will initialize it."
                st.rerun()

    reviewed_dates = pd.to_datetime(kb.get("last_reviewed", pd.Series(dtype=str)), errors="coerce")
    latest_review = reviewed_dates.max()
    k1, k2, k3, k4 = st.columns(4)
    k1.metric("Evidence records", len(kb), help="One strain–editing-route–source combination per record.")
    k2.metric("Microorganisms", kb["organism"].nunique() if not kb.empty else 0)
    k3.metric("Strains", kb[["organism", "strain"]].drop_duplicates().shape[0] if not kb.empty else 0)
    k4.metric("Last reviewed", latest_review.strftime("%Y-%m-%d") if pd.notna(latest_review) else "Not recorded")

    st.subheader("Search and filter the evidence")
    f1, f2, f3 = st.columns(3)
    search_kb = f1.text_input(
        "Search organism, strain, role or editing route",
        key="agricultural_kb_search",
        placeholder="e.g. Bacillus, biocontrol, Cas9",
    )
    evidence_options = sorted(kb["evidence_level"].dropna().unique()) if not kb.empty else []
    evidence_filter = f2.multiselect("Evidence level", evidence_options, default=evidence_options)
    source_options = sorted(kb["source_type"].dropna().unique()) if not kb.empty else []
    source_filter = f3.multiselect("Source type", source_options, default=source_options)
    f4, f5, f6 = st.columns(3)
    status_options = sorted(kb["biological_status"].dropna().unique()) if not kb.empty else []
    status_filter = f4.multiselect("Agricultural status", status_options, default=status_options)
    route_options = sorted(kb["editing_route"].dropna().unique()) if not kb.empty else []
    route_filter = f5.multiselect("Editing route", route_options, default=route_options)
    pam_filter = f6.selectbox("PAM evidence", ["All", "Reported or tool-defined", "Unknown/not inferred"])

    visible_kb = kb.copy()
    if not visible_kb.empty:
        visible_kb = visible_kb[
            visible_kb["evidence_level"].isin(evidence_filter)
            & visible_kb["source_type"].isin(source_filter)
            & visible_kb["biological_status"].isin(status_filter)
            & visible_kb["editing_route"].isin(route_filter)
        ]
        if search_kb.strip():
            searchable = visible_kb[["organism", "strain", "agricultural_role", "editing_route"]].astype(str).agg(" ".join, axis=1)
            visible_kb = visible_kb[searchable.str.contains(re.escape(search_kb.strip()), case=False, na=False)]
        if pam_filter == "Reported or tool-defined":
            visible_kb = visible_kb[~visible_kb["pam_or_target_requirement"].str.lower().isin(["", "unknown", "not inferred"])]
        elif pam_filter == "Unknown/not inferred":
            visible_kb = visible_kb[visible_kb["pam_or_target_requirement"].str.lower().isin(["", "unknown", "not inferred"])]

    st.dataframe(
        visible_kb,
        hide_index=True,
        width="stretch",
        column_config={
            "organism": st.column_config.TextColumn("Organism", pinned=True),
            "strain": st.column_config.TextColumn("Strain", pinned=True),
            "source_url": st.column_config.LinkColumn("Open evidence", display_text="Open source"),
            "pam_or_target_requirement": st.column_config.TextColumn("PAM / target requirement"),
            "pam_evidence": st.column_config.TextColumn("Meaning of PAM entry", width="large"),
            "limitations": st.column_config.TextColumn("Limitations", width="large"),
            "recommended_validation": st.column_config.TextColumn("Required validation", width="large"),
        },
    )
    st.caption(f"Showing {len(visible_kb)} of {len(kb)} evidence records. Every row retains its source and review date.")
    st.download_button(
        "Download filtered knowledgebase",
        visible_kb.to_csv(sep="\t", index=False).encode("utf-8"),
        "agricultural_microorganism_editing_knowledgebase.tsv",
        mime="text/tab-separated-values",
        icon=":material/download:",
    )

    with st.expander("Evidence interpretation and update policy", expanded=True):
        st.markdown(
            """
            - **Experimentally demonstrated:** the cited strain was edited or regulated using the stated route.
            - **Computationally inferred:** genome evidence supports a hypothesis, but editing has not been demonstrated.
            - **Transferred evidence:** a related strain or taxon supports feasibility; the target strain still requires testing.
            - **Unknown:** no defensible conclusion is available from the reviewed sources.

            The bundled release is a reviewed snapshot, not a universal census. New papers, patents and NCBI records must
            pass curator review before they become evidence records. This prevents automatic searches from being presented
            as verified biological facts. The update date records curation, not publication or database-release dates.
            """
        )

    with st.expander("Target-strain editing protocol planner", expanded=True):
        st.caption(
            "Select a literature-supported strain and route to generate a stage-gated experimental plan. "
            "Exact culture, transformation and selection parameters must come from the cited strain-specific method "
            "and your institutionally approved SOP."
        )
        if kb.empty:
            st.info("No knowledgebase record is available for protocol planning.")
        else:
            protocol_options = kb.apply(
                lambda row: f"{row['organism']} — {row['strain']} — {row['editing_route']}", axis=1
            ).tolist()
            selected_protocol_label = st.selectbox(
                "Target organism and editing route",
                protocol_options,
                key="agricultural_protocol_record",
            )
            selected_protocol = kb.iloc[protocol_options.index(selected_protocol_label)]
            protocol_steps = pd.DataFrame([
                {"Step": 1, "Stage": "Define the biological objective", "Required action": "Specify the intended edit, measurable agricultural phenotype and functions that must be preserved.", "Output / decision gate": "Approved target-product profile and stop criteria."},
                {"Step": 2, "Stage": "Confirm strain identity and genome", "Required action": "Verify strain provenance, assembly quality, taxonomic identity and the exact sequence of the intended locus.", "Output / decision gate": "Versioned genome and confirmed target locus."},
                {"Step": 3, "Stage": "Review biosafety and permissions", "Required action": "Determine containment, institutional approvals, antimicrobial-marker restrictions and whether the organism is beneficial, opportunistic or pathogenic.", "Output / decision gate": "Written authorization and approved laboratory SOP."},
                {"Step": 4, "Stage": "Choose the editing route", "Required action": f"Evaluate the selected route: {selected_protocol['editing_route']}. Check component completeness, PAM/target requirements, repair compatibility and evidence transferability.", "Output / decision gate": "Documented go/no-go rationale and backup route."},
                {"Step": 5, "Stage": "Design and rank candidates", "Required action": "Generate at least two top-ranked candidates and one lower-ranked comparator; evaluate uniqueness, gene consequence, defence barriers, agronomic value and preservation risk.", "Output / decision gate": "Auditable candidate table and Pareto ranking."},
                {"Step": 6, "Stage": "Design controls and readouts", "Required action": "Include parental, delivery-only or no-guide controls, a lower-ranked comparator, a molecular confirmation assay and quantitative phenotype readouts.", "Output / decision gate": "Pre-registered control matrix and acceptance thresholds."},
                {"Step": 7, "Stage": "Prepare the editing construct", "Required action": "Build or obtain the nuclease/regulatory component, guide or targeting element and repair template where required, following the cited method and approved SOP.", "Output / decision gate": "Sequence-verified construct or validated delivery material."},
                {"Step": 8, "Stage": "Establish delivery conditions", "Required action": "Use a validated strain-specific transformation or delivery procedure; first quantify delivery and survival with a non-editing control.", "Output / decision gate": "Acceptable delivery, viability and containment performance."},
                {"Step": 9, "Stage": "Generate and isolate candidates", "Required action": "Perform the approved editing workflow, recover independent candidate clones and retain complete sample provenance.", "Output / decision gate": "Traceable candidate isolates and matched controls."},
                {"Step": 10, "Stage": "Confirm genotype", "Required action": "Screen the intended locus and confirm edit junctions and sequence; check construct loss or persistence as required by the design.", "Output / decision gate": "At least one sequence-confirmed edit and documented negative candidates."},
                {"Step": 11, "Stage": "Measure phenotype and preservation", "Required action": "Quantify the intended reporter or agricultural phenotype, growth/fitness, one preserved nutrient-related function and a relevant non-target interaction.", "Output / decision gate": "Raw measurements linked to genotype and controls."},
                {"Step": 12, "Stage": "Analyse and feed back", "Required action": "Apply the pre-specified statistics, compare predicted and observed rankings, document failures and import the results into the digital strain-to-experiment ledger.", "Output / decision gate": "Reproducible result package and updated evidence status."},
            ])
            st.dataframe(protocol_steps, hide_index=True, width="stretch")
            p1, p2, p3 = st.columns(3)
            p1.metric("PAM / target requirement", str(selected_protocol["pam_or_target_requirement"]))
            p2.metric("Evidence", str(selected_protocol["evidence_level"]))
            p3.metric("Source", str(selected_protocol["source_id"]))
            st.info(f"Record-specific validation recommendation: {selected_protocol['recommended_validation']}")
            st.link_button("Open the supporting source", str(selected_protocol["source_url"]), icon=":material/open_in_new:")
            protocol_markdown = "\n".join([
                f"# AgriPAM-AI editing plan: {selected_protocol['organism']} {selected_protocol['strain']}",
                "",
                f"- Proposed route: {selected_protocol['editing_route']}",
                f"- PAM or target requirement: {selected_protocol['pam_or_target_requirement']}",
                f"- Evidence level: {selected_protocol['evidence_level']}",
                f"- Supporting source: {selected_protocol['source_id']} — {selected_protocol['source_url']}",
                f"- Known limitation: {selected_protocol['limitations']}",
                f"- Record-specific validation: {selected_protocol['recommended_validation']}",
                "",
                "## Stage-gated workflow",
                "",
                *[
                    f"{int(row['Step'])}. **{row['Stage']}** — {row['Required action']}  \n   Gate: {row['Output / decision gate']}"
                    for _, row in protocol_steps.iterrows()
                ],
                "",
                "## Interpretation boundary",
                "",
                "This plan is a decision and documentation framework. Exact culture conditions, reagent quantities, "
                "delivery settings and selection parameters must be taken from the cited strain-specific method and an "
                "institutionally approved SOP. A computationally compatible design is not evidence of successful editing.",
            ])
            st.download_button(
                "Download selected step-by-step plan",
                protocol_markdown.encode("utf-8"),
                file_name="rhizoforge_target_strain_editing_plan.md",
                mime="text/markdown",
                icon=":material/download:",
            )


with tabs[4], st.expander("Editing-system modules and capabilities", expanded=True):
    st.header("Integrative editing-systems atlas")
    st.info("This page expands beyond the demonstrated TTC workflow. Each module is an evidence screen; a detected component is not proof of editing activity.", icon=":material/info:")
    systems = pd.DataFrame([
        {"Module": "CRISPR-Cas and PAMs", "Evidence returned": "Cas genes, arrays, subtype, PAMs and coordinates", "Application": "RNA-guided targeting", "Confidence": "Database/tool-supported", "Validation": "System-specific reporter or target assay"},
        {"Module": "Recombinases and integrases", "Evidence returned": "Candidate genes, locus and neighboring attachment sites", "Application": "Site-specific integration or rearrangement", "Confidence": "Exploratory unless complete locus", "Validation": "Integration-junction assay"},
        {"Module": "Retrons / reverse transcriptases", "Evidence returned": "Reverse-transcriptase and ncRNA neighborhood", "Application": "Template-mediated editing hypothesis", "Confidence": "Exploratory", "Validation": "Targeted editing/readout assay"},
        {"Module": "CRISPR-associated transposases", "Evidence returned": "Cas/transposase components and target-adjacent loci", "Application": "RNA-guided insertion hypothesis", "Confidence": "Exploratory", "Validation": "Insertion mapping assay"},
        {"Module": "Homology-directed repair", "Evidence returned": "RecA/RecFOR and repair-module indicators", "Application": "Template-based repair compatibility", "Confidence": "Context-dependent", "Validation": "Sequence-confirmed repair assay"},
        {"Module": "Non-Type-I-C nucleases", "Evidence returned": "Other nuclease families and guide requirements", "Application": "Alternative RNA-guided targeting", "Confidence": "Database-supported where present", "Validation": "Nuclease activity assay"},
        {"Module": "Delivery and restriction barriers", "Evidence returned": "Restriction-modification, defence and mobility context", "Application": "Construct delivery/stability risk", "Confidence": "Screening evidence", "Validation": "Transformation and stability comparison"},
    ])
    st.dataframe(systems, hide_index=True, width="stretch")
    full_atlas = st.session_state.get("full_workflow")
    if full_atlas:
        consequences = pd.DataFrame(full_atlas["tables"].get("pam_gene_context.tsv", []))
        if not consequences.empty:
            st.subheader("Potentially interrupted genes at supported target sites")
            cols = [c for c in ["contig", "start", "end", "pam", "gene_relation", "selected_locus_tag", "product", "potential_consequence", "growth_relevance_screen", "agronomic_category_screen"] if c in consequences.columns]
            st.dataframe(consequences[cols], hide_index=True, width="stretch")
            st.caption("Gene impacts are computational context predictions and must be confirmed experimentally.")
    else:
        st.warning("Run a genome workflow first to populate system results and site-to-gene consequences.")

with knowledgebase_bank_entry, st.expander("Upload your microbial bank and calculate the SynCom ranking", expanded=True):
    st.write(
        "Enter the tested strains in the Excel template: one row per bacterium or fungus, agricultural functions "
        "scored from 0–5, blank values for measurements not performed, and the available pairwise interaction data. "
        "After upload, the software can assemble a compatible community, identify missing functions and rank which "
        "member is the most appropriate genome-editing chassis. The ranking separates measured evidence from missing "
        "data and applies the project biosafety exclusions. Nothing is transmitted elsewhere; the workbook is analyzed "
        "only in this application session."
    )
    syncom_dir_in = ROOT / "data" / "syncom"
    d1, d2 = st.columns(2)
    for column, filename, label in (
        (d1, "SynCom_input_template.xlsx", "Download the blank template (.xlsx)"),
        (d2, "SynCom_example_filled.xlsx", "Download a filled synthetic example (.xlsx)"),
    ):
        file_path = syncom_dir_in / filename
        if file_path.exists():
            column.download_button(label, file_path.read_bytes(), filename,
                                   mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                                   key=f"dl_{filename}", width="stretch")
    uploaded_bank = st.file_uploader("Upload your completed workbook", type=["xlsx"], key="syncom_upload")
    use_example = st.checkbox("Try the synthetic example instead (invented strains S1–S6, not real data)",
                              key="syncom_use_example")
    bank_source = uploaded_bank if uploaded_bank is not None else (
        io.BytesIO((syncom_dir_in / "SynCom_example_filled.xlsx").read_bytes())
        if use_example and (syncom_dir_in / "SynCom_example_filled.xlsx").exists() else None)
    if bank_source is None:
        st.info("Upload a workbook, or tick the box above to try the example.")
    else:
        try:
            from agripam import syncom_io
            sheets_in = syncom_io.read_workbook(bank_source)
            user_bank, bank_issues = syncom_io.build_bank(sheets_in)
        except ImportError:
            user_bank, bank_issues, sheets_in = None, [], {}
            st.error("Reading Excel files needs the openpyxl package: run `pip install openpyxl` and restart the app.")
        except Exception as error:  # unreadable or corrupt workbook
            user_bank, bank_issues, sheets_in = None, [], {}
            st.error(f"The workbook could not be read: {error}")
        for issue in bank_issues:
            (st.error if issue["level"] == "error" else st.warning)(f"{issue['where']}: {issue['message']}")
        if user_bank is not None:
            strains_in = list(user_bank["traits"])
            fungi_seen = sorted({fn for per in user_bank["fungi"].values() for fn in per}
                                | {s for s, k in user_bank.get("kind", {}).items() if k == "fungus"})
            sheet_members, sheet_fungi = syncom_io.community_from_sheet(sheets_in, user_bank)
            st.success(f"Workbook read: {len(strains_in)} strains, {len(user_bank['functions'])} functions, "
                       f"{len(fungi_seen)} fungal partner(s) with compatibility data.")
            s1, s2, s3 = st.columns(3)
            partner_fungi = s1.multiselect("Fungal anchor(s) of the community", fungi_seen,
                                           default=[f for f in sheet_fungi if f in fungi_seen], key="syncom_fungi")
            mode = s2.radio("Community", ["Assemble one for me", "I choose the members"],
                            index=1 if sheet_members else 0, key="syncom_mode")
            threshold_in = s3.slider("A function counts as delivered at score ≥", 1, 5, 2, key="syncom_threshold")
            chosen_members, size_in = [], 4
            if mode == "I choose the members":
                chosen_members = st.multiselect("Community members", strains_in,
                                                default=[m for m in sheet_members if m in strains_in], key="syncom_members")
            else:
                size_in = st.slider("Maximum community size", 2, 8, 4, key="syncom_size")
            if st.button("Run SynCom analysis", type="primary", key="syncom_run"):
                if mode == "I choose the members" and not chosen_members:
                    st.warning("Choose at least one member, or switch to “Assemble one for me”.")
                else:
                    st.session_state["syncom_result"] = (
                        syncom_io.analyze(user_bank, chosen_members or None, partner_fungi, size_in, threshold_in),
                        bank_issues,
                    )
            if "syncom_result" in st.session_state:
                result_in, issues_in = st.session_state["syncom_result"]
                st.subheader("Result")
                if not result_in["members"]:
                    st.warning("No eligible community could be assembled. Check fungal compatibility, biosafety flags "
                               "and the delivery threshold.")
                else:
                    r1, r2, r3 = st.columns(3)
                    r1.metric("Community", ", ".join(result_in["members"]))
                    top_row = next((r for r in result_in["ranking"] if r["rank"] == 1), None)
                    r2.metric("Top-ranked chassis", top_row["strain"] if top_row else "none eligible")
                    r3.metric("Functions missing", len(result_in["missing"]))
                    for warning_text in result_in.get("warnings", []):
                        st.warning(warning_text)
                    excluded = [r["strain"] for r in result_in["ranking"] if not r["safety_eligible"]]
                    if excluded:
                        st.warning("Not ranked as chassis because of a biosafety flag: " + ", ".join(excluded))
                    st.markdown("**Chassis ranking**")
                    user_ranking = pd.DataFrame(result_in["ranking"])
                    st.dataframe(user_ranking, hide_index=True, width="stretch")
                    st.markdown("**Function coverage (best score in the community)**")
                    st.dataframe(pd.DataFrame({"function": list(result_in["coverage"]),
                                               "best score": list(result_in["coverage"].values())}),
                                 hide_index=True, width="stretch")
                    if result_in["gaps"]:
                        st.markdown("**Functional gaps and possible donors**")
                        st.dataframe(pd.DataFrame(result_in["gaps"]), hide_index=True, width="stretch")
                    st.markdown("**Send a bank decision to genome design**")
                    eligible_user_chassis = [str(row["strain"]) for row in result_in["ranking"] if bool(row.get("safety_eligible", False))]
                    if eligible_user_chassis:
                        chosen_user_chassis = st.selectbox("Selected chassis from this ranking", eligible_user_chassis, key="user_bank_selected_chassis")
                        gap_names = [str(row.get("missing_function", "")) for row in result_in["gaps"] if row.get("missing_function")]
                        chosen_user_objective = st.selectbox(
                            "Agricultural objective",
                            gap_names + ["Add another agricultural benefit — specify", "Preserve or regulate an existing function"],
                            key="user_bank_selected_objective",
                        )
                        if chosen_user_objective == "Add another agricultural benefit — specify":
                            chosen_user_objective = st.text_input("Benefit to add", key="user_bank_custom_objective", placeholder="e.g. beta-glucanase or phosphate solubilization").strip()
                        protected_user_functions = st.multiselect(
                            "Functions that must not be disrupted", list(result_in["coverage"]),
                            default=[name for name, score in result_in["coverage"].items() if float(score) >= threshold_in],
                            key="user_bank_protected_functions",
                        )
                        if st.button("Use this SynCom decision in Genome evaluation", type="primary", key="user_bank_handoff"):
                            if not chosen_user_objective:
                                st.warning("Define the agricultural objective before continuing.")
                            else:
                                st.session_state["design_handoff"] = {
                                    "source": "uploaded microbial bank", "community": list(result_in["members"]),
                                    "chassis": chosen_user_chassis, "objective": chosen_user_objective,
                                    "protected_functions": protected_user_functions,
                                }
                                st.success("Design saved. Open Genome evaluation and provide the exact selected-isolate genome.")
                    st.caption("Top choice under 2,000 random weightings (seed 42): "
                               + ", ".join(f"{k} {v:.0%}" for k, v in result_in["sensitivity"].items()))
                    st.download_button("Download the results (.xlsx)", syncom_io.results_to_xlsx(result_in, issues_in),
                                       "SynCom_chassis_ranking.xlsx",
                                       mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                                       key="syncom_dl_results")
                st.warning("Decision support only. Scores are not measured editing efficiency, safety clearance or "
                           "field performance; compatibility rests on the tests you entered.")

with tabs[2], st.expander("Community-aware chassis selection: which SynCom member should be edited?", expanded=True):
    st.write(
        "Starting from strains whose agricultural effects were already tested in the laboratory, this module asks "
        "which member of a compatible community is the best genome-editing chassis. It favours a member whose own "
        "functions the rest of the community already covers, so the edit can add a function without removing one."
    )
    syncom_dir = ROOT / "results" / "syncom_selector"
    if not (syncom_dir / "summary.json").exists():
        st.info("Run `python scripts/run_syncom_selector.py` to generate the SynCom selector results.")
    else:
        syncom_summary = json.loads((syncom_dir / "summary.json").read_text())
        scenario = st.selectbox("Community", list(syncom_summary), key="syncom_scenario",
                                format_func=lambda key: key.replace("_", " "))
        info = syncom_summary[scenario]
        c1, c2, c3 = st.columns(3)
        c1.metric("Bacterial members", ", ".join(info["bacterial_members"]))
        c2.metric("Top-ranked chassis", info["top_chassis"] or "none eligible")
        if info.get("top_bacterial_chassis") and info["top_bacterial_chassis"] != info["top_chassis"]:
            st.caption(f"Top bacterium: {info['top_bacterial_chassis']}. A fungal anchor can tie in the ranking; editing a fungus needs its own biosafety and precedent review.")
        c3.metric("Functions missing", len(info["functions_missing"]))
        st.caption("Fungal partner used for compatibility: " + ", ".join(info["fungal_members"]))
        if info["chassis_excluded_for_safety"]:
            st.warning(
                "Excluded from chassis ranking by the project's biosafety flags: "
                + ", ".join(info["chassis_excluded_for_safety"])
                + ". Their inclusion in any field community also needs biosafety review."
            )
        st.subheader("Chassis ranking")
        ranking = pd.read_csv(syncom_dir / f"{scenario}_chassis_ranking.tsv", sep="\t")
        st.dataframe(ranking, hide_index=True, width="stretch")
        gaps_path = syncom_dir / f"{scenario}_gap_analysis.tsv"
        st.subheader("Functional gaps")
        scenario_gaps = pd.DataFrame()
        if gaps_path.exists():
            scenario_gaps = pd.read_csv(gaps_path, sep="\t")
            st.dataframe(scenario_gaps, hide_index=True, width="stretch")
        else:
            st.success("No function tested in the bank is missing from this community.")
        st.caption(
            "Top-choice frequency under 2,000 random weightings (seed 42): "
            + ", ".join(f"{k} {v:.0%}" for k, v in info["top_choice_frequency_under_random_weights"].items())
        )
        st.warning(
            "Decision support only. Compatibility rests on few plate tests (see fit_evidence_n); editing precedent is "
            "a genus/species lookup, and genome-based editability is not yet computed. Genome accessions in "
            "data/syncom/strain_metadata.tsv are same-species references, not the isolates."
        )
        st.markdown("### Send this SynCom decision to genome design")
        safety_values = ranking["safety_eligible"].astype(str).str.lower().isin(["true", "1", "yes"])
        eligible_scenario_chassis = ranking.loc[safety_values, "strain"].astype(str).tolist()
        if eligible_scenario_chassis:
            default_chassis = info.get("top_chassis") if info.get("top_chassis") in eligible_scenario_chassis else eligible_scenario_chassis[0]
            chosen_scenario_chassis = st.selectbox("Selected chassis", eligible_scenario_chassis, index=eligible_scenario_chassis.index(default_chassis), key="scenario_selected_chassis")
            gap_column = "missing_function" if "missing_function" in scenario_gaps.columns else None
            scenario_gap_names = scenario_gaps[gap_column].dropna().astype(str).tolist() if gap_column else list(info.get("functions_missing", []))
            chosen_scenario_objective = st.selectbox(
                "Agricultural objective",
                scenario_gap_names + ["Add another agricultural benefit — specify", "Preserve or regulate an existing function"],
                key="scenario_selected_objective",
            )
            if chosen_scenario_objective == "Add another agricultural benefit — specify":
                chosen_scenario_objective = st.text_input("Benefit to add", key="scenario_custom_objective", placeholder="e.g. beta-glucanase or phosphate solubilization").strip()
            protected_options = sorted(set(info.get("functions_covered", [])))
            protected_scenario_functions = st.multiselect("Functions that must not be disrupted", protected_options, default=protected_options, key="scenario_protected_functions")
            if st.button("Use this decision in Genome evaluation", type="primary", key="scenario_handoff"):
                if not chosen_scenario_objective:
                    st.warning("Define the agricultural objective before continuing.")
                else:
                    st.session_state["design_handoff"] = {
                        "source": scenario.replace("_", " "),
                        "community": list(info.get("bacterial_members", [])) + list(info.get("fungal_members", [])),
                        "chassis": chosen_scenario_chassis, "objective": chosen_scenario_objective,
                        "protected_functions": protected_scenario_functions,
                    }
                    st.success("Design saved. Open Genome evaluation and provide the exact selected-isolate genome.")

with tabs[2], st.expander("Native-bank evidence and candidate selection", expanded=True):
    st.header("Laboratory bank: candidate SynCom (worked example)")
    st.write(
        "This page evaluates candidate members of a native synthetic microbial community designed to protect crops "
        "and support yield. Selection is multi-objective: pathogen antagonism, nutrient or plant-growth support, "
        "compatibility with other beneficial isolates, biosafety, and the feasibility of strain-specific engineering."
    )
    st.info(
        "CRISPR-Cas evidence is one enabling trait—not the definition of a good SynCom member. A strong organism may "
        "remain valuable without a native CRISPR system because it can be used unmodified or evaluated with an "
        "introduced editing tool. Conversely, a native system does not compensate for poor compatibility or safety.",
        icon=":material/ecg_heart:",
    )
    if not native_cohort.empty:
        st.subheader("Native-bank phenotypes and proposed SynCom roles")
        syncom_view = native_cohort.copy()
        syncom_view["Proposed SynCom contribution"] = syncom_view.apply(
            lambda row: (
                "Antagonism + community compatibility"
                if str(row.get("pathogen_score", "Not measured")) != "Not measured"
                and str(row.get("beneficial_compatibility", "Not measured")) != "Not measured"
                else "Evidence incomplete — phenotype or compatibility measurement required"
            ), axis=1)
        syncom_columns = [column for column in ["rank", "strain", "identification", "pathogen_score",
                          "beneficial_compatibility", "native_discovery_priority", "study_role",
                          "Proposed SynCom contribution", "required_next_step"] if column in syncom_view.columns]
        st.dataframe(
            syncom_view[syncom_columns], hide_index=True, width="stretch",
            column_config={
                "strain": st.column_config.TextColumn("Native-bank isolate", pinned=True),
                "identification": st.column_config.TextColumn("Provisional identification", pinned=True),
                "pathogen_score": st.column_config.TextColumn("Pathogen-antagonism evidence"),
                "beneficial_compatibility": st.column_config.TextColumn("Beneficial compatibility"),
                "native_discovery_priority": st.column_config.TextColumn("Native-system priority"),
                "required_next_step": st.column_config.TextColumn("Evidence needed next", width="large"),
            },
        )
        st.caption(
            "The phenotype scores are experimental inputs supplied for the native bank. Public genomes from related "
            "strains provide transferred evidence only; they do not replace sequencing or testing the submitted isolate."
        )

    st.subheader("Reference discovery across the candidate collection")
    st.subheader("Corrected public-reference screen: current sequencing candidates")
    st.caption("Results are grouped by submitted isolate and updated as checkpointed accession waves finish. Positive reference evidence prioritizes sequencing; it does not establish the submitted isolate's genotype.")
    if corrected_bank_results.empty:
        st.warning("The prioritized reference analysis has not started.")
    else:
        priority_view = corrected_bank_results.copy()
        priority_view["pipeline_validation"] = priority_view.get(
            "parser_revision", pd.Series("", index=priority_view.index)).map(
                lambda value: "Validated parser v2" if value == "confirmed_cas_operons_tab_v2"
                else "WITHDRAWN — rerun required")
        valid_priority = priority_view["pipeline_validation"] == "Validated parser v2"
        priority_view["has_complete_native_locus"] = (
            pd.to_numeric(priority_view.get("complete_interference_operons", 0), errors="coerce").fillna(0) > 0
        ) & valid_priority
        priority_view["interpretation"] = priority_view.apply(
            lambda row: ("Candidate native locus; proceed to locus review and spacer/mobile matching"
                         if int(row.get("predicted_cas_operons", 0) or 0) > 0 else
                         "No qualifying native Cas operon detected in this public reference"), axis=1)
        priority_columns = [column for column in ["isolate", "reference_species", "canonical_assembly",
                            "assembly_level", "crispr_arrays", "trusted_arrays", "predicted_cas_operons",
                            "complete_interference_operons", "predicted_subtypes", "analysis_status",
                            "pipeline_validation", "interpretation"] if column in priority_view.columns]
        p1, p2, p3 = st.columns(3)
        p1.metric("Validated priority assemblies", int(valid_priority.sum()))
        p2.metric("Validated native Cas candidates", int(((pd.to_numeric(priority_view.get("predicted_cas_operons", 0), errors="coerce").fillna(0) > 0) & valid_priority).sum()))
        p3.metric("Priority assemblies requiring rerun", int((~valid_priority).sum()))
        if (~valid_priority).any():
            st.error("The earlier priority negatives are withdrawn: the batch parser read CRISPRCasTyper's rejected-candidate table instead of its confirmed-system table. Rerun with parser v2 is required.")
        summary_rows = []
        for isolate, group in priority_view[valid_priority].groupby("isolate", sort=False):
            positive = group[group["has_complete_native_locus"]]
            subtypes = sorted({str(value) for value in positive.get("predicted_subtypes", []) if str(value) and str(value) != "Not detected"})
            summary_rows.append({
                "Target isolate": isolate,
                "Reference species": ", ".join(sorted(set(group["reference_species"].astype(str)))),
                "Public references screened": len(group),
                "References with complete native loci": len(positive),
                "Detected subtypes": ", ".join(subtypes) if subtypes else "None detected",
                "What this permits us to say": (
                    "Strong transferred reference evidence; sequence the submitted isolate before claiming presence"
                    if len(positive) else
                    "No locus detected in the screened references; this does not establish absence from the submitted isolate"
                ),
            })
        if summary_rows:
            st.markdown("#### Result by submitted isolate")
            st.dataframe(pd.DataFrame(summary_rows), hide_index=True, width="stretch")
            b39_positive = priority_view[(priority_view["isolate"] == "B39") & priority_view["has_complete_native_locus"]]
            if not b39_positive.empty:
                subtype_counts = b39_positive["predicted_subtypes"].value_counts().to_dict()
                subtype_text = ", ".join(f"{count} Type {subtype}" for subtype, count in subtype_counts.items())
                st.success(
                    f"Current discovery lead (updated from analyzed evidence): {len(b39_positive)} of "
                    f"{int((priority_view['isolate'] == 'B39').sum())} screened "
                    f"*E. mendocina* references contain a complete computationally detected native locus "
                    f"({subtype_text}). This prioritizes sequencing B39; it is not evidence that B39 itself carries the locus."
                )
                st.warning(
                    "Novelty status: genomic occurrence is already database-detectable. A defensible publication claim would require "
                    "a systematic literature/patent review plus strain-specific confirmation, expression/activity evidence, PAM inference, "
                    "and experimental repurposing. Do not describe this as a first discovery yet."
                )
        st.markdown("#### Assembly-level evidence")
        st.dataframe(priority_view[priority_columns], hide_index=True, width="stretch")
        st.info("FI20 has no complete/chromosome *T. yunnanense* assembly in NCBI. Analyze the sequenced FI20 genome for fungal RNAi and repair machinery plus compatibility with an introduced editor.")
        st.download_button("Download priority reference results", priority_view.to_csv(sep="\t", index=False).encode("utf-8"),
                           "priority_chassis_crispr_cas_results.tsv", mime="text/tab-separated-values")
    st.subheader("FI20 fungal editing-readiness reference")
    st.info(
        "**Fungal native-machinery legend:** for FI20 (*T. yunnanense*), native editing readiness does not mean a bacterial CRISPR-Cas locus. It refers to endogenous RNAi machinery (such as Argonaute, Dicer and RdRP), homologous-recombination and non-homologous-end-joining repair capacity, and other features that may support or constrain an introduced editor. Any editor, target requirement and edit outcome must still be selected and validated in the FI20 isolate.",
        icon=":material/biotech:",
    )
    if fungal_reference_screen.empty:
        st.warning("No exact-species fungal reference has been screened.")
    else:
        fungal_view = fungal_reference_screen.copy()
        fungal_view["interpretation"] = fungal_view.apply(
            lambda row: (f"{int(row['annotation_hits']):,} sequence motifs in the draft reference; guide uniqueness and activity remain untested"
                         if str(row["evidence_category"]).startswith("Introduced ") else
                         (f"{int(row['annotation_hits'])} annotation matches; confirm protein domains and locus integrity"
                          if int(row["annotation_hits"]) else
                          "Not identified by annotation keywords; this is not evidence of biological absence")), axis=1)
        st.dataframe(fungal_view[["submitted_isolate", "reference_species", "accession", "evidence_category",
                                  "annotation_hits", "interpretation", "evidence_level", "boundary"]],
                     hide_index=True, width="stretch")
        st.caption("Fungal RNAi and DNA-repair systems are not native CRISPR-Cas. FI20 must be sequenced and analyzed directly before selecting an introduced editor or making a novelty claim.")
        st.download_button("Download FI20 fungal reference screen", fungal_view.to_csv(sep="\t", index=False).encode("utf-8"),
                           "FI20_fungal_reference_screen.tsv", mime="text/tab-separated-values")
    st.info(
        "The accession inventory covers every current exact-species assembly returned by NCBI for the identified "
        "microbes in the isolate collection. Sequence-level analysis is checkpointed and may still be pending. "
        "The smaller table below is the completed one-reference-per-species pilot, not the final all-accession result.",
        icon=":material/database:",
    )
    if not all_ncbi_inventory.empty:
        inv1, inv2, inv3, inv4 = st.columns(4)
        inv1.metric("Target species", all_ncbi_inventory["reference_species"].nunique())
        inv2.metric("NCBI accession records", len(all_ncbi_inventory))
        inv3.metric("Unique paired assemblies", all_ncbi_inventory["canonical_assembly"].nunique())
        completed_accessions = int(((all_accession_results.get("analysis_status", pd.Series(dtype=str)) == "Analyzed") &
                                    (all_accession_results.get("pipeline_validation", pd.Series(dtype=str)) == "Validated parser v2")).sum())
        inv4.metric("Assemblies fully analyzed", completed_accessions)
        st.progress(completed_accessions / max(all_ncbi_inventory["canonical_assembly"].nunique(), 1),
                    text=f"Sequence-level all-subtype survey: {completed_accessions:,} completed")
        st.caption(
            f"NCBI inventory updated {all_ncbi_inventory_manifest.get('created_utc', 'date unavailable')}. "
            "GCA/GCF pairs are retained as accession records but counted once as a biological assembly."
        )
        st.subheader("Live all-accession CRISPR-Cas results")
        if all_accession_results.empty:
            st.warning("No accession has completed sequence-level analysis yet.")
        else:
            visible_result_columns = [column for column in [
                "related_submitted_isolates", "reference_species", "canonical_assembly",
                "assembly_level", "analysis_status", "crispr_arrays", "trusted_arrays",
                "predicted_cas_operons", "complete_interference_operons", "predicted_subtypes",
                "mobile_element_status", "pam_status", "runtime_seconds", "analysis_utc", "error",
            ] if column in all_accession_results.columns]
            st.dataframe(all_accession_results[visible_result_columns], hide_index=True, width="stretch")
            st.caption(
                "These are completed all-subtype screens. 'Not detected' means no qualifying Cas operon was found "
                "in that accession; it does not mean the submitted laboratory isolate is negative."
            )
            st.download_button("Download completed all-subtype results",
                               all_accession_results.to_csv(sep="\t", index=False).encode("utf-8"),
                               "all_accession_crispr_cas_results.tsv", mime="text/tab-separated-values")
        with st.expander("All NCBI accessions: scope, status and per-species counts", expanded=False):
            status_by_species = (all_ncbi_inventory.groupby("reference_species", as_index=False)
                                 .agg(accession_records=("accession", "size"),
                                      unique_assemblies=("canonical_assembly", "nunique")))
            complete_by_species = (all_accession_results.loc[all_accession_results.get("analysis_status", "") == "Analyzed"]
                                   .groupby("reference_species")["canonical_assembly"].nunique()
                                   if not all_accession_results.empty else pd.Series(dtype=int))
            status_by_species["sequence_analyses_complete"] = status_by_species["reference_species"].map(complete_by_species).fillna(0).astype(int)
            st.dataframe(status_by_species, hide_index=True, width="stretch")
            st.download_button("Download complete NCBI accession inventory",
                               all_ncbi_inventory.to_csv(sep="\t", index=False).encode("utf-8"),
                               "all_ncbi_accession_inventory.tsv", mime="text/tab-separated-values")
        st.warning(
            "Inventory is not discovery evidence. A publishable native-system candidate requires a coherent Cas locus "
            "and compatible array, subtype assignment, spacer evidence against an independently sourced mobile-element "
            "sequence, PAM inference, novelty checking, and strain-specific laboratory validation."
        )
    else:
        st.warning("The all-accession NCBI inventory has not been generated yet.")

    st.subheader("Completed pilot: one selected reference per species")
    if public_native_survey.empty:
        st.warning("The public reference survey has not been generated yet.")
    else:
        survey = public_native_survey.copy()
        isolate_map = native_cohort.groupby("identification")["strain"].apply(lambda values: ", ".join(values)).to_dict()
        survey.insert(0, "submitted_isolates", survey["reference_species"].map(isolate_map).fillna(""))
        survey["ncbi_assembly"] = survey["accession"].map(lambda accession: f"https://www.ncbi.nlm.nih.gov/datasets/genome/{accession}/" if accession else "")
        survey["discovery_interpretation"] = survey.apply(
            lambda row: ("Candidate native Cas system in this public reference; prioritize isolate sequencing"
                         if int(row.get("predicted_cas_operons", 0) or 0) > 0 else
                         "Orphan CRISPR array(s) in this public reference; Cas machinery not detected"
                         if int(row.get("trusted_arrays", 0) or 0) > 0 else
                         "No system detected in this selected reference; other strains may differ"), axis=1,
        )
        s1, s2, s3 = st.columns(3)
        s1.metric("Species references analyzed", len(survey))
        s2.metric("References with trusted arrays", int((pd.to_numeric(survey["trusted_arrays"], errors="coerce").fillna(0) > 0).sum()))
        s3.metric("References with predicted Cas operons", int((pd.to_numeric(survey["predicted_cas_operons"], errors="coerce").fillna(0) > 0).sum()))
        st.dataframe(
            survey,
            hide_index=True,
            width="stretch",
            column_config={
                "submitted_isolates": st.column_config.TextColumn("Related submitted isolate(s)", pinned=True),
                "reference_species": st.column_config.TextColumn("Reference species", pinned=True),
                "ncbi_assembly": st.column_config.LinkColumn("NCBI assembly", display_text="Open NCBI"),
                "discovery_interpretation": st.column_config.TextColumn("Permitted discovery interpretation", width="large"),
                "isolate_interpretation": st.column_config.TextColumn("Boundary", width="large"),
            },
        )
        st.download_button(
            "Download complete public native-system survey",
            survey.to_csv(sep="\t", index=False).encode("utf-8"),
            "public_native_system_survey.tsv",
            mime="text/tab-separated-values",
        )
        st.warning(
            "A CRISPR array alone is not an editing system. A complete Cas locus, array–locus compatibility, PAM "
            "evidence, expression and functional validation are still required. A negative reference result is not "
            "a negative result for the corresponding submitted isolate."
        )

    st.divider()
    st.subheader("Which native-bank isolates should be sequenced and validated?")
    st.caption(
        "This decision table joins the bank phenotype/compatibility measurements with current public-reference "
        "CRISPR-Cas evidence. It updates as corrected accession analyses are added. Reference evidence never proves "
        "that the corresponding submitted isolate carries the same system."
    )
    sequencing_rows = []
    if not native_cohort.empty:
        validation_by_isolate = (native_validation.set_index("strain").to_dict("index")
                                 if not native_validation.empty and "strain" in native_validation.columns else {})
        valid_results = corrected_bank_results[
            corrected_bank_results.get("parser_revision", pd.Series("", index=corrected_bank_results.index)).eq(
                "confirmed_cas_operons_tab_v2")
        ].copy() if not corrected_bank_results.empty else pd.DataFrame()
        for _, isolate_row in native_cohort.iterrows():
            isolate = str(isolate_row["strain"])
            organism = str(isolate_row["identification"])
            matches = valid_results[valid_results["isolate"].eq(isolate)] if not valid_results.empty else pd.DataFrame()
            screened = int(matches["canonical_assembly"].nunique()) if not matches.empty else 0
            if not matches.empty:
                complete_values = pd.to_numeric(matches["complete_interference_operons"], errors="coerce").fillna(0)
                trusted_values = pd.to_numeric(matches["trusted_arrays"], errors="coerce").fillna(0)
                positives = matches[(complete_values > 0) & (trusted_values > 0)]
            else:
                positives = pd.DataFrame()
            positive_count = len(positives)
            subtypes = sorted(set(positives["predicted_subtypes"].astype(str))) if positive_count else []
            pam_evidence = "Pending spacer/mobile matching" if positive_count else "Not inferred"
            evidence = (f"{positive_count}/{screened} complete references; Type {', Type '.join(subtypes)}"
                        if positive_count else
                        (f"0/{screened} complete references" if screened else "Public-reference screen pending"))
            recommendation = "Screen public references before selecting for native-system validation"
            if isolate_row.get("native_discovery_priority") == "Hold":
                recommendation = "HOLD: resolve identity and biosafety before sequencing for engineering"
            elif isolate == "B39" and positive_count:
                recommendation = "SEQUENCE NOW: leading strain-specific native-system validation candidate"
            elif positive_count:
                recommendation = "SEQUENCE HIGH PRIORITY: confirm the native locus in the submitted isolate and integrate SynCom phenotype evidence"
            elif isolate == "B34":
                evidence = "Reference evidence: Type I-C 7/24; Type III-B 5/24"
                subtypes = ["I-C", "III-B"]
                pam_evidence = "TTC inferred for Type I-C; Type III-B not assigned"
                recommendation = "Resolve compatibility before considering sequencing for SynCom validation"
            elif screened and not positive_count:
                recommendation = "Analyze remaining drafts; sequence isolate if SynCom phenotype justifies discovery effort"
            validation_record = validation_by_isolate.get(isolate, {})
            validation_status = validation_record.get("overall_validation", "NOT STARTED")
            if validation_status == "NOT STARTED" and (positive_count or isolate == "B34"):
                validation_status = "REFERENCE EVIDENCE ONLY"
            sequencing_rows.append({
                "Bank isolate": isolate,
                "Provisional organism": organism,
                "Antagonism input": isolate_row.get("pathogen_score", "Not measured"),
                "Beneficial compatibility": isolate_row.get("beneficial_compatibility", "Not measured"),
                "Native-system reference evidence": evidence,
                "Native system type(s)": ", ".join(f"Type {value}" for value in subtypes) if subtypes else "Not detected / pending",
                "PAM evidence": pam_evidence,
                "Submitted-isolate validation": validation_status,
                "Validation evidence record": validation_record.get("evidence_reference", "None supplied"),
                "Decision": recommendation,
            })
    sequencing_decisions = pd.DataFrame(sequencing_rows)
    if sequencing_decisions.empty:
        st.info("The strain-level sequencing decision table will appear when the native-bank registry is available.")
    else:
        decision_order = {"SEQUENCE NOW": 0, "Analyze": 1, "Resolve": 2, "Screen": 3, "HOLD": 4}
        sequencing_decisions["_order"] = sequencing_decisions["Decision"].map(
            lambda value: next((rank for prefix, rank in decision_order.items() if str(value).startswith(prefix)), 3)
        )
        sequencing_decisions = sequencing_decisions.sort_values(["_order", "Bank isolate"]).drop(columns="_order")
        evidence_filter = st.selectbox(
            "Show strains",
            ["Native-system candidates", "All native-bank strains", "Experimentally validated"],
            help="A candidate has a complete native system in at least one related public reference. Validated means the submitted isolate itself has passed sequence and laboratory confirmation.",
        )
        if evidence_filter == "Native-system candidates":
            displayed_decisions = sequencing_decisions[
                sequencing_decisions["Native-system reference evidence"].str.contains("Type I", regex=False)
                | sequencing_decisions["Native-system reference evidence"].str.contains("complete references", regex=False)
            ]
            displayed_decisions = displayed_decisions[
                ~displayed_decisions["Native-system reference evidence"].str.startswith("0/")
            ]
        elif evidence_filter == "Experimentally validated":
            displayed_decisions = sequencing_decisions[
                sequencing_decisions["Submitted-isolate validation"].str.startswith("VALIDATED")
            ]
        else:
            displayed_decisions = sequencing_decisions
        st.dataframe(displayed_decisions, hide_index=True, width="stretch",
                     column_config={
                         "Decision": st.column_config.TextColumn("Sequencing/validation decision", width="large"),
                         "Submitted-isolate validation": st.column_config.TextColumn("Validation status", width="large"),
                     })
        if displayed_decisions.empty:
            st.info("No submitted isolate has reached this evidence level yet.")
        st.caption(
            "A strain enters the validated view only after its persistent evidence record is marked VALIDATED. "
            "Required evidence: isolate genome and locus confirmation, expression/activity, PAM or targeting "
            "requirement, editing outcome, sequence confirmation and phenotype/safety assessment."
        )
        st.download_button("Download strain sequencing decisions",
                           sequencing_decisions.to_csv(sep="\t", index=False).encode("utf-8"),
                           "native_bank_sequencing_validation_decisions.tsv", mime="text/tab-separated-values")
    st.info(
        "A strain can still be edited with an introduced system. For example, the Agricultural editing knowledgebase "
        "documents experimental SpCas9 editing in B. subtilis and B. velezensis strains. Its NGG requirement belongs "
        "to the introduced SpCas9 tool, not to a native Bacillus CRISPR system."
    )

with tabs[2], st.expander("Reference 64-triplet PAM model"):
    st.subheader("Reference 64-triplet PAM model")
    st.caption(
        "Trained from the comparative P. polymyxa Type I-C evidence. It is a reference "
        "prioritization model, not a strain-specific result for the uploaded genome."
    )
    pam = st.text_input("Candidate 3-nt PAM", "TTC", max_chars=3).upper()
    mobile = st.toggle("Target lies near mobile-element features", value=True)
    try:
        result = score_pam(pam, mobile, score_rows)
        p1, p2, p3 = st.columns(3)
        p1.metric("Model probability", f"{result['percent']:.1f}%")
        p2.metric("Rank in context", f"{result['rank']} / 64")
        p3.metric("Mobile-context feature", "Yes" if mobile else "No")
        st.info(
            "This score ranks sequence contexts; it is not a measured cleavage "
            "or editing efficiency."
        )
    except (ValueError, KeyError) as error:
        st.error(str(error))

    context_value = 1 if mobile else 0
    ranking = pd.DataFrame(score_rows)
    ranking = ranking[ranking["mobile_nearby"].astype(int) == context_value].copy()
    ranking["predicted_percent"] = ranking["predicted_percent"].astype(float)
    ranking = ranking.sort_values("rank_within_context").head(10)
    st.dataframe(
        ranking[["rank_within_context", "pam3", "predicted_percent"]],
        hide_index=True,
        width="stretch",
    )

with tabs[2], st.expander("Editor-specific target designer: TTC, NGG and TTTV", expanded=True):
    st.subheader("Editor-specific target designer")
    st.write(
        "Choose the editing route first. AgriPAM-AI then applies that system's PAM orientation and guide length "
        "to both DNA strands. TTC remains a reference hypothesis for native Type I-C; NGG and TTTV belong to "
        "introduced editors and do not imply that the organism carries those systems naturally."
    )
    designer_routes = {
        "Native Type I-C candidate (reference TTC hypothesis)": {
            "editor": "Native Type I-C candidate", "pam_pattern": "TTC", "pam_side": "5prime",
            "protospacer_length": 35, "action": "Cascade recognition followed by Cas3 interference",
            "evidence": "Reference-supported hypothesis; strain-specific activity requires validation",
        },
        "SpCas9 (introduced)": {**EDITOR_PRESETS["SpCas9"], "editor": "SpCas9", "action": "DNA cleavage", "evidence": "Introduced-editor targetability"},
        "dCas9 / CRISPRi (introduced)": {**EDITOR_PRESETS["dCas9/CRISPRi"], "editor": "dCas9/CRISPRi", "action": "Transcriptional repression without DNA cleavage", "evidence": "Introduced-editor targetability"},
        "Cas12a / Cpf1 (introduced)": {**EDITOR_PRESETS["Cas12a/Cpf1"], "editor": "Cas12a/Cpf1", "action": "Staggered DNA cleavage", "evidence": "Introduced-editor targetability"},
    }
    selected_route_name = st.selectbox("Editing route", list(designer_routes), key="knowledgebase_designer_route")
    selected_route = designer_routes[selected_route_name]
    r1, r2, r3 = st.columns(3)
    r1.metric("PAM rule", selected_route["pam_pattern"])
    r2.metric("PAM position", "5′ of target" if selected_route["pam_side"] == "5prime" else "3′ of target")
    r3.metric("Default guide length", f"{selected_route['protospacer_length']} nt")
    st.caption(f"Action: {selected_route['action']}. Evidence: {selected_route['evidence']}.")
    demo = "GCGATATTC" + "AACTAAATAAACAACAAAGGACTCCATACTGGTT" + "GCGATCGATC"
    uploaded = st.file_uploader("Upload a FASTA file", type=["fa", "fasta", "fna"], key="knowledgebase_designer_fasta")
    pasted = st.text_area("Or paste a DNA sequence", demo, height=130, key="knowledgebase_designer_sequence")
    length = st.slider("Guide or protospacer length", min_value=18, max_value=40,
                       value=int(selected_route["protospacer_length"]), key="knowledgebase_designer_length")

    raw = uploaded.getvalue().decode() if uploaded else pasted
    if st.button("Find targets for the selected route", type="primary", key="knowledgebase_find_editor_targets"):
        try:
            records = parse_fasta(raw if raw.lstrip().startswith(">") else f">pasted_sequence\n{raw}")
            candidates = scan_editor_targets(records, selected_route["editor"], selected_route["pam_pattern"], selected_route["pam_side"], length)
            table = pd.DataFrame(candidates)
            st.metric(f"{selected_route['pam_pattern']}-compatible candidates", len(table))
            if table.empty:
                st.info(f"No {selected_route['pam_pattern']}-compatible candidates were found at this guide length. Try a different route, a regulatory region, another documented nuclease or a validated neutral insertion site.")
            else:
                table.insert(0, "candidate_id", [f"AGP{i:04d}" for i in range(1, len(table)+1)])
                table["evidence_scope"] = selected_route["evidence"]
                table["model_status"] = "PAM-compatible; gene consequence, delivery, repair and off-target validation required"
                st.dataframe(table, hide_index=True, width="stretch")
                st.download_button(
                    "Download candidate TSV",
                    table.to_csv(sep="\t", index=False),
                    file_name="agripam_editor_specific_candidates.tsv",
                    mime="text/tab-separated-values",
                )
                report = {
                    "tool": "AgriPAM-AI",
                    "route": selected_route_name,
                    "pam_pattern": selected_route["pam_pattern"],
                    "pam_side": selected_route["pam_side"],
                    "action": selected_route["action"],
                    "evidence_scope": selected_route["evidence"],
                    "functional_validation": "pending",
                    "candidates": table.to_dict(orient="records"),
                }
                st.download_button(
                    "Download machine-readable JSON",
                    json.dumps(report, indent=2),
                    file_name="agripam_editor_specific_candidates.json",
                    mime="application/json",
                )
        except ValueError as error:
            st.error(str(error))

with tabs[2], st.expander("Auditable model and control evidence"):
    st.subheader("Auditable reference-study evidence")
    st.caption("Precomputed comparative P. polymyxa evidence; independent of the currently uploaded genome.")
    st.write("Context-level Type I-C evidence")
    st.dataframe(contexts, hide_index=True, width="stretch")
    st.write("Grouped cross-validation")
    st.dataframe(models, hide_index=True, width="stretch")
    st.write("Matched controls and Wilson intervals")
    st.dataframe(summary, hide_index=True, width="stretch")

with tabs[2], st.expander("Native-bank defence and editing-system atlas"):
    st.subheader("Native-bank isolates with native-system reference evidence")
    st.caption(
        "Only corrected complete-locus results are included. The counts describe related public genomes, not the "
        "submitted isolate. An isolate is not labelled validated until its own sequence and "
        "laboratory evidence are recorded."
    )
    atlas_rows = []
    if not corrected_bank_results.empty:
        corrected = corrected_bank_results[
            corrected_bank_results.get("parser_revision", pd.Series("", index=corrected_bank_results.index)).eq(
                "confirmed_cas_operons_tab_v2")
        ].copy()
        corrected["complete_interference_operons"] = pd.to_numeric(
            corrected["complete_interference_operons"], errors="coerce").fillna(0)
        corrected["trusted_arrays"] = pd.to_numeric(corrected["trusted_arrays"], errors="coerce").fillna(0)
        corrected = corrected[(corrected["complete_interference_operons"] > 0) & (corrected["trusted_arrays"] > 0)]
        for (isolate, organism), group in corrected.groupby(["isolate", "reference_species"]):
            validation_record = (native_validation[native_validation["strain"].eq(isolate)].iloc[0].to_dict()
                                 if not native_validation.empty and native_validation["strain"].eq(isolate).any() else {})
            atlas_rows.append({
                "Bank isolate": isolate,
                "Provisional organism": organism,
                "Native system type(s)": ", ".join(f"Type {value}" for value in sorted(set(group["predicted_subtypes"].astype(str)))),
                "Positive public-reference evidence": f"{int(group['canonical_assembly'].nunique())} complete genome(s)",
                "Trusted arrays": int(pd.to_numeric(group["trusted_arrays"], errors="coerce").fillna(0).sum()),
                "PAM status": "Not yet inferred; spacer/mobile matching required",
                "Submitted-isolate status": validation_record.get("overall_validation", "REFERENCE EVIDENCE ONLY"),
                "Next decision": "Prioritize isolate sequencing, locus confirmation and expression/activity testing",
            })
    # The established P. polymyxa panel remains one bank-wide candidate signal, not a privileged benchmark.
    if not native_cohort.empty and native_cohort["strain"].eq("B34").any():
        b34_validation = (native_validation[native_validation["strain"].eq("B34")].iloc[0].to_dict()
                          if not native_validation.empty and native_validation["strain"].eq("B34").any() else {})
        atlas_rows.append({
            "Bank isolate": "B34",
            "Provisional organism": "Paenibacillus polymyxa",
            "Native system type(s)": "Type I-C; Type III-B",
            "Positive public-reference evidence": "Type I-C: 7/24; Type III-B: 5/24",
            "Trusted arrays": "Panel evidence",
            "PAM status": "TTC inferred for Type I-C; Type III-B not assigned",
            "Submitted-isolate status": b34_validation.get("overall_validation", "REFERENCE EVIDENCE ONLY"),
            "Next decision": "Resolve poor measured community compatibility before prioritizing SynCom validation",
        })
    native_atlas = pd.DataFrame(atlas_rows)
    if not native_atlas.empty:
        # Keep exported/rendered columns Arrow-compatible while preserving explicit evidence wording.
        native_atlas = native_atlas.astype(str)
    if native_atlas.empty:
        st.info("No corrected complete native-system reference evidence is currently available for the bank.")
    else:
        st.dataframe(native_atlas, hide_index=True, width="stretch",
                     column_config={"Next decision": st.column_config.TextColumn(width="large")})
        st.download_button("Download native-bank system atlas",
                           native_atlas.to_csv(sep="\t", index=False).encode("utf-8"),
                           "native_bank_defence_editing_system_atlas.tsv", mime="text/tab-separated-values")
    st.warning(
        "REFERENCE EVIDENCE ONLY is not a native-system detection in the submitted isolate. The VALIDATED label "
        "requires isolate sequence, complete locus and array confirmation, activity/expression, targeting requirement "
        "or PAM evidence, editing outcome, and phenotype/safety confirmation."
    )

with tabs[4], st.expander("Reproducibility contract", expanded=True):
    st.subheader("Reproducibility contract")
    st.markdown(
        """
        - Exact assembly accessions and versions are retained.
        - Fixed random seeds are recorded in `config/project.yaml`.
        - Close homologues are grouped during model cross-validation.
        - Raw predictions are separated from curated biological conclusions.
        - All inferred PAMs remain predictions until experimental validation.
        - Tables used by this application are downloadable, plain-text TSV files.
        """
    )
    st.code("python -m unittest discover -s tests -v", language="bash")

with tabs[2], st.expander("Type I-C TTC benchmark genome scan"):
    st.subheader("Type I-C TTC benchmark genome scan")
    st.write(
        "Upload a genome FASTA or enter a versioned NCBI assembly accession. This quick demonstration "
        "AgriPAM-AI validates the sequence, summarizes the assembly, finds "
        "5′-TTC-compatible Type I-C targets on both strands, applies sequence "
        "quality checks and exports an auditable candidate table."
    )
    st.warning(
        "This applies the P. polymyxa reference TTC hypothesis and does not establish a strain-specific PAM. "
        "It ranks candidate sequences but does not run gene "
        "annotation, DefenseFinder, approximate off-target alignment or wet-lab validation."
    )

    with st.form("genome_analysis_form"):
        source = st.selectbox(
            "Genome source",
            ["Upload FASTA", "NCBI assembly accession", "Built-in demonstration"],
        )
        uploaded_genome = st.file_uploader(
            "Genome FASTA",
            type=["fa", "fasta", "fna"],
            disabled=source != "Upload FASTA",
        )
        accession = st.text_input(
            "Versioned NCBI assembly accession",
            placeholder="GCF_000597985.1",
            disabled=source != "NCBI assembly accession",
        )
        protospacer_length = st.slider("Protospacer length", 30, 40, 35)
        mobile_assumption = st.toggle(
            "Apply the mobile-associated context score",
            value=False,
            help="Use only when the target region has independent mobile-element evidence.",
        )
        maximum_rows = st.number_input(
            "Maximum candidates displayed",
            min_value=10,
            max_value=5000,
            value=250,
            step=10,
        )
        analyze = st.form_submit_button(
            "Analyze genome",
            type="primary",
            icon=":material/biotech:",
        )

    if analyze:
        try:
            with st.status("Preparing genome", expanded=True) as status:
                if source == "Upload FASTA":
                    if uploaded_genome is None:
                        raise ValueError("Choose a FASTA file before starting the analysis")
                    fasta_text = uploaded_genome.getvalue().decode("utf-8")
                    source_label = uploaded_genome.name
                elif source == "NCBI assembly accession":
                    clean_accession = validate_assembly_accession(accession)
                    st.write(f"Downloading {clean_accession} from NCBI Datasets")
                    fasta_text = download_genome_fasta(clean_accession, ROOT)
                    source_label = clean_accession
                else:
                    fasta_text = ">demo_contig\n" + ("ACGT" * 100) + "TTC" + ("GCGT" * 9) + ("ACGT" * 100)
                    source_label = "AgriPAM-AI demonstration sequence"
                records = parse_fasta(fasta_text)
                st.write(f"Validated {len(records):,} FASTA record(s)")
                genome_stats, genome_candidates = analyze_ttc_genome(
                    records,
                    protospacer_length,
                    score_rows,
                    mobile_assumption,
                )
                status.update(label="Genome analysis complete", state="complete", expanded=False)
            st.session_state["genome_analysis"] = {
                "source": source_label,
                "summary": genome_stats,
                "candidates": genome_candidates,
                "protospacer_length": protospacer_length,
                "records": records,
            }
        except (UnicodeDecodeError, ValueError, RuntimeError, TimeoutError) as error:
            st.error(str(error))

    analysis = st.session_state.get("genome_analysis")
    if analysis:
        genome_stats = analysis["summary"]
        genome_candidates = analysis["candidates"]
        st.caption(f"Current result: {analysis['source']}")
        c1, c2, c3, c4 = st.columns(4)
        c1.metric("Genome size", f"{int(genome_stats['total_length_bp']):,} bp")
        c2.metric("FASTA records", f"{int(genome_stats['contigs']):,}")
        c3.metric("GC content", f"{float(genome_stats['gc_percent']):.2f}%")
        c4.metric("TTC candidates", f"{int(genome_stats['ttc_candidates']):,}")
        st.caption(
            f"Ambiguous bases: {int(genome_stats['ambiguous_bases']):,} • "
            f"Sequence SHA-256: {genome_stats['sha256']}"
        )
        candidate_table = pd.DataFrame(genome_candidates)
        if candidate_table.empty:
            st.info("No unambiguous 5′-TTC target was found at the selected length.")
        else:
            shown = candidate_table.head(int(maximum_rows))
            st.dataframe(
                shown,
                hide_index=True,
                width="stretch",
                column_config={
                    "design_score": st.column_config.ProgressColumn(
                        "Design score", min_value=0, max_value=100, format="%.1f"
                    ),
                    "gc_percent": st.column_config.NumberColumn("GC (%)", format="%.2f"),
                    "pam_model_percent": st.column_config.NumberColumn("PAM model (%)", format="%.2f"),
                },
            )
            st.caption(
                f"Showing {len(shown):,} of {len(candidate_table):,} candidates. "
                "Exact genome copies are reported; approximate off-target search remains required."
            )
            # Interactive genomic target explorer
            if analysis.get("records"):
                st.markdown("### Genomic Target Explorer")
                st.caption(
                    "Select a ranked candidate to inspect its real genomic coordinates, "
                    "local target region, strand, PAM and protospacer sequence."
                )

                candidate_ids = shown["candidate_id"].astype(str).tolist()
                selected_id = st.selectbox(
                    "Candidate to visualize",
                    candidate_ids,
                    key="genome_target_explorer_candidate",
                )

                selected_candidate = next(
                    row
                    for row in genome_candidates
                    if str(row["candidate_id"]) == selected_id
                )

                explorer_html = build_target_explorer(
                    records=analysis["records"],
                    candidates=genome_candidates,
                    selected_candidate=selected_candidate,
                )

                st.components.v1.html(
                    explorer_html,
                    height=520,
                    scrolling=False,
                )

            result_report = {
                "tool": "AgriPAM-AI",
                "product": "AgriPAM-AI",
                "analysis": "Type I-C 5'-TTC quick genome scan",
                "source": analysis["source"],
                "protospacer_length": analysis["protospacer_length"],
                "interpretation": "prioritization score; not measured editing efficiency",
                "summary": genome_stats,
                "candidates": genome_candidates,
            }
            manifest = {
                "product": "AgriPAM-AI",
                "engine": "AgriPAM-AI",
                "run_timestamp_utc": datetime.now(timezone.utc).isoformat(),
                "source": analysis["source"],
                "parameters": {
                    "protospacer_length": analysis["protospacer_length"],
                    "pam": "5'-TTC",
                    "mobile_context_assumption": int(genome_stats["mobile_context_assumption"]),
                },
                "environment": {
                    "python": platform.python_version(),
                    "platform": platform.platform(),
                    "streamlit": getattr(st, "__version__", "unknown"),
                },
                "outputs": {
                    "summary": genome_stats,
                    "candidate_count": len(genome_candidates),
                    "tsv_filename": "agripam_genome_candidates.tsv",
                    "json_filename": "agripam_genome_analysis.json",
                },
                "interpretation": "Design-prioritization evidence only; verify candidate activity experimentally.",
            }
            with st.container(horizontal=True):
                st.download_button(
                    "Download all candidates (TSV)",
                    candidate_table.to_csv(sep="\t", index=False),
                    file_name="agripam_genome_candidates.tsv",
                    mime="text/tab-separated-values",
                    icon=":material/download:",
                )
                st.download_button(
                    "Download analysis report (JSON)",
                    json.dumps(result_report, indent=2),
                    file_name="agripam_genome_analysis.json",
                    mime="application/json",
                    icon=":material/data_object:",
                )
                st.download_button(
                    "Download run manifest",
                    json.dumps(manifest, indent=2),
                    file_name="rhizoforge_select_run_manifest.json",
                    mime="application/json",
                    icon=":material/description:",
                )

with tabs[1], st.container(border=True):
    st.header("Multi-objective design portfolio")
    st.write(
        "When supported edit sites are available, this section compares alternatives across four objectives. "
        "The values support decisions; they are not probabilities of successful editing."
    )
    st.dataframe(pd.DataFrame([
        {"Objective": "Editability", "Evidence": "PAM support, guide uniqueness and repair compatibility"},
        {"Objective": "Deliverability", "Evidence": "Restriction–modification and defence-system barriers"},
        {"Objective": "Agronomic value", "Evidence": "Antifungal, nutrient-processing and colonization relevance"},
        {"Objective": "Preservation and safety", "Evidence": "Growth flags, beneficial functions, repeated targets and coding context"},
    ]), hide_index=True, width="stretch")
    portfolio_run = st.session_state.get("full_workflow")
    if not portfolio_run:
        st.info("Run the Genome workspace first. A portfolio requires a completed strain analysis.")
    else:
        portfolio = pd.DataFrame(portfolio_run["tables"].get("multiobjective_design_portfolio.tsv", []))
        systems_screen = pd.DataFrame(portfolio_run["tables"].get("editing_systems_screen.tsv", []))
        system_assessments = pd.DataFrame(portfolio_run["tables"].get("editing_system_assessments.tsv", []))
        if not system_assessments.empty:
            st.subheader("Editing-system completeness and readiness assessment")
            st.caption(
                "System-level interpretation of required components, locus organization, integrity flags, targeting "
                "requirements and delivery barriers. Expression and activity cannot be concluded from a genome alone."
            )
            st.dataframe(system_assessments, hide_index=True, width="stretch")
            st.download_button(
                "Download editing-system assessment",
                system_assessments.to_csv(sep="\t", index=False).encode(),
                "editing_system_assessments.tsv", mime="text/tab-separated-values",
            )
        else:
            st.info(
                "This saved result predates the system-level assessment. Rerun the genome workflow to calculate "
                "component completeness, locus organization, integrity flags, targeting requirements and delivery barriers."
            )
        if not systems_screen.empty:
            with st.expander("Raw candidate-component matches — technical detail"):
                st.caption(
                    "Genes or proteins whose annotations match editing-related terms. This diagnostic table supports "
                    "the system-level assessment; individual rows are not complete or active editing systems."
                )
                st.dataframe(systems_screen, hide_index=True, width="stretch")
        if portfolio.empty:
            st.warning("No defensible CRISPR portfolio was generated. This usually means that no strain-specific PAM passed the evidence threshold. It is insufficient evidence—not proof that the strain cannot be edited.")
        else:
            st.subheader("Pareto-ranked portfolio")
            highlighted = portfolio[portfolio["portfolio_role"] != "retained candidate"]
            st.dataframe(highlighted, hide_index=True, width="stretch")
            with st.expander("Show all retained and rejected candidates"):
                st.dataframe(portfolio, hide_index=True, width="stretch")
            st.download_button("Download complete design portfolio", portfolio.to_csv(sep="\t", index=False).encode(), "multiobjective_design_portfolio.tsv", mime="text/tab-separated-values")

        st.subheader("Digital strain-to-experiment loop")
        st.markdown("Genome or reads → chassis identity → compatible editing routes → PAM evidence → delivery barriers → function preservation → candidate ranking → experimental controls → measured results → model update")
        template = pd.DataFrame(portfolio_run["tables"].get("experimental_results_template.tsv", []))
        if not template.empty:
            st.download_button("Download experimental-results template", template.to_csv(sep="\t", index=False).encode(), "experimental_results_template.tsv", mime="text/tab-separated-values")
        measured = st.file_uploader("Upload completed experimental-results table", type=["tsv", "csv"], help="Use the exported template. Results are joined to predictions by candidate_id.")
        if measured is not None:
            observed = pd.read_csv(measured, sep="\t" if measured.name.lower().endswith(".tsv") else ",")
            if "candidate_id" not in observed.columns:
                st.error("The results table must contain candidate_id.")
            elif portfolio.empty:
                st.error("No current design portfolio is available for matching these results.")
            else:
                joined = portfolio.merge(observed, on="candidate_id", how="left", suffixes=("_prediction", "_experiment"))
                st.success("Predictions and experimental results linked by candidate ID.")
                st.dataframe(joined, hide_index=True, width="stretch")
                st.download_button("Download prediction–experiment ledger", joined.to_csv(sep="\t", index=False).encode(), "prediction_experiment_ledger.tsv", mime="text/tab-separated-values")
        st.caption("This release records experimental feedback for prospective model updates. Automatic retraining is intentionally withheld until measurements pass quality control.")

with tabs[5]:
    st.header("Independent retrospective validation")
    st.info(
        "This page tests frozen predictions against deposited outcomes that were not used to generate them. "
        "Registered does not mean completed: metrics appear only after a real normalized outcome table is supplied.",
        icon=":material/science:",
    )
    if not external_validation.empty:
        with st.expander("Dataset provenance and eligibility", expanded=True):
            st.dataframe(external_validation, hide_index=True, width="stretch")
            primary = external_validation[external_validation["dataset_id"] == "YU2024_GSE196911"]
            if not primary.empty:
                st.warning(str(primary.iloc[0]["local_status"]))

    if gse196911_manifest:
        st.subheader("Deposited held-out outcome package")
        e1, e2, e3 = st.columns(3)
        e1.metric("Deposited guides", int(gse196911_manifest["n_guides"]))
        e2.metric("Target genes", int(gse196911_manifest["n_genes"]))
        e3.metric("Predictions", "Frozen" if gse196911_validation_metrics else "Withheld")
        if gse196911_validation_metrics:
            observed = gse196911_validation_metrics["analysis"]
            control = gse196911_validation_metrics["negative_control"]
            st.success("The released Yu et al. forest was run on the outcome-blind input, predictions were frozen, and only then were outcomes joined.")
            v1, v2, v3 = st.columns(3)
            v1.metric("Pooled Spearman ρ", f"{observed['spearman_rho']:.3f}")
            v2.metric("Median within-gene ρ", f"{observed['median_within_gene_spearman']:.3f}")
            v3.metric("Permutation-control ρ", f"{control['spearman_rho']:.3f}")
            st.caption("Primary endpoint: negative OD1 log2 fold-change (higher = stronger depletion). Nine genes, one E. coli genome; this is not a held-out-genome test.")
            bootstrap = gse196911_validation_metrics.get("cluster_bootstrap", {})
            null = gse196911_validation_metrics.get("permutation_null", {})
            if bootstrap.get("status") == "complete":
                pooled_ci = bootstrap["pooled_spearman"]
                within_ci = bootstrap["median_within_cluster_spearman"]
                with st.expander("Gene-cluster uncertainty", expanded=True):
                    st.write(
                        f"2,000 gene-cluster bootstrap resamples: pooled ρ 95% interval "
                        f"{pooled_ci['ci95_low']:.3f}–{pooled_ci['ci95_high']:.3f}; "
                        f"median within-gene ρ 95% interval {within_ci['ci95_low']:.3f}–{within_ci['ci95_high']:.3f}."
                    )

            if null.get("status") == "complete":
                pooled_null = null["pooled_spearman"]
                within_null = null["median_within_cluster_spearman"]
                st.write(
                    f"Within-gene permutation null (2,000 resamples): empirical two-sided p = "
                    f"{pooled_null['empirical_two_sided_p']:.4g} for pooled ρ and "
                    f"{within_null['empirical_two_sided_p']:.4g} for median within-gene ρ."
                )
                st.caption("The minimum attainable empirical p-value is 1/2,001; this is a computational randomization test, not evidence of cross-genome generalization.")
            st.dataframe(pd.DataFrame(observed["per_gene_spearman"]), hide_index=True, width="stretch")
            influence = gse196911_validation_metrics.get("leave_one_group_out", {})
            meta = gse196911_validation_metrics.get("random_effects_meta_analysis", {})
            if influence.get("status") == "complete" and meta.get("status") == "complete":
                with st.expander("Robustness and heterogeneity", expanded=True):
                    r1, r2, r3 = st.columns(3)
                    r1.metric("Leave-one-gene-out minimum ρ", f"{influence['minimum_leave_one_out_spearman']:.3f}")
                    r2.metric("Random-effects pooled ρ", f"{meta['pooled_spearman']:.3f}")
                    r3.metric("Between-gene I²", f"{meta['i_squared_percent']:.1f}%")
                    st.write(
                        f"Removing any one gene leaves pooled ρ between "
                        f"{influence['minimum_leave_one_out_spearman']:.3f} and "
                        f"{influence['maximum_leave_one_out_spearman']:.3f}. The exploratory random-effects "
                        f"estimate is {meta['pooled_spearman']:.3f} (95% CI {meta['ci95_low']:.3f}–{meta['ci95_high']:.3f})."
                    )
                    st.dataframe(pd.DataFrame(influence["rows"]), hide_index=True, width="stretch")
                    st.caption("I² shows substantial between-gene heterogeneity; the model is useful on average but not uniformly successful. The Fisher-z meta-analysis is explicitly exploratory.")

            with st.expander("Ranking utility", expanded=True):
                st.dataframe(pd.DataFrame(observed["continuous_top_k"]), hide_index=True, width="stretch")
                top_k_test = gse196911_validation_metrics.get("top_k_permutation", {})
                if top_k_test.get("status") == "complete":
                    st.markdown("**Top-k randomization test**")
                    st.dataframe(pd.DataFrame(top_k_test["rows"]), hide_index=True, width="stretch")
                    st.caption("Outcomes were permuted within genes, preserving gene-level distributions. This tests prioritization utility without inventing an activity threshold.")
                st.caption("ROC/PR and calibration are intentionally absent: no binary activity threshold was pre-specified, and the forest emits a continuous score rather than a probability.")
                sensitivity = pd.DataFrame(gse196911_validation_metrics.get("feature_group_sensitivity", []))
                if not sensitivity.empty:
                    st.subheader("Feature-group reliance sensitivity")
                    st.dataframe(sensitivity, hide_index=True, width="stretch")
                    st.caption(
                        "Each feature group was jointly permuted with a fixed seed before outcomes were joined. "
                        "The model was not retrained, so this measures reliance under perturbation—not a retrained ablation or causal feature effect."
                    )
                st.download_button("Download reproduced validation table", gse196911_validation_table.to_csv(sep="\t", index=False).encode(),
                                   "GSE196911_frozen_predictions_and_outcomes.tsv", mime="text/tab-separated-values")
                st.download_button("Download reproduced metrics", json.dumps(gse196911_validation_metrics, indent=2, allow_nan=True),
                                   "GSE196911_validation_metrics.json", mime="application/json")

        else:
            st.success("GSE196911 outcomes were downloaded, normalized and checksummed. No performance is shown because independent predictions have not yet been generated.")
        st.caption(f"Source SHA-256: {gse196911_manifest['source_sha256']} • normalized SHA-256: {gse196911_manifest['normalized_sha256']}")
        if not gse196911_outcomes.empty:
            gene_summary = (gse196911_outcomes.groupby("gene_id", as_index=False)
                            .agg(guides=("guide_id", "size"), median_depletion=("outcome", "median")))
            st.dataframe(gene_summary, hide_index=True, width="stretch")
            st.bar_chart(gene_summary, x="gene_id", y="median_depletion")
            st.download_button("Download normalized held-out outcomes", gse196911_outcomes.to_csv(sep="\t", index=False).encode(),
                               "GSE196911_heldout_outcomes.tsv", mime="text/tab-separated-values")
        if gse196911_mapping_manifest:
            st.success(
                f"Sequence mapping complete: {gse196911_mapping_manifest['matched_guides']} deposited outcomes matched "
                f"{gse196911_mapping_manifest['unique_sequences']} unique supplementary guide sequences. "
                "The prediction-input file contains no logFC or p-value columns."
            )
        if not published_validation_summary.empty:
            with st.expander("Published external-validation reference results", expanded=True):
                st.warning("These headline values are transcribed from Supplementary Table S14. They use the authors' aggregate estimand (median excluding purE/purK across time points), not the primary OD1/all-nine-gene estimand shown above.")
                st.dataframe(published_validation_summary, hide_index=True, width="stretch")
    if crisprhal_metrics:
        st.subheader("Predefined bacterial Cas9 released holdout")
        st.write(
            "A transparent sequence model was fitted on the released crisprHAL TevSpCas9 training table only. "
            "The 5,049 test sequences were separated from their activities, predictions were frozen and "
            "checksummed, and the outcomes were joined afterward."
        )
        primary = crisprhal_metrics["primary"]
        baseline = crisprhal_metrics["gc_baseline"]
        uncertainty = crisprhal_metrics["uncertainty_and_null"]
        c1, c2, c3, c4 = st.columns(4)
        c1.metric("Released test guides", f"{crisprhal_metrics['n_testing']:,}")
        c2.metric("Sequence-model ρ", f"{primary['spearman_rho']:.3f}")
        c3.metric("GC-only baseline ρ", f"{baseline['spearman_rho']:.3f}")
        c4.metric("Δρ over GC baseline", f"{crisprhal_metrics['delta_spearman_over_gc_baseline']:.3f}")
        st.caption(
            f"Guide-bootstrap 95% interval {uncertainty['bootstrap_ci95_low']:.3f}–"
            f"{uncertainty['bootstrap_ci95_high']:.3f}; permutation p="
            f"{uncertainty['permutation_two_sided_p']:.4g}; exact train/test sequence overlap: "
            f"{crisprhal_metrics['exact_sequence_overlap_train_test']}."
        )
        similarity = crisprhal_metrics.get("sequence_similarity_audit", {})
        if similarity:
            st.markdown("**Sequence-similarity leakage audit**")
            s1, s2, s3 = st.columns(3)
            s1.metric("Near-identical guides (≥18/20)", similarity["test_guides_with_guide_identity_at_least_18_of_20"])
            s2.metric("Distant test subset", f"{similarity['distant_subset_n']:,}")
            s3.metric("Distant-subset ρ", f"{similarity['distant_subset_metrics']['spearman_rho']:.3f}")
            st.caption(
                "For every test guide, maximum positional identity to all 20,195 training guides was computed. "
                "After excluding every test guide with a nearest training guide above 17/20 identity, the "
                "association is unchanged; this reduces concern that a few near duplicates drive performance."
            )
        st.dataframe(pd.DataFrame(crisprhal_metrics["top_k_permutation"]["rows"]), hide_index=True, width="stretch")
        st.warning(
            "This is a predefined guide-level holdout from one study and one chassis—not external-study or "
            "held-out-genome validation. It tests sequence-learning and outcome separation, not the complete "
            "AgriPAM chassis model."
        )
        with st.expander("Model, provenance and limitations"):
            st.write(crisprhal_metrics["benchmark_role"])
            st.write(f"Ridge alpha selected by training-only five-fold CV: {crisprhal_metrics['ridge_alpha_selected_by_training_only_5fold_cv']:.6g}")
            for limitation in crisprhal_metrics["limitations"]:
                st.markdown(f"- {limitation}")
        if not crisprhal_table.empty:
            st.download_button("Download crisprHAL released-holdout table", crisprhal_table.to_csv(sep="\t", index=False).encode(),
                               "crisprHAL_Tev_released_holdout.tsv", mime="text/tab-separated-values")
            st.download_button("Download crisprHAL released-holdout metrics", json.dumps(crisprhal_metrics, indent=2),
                               "crisprHAL_Tev_released_holdout_metrics.json", mime="application/json")
    if hawkins_metrics:
        st.subheader("Cross-chassis transfer benchmark")
        st.write(
            "Hawkins et al. measured the same 33 fully matched GFP guides in both *E. coli* and "
            "*B. subtilis*. The mean *E. coli* activity was frozen as a sequence-matched transfer "
            "score before the separately stored *B. subtilis* outcomes were joined."
        )
        h1, h2, h3 = st.columns(3)
        h1.metric("Matched guides", int(hawkins_metrics["n_guides"]))
        h2.metric("Cross-chassis Spearman ρ", f"{hawkins_metrics['spearman_rho']:.3f}")
        h3.metric("Permutation p", f"{hawkins_metrics['permutation']['empirical_two_sided_p']:.4f}")
        ci = hawkins_metrics["bootstrap"]
        reliability = hawkins_metrics["replicate_reliability"]
        st.caption(
            f"Guide-bootstrap 95% interval: {ci['ci95_low']:.3f}–{ci['ci95_high']:.3f}. "
            f"Replicate reliability: E. coli ρ={reliability['ecoli_spearman']:.3f}; "
            f"B. subtilis ρ={reliability['bsubtilis_spearman']:.3f}."
        )
        st.dataframe(pd.DataFrame(hawkins_metrics["continuous_top_k"]), hide_index=True, width="stretch")
        st.warning(
            "This is evidence that guide rank transfers only partially between two chassis for one GFP reporter. "
            "The E. coli measurement is an experimental transfer baseline, not an AgriPAM-trained prediction, "
            "and the result is not genome-wide validation or evidence for B26."
        )
        with st.expander("Provenance, separation and limitations"):
            st.markdown(
                "**Source:** Hawkins et al. (2020), DOI 10.1016/j.cels.2020.09.009; authors' public "
                "`mismatch_crispri/gfpdata` release.  \n"
                "**Separation:** E. coli replicate means are written and checksummed separately from "
                "the held-out B. subtilis replicate means; joining occurs by guide ID afterward.  \n"
                "**Classification metrics:** omitted because no binary biological threshold or calibrated "
                "probability was pre-specified."
            )
            for limitation in hawkins_metrics["limitations"]:
                st.markdown(f"- {limitation}")
        if not hawkins_table.empty:
            st.download_button("Download cross-chassis validation table", hawkins_table.to_csv(sep="\t", index=False).encode(),
                               "Hawkins2020_cross_chassis_validation.tsv", mime="text/tab-separated-values")
            st.download_button("Download cross-chassis metrics", json.dumps(hawkins_metrics, indent=2),
                               "Hawkins2020_cross_chassis_metrics.json", mime="application/json")
    if challenge_manifest:
        with st.expander("Frozen unseen-genome challenge", expanded=True):
            st.write(
                "A deposited genome that was not used to fit the reference PAM model is analyzed without any "
                "experimental outcome. The challenge tests whether the platform reports a useful design package "
                "and abstains when the biological scope does not support a strain-specific claim."
            )
            u1, u2, u3 = st.columns(3)
            u1.metric("Genome", challenge_manifest["accession"])
            u2.metric("Reference TTC candidates", f"{challenge_manifest['summary']['ttc_candidates']:,}")
            u3.metric("Decision", "ABSTAIN")
            st.warning(f"{challenge_manifest['transfer_decision']}: {challenge_manifest['reason']}")
            st.dataframe(pd.DataFrame(challenge_manifest["qc"]), hide_index=True, width="stretch")
            st.markdown("**Uncertainty**")
            for item in challenge_manifest["uncertainty"]:
                st.markdown(f"- {item}")
            st.markdown("**Next physical verification**")
            for item in challenge_manifest["next_experiment"]:
                st.markdown(f"- {item}")
            st.caption(
                f"Frozen {challenge_manifest['frozen_at_utc']} • input SHA-256 "
                f"{challenge_manifest['input_fasta_sha256']} • outcomes available to analysis: no"
            )
            st.download_button("Download challenge manifest", json.dumps(challenge_manifest, indent=2),
                               "unseen_genome_challenge_manifest.json", mime="application/json")
            if not challenge_candidates.empty:
                st.download_button("Download top frozen candidates", challenge_candidates.to_csv(sep="\t", index=False).encode(),
                                   "unseen_genome_top_candidates.tsv", mime="text/tab-separated-values")

    with st.expander("Leakage-control contract", expanded=True):
        st.markdown(
            "1. Register and checksum the deposited release.  \n"
            "2. Hold out complete genomes/studies—not random guides alone.  \n"
            "3. Fit preprocessing and models on training groups only.  \n"
            "4. Freeze held-out predictions before joining experimental outcomes.  \n"
            "5. Apply the same split to negative controls and feature ablations."
        )
        st.caption("Full protocol: docs/EXTERNAL_VALIDATION.md")

    with st.expander("Evaluate a frozen held-out table", expanded=False):
        st.write(
            "Upload a normalized TSV/CSV with dataset_id, organism, genome_id, guide_id, prediction and outcome. "
            "Optional columns named prediction_ablation_* are compared on exactly the same rows."
        )
        validation_upload = st.file_uploader("Frozen predictions + subsequently joined outcomes", type=["tsv", "csv"], key="external_validation_upload")
        use_threshold = st.toggle("A biological active/inactive threshold was pre-specified", value=False)
        threshold = st.number_input("Activity threshold", value=0.5, disabled=not use_threshold)
        if validation_upload is not None:
            try:
                held_out = pd.read_csv(validation_upload, sep="\t" if validation_upload.name.lower().endswith(".tsv") else ",")
                config = ValidationConfig(positive_threshold=float(threshold) if use_threshold else None)
                metrics = evaluate_predictions(held_out, config)
                control = permuted_outcome_control(held_out, config)
                m1, m2, m3 = st.columns(3)
                m1.metric("Held-out guides", int(metrics["n_guides"]))
                m2.metric("Held-out genomes", int(metrics["n_genomes"]))
                m3.metric("Spearman ρ", f"{metrics['spearman_rho']:.3f}")
                if "median_within_gene_spearman" in metrics:
                    st.metric("Median within-gene Spearman ρ", f"{metrics['median_within_gene_spearman']:.3f}")
                    st.dataframe(pd.DataFrame(metrics["per_gene_spearman"]), hide_index=True, width="stretch")
                if use_threshold:
                    c1, c2, c3 = st.columns(3)
                    c1.metric("ROC AUC", f"{metrics['roc_auc']:.3f}")
                    c2.metric("Average precision", f"{metrics['average_precision']:.3f}")
                    c3.metric("Positive prevalence", f"{metrics['positive_prevalence']:.1%}")
                    st.dataframe(pd.DataFrame(metrics["top_k_enrichment"]), hide_index=True, width="stretch")
                    if "calibration" in metrics:
                        st.subheader("Calibration")
                        st.dataframe(pd.DataFrame(metrics["calibration"]), hide_index=True, width="stretch")
                        st.caption(f"Brier score: {metrics['brier_score']:.4f}")
                    else:
                        st.info(metrics.get("calibration_note", "Calibration not available."))
                else:
                    st.info(metrics["classification_note"])
                st.subheader("Negative control")
                st.write(f"Outcome-permutation Spearman ρ (fixed seed 42): {control['spearman_rho']:.3f}")
                ablations = [column for column in held_out.columns if column.startswith("prediction_ablation_")]
                if ablations:
                    st.subheader("Feature ablation")
                    st.dataframe(pd.DataFrame(evaluate_ablations(held_out, ["prediction", *ablations], config)), hide_index=True, width="stretch")
                report = {"dataset_provenance": validation_upload.name, "metrics": metrics,
                          "negative_control": control, "limitations": "Retrospective computational validation; not wet-lab validation of B26 or the native Type I-C system."}
                st.download_button("Download validation report", json.dumps(report, indent=2, allow_nan=True),
                                   "external_validation_report.json", mime="application/json")
            except (ValueError, pd.errors.ParserError) as error:
                st.error(str(error))

    with st.expander("B26 prospective case", expanded=False):
        st.warning(
            "B26 (*L. fusiformis*) WGS is pending. No B26 genome-specific PAM, guide, editing outcome "
            "or wet-lab validation is claimed. When WGS arrives, the design will be locked before experimental outcomes are reviewed."
        )

with tabs[6]:
    from agripam.constructs_ui import render as render_constructs_tab
    render_constructs_tab(ROOT)
