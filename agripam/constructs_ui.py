"""Streamlit tab: BioBrick parts, construct design for the chosen edit, and an experiment tracker."""
from __future__ import annotations

import io
from pathlib import Path

import pandas as pd
import streamlit as st

from . import biobricks as bb
from . import parts_io
from .core import parse_fasta

LEVEL_ICON = {"ok": "✅", "warning": "⚠️", "error": "⛔"}


def _library() -> dict:
    return st.session_state.setdefault("parts_library", {})


def _genome_records(root: Path, uploaded) -> tuple[list[tuple[str, str]], str, Path | None]:
    """Records from the analysed genome (Genome evaluation tab) or from an uploaded FASTA."""
    full = st.session_state.get("full_workflow")
    result_dir = Path(full["result_dir"]) if full and full.get("result_dir") else None
    if uploaded is not None:
        return parse_fasta(uploaded.getvalue().decode("utf-8", "replace")), f"uploaded file {uploaded.name}", result_dir
    if result_dir:
        for name in ("genome.fna", "normalized_genome.fna"):
            path = result_dir / name
            if path.exists():
                return parse_fasta(path.read_text()), f"genome analysed for {full.get('source', 'the current strain')}", result_dir
    return [], "", result_dir


def _read_table(result_dir: Path | None, name: str) -> pd.DataFrame:
    if result_dir and (result_dir / name).exists():
        return pd.read_csv(result_dir / name, sep="\t")
    return pd.DataFrame()


def _pick(label: str, parts: dict, roles: list[str], key: str, optional: bool = False):
    options = [pid for pid, p in parts.items() if p["role"] in roles]
    if optional:
        options = ["(none)"] + options
    if not options:
        st.warning(f"No part with role {' / '.join(roles)} in the library yet.")
        return None
    choice = st.selectbox(label, options, key=key,
                          format_func=lambda pid: pid if pid == "(none)" else f"{pid}: {parts[pid]['name'][:48]}")
    return None if choice == "(none)" else parts[choice]


def render(root: Path) -> None:
    st.header("Parts and constructs: from a chosen edit to a followed experiment")
    st.info(
        "Choose BioBrick-style parts, assemble the construct for the edit you selected, screen it against standard-assembly "
        "rules and the chassis's restriction motifs, then track what happens at the bench. A compatible part is not a "
        "working part: activity in your chassis still has to be measured.",
        icon=":material/build:")
    library = _library()

    # ---------------------------------------------------------------- 1. library
    with st.expander("1. Parts library: starter set, Excel sheet or iGEM Registry", expanded=True):
        c1, c2, c3 = st.columns(3)
        starter = root / "data" / "biobricks" / "starter_library.tsv"
        if c1.button("Load the starter library", key="bb_starter", disabled=not starter.exists(), width="stretch"):
            for part in parts_io.read_parts_tsv(starter):
                library[part["part_id"]] = part
        template = root / "data" / "biobricks" / "Parts_input_template.xlsx"
        if template.exists():
            c2.download_button("Download the parts template (.xlsx)", template.read_bytes(), template.name,
                               mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet", key="bb_tpl", width="stretch")
        uploaded_parts = c3.file_uploader("Upload your parts workbook", type=["xlsx"], key="bb_upload", label_visibility="collapsed")
        if uploaded_parts is not None and st.session_state.get("bb_last_upload") != uploaded_parts.file_id:
            try:
                parts, issues = parts_io.read_parts_workbook(uploaded_parts)
                for issue in issues:
                    (st.error if issue["level"] == "error" else st.warning)(f"{issue['where']}: {issue['message']}")
                for part in parts:
                    library[part["part_id"]] = part
                st.session_state["bb_last_upload"] = uploaded_parts.file_id
                if parts:
                    st.success(f"Added {len(parts)} part(s) from the workbook.")
            except Exception as error:
                st.error(f"The workbook could not be read: {error}")
        f1, f2 = st.columns([3, 1], vertical_alignment="bottom")
        ids = f1.text_input("Fetch from the iGEM Registry (IDs separated by commas)", "BBa_J23100, BBa_B0034, BBa_E0040, BBa_B0015",
                            key="bb_ids")
        if f2.button("Fetch", key="bb_fetch", width="stretch"):
            for pid in [i.strip() for i in ids.split(",") if i.strip()]:
                try:
                    part = bb.fetch_registry_part(pid)
                    part["hosts"] = library.get(pid, {}).get("hosts", "")
                    library[part["part_id"]] = part
                    st.success(f"{pid}: {len(part['sequence'])} bp, role set to {part['role']} (edit below if needed).")
                except (ValueError, RuntimeError) as error:
                    st.error(str(error))
        if library:
            frame = pd.DataFrame([{"part_id": p["part_id"], "name": p["name"], "role": p["role"], "hosts": p.get("hosts", ""),
                                   "length_bp": len(p["sequence"]), "GC %": bb.gc_percent(p["sequence"]),
                                   "source": p.get("source", "")} for p in library.values()])
            st.caption("Edit the role or the declared host range (semicolon-separated, for example “E. coli; B. subtilis”). "
                       "Host ranges are used to warn when a part has not been declared for your chassis.")
            edited = st.data_editor(frame, hide_index=True, width="stretch", key="bb_editor",
                                    disabled=["part_id", "name", "length_bp", "GC %", "source"],
                                    column_config={"role": st.column_config.SelectboxColumn("role", options=bb.ROLES)})
            for _, row in edited.iterrows():
                library[row["part_id"]]["role"] = row["role"]
                library[row["part_id"]]["hosts"] = row["hosts"] or ""
            if st.button("Clear the library", key="bb_clear"):
                library.clear()
                st.rerun()
        else:
            st.info("The library is empty. Load the starter set, upload a workbook, or fetch parts by ID.")

    # ---------------------------------------------------------------- 2. design
    with st.expander("2. Choose the edit and build the construct", expanded=True):
        full = st.session_state.get("full_workflow")
        default_chassis = (full or {}).get("identity", {}).get("organism", "") if full else ""
        s1, s2 = st.columns(2)
        chassis = s1.text_input("Chassis (species)", default_chassis, placeholder="for example B. subtilis", key="bb_chassis")
        standards = s2.multiselect("Assembly standards to screen", list(bb.STANDARDS), default=["RFC10 (BioBrick)"], key="bb_std")
        motif_text = st.text_area(
            "Chassis restriction motifs (optional, one per line as name = sequence, IUPAC allowed)", "", height=80, key="bb_motifs",
            placeholder="HindIII-like = AAGCTT\nType I site = GCANNNNNNNTGC")
        motifs = {}
        for line in motif_text.splitlines():
            if "=" in line:
                label, motif = (x.strip() for x in line.split("=", 1))
                motifs[label] = motif
        scar_label = st.selectbox("Junction between parts", list(bb.SCARS), key="bb_scar")
        edit_type = st.radio("Edit type", ["Reporter insertion at a neutral site", "CRISPRi guide cassette"], horizontal=True, key="bb_edit")
        genome_upload = st.file_uploader("Genome FASTA (optional if you analysed a genome in the Genome evaluation tab)",
                                         type=["fa", "fasta", "fna"], key="bb_genome")
        records, genome_label, result_dir = _genome_records(root, genome_upload)
        if records:
            st.caption(f"Using {genome_label}: {len(records)} contig(s), {sum(len(s) for _, s in records):,} bp.")
        built = None

        if edit_type == "Reporter insertion at a neutral site":
            regions = _read_table(result_dir, "candidate_neutral_regions.tsv")
            contig = position = None
            if not regions.empty:
                st.caption("Candidate neutral regions from the genome workflow (candidates only: essentiality, synteny and "
                           "stability still need checking).")
                st.dataframe(regions.head(20), hide_index=True, width="stretch")
            if records:
                default_contig = records[0][0]
                if not regions.empty:
                    default_contig = str(regions.iloc[0]["contig"])
                names = [c for c, _ in records]
                contig = st.selectbox("Contig", names, index=names.index(default_contig) if default_contig in names else 0, key="bb_contig")
                length = len(dict(records)[contig])
                mid = int((regions.iloc[0]["start"] + regions.iloc[0]["end"]) // 2) if (not regions.empty and str(regions.iloc[0]["contig"]) == contig) else length // 2
                position = st.number_input("Insert after base (1-based)", 1, max(1, length - 1), min(mid, max(1, length - 1)), key="bb_pos")
            else:
                st.info("Upload a genome FASTA, or analyse a genome first, to set the insertion site.")
            a1, a2 = st.columns(2)
            arm = a1.slider("Homology arm length (bp)", 200, 1500, 750, 50, key="bb_arm")
            orientation = a2.radio("Cassette orientation vs. the + strand", ["+", "-"], horizontal=True, key="bb_orient")
            p1, p2, p3, p4 = st.columns(4)
            with p1:
                promoter = _pick("Promoter", library, ["promoter"], "bb_pr")
            with p2:
                rbs = _pick("RBS (optional)", library, ["rbs"], "bb_rbs", optional=True)
            with p3:
                reporter = _pick("Reporter / coding part", library, ["reporter", "cds"], "bb_rep")
            with p4:
                terminator = _pick("Terminator (optional)", library, ["terminator"], "bb_term", optional=True)
            if st.button("Build the construct", type="primary", key="bb_build_rep", disabled=not (records and promoter and reporter)):
                try:
                    chosen = [x for x in (promoter, rbs, reporter, terminator) if x]
                    cassette = bb.assemble(chosen, scar=bb.SCARS[scar_label])
                    arms = bb.homology_arms(records, contig, int(position), int(arm))
                    construct = bb.integration_construct(cassette, arms, orientation)
                    built = {"edit_type": "reporter insertion", "target": f"{contig} after base {int(position)}",
                             "parts": " + ".join(x["part_id"] for x in chosen), "sequence": construct["sequence"],
                             "annotations": construct["annotations"], "pcr": bb.diagnostic_pcr(arms, len(cassette["sequence"])),
                             "part_list": chosen, "predicted_score": ""}
                except ValueError as error:
                    st.error(str(error))
        else:
            targets = _read_table(result_dir, "introduced_editor_targets.tsv")
            portfolio = _read_table(result_dir, "multiobjective_design_portfolio.tsv")
            spacer, chosen_label, score = "", "", ""
            if not targets.empty and "editor" in targets.columns:
                dcas = targets[targets["editor"].astype(str).str.contains("dCas9", na=False)].copy()
                if not dcas.empty:
                    dcas["label"] = dcas.apply(lambda r: f"{r.get('contig', '')}:{r.get('start', '')}{r.get('strand', '')} {r.get('protospacer', '')}", axis=1)
                    pick = st.selectbox("Guide from the genome workflow (dCas9/CRISPRi sites)", ["(paste my own)"] + dcas["label"].head(200).tolist(), key="bb_guide")
                    if pick != "(paste my own)":
                        row = dcas[dcas["label"] == pick].iloc[0]
                        spacer, chosen_label = str(row["protospacer"])[:20], pick
                        if not portfolio.empty and "protospacer" in portfolio.columns and "balanced_score" in portfolio.columns:
                            match = portfolio[portfolio["protospacer"] == row["protospacer"]]
                            score = float(match.iloc[0]["balanced_score"]) if not match.empty else ""
            if not spacer:
                spacer = st.text_input("20-nt spacer (the protospacer, without the PAM)", "", key="bb_spacer", max_chars=20)
                chosen_label = f"user spacer {spacer}"
            p1, p2, p3 = st.columns(3)
            with p1:
                promoter = _pick("sgRNA promoter", library, ["promoter"], "bb_pr_g")
            with p2:
                scaffold = _pick("sgRNA scaffold", library, ["scaffold"], "bb_scaf")
            with p3:
                terminator = _pick("Terminator (optional)", library, ["terminator"], "bb_term_g", optional=True)
            st.caption("dCas9 is assumed to be expressed already (from the chassis or an existing plasmid). "
                       "Add a dCas9 expression cassette by building a reporter-style construct with your dCas9 coding part.")
            if st.button("Build the guide cassette", type="primary", key="bb_build_g", disabled=not (promoter and scaffold and len(spacer) == 20)):
                try:
                    cassette = bb.sgrna_cassette(promoter, spacer, scaffold, terminator, scar=bb.SCARS[scar_label])
                    built = {"edit_type": "CRISPRi", "target": chosen_label, "parts": f"{promoter['part_id']} + spacer + {scaffold['part_id']}"
                             + (f" + {terminator['part_id']}" if terminator else ""), "sequence": cassette["sequence"],
                             "annotations": cassette["annotations"], "pcr": None,
                             "part_list": [x for x in (promoter, scaffold, terminator) if x], "predicted_score": score}
                except ValueError as error:
                    st.error(str(error))

        if built:
            checks = []
            for part in built["part_list"]:
                checks += bb.check_sequence(part["part_id"], part["sequence"], standards, motifs, chassis, part.get("hosts", ""), True)
            checks += bb.check_sequence("whole construct", built["sequence"], standards, motifs, chassis, "", False)
            constructs = st.session_state.setdefault("bb_constructs", [])
            built["construct_id"] = f"K{len(constructs) + 1:03d}"
            built["checks"] = checks
            built["chassis"] = chassis
            constructs.append(built)
            st.success(f"Built {built['construct_id']} ({len(built['sequence']):,} bp). See section 3.")

    # ---------------------------------------------------------------- 3. results
    constructs = st.session_state.get("bb_constructs", [])
    with st.expander("3. Constructs, checks and downloads", expanded=True):
        if not constructs:
            st.info("Build a construct in section 2 to see its checks and download it.")
        else:
            chosen_id = st.selectbox("Construct", [c["construct_id"] for c in constructs], index=len(constructs) - 1, key="bb_show")
            c = next(x for x in constructs if x["construct_id"] == chosen_id)
            m1, m2, m3, m4 = st.columns(4)
            m1.metric("Length", f"{len(c['sequence']):,} bp")
            m2.metric("GC", f"{bb.gc_percent(c['sequence'])}%")
            errors = sum(1 for f in c["checks"] if f["level"] == "error")
            warns = sum(1 for f in c["checks"] if f["level"] == "warning")
            m3.metric("Errors", errors)
            m4.metric("Warnings", warns)
            if errors:
                st.error("This construct contains sites forbidden by the selected standard. It may still work by other assembly "
                         "methods, but not by standard BioBrick cloning.")
            st.markdown("**Checks**")
            st.dataframe(pd.DataFrame([{"": LEVEL_ICON[f["level"]], **{k: f[k] for k in ("item", "check", "detail")}} for f in c["checks"]]),
                         hide_index=True, width="stretch")
            st.markdown("**Annotated map**")
            st.dataframe(pd.DataFrame(c["annotations"])[[k for k in ("label", "role", "start", "end", "strand") if k in pd.DataFrame(c["annotations"]).columns]],
                         hide_index=True, width="stretch")
            if c.get("pcr"):
                st.caption(f"Diagnostic PCR with primers {c['pcr']['outside_bp']} bp outside each arm: wild type {c['pcr']['wild_type_bp']} bp, "
                           f"edited {c['pcr']['edited_bp']} bp. A different size or no band means the edit needs sequencing.")
            d1, d2 = st.columns(2)
            d1.download_button("Download FASTA", bb.to_fasta(c["construct_id"], c["sequence"]), f"{c['construct_id']}.fasta", key=f"bb_fa_{chosen_id}", width="stretch")
            d2.download_button("Download GenBank", bb.to_genbank(c["construct_id"], c["sequence"], c["annotations"], c["parts"]),
                               f"{c['construct_id']}.gb", key=f"bb_gb_{chosen_id}", width="stretch")

    # ---------------------------------------------------------------- 4. tracker
    with st.expander("4. Follow the output: experiment tracker", expanded=True):
        st.write("Every construct gets a tracker row, plus the controls an experiment needs (parental, negative, and a low-ranked guide "
                 "for CRISPRi). Fill the stage columns and the measurements, download the file, and upload it again to see progress "
                 "and how the predicted ranking compares with what was measured.")
        group = st.text_input("Experiment name", "experiment 1", key="bb_group")
        t1, t2 = st.columns(2)
        if t1.button("Create tracker rows from the built constructs", key="bb_mk_tracker", disabled=not constructs, width="stretch"):
            st.session_state["bb_tracker"] = bb.tracker_rows(
                [{"construct_id": c["construct_id"], "edit_type": c["edit_type"], "target": c["target"], "parts": c["parts"],
                  "length_bp": len(c["sequence"]), "predicted_score": c.get("predicted_score", "")} for c in constructs], group)
        tracker_upload = t2.file_uploader("Or upload a filled tracker (.xlsx)", type=["xlsx"], key="bb_tracker_up", label_visibility="collapsed")
        if tracker_upload is not None and st.session_state.get("bb_last_tracker") != tracker_upload.file_id:
            try:
                st.session_state["bb_tracker"] = parts_io.read_tracker(tracker_upload)
                st.session_state["bb_last_tracker"] = tracker_upload.file_id
            except Exception as error:
                st.error(f"The tracker could not be read: {error}")
        rows = st.session_state.get("bb_tracker")
        if not rows:
            st.info("No tracker yet: build a construct and create the rows, or upload a tracker.")
        else:
            edited = st.data_editor(pd.DataFrame(rows, columns=bb.TRACKER_COLUMNS), hide_index=True, width="stretch",
                                    num_rows="dynamic", key="bb_tracker_editor")
            rows = edited.where(pd.notna(edited), "").to_dict(orient="records")
            st.session_state["bb_tracker"] = rows
            analysis = bb.analyze_tracker(rows)
            k1, k2, k3 = st.columns(3)
            k1.metric("Constructs designed", analysis["funnel"]["designed"])
            k2.metric("Edit confirmed", analysis["funnel"]["edit_confirmed"])
            k3.metric("Predicted vs measured ρ", "n/a" if analysis["spearman_rho"] is None else f"{analysis['spearman_rho']:.2f}",
                      help="Spearman correlation between the predicted score and the normalised effect; needs 5 or more constructs with both.")
            st.bar_chart(pd.DataFrame({"constructs": list(analysis["funnel"].values())}, index=[s.replace("_", " ") for s in analysis["funnel"]]))
            for warning in analysis["warnings"]:
                st.warning(warning)
            if analysis["spearman_rho"] is not None:
                st.caption(f"Based on {analysis['n_paired']} constructs with both a predicted score and a measured effect.")
            st.download_button("Download the tracker (.xlsx)", parts_io.tracker_to_xlsx(rows), "construct_tracker.xlsx",
                               mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet", key="bb_tracker_dl")
