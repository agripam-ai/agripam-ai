"""Reusable visualizations for RhizoForge-Select genomic target exploration."""

from __future__ import annotations

from html import escape


def _reverse_complement(sequence):
    """Return the reverse complement of a DNA sequence."""
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return sequence.translate(table)[::-1]


def _find_record(records, contig_name):
    """Return the DNA sequence corresponding to a candidate contig."""
    for name, sequence in records:
        if name == contig_name:
            return sequence
    raise ValueError(f"Contig not found: {contig_name}")


def build_target_explorer(
    records,
    candidates,
    selected_candidate,
    flank_bp=250,
):
    """Build an SVG showing genomic context and nucleotide-level guide detail."""

    contig = str(selected_candidate["contig"])
    sequence = _find_record(records, contig)

    start = int(selected_candidate["start_1based"])
    end = int(selected_candidate["end_1based"])
    strand = str(selected_candidate["strand"])
    pam = str(selected_candidate["pam"])
    protospacer = str(selected_candidate["protospacer"])
    candidate_id = str(selected_candidate["candidate_id"])

    region_start = max(1, start - flank_bp)
    region_end = min(len(sequence), end + flank_bp)

    local_candidates = [
        row
        for row in candidates
        if str(row["contig"]) == contig
        and int(row["end_1based"]) >= region_start
        and int(row["start_1based"]) <= region_end
    ]

    plot_left = 80
    plot_right = 1120
    plot_width = plot_right - plot_left

    def xpos(position):
        if region_end == region_start:
            return plot_left
        fraction = (position - region_start) / (region_end - region_start)
        return plot_left + fraction * plot_width

    markers = []

    for row in local_candidates:
        row_start = int(row["start_1based"])
        row_end = int(row["end_1based"])
        row_id = str(row["candidate_id"])
        row_strand = str(row["strand"])

        x1 = xpos(row_start)
        x2 = xpos(row_end)
        width = max(3, x2 - x1)

        selected = row_id == candidate_id
        fill = "#135EA8" if selected else "#8FA8C4"
        y = 155 if row_strand == "+" else 190

        markers.append(
            f'<rect x="{x1:.1f}" y="{y}" width="{width:.1f}" height="14" '
            f'rx="4" fill="{fill}" opacity="{"1" if selected else "0.65"}"/>'
        )

        if selected:
            markers.append(
                f'<text x="{(x1 + x2) / 2:.1f}" y="{y - 9}" '
                f'text-anchor="middle" class="selected-label">'
                f'{escape(row_id)}</text>'
            )

    context_start = max(0, start - 1 - 20)
    context_end = min(len(sequence), end + 20)
    sequence_context = sequence[context_start:context_end]

    upstream_length = (start - 1) - context_start
    target_length = end - start + 1

    upstream = sequence_context[:upstream_length]
    target = sequence_context[
        upstream_length:upstream_length + target_length
    ]
    downstream = sequence_context[
        upstream_length + target_length:
    ]

    if strand == "+":
        molecular_display = (
            f"{escape(upstream)} "
            f'<tspan class="pam">{escape(pam)}</tspan>'
            f'<tspan class="guide">{escape(protospacer)}</tspan> '
            f"{escape(downstream)}"
        )
    else:
        molecular_display = (
            f"{escape(upstream)} "
            f'<tspan class="guide">{escape(target)}</tspan> '
            f"{escape(downstream)}"
        )

    score = float(selected_candidate.get("design_score", 0))
    gc = float(selected_candidate.get("gc_percent", 0))
    copies = int(selected_candidate.get("exact_genome_copies", 0))
    pam_score = float(selected_candidate.get("pam_model_percent", 0))

    svg = f"""
    <style>
      html, body {{
        margin: 0;
        padding: 0;
        background: #F7F9FC;
        font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      }}
      .frame {{
        padding: 18px;
        box-sizing: border-box;
      }}
      svg {{
        width: 100%;
        height: auto;
        display: block;
      }}
      .title {{
        fill: #26334A;
        font-size: 22px;
        font-weight: 700;
      }}
      .subtitle {{
        fill: #53627A;
        font-size: 14px;
      }}
      .axis-label {{
        fill: #53627A;
        font-size: 12px;
      }}
      .selected-label {{
        fill: #135EA8;
        font-size: 12px;
        font-weight: 700;
      }}
      .sequence {{
        fill: #26334A;
        font-family: Menlo, Monaco, monospace;
        font-size: 15px;
      }}
      .pam {{
        fill: #C83E4D;
        font-weight: 800;
      }}
      .guide {{
        fill: #135EA8;
        font-weight: 700;
      }}
      .metric-title {{
        fill: #53627A;
        font-size: 12px;
      }}
      .metric-value {{
        fill: #26334A;
        font-size: 18px;
        font-weight: 700;
      }}
    </style>

    <div class="frame">
    <svg viewBox="0 0 1200 500"
         role="img"
         aria-label="RhizoForge genomic target explorer">

      <text x="60" y="45" class="title">
        Genomic Target Explorer
      </text>

      <text x="60" y="72" class="subtitle">
        {escape(contig)} · selected candidate {escape(candidate_id)}
      </text>

      <text x="60" y="115" class="subtitle">
        Local genomic context ({region_start:,}–{region_end:,} bp)
      </text>

      <line x1="{plot_left}" y1="175"
            x2="{plot_right}" y2="175"
            stroke="#65738A" stroke-width="3"/>

      <text x="{plot_left}" y="225"
            class="axis-label">{region_start:,}</text>

      <text x="{plot_right}" y="225"
            text-anchor="end"
            class="axis-label">{region_end:,}</text>

      <text x="60" y="158" class="axis-label">+ strand</text>
      <text x="60" y="203" class="axis-label">− strand</text>

      {''.join(markers)}

      <line x1="{xpos(start):.1f}" y1="130"
            x2="{xpos(start):.1f}" y2="215"
            stroke="#135EA8"
            stroke-width="2"
            stroke-dasharray="4 4"/>

      <text x="60" y="275" class="title">
        Selected guide
      </text>

      <text x="60" y="305" class="subtitle">
        Position {start:,}–{end:,} · strand {escape(strand)}
      </text>

      <rect x="60" y="330" width="1080" height="58"
            rx="8" fill="#FFFFFF" stroke="#D7DFEA"/>

      <text x="80" y="365" class="sequence">
        5′ {molecular_display} 3′
      </text>

      <text x="60" y="425" class="metric-title">DESIGN SCORE</text>
      <text x="60" y="452" class="metric-value">{score:.1f}</text>

      <text x="280" y="425" class="metric-title">GC CONTENT</text>
      <text x="280" y="452" class="metric-value">{gc:.1f}%</text>

      <text x="500" y="425" class="metric-title">PAM MODEL</text>
      <text x="500" y="452" class="metric-value">{pam_score:.1f}%</text>

      <text x="720" y="425" class="metric-title">EXACT COPIES</text>
      <text x="720" y="452" class="metric-value">{copies}</text>

      <text x="940" y="425" class="metric-title">STRAND</text>
      <text x="940" y="452" class="metric-value">{escape(strand)}</text>

    </svg>
    </div>
    """

    return svg



def build_strategy_map(
    native_systems,
    introduced_editors,
    retained_designs=0,
    pareto_designs=0,
):
    """Build an evidence-aware SVG overview of candidate editing strategies."""

    def native_status(row):
        completeness = str(row.get("completeness", "")).lower()
        if completeness == "minimum component model satisfied":
            return "SUPPORTED", "#2E7D5B"
        if completeness == "partial candidate":
            return "PARTIAL", "#B7791F"
        if completeness == "not detected":
            return "NOT DETECTED", "#7A8494"
        return "UNRESOLVED", "#7A8494"

    def introduced_status(row):
        sites = int(row.get("pam_compatible_sites", 0) or 0)
        tier = str(row.get("compatibility_tier", "")).lower()

        if sites == 0:
            return "NO TARGETS", "#7A8494"

        if "computationally compatible" in tier:
            return "COMPATIBLE", "#2E7D5B"

        if "conditional" in tier:
            return "CONDITIONAL", "#B7791F"

        return "TARGETABLE", "#356FA8"

    native_rows = list(native_systems or [])
    introduced_rows = list(introduced_editors or [])

    # Compact names for the visual map.
    native_labels = {
        "CRISPR interference machinery": "Native CRISPR",
        "Recombinase / integrase route": "Recombinase",
        "Retron / reverse-transcriptase route": "Retron",
        "CRISPR-associated transposase route": "CAST",
        "Homology-directed repair route": "HDR / repair",
        "Alternative RNA-guided nuclease route": "Alternative CRISPR",
    }

    native_cards = []
    for row in native_rows:
        name = native_labels.get(
            str(row.get("system_model", "")),
            str(row.get("system_model", "")),
        )
        status, color = native_status(row)
        native_cards.append((name, status, color))

    introduced_cards = []
    for row in introduced_rows:
        name = str(row.get("editor", ""))
        status, color = introduced_status(row)
        sites = int(row.get("pam_compatible_sites", 0) or 0)
        pam = str(row.get("pam_rule", ""))
        detail = f"{sites:,} sites · {pam}" if sites else pam
        introduced_cards.append((name, status, color, detail))

    def card(x, y, width, title, status, color, detail=""):
        safe_title = escape(title)
        safe_status = escape(status)
        safe_detail = escape(detail)

        detail_svg = (
            f'<text x="{x + 14}" y="{y + 68}" class="card-detail">'
            f'{safe_detail}</text>'
            if safe_detail
            else ""
        )

        return f"""
        <rect x="{x}" y="{y}" width="{width}" height="82"
              rx="10" fill="#FFFFFF" stroke="#D8E0EA"/>
        <rect x="{x}" y="{y}" width="6" height="82"
              rx="3" fill="{color}"/>
        <text x="{x + 14}" y="{y + 25}" class="card-title">
            {safe_title}
        </text>
        <text x="{x + 14}" y="{y + 48}"
              class="card-status" fill="{color}">
            {safe_status}
        </text>
        {detail_svg}
        """

    native_svg = []
    native_width = 245
    for index, (name, status, color) in enumerate(native_cards[:6]):
        col = index % 2
        row_index = index // 2
        x = 55 + col * 265
        y = 205 + row_index * 100
        native_svg.append(
            card(x, y, native_width, name, status, color)
        )

    introduced_svg = []
    introduced_width = 245
    for index, (name, status, color, detail) in enumerate(
        introduced_cards[:3]
    ):
        x = 655
        y = 205 + index * 100
        introduced_svg.append(
            card(x, y, introduced_width, name, status, color, detail)
        )

    retained = int(retained_designs or 0)
    pareto = int(pareto_designs or 0)

    svg = f"""
    <style>
      html, body {{
        margin: 0;
        padding: 0;
        background: #F7F9FC;
        font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      }}
      .frame {{
        padding: 10px;
        box-sizing: border-box;
      }}
      svg {{
        width: 100%;
        height: auto;
        display: block;
      }}
      .title {{
        fill: #26334A;
        font-size: 24px;
        font-weight: 700;
      }}
      .subtitle {{
        fill: #617086;
        font-size: 13px;
      }}
      .section {{
        fill: #26334A;
        font-size: 15px;
        font-weight: 700;
      }}
      .card-title {{
        fill: #26334A;
        font-size: 13px;
        font-weight: 700;
      }}
      .card-status {{
        font-size: 11px;
        font-weight: 800;
      }}
      .card-detail {{
        fill: #617086;
        font-size: 10px;
      }}
      .flow-title {{
        fill: #26334A;
        font-size: 14px;
        font-weight: 700;
      }}
      .flow-detail {{
        fill: #617086;
        font-size: 11px;
      }}
      .metric {{
        fill: #26334A;
        font-size: 20px;
        font-weight: 800;
      }}
    </style>

    <div class="frame">
    <svg viewBox="0 0 1200 620"
         role="img"
         aria-label="RhizoForge editing strategy map">

      <text x="50" y="42" class="title">
        Genome-informed Editing Strategy Map
      </text>

      <text x="50" y="67" class="subtitle">
        Strategy evidence is separated from predicted targetability and
        experimental editing efficiency.
      </text>

      <rect x="390" y="92" width="420" height="58"
            rx="12" fill="#FFFFFF" stroke="#C9D5E3"/>
      <text x="600" y="117" text-anchor="middle" class="flow-title">
        UPLOADED BACTERIAL GENOME
      </text>
      <text x="600" y="137" text-anchor="middle" class="flow-detail">
        annotation · defense · repair · CRISPR/PAM evidence
      </text>

      <line x1="600" y1="150" x2="600" y2="175"
            stroke="#8A98AA" stroke-width="2"/>
      <line x1="310" y1="175" x2="790" y2="175"
            stroke="#8A98AA" stroke-width="2"/>
      <line x1="310" y1="175" x2="310" y2="190"
            stroke="#8A98AA" stroke-width="2"/>
      <line x1="790" y1="175" x2="790" y2="190"
            stroke="#8A98AA" stroke-width="2"/>

      <text x="310" y="195" text-anchor="middle" class="section">
        GENOME-DERIVED ROUTES
      </text>

      <text x="790" y="195" text-anchor="middle" class="section">
        INTRODUCED EDITORS
      </text>

      {''.join(native_svg)}
      {''.join(introduced_svg)}

      <line x1="310" y1="505" x2="310" y2="530"
            stroke="#8A98AA" stroke-width="2"/>
      <line x1="790" y1="505" x2="790" y2="530"
            stroke="#8A98AA" stroke-width="2"/>
      <line x1="310" y1="530" x2="790" y2="530"
            stroke="#8A98AA" stroke-width="2"/>
      <line x1="550" y1="530" x2="550" y2="548"
            stroke="#8A98AA" stroke-width="2"/>

      <rect x="390" y="548" width="320" height="55"
            rx="10" fill="#FFFFFF" stroke="#C9D5E3"/>
      <text x="550" y="570" text-anchor="middle" class="flow-title">
        MULTI-OBJECTIVE DESIGN
      </text>
      <text x="550" y="590" text-anchor="middle" class="flow-detail">
        editability · deliverability · agronomic value · preservation
      </text>

      <rect x="735" y="548" width="180" height="55"
            rx="10" fill="#FFFFFF" stroke="#C9D5E3"/>
      <text x="825" y="570" text-anchor="middle" class="metric">
        {retained:,}
      </text>
      <text x="825" y="590" text-anchor="middle" class="flow-detail">
        retained designs
      </text>

      <rect x="935" y="548" width="180" height="55"
            rx="10" fill="#FFFFFF" stroke="#C9D5E3"/>
      <text x="1025" y="570" text-anchor="middle" class="metric">
        {pareto:,}
      </text>
      <text x="1025" y="590" text-anchor="middle" class="flow-detail">
        Pareto designs
      </text>

    </svg>
    </div>
    """

    return svg
def build_editor_target_explorer(
    candidate,
    features,
    flank_bp=2500,
):
    """Build an SVG explorer for an introduced-editor target and real GFF features."""

    contig = str(candidate.get("contig", ""))
    start = int(candidate.get("start", 0))
    end = int(candidate.get("end", 0))
    pam_start = int(candidate.get("pam_start", 0))
    pam_end = int(candidate.get("pam_end", 0))

    strand = str(candidate.get("strand", "+"))
    editor = str(candidate.get("editor", "Editor"))
    application = str(candidate.get("application", ""))
    pam_pattern = str(candidate.get("pam_pattern", ""))
    pam = str(candidate.get("pam", ""))
    protospacer = str(candidate.get("protospacer", ""))

    region_start = max(1, start - flank_bp)
    region_end = end + flank_bp
    region_length = max(1, region_end - region_start + 1)

    # Only real CDS features overlapping the displayed genomic window.
    local_features = [
        row
        for row in features
        if str(row.get("contig", "")) == contig
        and str(row.get("type", "")).upper() == "CDS"
        and int(row.get("end", 0)) >= region_start
        and int(row.get("start", 0)) <= region_end
    ]

    def xcoord(position):
        left = 70
        width = 1060
        fraction = (int(position) - region_start) / region_length
        return left + max(0.0, min(1.0, fraction)) * width

    target_x1 = xcoord(start)
    target_x2 = xcoord(end)
    target_width = max(4, target_x2 - target_x1)

    pam_x1 = xcoord(pam_start)
    pam_x2 = xcoord(pam_end)
    pam_width = max(4, pam_x2 - pam_x1)

    # Derive the PAM side from the explicit genomic coordinates.
    if strand == "+":
        pam_side = "5′" if pam_end < (start + len(protospacer)) else "3′"
    else:
        pam_side = "5′" if pam_start > (end - len(protospacer)) else "3′"

    # The stored protospacer/PAM sequences are already editor-oriented.
    if pam_side == "3′":
        oriented_target = f"{protospacer}  {pam}"
        sequence_labels = f"5′  {protospacer}  {pam}  3′"
    else:
        oriented_target = f"{pam}  {protospacer}"
        sequence_labels = f"5′  {pam}  {protospacer}  3′"

    gc = candidate.get("gc_percent", "")
    uniqueness = str(
        candidate.get(
            "exact_uniqueness",
            candidate.get("guide_uniqueness", "not reported"),
        )
    )
    score = candidate.get(
        "offtarget_adjusted_priority_score",
        candidate.get("guide_priority_score", ""),
    )
    off_status = str(candidate.get("off_target_screen_status", "not reported"))
    off_count = candidate.get("approx_offtargets_le2", "")

    if off_status == "screened":
        if str(off_count) not in ("", "nan", "None"):
            off_label = f"{off_count} alternative PAM-compatible sites ≤2 substitutions"
        else:
            off_label = "Completed substitution-only screen"
    elif off_status == "not_screened_reference_limit":
        off_label = "Not screened — PAM-compatible reference space exceeded limit"
    elif off_status == "not_prioritized":
        off_label = "Not screened — outside prioritized screening set"
    else:
        off_label = off_status.replace("_", " ")

    feature_svg = []
    for idx, feature in enumerate(local_features):
        fs = max(region_start, int(feature["start"]))
        fe = min(region_end, int(feature["end"]))
        x1 = xcoord(fs)
        x2 = xcoord(fe)
        width = max(5, x2 - x1)

        fstrand = str(feature.get("strand", "."))
        gene = str(feature.get("gene", "") or "")
        fid = str(feature.get("id", "") or "")
        product = str(feature.get("product", "") or "")

        label = gene or fid or "CDS"
        if product and product != label:
            tooltip = f"{label}: {product}"
        else:
            tooltip = label

        y = 222 if idx % 2 == 0 else 262

        if fstrand == "-":
            points = (
                f"{x1 + 10},{y - 13} "
                f"{x2},{y - 13} "
                f"{x2},{y + 13} "
                f"{x1 + 10},{y + 13} "
                f"{x1},{y}"
            )
        else:
            points = (
                f"{x1},{y - 13} "
                f"{max(x1, x2 - 10)},{y - 13} "
                f"{x2},{y} "
                f"{max(x1, x2 - 10)},{y + 13} "
                f"{x1},{y + 13}"
            )

        feature_svg.append(
            f"""
            <g>
              <title>{escape(tooltip)}</title>
              <polygon points="{points}" class="cds"/>
              <text x="{(x1 + x2) / 2:.1f}" y="{y + 31}"
                    text-anchor="middle" class="gene-label">
                {escape(label[:24])}
              </text>
            </g>
            """
        )

    score_text = (
        f"{float(score):.1f}"
        if str(score) not in ("", "nan", "None")
        else "not reported"
    )
    gc_text = (
        f"{float(gc):.1f}%"
        if str(gc) not in ("", "nan", "None")
        else "not reported"
    )

    svg = f"""
    <div style="font-family:Inter,Arial,sans-serif;width:100%;">
    <svg viewBox="0 0 1200 650"
         xmlns="http://www.w3.org/2000/svg"
         style="width:100%;height:auto;background:#ffffff;border-radius:14px;">

      <style>
        .title {{ font-size:25px;font-weight:700;fill:#14213d; }}
        .subtitle {{ font-size:14px;fill:#536273; }}
        .section {{ font-size:17px;font-weight:700;fill:#14213d; }}
        .small {{ font-size:12px;fill:#536273; }}
        .metric {{ font-size:14px;font-weight:600;fill:#14213d; }}
        .gene-label {{ font-size:10px;fill:#34495e; }}
        .cds {{ fill:#9bb8d3;stroke:#557a9e;stroke-width:1; }}
        .target {{ fill:#f6bd60;stroke:#d58a00;stroke-width:2; }}
        .pam {{ fill:#e76f51;stroke:#b9432d;stroke-width:2; }}
        .axis {{ stroke:#64748b;stroke-width:2; }}
        .guide {{ font-family:monospace;font-size:20px;font-weight:700;fill:#14213d; }}
        .box {{ fill:#f8fafc;stroke:#d7dee7;stroke-width:1; }}
      </style>

      <text x="55" y="48" class="title">
        {escape(editor)} target explorer
      </text>
      <text x="55" y="73" class="subtitle">
        {escape(application)} • {escape(contig)}
      </text>

      <rect x="55" y="95" width="1090" height="58" rx="10" class="box"/>
      <text x="75" y="119" class="metric">
        Target: {start:,}–{end:,} bp • strand {escape(strand)}
      </text>
      <text x="75" y="141" class="small">
        PAM rule {escape(pam_pattern)} • observed PAM {escape(pam)} • PAM coordinates {pam_start:,}–{pam_end:,}
      </text>

      <text x="55" y="190" class="section">Genomic context</text>
      <line x1="70" y1="242" x2="1130" y2="242" class="axis"/>

      {''.join(feature_svg)}

      <rect x="{target_x1:.1f}" y="205"
            width="{target_width:.1f}" height="74"
            rx="4" class="target"/>
      <rect x="{pam_x1:.1f}" y="197"
            width="{pam_width:.1f}" height="90"
            rx="3" class="pam"/>

      <text x="70" y="318" class="small">{region_start:,} bp</text>
      <text x="1130" y="318" text-anchor="end" class="small">{region_end:,} bp</text>

      <rect x="55" y="345" width="1090" height="116" rx="12" class="box"/>
      <text x="75" y="375" class="section">Selected guide</text>
      <text x="75" y="410" class="guide">{escape(sequence_labels)}</text>
      <text x="75" y="438" class="small">
        protospacer + PAM shown in editor-oriented 5′→3′ direction • PAM side {escape(pam_side)}
      </text>

      <rect x="55" y="485" width="255" height="105" rx="12" class="box"/>
      <text x="75" y="515" class="small">GC content</text>
      <text x="75" y="550" class="title">{escape(gc_text)}</text>

      <rect x="325" y="485" width="255" height="105" rx="12" class="box"/>
      <text x="345" y="515" class="small">Exact genomic uniqueness</text>
      <text x="345" y="550" class="metric">{escape(uniqueness)}</text>

      <rect x="595" y="485" width="255" height="105" rx="12" class="box"/>
      <text x="615" y="515" class="small">Priority score</text>
      <text x="615" y="550" class="title">{escape(score_text)}</text>

      <rect x="865" y="485" width="280" height="105" rx="12" class="box"/>
      <text x="885" y="515" class="small">Off-target evidence</text>
      <foreignObject x="885" y="530" width="240" height="50">
        <div xmlns="http://www.w3.org/1999/xhtml"
             style="font-size:12px;font-weight:600;color:#14213d;line-height:1.3;">
          {escape(off_label)}
        </div>
      </foreignObject>

      <text x="55" y="625" class="small">
        CDS arrows are derived from the workflow GFF annotation. Targetability and scores are computational decision-support evidence, not experimental editing efficiency.
      </text>
    </svg>
    </div>
    """

    return svg

def build_introduced_editor_funnel(rows):
    """Build an evidence-aware funnel for introduced-editor targetability."""

    cards = []

    for row in rows or []:
        editor = str(row.get("editor", "Editor"))
        pam_rule = str(row.get("pam_rule", ""))
        pam_sites = int(row.get("pam_compatible_sites", 0) or 0)
        unique_sites = int(row.get("exact_unique_sites", 0) or 0)

        status = str(row.get("off_target_screening_status", ""))
        screened = int(
            row.get("approximately_screened_top_sites", 0) or 0
        )
        clean = row.get(
            "screened_sites_without_1_or_2_mismatch_hit",
            None,
        )

        if status == "screened" and screened:
            if clean is not None:
                specificity = (
                    f"{int(clean):,} / {screened:,}<br>"
                    "<small>without ≤2-mismatch hits</small>"
                )
            else:
                specificity = (
                    f"{screened:,} screened<br>"
                    "<small>clean count unresolved</small>"
                )
            specificity_class = "screened"
        elif status == "not_screened_reference_limit":
            specificity = (
                "<b>SCREENING SKIPPED</b><br>"
                "<small>Reference space exceeded 500,000 sequences.<br>"
                "Specificity remains unresolved.</small>"
            )
            specificity_class = "unresolved"
        else:
            specificity = (
                "<b>NOT SCREENED</b><br>"
                "<small>Specificity remains unresolved.</small>"
            )
            specificity_class = "unresolved"

        cards.append(
            f"""
            <div class="rf-editor-card">
              <div class="rf-editor-title">{escape(editor)}</div>
              <div class="rf-editor-subtitle">
                PAM rule: {escape(pam_rule)}
              </div>

              <div class="rf-stage">
                <span>PAM-compatible sites</span>
                <strong>{pam_sites:,}</strong>
              </div>

              <div class="rf-arrow">↓</div>

              <div class="rf-stage">
                <span>Exact-unique sites</span>
                <strong>{unique_sites:,}</strong>
              </div>

              <div class="rf-arrow">↓</div>

              <div class="rf-stage">
                <span>Priority set</span>
                <strong>Top 2,000</strong>
              </div>

              <div class="rf-arrow">↓</div>

              <div class="rf-specificity {specificity_class}">
                <span>Bounded specificity evidence</span>
                <div>{specificity}</div>
              </div>
            </div>
            """
        )

    return f"""
    <div class="rf-funnel">
      <style>
        .rf-funnel {{
          font-family: Inter, Arial, sans-serif;
          width: 100%;
          box-sizing: border-box;
        }}

        .rf-funnel-grid {{
          display: grid;
          grid-template-columns: repeat({max(1, len(cards))}, 1fr);
          gap: 16px;
        }}

        .rf-editor-card {{
          border: 1px solid #d7dee7;
          border-radius: 14px;
          padding: 16px;
          background: #ffffff;
          box-sizing: border-box;
        }}

        .rf-editor-title {{
          font-size: 20px;
          font-weight: 700;
          color: #14213d;
        }}

        .rf-editor-subtitle {{
          font-size: 12px;
          color: #536273;
          margin: 4px 0 14px;
        }}

        .rf-stage {{
          border: 1px solid #d7dee7;
          border-radius: 10px;
          padding: 10px;
          text-align: center;
          background: #f8fafc;
        }}

        .rf-stage span {{
          display: block;
          font-size: 12px;
          color: #536273;
        }}

        .rf-stage strong {{
          display: block;
          margin-top: 3px;
          font-size: 19px;
          color: #14213d;
        }}

        .rf-arrow {{
          text-align: center;
          font-size: 18px;
          color: #64748b;
          padding: 4px;
        }}

        .rf-specificity {{
          border-radius: 10px;
          padding: 12px;
          text-align: center;
          font-size: 13px;
          line-height: 1.35;
        }}

        .rf-specificity span {{
          display: block;
          font-size: 12px;
          margin-bottom: 6px;
          color: #536273;
        }}

        .rf-specificity small {{
          font-size: 11px;
        }}

        .rf-specificity.screened {{
          background: #eef8f0;
          border: 1px solid #9bc8a1;
          color: #245b2b;
        }}

        .rf-specificity.unresolved {{
          background: #fff7e6;
          border: 1px solid #e3bd70;
          color: #76520b;
        }}
      </style>

      <div class="rf-funnel-grid">
        {''.join(cards)}
      </div>
    </div>
    """