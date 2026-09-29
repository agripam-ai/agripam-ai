"""Small layout kit: one slide description, two renderers (python-pptx file and PNG preview).

Text is measured with real Arial metrics and shrunk until it fits its box, so the
preview and the .pptx agree and nothing overflows. Requires python-pptx and Pillow.
"""
from __future__ import annotations

from pathlib import Path

from PIL import Image, ImageDraw, ImageFont
from pptx import Presentation
from pptx.chart.data import CategoryChartData
from pptx.dml.color import RGBColor
from pptx.enum.chart import XL_CHART_TYPE, XL_LABEL_POSITION, XL_TICK_LABEL_POSITION
from pptx.enum.shapes import MSO_SHAPE
from pptx.enum.text import MSO_ANCHOR, PP_ALIGN
from pptx.oxml.ns import qn
from pptx.util import Emu, Inches, Pt
from lxml import etree

FONT_DIR = Path("/System/Library/Fonts/Supplemental")
FONT = "Arial"
_cache: dict = {}


def font(size: float, bold: bool = False):
    key = (round(size * 4), bold)
    if key not in _cache:
        name = "Arial Bold.ttf" if bold else "Arial.ttf"
        _cache[key] = ImageFont.truetype(str(FONT_DIR / name), int(round(size * 4)))  # 4x for precision
    return _cache[key]


def measure(text: str, size: float, bold: bool) -> float:
    return font(size, bold).getlength(text) / 4.0  # points


def rgb(hexstr: str) -> RGBColor:
    return RGBColor.from_string(hexstr)


def norm_paragraphs(content, bold=False, color=None):
    """content: str | list of paragraphs; paragraph: str | list of runs; run: str | (text, {bold,color})."""
    if isinstance(content, str):
        content = [content]
    out = []
    for para in content:
        opts = {}
        if isinstance(para, dict):
            opts, para = para, para["runs"]
        if isinstance(para, str):
            para = [para]
        runs = []
        for run in para:
            text, style = (run, {}) if isinstance(run, str) else run
            runs.append((text, {"bold": style.get("bold", bold), "color": style.get("color", color)}))
        out.append({"runs": runs, "bullet": opts.get("bullet", False), "after": opts.get("after", None)})
    return out


def wrap(paras, size, width_pt, bullet_indent_pt):
    """Greedy word wrap. Returns list of lines; each line = (para_index, [(word, bold, color)], is_first)."""
    lines = []
    for pi, para in enumerate(paras):
        avail = width_pt - (bullet_indent_pt if para["bullet"] else 0)
        words = []
        for text, style in para["runs"]:
            for i, w in enumerate(text.split(" ")):
                if w == "" and i > 0:
                    continue
                words.append((w, style["bold"], style["color"]))
        line, used, first = [], 0.0, True
        space = measure(" ", size, False)
        for w, b, c in words:
            wl = measure(w, size, b)
            add = wl if not line else wl + space
            if line and used + add > avail:
                lines.append((pi, line, first))
                line, used, first = [], 0.0, False
                add = wl
            line.append((w, b, c))
            used += add
        lines.append((pi, line, first))
    return lines


def layout_text(el):
    """Return (size, lines) after shrinking to fit; sets el['overflow'] if it cannot fit."""
    paras = norm_paragraphs(el["text"], el.get("bold", False), el.get("color"))
    size = el.get("size", 14)
    minimum = el.get("min_size", max(8, size * 0.6))
    width_pt = el["w"] * 72 - 2 * el.get("pad", 0) * 72
    height_pt = el["h"] * 72 - 2 * el.get("pad", 0) * 72
    indent = 14.0
    ls = el.get("line_spacing", 1.18)
    while True:
        lines = wrap(paras, size, width_pt * 0.985, indent)
        after = sum((p["after"] if p["after"] is not None else el.get("para_after", 0)) for p in paras[:-1])
        total = len(lines) * size * ls + after
        widest = max((sum(measure(w, size, b) for w, b, _ in ln[1]) for ln in lines), default=0)
        if total <= height_pt and (widest <= width_pt or True):
            break
        if size - 0.5 < minimum:
            el["overflow"] = True
            break
        size -= 0.5
    return size, lines, paras, ls


class Slide:
    def __init__(self, width, height, bg="FFFFFF"):
        self.w, self.h, self.bg = width, height, bg
        self.els: list[dict] = []
        self.notes = ""

    def rect(self, x, y, w, h, fill=None, line=None, radius=None, lw=0.75):
        self.els.append(dict(kind="rect", x=x, y=y, w=w, h=h, fill=fill, line=line, radius=radius, lw=lw))

    def oval(self, x, y, w, h, fill=None, line=None):
        self.els.append(dict(kind="oval", x=x, y=y, w=w, h=h, fill=fill, line=line, lw=0.75))

    def arrow(self, x, y, w, h, fill):
        self.els.append(dict(kind="arrow", x=x, y=y, w=w, h=h, fill=fill))

    def text(self, x, y, w, h, text, size=14, color="17322E", bold=False, align="l", valign="t", **kw):
        self.els.append(dict(kind="text", x=x, y=y, w=w, h=h, text=text, size=size, color=color, bold=bold,
                             align=align, valign=valign, **kw))

    def image(self, x, y, w, h, path, align="c"):
        self.els.append(dict(kind="image", x=x, y=y, w=w, h=h, path=str(path), align=align))

    def bar_chart(self, x, y, w, h, categories, values, color, title=None, number_format="0.00", vmin=None,
                  vmax=None, horizontal=False, point_colors=None):
        self.els.append(dict(kind="chart", x=x, y=y, w=w, h=h, categories=categories, values=values,
                             color=color, title=title, fmt=number_format, vmin=vmin, vmax=vmax,
                             horizontal=horizontal, point_colors=point_colors))


def image_box(el):
    with Image.open(el["path"]) as im:
        iw, ih = im.size
    scale = min(el["w"] / iw, el["h"] / ih)
    w, h = iw * scale, ih * scale
    x = el["x"] + (el["w"] - w) / 2 if el["align"] == "c" else el["x"]
    y = el["y"] + (el["h"] - h) / 2
    return x, y, w, h


# ---------------------------------------------------------------- PPTX
def _set_bullet(paragraph, on: bool):
    pPr = paragraph._p.get_or_add_pPr()
    for tag in ("a:buNone", "a:buChar", "a:buAutoNum"):
        for child in pPr.findall(qn(tag)):
            pPr.remove(child)
    if on:
        pPr.set("marL", str(int(14 * 12700)))
        pPr.set("indent", str(int(-14 * 12700)))
        bu = etree.SubElement(pPr, qn("a:buChar"))
        bu.set("char", "•")
    else:
        pPr.set("marL", "0")
        pPr.set("indent", "0")
        etree.SubElement(pPr, qn("a:buNone"))


def build_pptx(slides: list[Slide], path: Path, width_in, height_in):
    prs = Presentation()
    prs.slide_width, prs.slide_height = Inches(width_in), Inches(height_in)
    blank = prs.slide_layouts[6]
    for s in slides:
        slide = prs.slides.add_slide(blank)
        bg = slide.background.fill
        bg.solid()
        bg.fore_color.rgb = rgb(s.bg)
        for el in s.els:
            k = el["kind"]
            if k in ("rect", "oval", "arrow"):
                shape_type = (MSO_SHAPE.ROUNDED_RECTANGLE if k == "rect" and el.get("radius") else
                              MSO_SHAPE.RECTANGLE if k == "rect" else MSO_SHAPE.OVAL if k == "oval" else
                              MSO_SHAPE.RIGHT_ARROW)
                sh = slide.shapes.add_shape(shape_type, Inches(el["x"]), Inches(el["y"]), Inches(el["w"]), Inches(el["h"]))
                if k == "rect" and el.get("radius"):
                    sh.adjustments[0] = min(0.5, el["radius"] / min(el["w"], el["h"]))
                if el.get("fill"):
                    sh.fill.solid()
                    sh.fill.fore_color.rgb = rgb(el["fill"])
                else:
                    sh.fill.background()
                if el.get("line"):
                    sh.line.color.rgb = rgb(el["line"])
                    sh.line.width = Pt(el.get("lw", 0.75))
                else:
                    sh.line.fill.background()
                sh.shadow.inherit = False
            elif k == "text":
                size, lines, paras, ls = layout_text(el)
                tb = slide.shapes.add_textbox(Inches(el["x"]), Inches(el["y"]), Inches(el["w"]), Inches(el["h"]))
                tf = tb.text_frame
                tf.word_wrap = True
                tf.auto_size = None
                pad = Inches(el.get("pad", 0))
                tf.margin_left = tf.margin_right = tf.margin_top = tf.margin_bottom = pad
                tf.vertical_anchor = {"t": MSO_ANCHOR.TOP, "m": MSO_ANCHOR.MIDDLE, "b": MSO_ANCHOR.BOTTOM}[el["valign"]]
                for i, para in enumerate(paras):
                    p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
                    p.alignment = {"l": PP_ALIGN.LEFT, "c": PP_ALIGN.CENTER, "r": PP_ALIGN.RIGHT}[el["align"]]
                    p.line_spacing = ls / 1.2
                    after = para["after"] if para["after"] is not None else el.get("para_after", 0)
                    p.space_after = Pt(after)
                    _set_bullet(p, para["bullet"])
                    for text, style in para["runs"]:
                        r = p.add_run()
                        r.text = text
                        r.font.name = FONT
                        r.font.size = Pt(size)
                        r.font.bold = bool(style["bold"])
                        r.font.color.rgb = rgb(style["color"] or el["color"])
            elif k == "image":
                x, y, w, h = image_box(el)
                slide.shapes.add_picture(el["path"], Inches(x), Inches(y), Inches(w), Inches(h))
            elif k == "chart":
                cd = CategoryChartData()
                cd.categories = el["categories"]
                cd.add_series("value", el["values"])
                kind = XL_CHART_TYPE.BAR_CLUSTERED if el["horizontal"] else XL_CHART_TYPE.COLUMN_CLUSTERED
                gf = slide.shapes.add_chart(kind, Inches(el["x"]), Inches(el["y"]), Inches(el["w"]), Inches(el["h"]), cd)
                ch = gf.chart
                ch.has_legend = False
                ch.has_title = bool(el["title"])
                if el["title"]:
                    ch.chart_title.text_frame.text = el["title"]
                    run = ch.chart_title.text_frame.paragraphs[0].runs[0]
                    run.font.size, run.font.bold, run.font.name = Pt(12), True, FONT
                plot = ch.plots[0]
                plot.gap_width = 60
                plot.has_data_labels = True
                dl = plot.data_labels
                dl.number_format, dl.number_format_is_linked = el["fmt"], False
                dl.font.size, dl.font.name = Pt(12), FONT
                dl.position = XL_LABEL_POSITION.OUTSIDE_END
                series = plot.series[0]
                series.format.fill.solid()
                series.format.fill.fore_color.rgb = rgb(el["color"])
                if el["point_colors"]:
                    for idx, col in enumerate(el["point_colors"]):
                        pt = series.points[idx]
                        pt.format.fill.solid()
                        pt.format.fill.fore_color.rgb = rgb(col)
                va = ch.value_axis
                va.has_major_gridlines = True
                va.major_gridlines.format.line.color.rgb = rgb("DDE7E3")
                va.tick_labels.font.size, va.tick_labels.font.name = Pt(11), FONT
                if el["vmin"] is not None:
                    va.minimum_scale = el["vmin"]
                if el["vmax"] is not None:
                    va.maximum_scale = el["vmax"]
                ca = ch.category_axis
                ca.tick_labels.font.size, ca.tick_labels.font.name = Pt(12), FONT
                ca.tick_label_position = XL_TICK_LABEL_POSITION.LOW  # keep category names below negative bars and their labels
        if s.notes:
            slide.notes_slide.notes_text_frame.text = s.notes
    props = prs.core_properties
    props.author, props.last_modified_by, props.title = "", "", "AgriPAM-AI"
    props.comments, props.keywords = "", ""
    prs.save(str(path))


# ---------------------------------------------------------------- PNG preview
def render_preview(s: Slide, path: Path, px_per_in=110):
    W, H = int(s.w * px_per_in), int(s.h * px_per_in)
    img = Image.new("RGB", (W, H), "#" + s.bg)
    d = ImageDraw.Draw(img)
    k = px_per_in

    def col(c):
        return "#" + c if c else None

    for el in s.els:
        kind = el["kind"]
        box = [el["x"] * k, el["y"] * k, (el["x"] + el["w"]) * k, (el["y"] + el["h"]) * k]
        if kind == "rect":
            if el.get("radius"):
                d.rounded_rectangle(box, radius=el["radius"] * k, fill=col(el["fill"]), outline=col(el["line"]), width=1)
            else:
                d.rectangle(box, fill=col(el["fill"]), outline=col(el["line"]))
        elif kind == "oval":
            d.ellipse(box, fill=col(el["fill"]), outline=col(el["line"]))
        elif kind == "arrow":
            x0, y0, x1, y1 = box
            mid = (y0 + y1) / 2
            d.polygon([(x0, y0 + (y1 - y0) * .25), (x0 + (x1 - x0) * .6, y0 + (y1 - y0) * .25), (x0 + (x1 - x0) * .6, y0),
                       (x1, mid), (x0 + (x1 - x0) * .6, y1), (x0 + (x1 - x0) * .6, y0 + (y1 - y0) * .75),
                       (x0, y0 + (y1 - y0) * .75)], fill=col(el["fill"]))
        elif kind == "image":
            x, y, w, h = image_box(el)
            with Image.open(el["path"]) as im:
                im = im.convert("RGBA").resize((max(1, int(w * k)), max(1, int(h * k))))
                img.paste(im, (int(x * k), int(y * k)), im)
        elif kind == "chart":
            vals = el["values"]
            lo = el["vmin"] if el["vmin"] is not None else min(0, min(vals))
            hi = el["vmax"] if el["vmax"] is not None else max(vals) * 1.15
            x0, y0, x1, y1 = box
            top, bottom, left = y0 + 30, y1 - 40, x0 + 45
            d.rectangle(box, outline="#DDE7E3")
            n = len(vals)
            bw = (x1 - left - 20) / n
            zero_y = bottom - (0 - lo) / (hi - lo) * (bottom - top)
            d.line([(left, zero_y), (x1 - 10, zero_y)], fill="#999999")
            for i, v in enumerate(vals):
                bx = left + i * bw + bw * 0.2
                vy = bottom - (v - lo) / (hi - lo) * (bottom - top)
                colr = (el["point_colors"][i] if el["point_colors"] else el["color"])
                d.rectangle([bx, min(vy, zero_y), bx + bw * 0.6, max(vy, zero_y)], fill="#" + colr)
                label = format(v, ".3f" if el["fmt"] == "0.000" else ".2f" if el["fmt"] == "0.00" else ".1f")
                fpx = ImageFont.truetype(str(FONT_DIR / "Arial.ttf"), max(8, int(10 * k / 72 * 1.1)))
                d.text((bx, min(vy, zero_y) - 16), label, fill="#17322E", font=fpx)
                for li, part in enumerate(el["categories"][i].split("\n")):
                    d.text((bx - 4, bottom + 6 + li * 14), part[:22], fill="#17322E", font=fpx)
            if el["title"]:
                d.text((x0 + 8, y0 + 6), el["title"], fill="#17322E", font=ImageFont.truetype(str(FONT_DIR / "Arial Bold.ttf"), max(8, int(12 * k / 72 * 1.1))))
        elif kind == "text":
            size, lines, paras, ls = layout_text(el)
            f_px = size * k / 72.0
            line_h = size * ls * k / 72.0
            total = len(lines) * line_h + sum((p["after"] if p["after"] is not None else el.get("para_after", 0)) for p in paras[:-1]) * k / 72.0
            pad = el.get("pad", 0) * k
            avail_h = el["h"] * k - 2 * pad
            y = el["y"] * k + pad + {"t": 0, "m": (avail_h - total) / 2, "b": avail_h - total}[el["valign"]]
            width = el["w"] * k - 2 * pad
            prev_pi = None
            for pi, words, first in lines:
                if prev_pi is not None and pi != prev_pi:
                    y += (paras[prev_pi]["after"] if paras[prev_pi]["after"] is not None else el.get("para_after", 0)) * k / 72.0
                prev_pi = pi
                bullet = paras[pi]["bullet"]
                indent = 14 * k / 72.0 if bullet else 0
                text_w = sum(measure(w, size, b) for w, b, _ in words) * k / 72.0 + measure(" ", size, False) * k / 72.0 * (len(words) - 1)
                x = el["x"] * k + pad + indent
                if el["align"] == "c":
                    x = el["x"] * k + pad + (width - text_w) / 2
                elif el["align"] == "r":
                    x = el["x"] * k + pad + width - text_w
                if bullet and first:
                    d.text((el["x"] * k + pad, y), "•", fill="#" + (el["color"]), font=ImageFont.truetype(str(FONT_DIR / "Arial.ttf"), max(6, int(f_px))))
                for w, b, c in words:
                    fnt = ImageFont.truetype(str(FONT_DIR / ("Arial Bold.ttf" if b else "Arial.ttf")), max(6, int(round(f_px))))
                    d.text((x, y), w, fill="#" + (c or el["color"]), font=fnt)
                    x += fnt.getlength(w) + measure(" ", size, False) * k / 72.0
                y += line_h
    img.save(path)
    return [e for e in s.els if e.get("overflow")]
