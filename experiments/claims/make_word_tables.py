"""Write the two main-text tables (spatial claims; C4 per marked event, Form B) as Word tables.

Reads out/claims_table.csv and out/table_C4_events.csv (written by claims_test.ipynb) and writes
out/tables_main.docx. Run with the base miniconda python (python-docx is not in gmtenv):
    ~/miniconda3/bin/python3 make_word_tables.py
Plain tables: black text, no colours, bold only in the header row. Numbers come from the CSVs only.
"""
from pathlib import Path

import pandas as pd
from docx import Document
from docx.enum.table import WD_TABLE_ALIGNMENT, WD_CELL_VERTICAL_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Pt, RGBColor, Cm

HERE = Path(__file__).resolve().parent
OUT = HERE / 'out'

PLAIN = True   # plain table (Yuxin, 2026-10-30): black text, no colour, bold only in the header row
BLUE = RED = GREY = INK = None
FONT, FS, FS_SMALL, FS_NOTE = 'Times New Roman', 9, 7.5, 8

REGIONS = ['LMS', 'XSH', 'WYN', 'XJ', 'NCP', 'NW']
REGION_NAMES = {'LMS': 'Longmenshan', 'XSH': 'Xianshuihe', 'WYN': 'Western Yunnan', 'XJ': 'Xiaojiang',
                'NCP': 'North China Plain', 'NW': 'Northwest China'}
TABLE_ROWS = [  # (key, claim id, text, reference text)
    ('C12_rate', 'C1–C2', 'Low b with high strain rate, high b with low strain rate', 'chance 50%'),
    ('C1_regime', 'C1', 'Contractional nodes that are low\u2011b', 'share of low\u2011b nodes'),      # non-breaking hyphens
    ('C2_regime', 'C2', 'Extensional nodes that are high\u2011b', 'share of high\u2011b nodes'),
    ('C3_hf', 'C3', 'High\u2011heat\u2011flow nodes that are high\u2011b', 'share of high\u2011b nodes'),
]
WINDOW_NAMES = {'LMS_2006-2010': 'Longmenshan 2006–2010', 'LMS_2011-2016': 'Longmenshan 2011–2016',
                'XSH_2022-2023': 'Xianshuihe 2022–2023', 'WYN_2020-2023': 'Western Yunnan 2020–2023',
                'XJ_2016-2023': 'Xiaojiang 2016–2023', 'NCP_2008-2023': 'North China Plain 2008–2023',
                'NW_2006-2023': 'Northwest China 2006–2023'}
FILLED, OPEN, DOT = '●', '○', '·'


# ------------------------------------------------------------------ helpers
def run(par, text, size=FS, bold=False, color=None, italic=False, font=FONT):
    r = par.add_run(text)
    r.font.name, r.font.size, r.bold, r.italic = font, Pt(size), bold, italic
    r._element.rPr.rFonts.set(qn('w:eastAsia'), font)
    if color is not None:
        r.font.color.rgb = color
    return r


def set_cell_borders(cell, **edges):
    """edges: top/bottom/left/right = (sz in eighths of a point, 'single'|'nil')."""
    tcPr = cell._tc.get_or_add_tcPr()
    borders = tcPr.find(qn('w:tcBorders'))
    if borders is None:
        borders = OxmlElement('w:tcBorders')
        tcPr.append(borders)
    for edge, (sz, val) in edges.items():
        el = borders.find(qn(f'w:{edge}'))
        if el is None:
            el = OxmlElement(f'w:{edge}')
            borders.append(el)
        el.set(qn('w:val'), val)
        el.set(qn('w:sz'), str(sz))
        el.set(qn('w:space'), '0')
        el.set(qn('w:color'), '000000')


def no_borders(table):
    tblPr = table._tbl.tblPr
    b = OxmlElement('w:tblBorders')
    for edge in ('top', 'left', 'bottom', 'right', 'insideH', 'insideV'):
        el = OxmlElement(f'w:{edge}')
        el.set(qn('w:val'), 'nil')
        b.append(el)
    tblPr.append(b)


def book_rules(table, header_rows=1):
    """Three-line (booktabs) look: rule above the header, below the header, below the last row."""
    n = len(table.rows)
    for c in table.rows[0].cells:
        set_cell_borders(c, top=(8, 'single'))
    for c in table.rows[header_rows - 1].cells:
        set_cell_borders(c, bottom=(4, 'single'))
    for c in table.rows[n - 1].cells:
        set_cell_borders(c, bottom=(8, 'single'))


def cell_par(cell, align=WD_ALIGN_PARAGRAPH.CENTER, first=True):
    p = cell.paragraphs[0] if first else cell.add_paragraph()
    p.alignment = align
    p.paragraph_format.space_before = Pt(0)
    p.paragraph_format.space_after = Pt(0)
    p.paragraph_format.line_spacing = 1.0
    return p


def tight_cells(table, pad_pt=1.5):
    for row in table.rows:
        for c in row.cells:
            c.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER
            tcPr = c._tc.get_or_add_tcPr()
            mar = OxmlElement('w:tcMar')
            for side in ('top', 'left', 'bottom', 'right'):      # schema order
                el = OxmlElement(f'w:{side}')
                el.set(qn('w:w'), str(int(pad_pt * 20)) if side in ('top', 'bottom') else '25')
                el.set(qn('w:type'), 'dxa')
                mar.append(el)
            tcPr.append(mar)


def set_col_widths(table, widths_cm):
    """Fixed layout with the widths in the table grid (Word and Quick Look both honour tblGrid; cell.width alone is ignored)."""
    table.autofit = False
    tblPr = table._tbl.tblPr
    if tblPr.find(qn('w:tblLayout')) is None:      # autofit=False already adds one
        lay = OxmlElement('w:tblLayout')
        lay.set(qn('w:type'), 'fixed')
        tblPr.append(lay)
    tw = OxmlElement('w:tblW')
    tw.set(qn('w:w'), str(int(sum(widths_cm) / 2.54 * 1440)))
    tw.set(qn('w:type'), 'dxa')
    for old in tblPr.findall(qn('w:tblW')):
        tblPr.remove(old)
    tblPr.append(tw)
    grid = table._tbl.tblGrid
    for gc, w in zip(grid.findall(qn('w:gridCol')), widths_cm):
        gc.set(qn('w:w'), str(int(w / 2.54 * 1440)))
    for row in table.rows:
        for c, w in zip(row.cells, widths_cm):
            c.width = Cm(w)


def caption(doc, label, text):
    p = doc.add_paragraph()
    p.paragraph_format.space_after = Pt(3)
    run(p, label + ' ', bold=True)
    run(p, text)
    return p


def note(doc, parts):
    """parts: list of (text, style) with style in {'', 'b', 'i', 'blue', 'red'}."""
    p = doc.add_paragraph()
    p.paragraph_format.space_before = Pt(2)
    p.paragraph_format.space_after = Pt(10)
    for text, style in parts:
        run(p, text, size=FS_NOTE, bold=style == 'b', italic=style == 'i',
            color=None)
    return p


def pct_pair(par, pc, pp, size=FS, bold=False):
    """'35 | 46' (classical | positive)."""
    run(par, '—' if pd.isna(pc) else f'{pc:.0f}%', size=size, bold=bold)
    run(par, ' | ', size=size)
    run(par, '—' if pd.isna(pp) else f'{pp:.0f}%', size=size, bold=bold)


# ------------------------------------------------------------------ Table 1: spatial claims
def table_spatial(doc, tab):
    tab = tab[tab.row.isin([k for k, *_ in TABLE_ROWS])].copy()
    caption(doc, 'Table X.',
            'Share of grid nodes supporting the spatial claims C1–C3 in the 2008–2023 b-value maps of the six sub-regions. '
            'Each cell gives the percentage for b-classical | b-positive and the number of nodes tested (n).')
    t = doc.add_table(rows=1 + len(TABLE_ROWS), cols=2 + len(REGIONS))
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    no_borders(t)
    hdr = t.rows[0].cells
    run(cell_par(hdr[0], WD_ALIGN_PARAGRAPH.LEFT), 'Claim', bold=True)
    run(cell_par(hdr[1], WD_ALIGN_PARAGRAPH.LEFT), 'Nodes counted (% of nodes tested)', bold=True)
    for j, reg in enumerate(REGIONS):
        run(cell_par(hdr[2 + j]), REGION_NAMES[reg], bold=True)
    for i, (key, cid, text, reftxt) in enumerate(TABLE_ROWS):
        cells = t.rows[1 + i].cells
        run(cell_par(cells[0], WD_ALIGN_PARAGRAPH.LEFT), cid)
        p = cell_par(cells[1], WD_ALIGN_PARAGRAPH.LEFT)
        run(p, text)
        for j, reg in enumerate(REGIONS):
            c = tab[(tab.row == key) & (tab.region == reg) & (tab.estimator == 'classical')]
            pq = tab[(tab.row == key) & (tab.region == reg) & (tab.estimator == 'positive')]
            par = cell_par(cells[2 + j])
            if c.empty or pd.isna(c.pct.iloc[0]):
                run(par, 'n/a')
                continue
            c, pq = c.iloc[0], pq.iloc[0]
            pct_pair(par, c.pct, pq.pct)
            run(cell_par(cells[2 + j], first=False), f'n {int(c.n):,}', size=FS_SMALL)
    tight_cells(t)
    book_rules(t)
    set_col_widths(t, [1.1, 4.3] + [2.05] * len(REGIONS))
    shares = ';  '.join(f"{REGION_NAMES[r]} {tab[(tab.row=='C1_regime')&(tab.region==r)&(tab.estimator=='classical')].share_low.iloc[0]:.0f}/"
                        f"{tab[(tab.row=='C1_regime')&(tab.region==r)&(tab.estimator=='classical')].share_high.iloc[0]:.0f} | "
                        f"{tab[(tab.row=='C1_regime')&(tab.region==r)&(tab.estimator=='positive')].share_low.iloc[0]:.0f}/"
                        f"{tab[(tab.row=='C1_regime')&(tab.region==r)&(tab.estimator=='positive')].share_high.iloc[0]:.0f}"
                        for r in REGIONS)
    note(doc, [
        ('Note. ', 'b'),
        ('Low-b / high-b node = b below / above the whole-period b of the same region and estimator (classes relative to the '
         'regional average, Schorlemmer et al., 2005). High strain rate = second invariant above the region’s node median; '
         'contractional / extensional = lower / upper tercile of the dilatation rate; high heat flow = above the median of the '
         'covered nodes. For C1–C2, n counts all nodes with strain data and 50% is the share expected with no relation because the '
         'strain-rate split is at the median; for C1, C2 and C3 the share expected with no relation is the share of low-b or '
         'high-b nodes in the region (low/high, classical | positive: ' + shares + '). Percentages are area fractions; '
         'neighbourhoods overlap. '
         'Xianshuihe has no heat-flow site (n/a).', ''),
    ])


# ------------------------------------------------------------------ Table 2: C4 per marked event (Form B)
def table_c4(doc, long):
    wins = list(WINDOW_NAMES)
    n_ev = int(long.groupby('window').size().max())
    caption(doc, 'Table Y.',
            'Temporal claim C4: b-value drop before each marked event of Figs 4–10. One row per window; events in time '
            'order (magnitude; * = aftershock of an earlier marked event). Symbols: b-classical | b-positive; '
            '● = last window before the event at least 10% below the pre-sequence reference, ○ = no drop, '
            '· = no valid window. Last column: events with a drop / events tested.')
    t = doc.add_table(rows=1 + len(wins), cols=2 + n_ev)
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    no_borders(t)
    hdr = t.rows[0].cells
    run(cell_par(hdr[0], WD_ALIGN_PARAGRAPH.LEFT), 'Window', bold=True)
    for k in range(n_ev):
        run(cell_par(hdr[1 + k]), str(k + 1), bold=True)
    p = cell_par(hdr[1 + n_ev])
    run(p, 'Drop c | p', bold=True)
    for i, lab in enumerate(wins):
        g = long[long.window == lab].sort_values('idx')
        cells = t.rows[1 + i].cells
        p = cell_par(cells[0], WD_ALIGN_PARAGRAPH.LEFT)
        run(p, f'{g.fig.iloc[0]}  {WINDOW_NAMES[lab]}')
        for k, (_, e) in enumerate(g.iterrows()):
            par = cell_par(cells[1 + k])
            run(par, f'{e.mag:.1f}' + ('' if e.initiator else '*'), size=FS)
            p2 = cell_par(cells[1 + k], first=False)
            for v in (e.drop_c, e.drop_p):
                run(p2, DOT if pd.isna(v) else (FILLED if bool(v) else OPEN), size=FS + 1)
        par = cell_par(cells[1 + n_ev])
        kc, nc = int(g.drop_c.dropna().astype(bool).sum()), int(g.drop_c.notna().sum())
        kp, np_ = int(g.drop_p.dropna().astype(bool).sum()), int(g.drop_p.notna().sum())
        run(par, f'{kc}/{nc} | {kp}/{np_}')
    tight_cells(t)
    book_rules(t)
    set_col_widths(t, [4.9] + [1.05] * n_ev + [1.7])
    rise = ';  '.join(f"{long[long.window==l].fig.iloc[0]} "
                      f"{int(long[long.window==l].rise_c.dropna().astype(bool).sum())}/{int(long[long.window==l].rise_c.notna().sum())} | "
                      f"{int(long[long.window==l].rise_p.dropna().astype(bool).sum())}/{int(long[long.window==l].rise_p.notna().sum())}"
                      for l in wins)
    note(doc, [
        ('Note. ', 'b'),
        ('Marked events = all events with M ≥ the fifth-largest magnitude of the window (Section 2.2). Drop = last fixed-N window '
         'ending before the event ≥ 10% below the reference, the median of the earlier windows of the panel that share no event '
         'with it and end outside the Gardner–Knopoff window of an earlier marked event (Gulia & Wiemer, 2019); windows and Mc as '
         'in Figs 4–10, classical windows re-cut at the per-window MAXC + 0.2. Rise after the event (first post-event window ≥ 10% '
         'above the same reference), classical | positive: ' + rise + '.', ''),
    ])


def main():
    tab = pd.read_csv(OUT / 'claims_table.csv')
    long = pd.read_csv(OUT / 'table_C4_events.csv')
    doc = Document()
    sec = doc.sections[0]
    sec.page_width, sec.page_height = Cm(21.0), Cm(29.7)
    sec.left_margin = sec.right_margin = Cm(1.6)
    sec.top_margin = sec.bottom_margin = Cm(2.0)
    st = doc.styles['Normal']
    st.font.name, st.font.size = FONT, Pt(FS)
    st.element.rPr.rFonts.set(qn('w:eastAsia'), FONT)
    table_spatial(doc, tab)
    table_c4(doc, long)
    path = OUT / 'tables_main.docx'
    doc.save(path)
    print('wrote', path)


if __name__ == '__main__':
    main()
